from Bio.PDB import MMCIFParser
import numpy as np
from scipy.spatial.transform import Rotation as R
import sys
import json
import os
from tqdm import tqdm
import pandas as pd

# TODO weighted avg

def get_atom_coords_by_residue_from_cif(cif_file):
    """Extracts the coordinates of all atoms, grouped by residue, from a .cif file."""
    parser = MMCIFParser()
    structure = parser.get_structure('protein', cif_file)
    
    residue_atom_coords = {}
    
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_id = residue.id[1]  # Use chain ID and residue number as the key
                if residue_id not in residue_atom_coords:
                    residue_atom_coords[residue_id] = []
                
                for atom in residue:
                    residue_atom_coords[residue_id].append(list(atom.get_coord()))
    

    return residue_atom_coords

# folder1 = "fold_prion_tv3"
# cif_file1 = f'{folder1}/{folder1}_model_0.cif'
# json_file1 = f'{folder1}/{folder1}_full_data_0.json'  # Replace with your first JSON file


def extract_plddt_scores(json_file):
    """Extracts pLDDT scores for each atom from the AlphaFold JSON file."""
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    return data.get('atom_plddts', [])

def extract_plddt_and_coords(cif_file, json_file):
    residues2cords = get_atom_coords_by_residue_from_cif(cif_file)
    plddt_list = extract_plddt_scores(json_file)
    # print(residues2cords[1], plddt_list[1])
    # print(len(plddt_list), len(residues2cords))
    
    # print(type(plddt_list), type(residues2cords))
    strong_cords = []
    residues2plddt = {}
    ind = 0
    for k, v in residues2cords.items():
        residues2plddt[k] = []

        for _ in v:
            # if plddt_list[ind] > thresh:
            #     strong_cords.append(coord)
            # else:
            #     strong_cords.append([0, 0, 0])
            residues2plddt[k].append(plddt_list[ind])
            
            ind += 1
            
    # return np.array(strong_cords)
    return residues2cords, residues2plddt


def kabsch_numpy(P, Q):
    """
    Computes the optimal rotation and translation to align two sets of points (P -> Q),
    and their RMSD.

    :param P: A Nx3 matrix of points
    :param Q: A Nx3 matrix of points
    :return: A tuple containing the optimal rotation matrix, the optimal
             translation vector, and the RMSD.
    """
    
    # print(P.shape, Q.shape)
    assert P.shape == Q.shape, "Matrix dimensions must match"

    # Compute centroids
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    # Optimal translation
    t = centroid_Q - centroid_P

    # Center the points
    p = P - centroid_P
    q = Q - centroid_Q

    # Compute the covariance matrix
    H = np.dot(p.T, q)

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Validate right-handed coordinate system
    if np.linalg.det(np.dot(Vt.T, U.T)) < 0.0:
        Vt[-1, :] *= -1.0

    # Optimal rotation
    R = np.dot(Vt.T, U.T)

    # RMSD
    rmsd = np.sqrt(np.sum(np.square(np.dot(p, R.T) - q)) / P.shape[0])

    return R, t, rmsd

def weighted_kabsch_numpy(P, Q, avg_conf):
    """
    Computes the optimal rotation and translation to align two sets of points (P -> Q),
    and their RMSD.

    :param P: A Nx3 matrix of points
    :param Q: A Nx3 matrix of points
    :return: A tuple containing the optimal rotation matrix, the optimal
             translation vector, and the RMSD.
    """
    
    # print(P.shape, Q.shape)
    assert P.shape == Q.shape, "Matrix dimensions must match"

    # Compute centroids
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)


    #FIXME TODO we don't use the translation -- is this messing us up?
    # it may be ok because we center the points
    
    # Optimal translation

    t = centroid_Q - centroid_P

    # Center the points
    p = P - centroid_P
    q = Q - centroid_Q

    # Compute the covariance matrix
    H = np.dot(p.T, q)

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Validate right-handed coordinate system
    if np.linalg.det(np.dot(Vt.T, U.T)) < 0.0:
        Vt[-1, :] *= -1.0

    # Optimal rotation
    R = np.dot(Vt.T, U.T)

    # TODO Should we add translation here
    alignments = np.square(np.dot(p, R.T) - q)
    # print(alignments.shape, avg_conf.shape)

    weighted_alignments = alignments * avg_conf[:, np.newaxis]
    # RMSD
    rmsd = np.sqrt(np.sum(weighted_alignments) / P.shape[0])

    return R, t, rmsd

def findPandQ(cif1, cif2, json1, json2, thresh=70.0):
    """Compares two structures by calculating distances and extracting pLDDT scores."""
    residues2cords1, residues2plddt1 = extract_plddt_and_coords(cif1, json1)
    residues2cords2, residues2plddt2 = extract_plddt_and_coords(cif2, json2)
    
    if len(residues2cords1) != len(residues2cords2):
        return None

    # print(plddt1.keys())
    # conf_dists = []
    P, Q, avg_conf = [], [], []
    # print(len(coords1), len(coords2), len(plddt1), len(plddt2))
    ctr = 0
    
    for residue_id in residues2cords1.keys():
        if residue_id in residues2cords2:
            # print(residue_id)
            residue_coords1 = np.array(residues2cords1[residue_id])
            residue_coords2 = np.array(residues2cords2[residue_id])
            residue_plddt1 = np.array(residues2plddt1[residue_id])
            residue_plddt2 = np.array(residues2plddt2[residue_id])
            
            # print(atom_plddt1)
            
            # print("atom coords", atom_coords1, atom_coords2)
            # Check if both residues have the same number of atoms
            if len(residue_coords1) == len(residue_coords2):
                for c1, c2, p1, p2 in zip(residue_coords1, residue_coords2, residue_plddt1, residue_plddt2):
                    # print(ctr, (p1 + p2) / 2)
                    
                    # NOTE: CHANGED
                    # if (p1 + p2) / 2 > thresh:
                    P.append(c1)
                    Q.append(c2)
                    avg_conf.append((p1 + p2) / 200)   

                    # else:
                    #     # hmm, is this okay?
                    #     P.append([0,0,0])
                    #     Q.append([0,0,0])
                        
                        
                # if residue_plddt1[residue_id] > thresh
                # distances = calculate_distances(residue_coords1, residue_coords2)
                # print(distances)
                
            else:
                pass
                # print(f"Warning: Residue {residue_id} has different numbers of atoms between the two structures.")

            ctr += 1
            
        else:
            print("Error: we shouldnt be here")
            
    return np.array(P), np.array(Q), np.array(avg_conf)

def dual_weighted_rmsd_align(protein1, protein2, json1, json2, thresh=70.0):
    # P = get_alpha_carbons_from_cif(protein1_file)
    # P = extract_confident_atoms(protein1, json1)
    # Q = extract_confident_atoms(protein2, json2)
    
    # Q = get_alpha_carbons_from_cif(protein2_file)
    res = findPandQ(protein1, protein2, json1, json2, thresh)
    
    if res is None: 
        return None
    P, Q, avg_conf = res
    #FINDME CHanged
    rmsd1 = weighted_kabsch_numpy(P, Q, avg_conf)
    # rmsd1 = kabsch_numpy(P, Q, avg_conf)
    
    
    return rmsd1[2]


def align_all(ref_id, seq_ids, protein_name, alignment_function, log_fn):
    
    out_root = "/shared/25mdl4/af_output/"
    
    # ID,Old,White,Unhealthy,Align Score,Cancer,Sequence,is_pathogenic
    # ids = data["ID"]
    
    scores = []
    
    # name = name.replace(":", "_").replace(">", "_").lower()
    
    for i, idx in enumerate(seq_ids):
        
        if i % (len(seq_ids) // 10) == 0:
            log_fn(f"{(i / len(seq_ids)) * 100:.0f}% complete")
            
        folder1 = os.path.join(out_root, f"{protein_name}_{idx}")
        folder2 = os.path.join(out_root, f"{protein_name}_{ref_id}")
        
        if not os.path.isdir(folder1):
            log_fn(f"WARNING: {idx} invalid path -- skipping" )
            continue
        
        cif_file1 = f'{folder1}/seed-2_sample-0/model.cif'
        cif_file2 = f'{folder2}/seed-2_sample-0/model.cif'
        json_file1 = f'{folder1}/seed-2_sample-0/confidences.json'  # Replace with your first JSON file
        json_file2 = f'{folder2}/seed-2_sample-0/confidences.json'
        align_score = alignment_function(cif_file1, cif_file2, json_file1, json_file2)
        
        if align_score is None:
            log_fn(f"WARNING: nonsense mutation -- skipping")
            continue
        
        scores.append(align_score)
    
    return scores