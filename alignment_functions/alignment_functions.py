from Bio.PDB import MMCIFParser
import numpy as np
from scipy.spatial.transform import Rotation as R
import sys
import json
import os
from tqdm import tqdm
import pandas as pd
from tmtools import tm_align
from Bio.Data import IUPACData


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




def get_ca_coords_and_plddt_from_cif(cif_file):
    """
    Extracts Cα atom coordinates, pLDDT scores, and sequence from a .cif file.
    
    Returns:
        coords (np.ndarray): N x 3 array of Cα coordinates.
        seq (str): Amino acid sequence (1-letter codes) matching the coordinates.
        plddts (np.ndarray): N-length array of pLDDT scores for Cα atoms.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", cif_file)

    coords = []
    seq = []
    plddts = []

    # Map 3-letter residue names to 1-letter codes
    res3_to1 = {k.upper(): v for k, v in IUPACData.protein_letters_3to1.items()}

    atom_idx = 0  # Index in the mmCIF atom list (used to get pLDDT from atom records)
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().upper()
                if resname not in res3_to1:
                    continue  # skip unknown residues

                for atom in residue:
                    if atom.get_id() == "CA":
                        coords.append(atom.get_coord())
                        seq.append(res3_to1[resname])

                        # Try to get B-factor as pLDDT (AlphaFold stores pLDDT in B-factor field)
                        plddt = atom.get_bfactor()
                        plddts.append(plddt)
                        break

    return np.array(coords), ''.join(seq), np.array(plddts)

def apply_alignment(coords, R, t):
    """Apply TM-align's optimal rotation and translation."""
    return coords @ R.T + t

def weighted_rmsd(P, Q, weights):
    """
    Computes confidence-weighted RMSD between two aligned point sets P and Q.

    Args:
        P (np.ndarray): Nx3 array of coordinates (e.g., aligned Cα from structure 1)
        Q (np.ndarray): Nx3 array of coordinates (e.g., reference Cα from structure 2)
        weights (np.ndarray): N-length array of weights (e.g., average pLDDT per residue, scaled to [0, 1])

    Returns:
        float: Weighted RMSD
    """
    assert P.shape == Q.shape
    assert P.shape[0] == weights.shape[0]

    diffs = P - Q
    squared_diffs = np.sum(diffs ** 2, axis=1)
    weighted_sum = np.sum(weights * squared_diffs)
    norm = np.sum(weights)

    return np.sqrt(weighted_sum / norm)


# coords1, seq1, confs1 = get_ca_coords_and_plddt_from_cif("/shared/25mdl4/af_output/prion_ref/seed-2_sample-0/model.cif")
# coords2, seq2, confs2 = get_ca_coords_and_plddt_from_cif("/shared/25mdl4/af_output/prion_m129v/seed-2_sample-0/model.cif")

# # # print(coords.shape, seq)  
# res = tm_align(coords1, coords2, seq1, seq2)


# coords1_aligned_to_2 = apply_alignment(coords1, res.R, res.t)


# print(res.distances)
# print(res.tm_norm_chain1, res.tm_norm_chain2, res.rmsd)
# print(confs1.shape)
# print(confs1)



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

# note we should just reweight by 1 if we dont want to weight
def weighted_kabsch_numpy(P, Q, conf_weights):
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
    
    alignments = np.square(np.dot(p, R.T) - q)
    weighted_alignments = alignments * conf_weights[:, np.newaxis]
    # RMSD
    rmsd = np.sqrt(np.sum(weighted_alignments) / P.shape[0])

    return rmsd

def findPandQ_averaged(cif1, cif2, json1, json2):
    """Compares two structures by calculating distances and extracting pLDDT scores."""
    residues2cords1, residues2plddt1 = extract_plddt_and_coords(cif1, json1)
    residues2cords2, residues2plddt2 = extract_plddt_and_coords(cif2, json2)
    
    if len(residues2cords1) != len(residues2cords2):
        return None

    P, Q, avg_conf = [], [], []
    ctr = 0
    
    for residue_id in residues2cords1.keys():
        if residue_id in residues2cords2:
            # print(residue_id)
            residue_coords1 = np.array(residues2cords1[residue_id])
            residue_coords2 = np.array(residues2cords2[residue_id])
            residue_plddt1 = np.array(residues2plddt1[residue_id])
            residue_plddt2 = np.array(residues2plddt2[residue_id])
            
            # Check if both residues have the same number of atoms
            if len(residue_coords1) == len(residue_coords2):
                for c1, c2, p1, p2 in zip(residue_coords1, residue_coords2, residue_plddt1, residue_plddt2):

                    P.append(c1)
                    Q.append(c2)
                    avg_conf.append((p1 + p2) / 200)   
            ctr += 1
            
    return np.array(P), np.array(Q), np.array(avg_conf)

def findPandQ_refonly(cif1, cif2, json1, json2):
    """Compares two structures by calculating distances and extracting pLDDT scores."""
    residues2cords1, residues2plddt1 = extract_plddt_and_coords(cif1, json1)
    residues2cords2, residues2plddt2 = extract_plddt_and_coords(cif2, json2)
    
    if len(residues2cords1) != len(residues2cords2):
        return None

    P, Q, avg_conf = [], [], []
    ctr = 0
    
    for residue_id in residues2cords1.keys():
        if residue_id in residues2cords2:
            # print(residue_id)
            residue_coords1 = np.array(residues2cords1[residue_id])
            residue_coords2 = np.array(residues2cords2[residue_id])
            residue_plddt1 = np.array(residues2plddt1[residue_id])
            residue_plddt2 = np.array(residues2plddt2[residue_id])
            
            # Check if both residues have the same number of atoms
            if len(residue_coords1) == len(residue_coords2):
                for c1, c2, p1, p2 in zip(residue_coords1, residue_coords2, residue_plddt1, residue_plddt2):

                    P.append(c1)
                    Q.append(c2)
                    avg_conf.append(p1 / 100)   
            ctr += 1
            
    return np.array(P), np.array(Q), np.array(avg_conf)

#make sure first is ref
def rmsd_align(protein1, protein2, json1, json2, PandQfinder, confidence_weighted=True, reweighted=False):

    res = PandQfinder(protein1, protein2, json1, json2)
    
    if res is None: 
        return None
    P, Q, avg_conf = res
    print(P.shape, Q.shape)

    if confidence_weighted:
        conf_arr = avg_conf
    else:
        conf_arr = np.ones_like(avg_conf)
        
    rmsd = weighted_kabsch_numpy(P, Q, conf_arr)
    
    
    if reweighted:
        rmsd = rmsd * (1 / np.mean(conf_arr))

    return rmsd

# currently looks only at CA plddts
# need to play with reweighting
# should try tm_score but idk how reweighting works
def tm_align_wrapper(cif1, cif2):
    
    coords1, seq1, confs1 = get_ca_coords_and_plddt_from_cif(cif1)
    coords2, seq2, confs2 = get_ca_coords_and_plddt_from_cif(cif2)

    res = tm_align(coords1, coords2, seq1, seq2)


    coords1_aligned_to_2 = apply_alignment(coords1, res.u, res.t)
    
    score = weighted_rmsd(coords1_aligned_to_2, coords2, (confs1 + confs2) / 2)

    return score


# this wrapper sucks 
# def tm_align_wrapper(cif1, cif2):
    
#     coords1, seq1, confs1 = get_ca_coords_and_plddt_from_cif(cif1)
#     coords2, seq2, confs2 = get_ca_coords_and_plddt_from_cif(cif2)

#     res = tm_align(coords1, coords2, seq1, seq2)

#     score = 1 - (res.tm_norm_chain1 + res.tm_norm_chain2) / 2    
    
#     # coords1_aligned_to_2 = apply_alignment(coords1, res.u, res.t)
    
#     # score = weighted_rmsd(coords1_aligned_to_2, coords2, (confs1 + confs2) / 2)

#     return score
    
    

# coords1, seq1, confs1 = get_ca_coords_and_plddt_from_cif("/shared/25mdl4/af_output/prion_ref/seed-2_sample-0/model.cif")
# coords2, seq2, confs2 = get_ca_coords_and_plddt_from_cif("/shared/25mdl4/af_output/prion_m129v/seed-2_sample-0/model.cif")
# tm_align_wrapper("/shared/25mdl4/af_output/prion_ref/seed-2_sample-0/model.cif", "/shared/25mdl4/af_output/prion_a133v/seed-2_sample-0/model.cif")


def align_all(ref_id, seq_ids, protein_name, alignment_function, PandQfinder, confidence_weighted=True, reweighted=False, log_fn=print):
    
    out_root = "/shared/25mdl4/af_output/"
    
    # ID,Old,White,Unhealthy,Align Score,Cancer,Sequence,is_pathogenic
    # ids = data["ID"]
    
    scores = []
    
    # name = name.replace(":", "_").replace(">", "_").lower()
    
    for i, idx in enumerate(seq_ids):
        
        if i % (len(seq_ids) // 10) == 0:
            log_fn(f"{(i / len(seq_ids)) * 100:.0f}% complete")
            
        folder1 = os.path.join(out_root, f"{protein_name}_{ref_id}")
        folder2 = os.path.join(out_root, f"{protein_name}_{idx}")
        
        if not os.path.isdir(folder2):
            log_fn(f"WARNING: {idx} invalid path -- skipping" )
            continue
        
        cif_file1 = f'{folder1}/seed-2_sample-0/model.cif'
        cif_file2 = f'{folder2}/seed-2_sample-0/model.cif'
        json_file1 = f'{folder1}/seed-2_sample-0/confidences.json'
        json_file2 = f'{folder2}/seed-2_sample-0/confidences.json'
        align_score = alignment_function(cif_file1, cif_file2, json_file1, json_file2, PandQfinder, confidence_weighted, reweighted)
        
        if align_score is None:
            log_fn(f"WARNING: nonsense mutation -- skipping")
            continue
        
        scores.append(align_score)
    
    return scores

# this currently still uses rmsd but with tm align
def tm_align_all(ref_id, seq_ids, protein_name, log_fn=print):
    
    out_root = "/shared/25mdl4/af_output/"
    
    # ID,Old,White,Unhealthy,Align Score,Cancer,Sequence,is_pathogenic
    # ids = data["ID"]
    
    scores = []
    
    # name = name.replace(":", "_").replace(">", "_").lower()
    
    for i, idx in enumerate(seq_ids):
        
        if i % (len(seq_ids) // 10) == 0:
            log_fn(f"{(i / len(seq_ids)) * 100:.0f}% complete")
            
        folder1 = os.path.join(out_root, f"{protein_name}_{ref_id}")
        folder2 = os.path.join(out_root, f"{protein_name}_{idx}")
        
        if not os.path.isdir(folder2):
            log_fn(f"WARNING: {idx} invalid path -- skipping" )
            continue
        
        cif_file1 = f'{folder1}/seed-2_sample-0/model.cif'
        cif_file2 = f'{folder2}/seed-2_sample-0/model.cif'
        align_score = tm_align_wrapper(cif_file1, cif_file2)
        
        if align_score is None:
            log_fn(f"WARNING: nonsense mutation -- skipping")
            continue
        
        scores.append(align_score)
    
    return scores