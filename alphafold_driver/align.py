from Bio.PDB import MMCIFParser
import numpy as np
from scipy.spatial.transform import Rotation as R
import sys
import json
import os
from tqdm import tqdm
import pandas as pd

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


def findPandQ(cif1, cif2, json1, json2, thresh=70.0):
    """Compares two structures by calculating distances and extracting pLDDT scores."""
    residues2cords1, residues2plddt1 = extract_plddt_and_coords(cif1, json1)
    residues2cords2, residues2plddt2 = extract_plddt_and_coords(cif2, json2)
    
    if len(residues2cords1) != len(residues2cords2):
        return None

    # print(plddt1.keys())
    # conf_dists = []
    P, Q = [], []
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
                    if (p1 + p2) / 2 > thresh:
                        P.append(c1)
                        Q.append(c2)
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
            
    return np.array(P), np.array(Q)


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


# Main alignment routine
def tm_align_rmsd(protein1, protein2, json1, json2, thresh=70.0):
    # P = get_alpha_carbons_from_cif(protein1_file)
    # P = extract_confident_atoms(protein1, json1)
    # Q = extract_confident_atoms(protein2, json2)
    
    # Q = get_alpha_carbons_from_cif(protein2_file)
    res = findPandQ(protein1, protein2, json1, json2, thresh)
    
    if res is None: 
        return None
    P, Q = res
    rmsd1 = kabsch_numpy(P, Q)
    
    return rmsd1[2]
    # Calculate TM-align based RMSD
    # best_rmsd = calculate_tm_rmsd(P, Q, iterations)
    
    # print(f"Best RMSD after alignment: {best_rmsd:.4f}")


def augment_data():
    
    out_root = "/shared/25mdl4/af_output/"
    
    data = pd.read_csv("nice_data.csv")
    
    # ID,Old,White,Unhealthy,Align Score,Cancer,Sequence,is_pathogenic
    # ids = data["ID"]
    
    valid_rows = []
    
    # name = name.replace(":", "_").replace(">", "_").lower()
    for i, row in tqdm(data.iterrows()):
        
        name = row["ID"].replace(":", "_").replace(">", "_").lower()
        
        folder1 = os.path.join(out_root, name)
        folder2 = os.path.join(out_root, "reference")
        
        
        if not os.path.isdir(folder1):
            print(f"WARNING: {name} invalid path -- skipping" )
            continue
        
        
        cif_file1 = f'{folder1}/seed-2_sample-0/model.cif'
        cif_file2 = f'{folder2}/seed-2_sample-0/model.cif'
        json_file1 = f'{folder1}/seed-2_sample-0/confidences.json'  # Replace with your first JSON file
        json_file2 = f'{folder2}/seed-2_sample-0/confidences.json'
        align_score = tm_align_rmsd(cif_file1, cif_file2, json_file1, json_file2)
        
        if align_score is None:
            print(f"WARNING: nonsense mutation -- skipping")
            continue
            
        row["Align Score"] = align_score
        
        valid_rows.append(row)
        
    updated_data = pd.DataFrame(valid_rows)
    
    updated_data.to_csv("align.csv", index=False)
    
    
        
augment_data()


# folder1 = out_root + sys.argv[1]
# folder2 = out_root + "reference"
# # Example usage
# cif_file1 = f'{folder1}/seed-2_sample-0/model.cif'
# cif_file2 = f'{folder2}/seed-2_sample-0/model.cif'
# json_file1 = f'{folder1}/seed-2_sample-0/confidences.json'  # Replace with your first JSON file
# json_file2 = f'{folder2}/seed-2_sample-0/confidences.json'

# # print(extract_confident_atoms(cif_file1, json_file1))
# # findPandQ(cif_file1, cif_file2, json_file1, json_file2)
# tm_align_rmsd(cif_file1, cif_file2, json_file1, json_file2)


