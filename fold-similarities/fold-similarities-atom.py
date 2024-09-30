from Bio.PDB import MMCIFParser
import numpy as np
import json

# def get_alpha_carbons_from_cif(cif_file):
#     """Extracts the alpha-carbon (CA) coordinates from a .cif file."""
#     parser = MMCIFParser()
#     structure = parser.get_structure('protein', cif_file)
    
#     alpha_carbons = []
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 if 'CA' in residue:
#                     alpha_carbons.append(residue['CA'].get_coord())
    
#     return np.array(alpha_carbons)

def get_all_atom_coords_from_cif(cif_file):
    """Extracts the coordinates of all atoms from a .cif file."""
    parser = MMCIFParser()
    structure = parser.get_structure('protein', cif_file)
    
    all_atom_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    
                    all_atom_coords.append(atom.get_coord())
    
    return np.array(all_atom_coords)

def extract_plddt_scores(json_file):
    """Extracts pLDDT scores for each residue from the AlphaFold JSON file."""
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    return data.get('atom_plddts', [])

def calculate_distances(coords1, coords2):
    """Calculate the Euclidean distance between corresponding alpha-carbon atoms."""
    distances = np.linalg.norm(coords1 - coords2, axis=1)
    return distances

def compare_structures_with_plddt(cif_file1, cif_file2, json_file1, json_file2):
    """Compares two structures by calculating distances and extracting pLDDT scores."""
    coords1 = get_all_atom_coords_from_cif(cif_file1)
    coords2 = get_all_atom_coords_from_cif(cif_file2)
    
    if len(coords1) != len(coords2):
        print("Warning: The structures have a different number of alpha-carbons.")
        print(len(coords1), len(coords2))
        # return

    
    # Extract pLDDT scores from the JSON files
    plddt1 = extract_plddt_scores(json_file1)
    plddt2 = extract_plddt_scores(json_file2)
    
    
    distances = calculate_distances(coords1, coords2)
    
    if len(distances) != len(plddt1):
        print("Warning: In theory i should write something descriptive.")
        return 
    
    # Print distances and pLDDT scores for each residue
    # print(distances, plddt1, plddt2)
    
    conf_dists = []
    for i, (dist, score1, score2) in enumerate(zip(distances, plddt1, plddt2)):
        print(score1, score2)
        avg_score = (score1 + score2) / 2
        if avg_score > 80.0:
            conf_dists.append(dist)
        print(f"Residue {i+1}: Distance = {dist:.2f} Å | pLDDT1 = {score1:.2f}, pLDDT2 = {score2:.2f}, pLDDTm = {avg_score:.2f}")
    
    conf_dists = np.array(conf_dists)
    # Calculate and print RMSD
    rmsd = np.sqrt(np.mean(distances**2))
    
    conf_rmsd = np.sqrt(np.mean(conf_dists**2))
    
    # print(len(plddt1), len(distances))
    
    print(f"\nRMSD: {rmsd:.2f} Å")
    print(f"\nConfident RMSD: {conf_rmsd:.2f} Å")
    
folder1 = "fold_prion_tv3"
# folder2 = "fold_prion_tv3m_silent"
folder2 = "fold_prion_tv3m_p105l"


# Example usage
cif_file1 = f'{folder1}/{folder1}_model_0.cif'
cif_file2 = f'{folder2}/{folder2}_model_0.cif'
json_file1 = f'{folder1}/{folder1}_full_data_0.json'  # Replace with your first JSON file
json_file2 = f'{folder2}/{folder2}_full_data_0.json'
compare_structures_with_plddt(cif_file1, cif_file2, json_file1, json_file2)
