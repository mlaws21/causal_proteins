from Bio.PDB import MMCIFParser
import numpy as np
import json

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

folder1 = "fold_prion_tv3"
cif_file1 = f'{folder1}/{folder1}_model_0.cif'
json_file1 = f'{folder1}/{folder1}_full_data_0.json'  # Replace with your first JSON file


def extract_plddt_scores(json_file):
    """Extracts pLDDT scores for each residue from the AlphaFold JSON file."""
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    return data.get('atom_plddts', [])



# print(get_atom_coords_by_residue_from_cif(cif_file1)[("A", 1)])
# residues2cords = get_atom_coords_by_residue_from_cif(cif_file1)


# print(extract_plddt_scores(json_file1))

def extract_plddt_and_coords(cif_file, json_file):
    residues2cords = get_atom_coords_by_residue_from_cif(cif_file)
    plddt_list = extract_plddt_scores(json_file)
    
    residues2plddt = {}
    ind = 0
    for k, v in residues2cords.items():
        residues2plddt[k] = []
        for _ in v:
            
            residues2plddt[k].append(plddt_list[ind])
            
            ind += 1
            
    return residues2cords, residues2plddt

residues2cords, residues2plddt = extract_plddt_and_coords(cif_file1, json_file1)
print(len(residues2cords[2]))
print(len(residues2plddt[2]))


def calculate_distances(coords1, coords2):
    """Calculate the Euclidean distance between corresponding atoms."""
    distances = np.linalg.norm(coords1 - coords2, axis=1)
    return distances

def compare_structures(cif_file1, cif_file2, json_file1, json_file2):
    """Compares two structures by calculating distances and extracting pLDDT scores."""
    residues2cords1, residues2plddt1 = extract_plddt_and_coords(cif_file1, json_file1)
    residues2cords2, residues2plddt2 = extract_plddt_and_coords(cif_file2, json_file2)
    


    # print(plddt1.keys())
    conf_dists = []
    
    # print(len(coords1), len(coords2), len(plddt1), len(plddt2))
    for residue_id in coords1.keys():
        if residue_id in coords2:
            # print(residue_id)
            atom_coords1 = np.array(coords1[residue_id])
            atom_coords2 = np.array(coords2[residue_id])
            atom_plddt1 = np.array(plddt1[residue_id[1]])
            atom_plddt2 = np.array(plddt2[residue_id[1]])
            
            print(atom_plddt1)
            
            # print("atom coords", atom_coords1, atom_coords2)
            # Check if both residues have the same number of atoms
            if len(atom_coords1) == len(atom_coords2):
                distances = calculate_distances(atom_coords1, atom_coords2)
                # print(distances)
                
                
                # Get pLDDT scores for the residue
                # print("1")
                score1 = np.mean(atom_plddt1) 
                # print("2")
                score2 = np.mean(atom_plddt2)
                
                avg_score = (score1 + score2) / 2
                if avg_score > 80.0:
                    conf_dists.append(np.mean(distances))
                
                print(f"Residue {residue_id}: Distance = {np.mean(distances):.2f} Å | pLDDT1 = {score1:.2f}, pLDDT2 = {score2:.2f}, pLDDTm = {avg_score:.2f}")
            else:
                print(f"Warning: Residue {residue_id} has different numbers of atoms between the two structures.")
    
    conf_dists = np.array(conf_dists)
    
    print(conf_dists)
    
    # Calculate and print RMSD
    # rmsd = np.sqrt(np.mean(conf_dists**2))
    rmsd = np.sqrt(np.mean(conf_dists**2))
    
    print(f"\nRMSD for confident regions: {rmsd:.2f} Å")

folder1 = "fold_prion_tv3"
folder2 = "fold_prion_tv3m_p105l"

# Example usage
cif_file1 = f'{folder1}/{folder1}_model_0.cif'
cif_file2 = f'{folder2}/{folder2}_model_0.cif'
json_file1 = f'{folder1}/{folder1}_full_data_0.json'  # Replace with your first JSON file
json_file2 = f'{folder2}/{folder2}_full_data_0.json'
# compare_structures_with_plddt(cif_file1, cif_file2, json_file1, json_file2)
