from Bio.PDB import MMCIFParser
import numpy as np

def get_alpha_carbons_from_cif(cif_file):
    """Extracts the alpha-carbon (CA) coordinates from a .cif file."""
    parser = MMCIFParser()
    structure = parser.get_structure('protein', cif_file)
    
    alpha_carbons = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    alpha_carbons.append(residue['CA'].get_coord())
    
    return np.array(alpha_carbons)

def calculate_distances(coords1, coords2):
    """Calculate the Euclidean distance between corresponding alpha-carbon atoms."""
    distances = np.linalg.norm(coords1 - coords2, axis=1)
    return distances

def compare_cif_files(cif_file1, cif_file2):
    """Compares two CIF files by calculating the distances between corresponding alpha-carbons."""
    coords1 = get_alpha_carbons_from_cif(cif_file1)
    coords2 = get_alpha_carbons_from_cif(cif_file2)
    
    if len(coords1) != len(coords2):
        print("Warning: The structures have a different number of alpha-carbons.")
        return
    
    distances = calculate_distances(coords1, coords2)
    
    for i, dist in enumerate(distances):
        print(f"Residue {i+1}: Distance = {dist:.2f} Å")
    
    # Calculate and print RMSD
    rmsd = np.sqrt(np.mean(distances**2))
    print(f"\nRMSD: {rmsd:.2f} Å")

# Example usage
cif_file1 = 'fold_prion_tv3/fold_prion_tv3_model_0.cif'  # Replace with your actual CIF file
cif_file2 = 'fold_prion_tv3m_p105l/fold_prion_tv3m_p105l_model_0.cif'  # Replace with your second CIF file
# cif_file2 = 'fold_prion_tv3m_silent/fold_prion_tv3m_silent_model_0.cif'  # Replace with your second CIF file

compare_cif_files(cif_file1, cif_file2)
