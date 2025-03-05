import pandas as pd

# Load the CSV file into a DataFrame
file_path = "prion_uniprot.csv"  # Replace with your file path
df = pd.read_csv(file_path)

# Calculate the proportion of "pathogenic" in Clinical_significance_ENIGMA
pathogenic_count = sum(df["is_pathogenic"])
total_count = len(df)
proportion_pathogenic = pathogenic_count / total_count

# Display the result
print(f"Proportion of 'pathogenic': {proportion_pathogenic:.2%}")
