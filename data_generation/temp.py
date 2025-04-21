import pandas as pd

# Read the CSV
df = pd.read_csv("prion_uniprot.csv")

# Keep only specified columns (drop all others)
df = df[["ID", "Sequence", "is_pathogenic"]]

# Drop duplicate rows
df = df.drop_duplicates()

df.to_csv('prion_input.csv', index=False)
