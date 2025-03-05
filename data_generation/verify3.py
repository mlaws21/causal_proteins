import pandas as pd
from collections import Counter
# Read the CSV file
df = pd.read_csv("prion_uniprot.csv")
print(len(df))

# Check if the column exists
if "ID" in df.columns:
    # Get unique values in the column
    unique_values = df["ID"].dropna().unique()
    
    c = Counter(df["ID"])
    
    print(c)
    # Print the unique values
    print("Unique values in 'Genomic_Coordinate_hg38':")
    print(len(unique_values))
else:
    print("Column 'Genomic_Coordinate_hg38' not found in the CSV file.")
