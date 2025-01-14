import pandas as pd

# Read the CSV file
df = pd.read_csv("brca_mutations.csv")
print(len(df))

# Check if the column exists
if "Genomic_Coordinate_hg38" in df.columns:
    # Get unique values in the column
    unique_values = df["Genomic_Coordinate_hg38"].dropna().unique()
    
    # Print the unique values
    print("Unique values in 'Genomic_Coordinate_hg38':")
    print(len(unique_values))
else:
    print("Column 'Genomic_Coordinate_hg38' not found in the CSV file.")
