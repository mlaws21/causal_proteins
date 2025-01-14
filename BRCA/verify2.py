import pandas as pd

# Load the CSV file into a DataFrame
file_path = "pop2.csv"  # Replace with your file path
df = pd.read_csv(file_path)

# Calculate the proportion of "pathogenic" in Clinical_significance_ENIGMA
pathogenic_count = df["Clinical_significance_ENIGMA"].str.lower().str.contains("pathogenic").sum()
total_count = len(df)
proportion_pathogenic = pathogenic_count / total_count

# Display the result
print(f"Proportion of 'pathogenic': {proportion_pathogenic:.2%}")
