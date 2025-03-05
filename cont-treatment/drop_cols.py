
import pandas as pd

# Load the DataFrame
df = pd.read_csv("prion_data.csv")

# Remove specific columns
df = df.drop(columns=["Sequence", "is_pathogenic", "ID"])  # Replace with actual column names
# df.loc[:, df.dtypes == bool] = df.loc[:, df.dtypes == bool].astype(int)
# Save the modified DataFrame back to a CSV file
df.to_csv("prion_data.csv", index=False)

print("Columns removed and DataFrame saved as 'final_updated.csv'.")

