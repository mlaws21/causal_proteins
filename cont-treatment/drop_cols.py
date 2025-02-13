
import pandas as pd

# Load the DataFrame
df = pd.read_csv("final_clean.csv")

# Remove specific columns
df = df.drop(columns=["Sequence", "is_pathogenic", "ID"])  # Replace with actual column names
# df.loc[:, df.dtypes == bool] = df.loc[:, df.dtypes == bool].astype(int)
# Save the modified DataFrame back to a CSV file
df.to_csv("final_clean.csv", index=False)

print("Columns removed and DataFrame saved as 'final_updated.csv'.")

