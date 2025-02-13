import pandas as pd

# Load the data from a CSV file
file_path = "data.csv"  # Replace with your file path
df = pd.read_csv(file_path)

# Convert 'Cancer' to binary (1 for True, 0 for False) if not already binary
df["Cancer"] = df["Cancer"].astype(int)

# Group by the relevant columns and calculate percentages
grouped = df.groupby(["Old", "White", "Unhealthy"])["Cancer"].agg(["sum", "count"]).reset_index()
grouped["Cancer_Percentage"] = (grouped["sum"] / grouped["count"]) * 100

# Calculate total cancer cases and overall percentage
total_cancer_cases = df["Cancer"].sum()
total_cases = len(df)
total_percentage = (total_cancer_cases / total_cases) * 100

# Add the total percentage as a new row
grouped.loc[len(grouped.index)] = {
    "Old": "Total",
    "White": "Total",
    "Unhealthy": "Total",
    "sum": total_cancer_cases,
    "count": total_cases,
    "Cancer_Percentage": total_percentage
}

# Rename columns for clarity
grouped = grouped.rename(columns={"sum": "Cancer_Cases", "count": "Total_Cases"})

# Display the result
print(grouped)