import pandas as pd

# Load the CSV file
df = pd.read_csv("prion_data.csv")


# print(df["Align_Score"])
print(df["Align_Score"].max())
print(df["Align_Score"].min())
print(df["Align_Score"].mean())
print(df["Align_Score"].std())

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the DataFrame
# df = pd.read_csv("final_full.csv")

# Choose the column to plot
column_name = "Align_Score"  # Replace with the actual column name

# Plot the distribution
plt.figure(figsize=(8, 5))
print(len(df["Sequence"].unique()))
print(len(df["ID"].unique()))

sns.histplot(df[column_name].unique(), bins=30, kde=True)  # kde=True adds a density curve
plt.xlabel(column_name)
plt.ylabel("Frequency")
plt.title(f"Distribution of {column_name}")
plt.savefig("test.png", format="png", dpi=300)



# Count the number of True values in 'Cancer' and 'is_pathogenic' columns
# num_cancer_true = df["Cancer"].sum()  # Assumes 'True' is stored as a boolean or 1
# num_pathogenic_true = df["is_pathogenic"].sum()  # Assumes 'True' is stored as a boolean or 1

# # Print results
# print(f"Number of True in Cancer: {num_cancer_true}")
# print(f"Number of True in is_pathogenic: {num_pathogenic_true}")


