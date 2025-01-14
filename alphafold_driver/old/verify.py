import pandas as pd

# Step 1: Read data from file
df = pd.read_csv("data.csv")

# Get unique elements in Column1
unique_elements = df["ID"].unique()
print("Unique elements in Column1:", unique_elements)
