import pandas as pd

# Load the CSV file (update 'your_file.csv' with the actual filename)
df = pd.read_csv("final.csv")

# Remove rows with duplicate values in the "Sequence" column (keeping only the first occurrence)
df = df.drop_duplicates(subset="Sequence", keep="first")

df.to_csv("cleaned_data.csv", index=False)
