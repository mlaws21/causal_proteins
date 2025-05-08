import os

# Path to the directory containing the files
directory = "/shared/25mdl4/af_output"  # Change this to your specific path if needed

for filename in os.listdir(directory):
    if filename.startswith("prion411_"):
        new_name = filename.replace("prion411_", "prion_", 1)
        old_path = os.path.join(directory, filename)
        new_path = os.path.join(directory, new_name)
        os.rename(old_path, new_path)
        print(f"Renamed: {filename} → {new_name}")
