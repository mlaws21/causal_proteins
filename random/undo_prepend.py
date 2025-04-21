import os

# Set the path to your target directory
directory = '/shared/25mdl4/af_output'


for filename in os.listdir(directory):
    if filename.startswith("prion411_"):
        old_path = os.path.join(directory, filename)
        new_filename = filename[len("prion411_"):]
        new_path = os.path.join(directory, new_filename)
        os.rename(old_path, new_path)
        print(f"Reverted {filename} to {new_filename}")
