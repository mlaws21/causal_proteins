import os

# Set the path to your target directory
directory = '/shared/25mdl4/af_output'

for filename in os.listdir(directory):
    if not filename.startswith("chr17"):
        old_path = os.path.join(directory, filename)
        new_filename = "prion411_" + filename
        new_path = os.path.join(directory, new_filename)
        os.rename(old_path, new_path)
        print(f"Renamed {filename} to {new_filename}")
