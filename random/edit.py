import os

folder_path = '/shared/25mdl4/af_output'


for filename in os.listdir(folder_path):
    if "m1m" in filename:
        new_name = filename.replace("m1m", "ref")
        old_path = os.path.join(folder_path, filename)
        new_path = os.path.join(folder_path, new_name)
        os.rename(old_path, new_path)
        print(f"Renamed: {filename} -> {new_name}")
