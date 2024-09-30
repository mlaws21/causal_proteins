import json

def print_json_field_titles(data, indent=0):
    if isinstance(data, dict):
        for key, value in data.items():
            print('  ' * indent + f'{key}')
            if isinstance(value, dict):  # Only recurse into dictionaries
                print_json_field_titles(value, indent + 1)
            elif isinstance(value, list) and value and isinstance(value[0], dict):
                # If it's a list of dictionaries, recurse into the first item
                print('  ' * (indent + 1) + "[...]")
                print_json_field_titles(value[0], indent + 1)

# Function to read JSON from a file and parse it
def parse_json_from_file(filename):
    with open(filename, 'r') as file:
        data = json.load(file)
    print_json_field_titles(data)
# Example usage
file_path = 'fold_prion_tv3/fold_prion_tv3_full_data_0.json'  # Replace with your file path
parse_json_from_file(file_path)
