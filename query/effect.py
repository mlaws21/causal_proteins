import glob

from predict import calc_point_and_conf


# Define the function to process each file
def process_files():
    # Find all CSV files matching the pattern
    csv_files = glob.glob("post_treatment*.csv")

    with open("prion_output.txt", "a") as file:
        for filename in csv_files:
            # Assuming `mut` and `spline_filename` are defined somewhere
            
            spl = filename.split("_")
            
            print(f'{spl[2]} [{spl[3]}]: {calc_point_and_conf(filename, "prion_uniprot.csv", "prion_spline.pkl")}', file=file)

# Call the function
process_files()
