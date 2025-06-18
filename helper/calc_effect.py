import pandas as pd
import numpy as np
import pickle
import sys
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt
import glob




def F(a, series):
    # emperical distribution function
    # this models the cumulative distribution function by convention
    
    n = len(series)
    edf = np.sum(series <= a) / n
    return edf

def calc_effect(intervention_data_filename, full_data_filename, spline_filename, treatment, main_treatment, num_boostraps=250):

    
    full_data = pd.read_csv(full_data_filename)
    
    intervention_data = pd.read_csv(intervention_data_filename)

    with open(f"{spline_filename}", "rb") as f:
        x_loaded, y_loaded = pickle.load(f)

    # Reconstruct the interpolator
    spline_loaded = PchipInterpolator(x_loaded, y_loaded)



    point_probs_cancer = np.array([spline_loaded(F(i, full_data[main_treatment]), 1) for i in intervention_data[treatment]])
    
    mean_point_probs_cancer = np.mean(point_probs_cancer)
    
    
    boot_probs = []
    for i in range(num_boostraps):
        bootstrap_data = intervention_data.sample(n=len(intervention_data), replace=True, random_state=i)
        
        boot_probs_cancer = np.array([spline_loaded(F(x, full_data[main_treatment]), 1) for x in bootstrap_data[treatment]])
        
        mean_boot_probs_cancer = np.mean(boot_probs_cancer)
        
        boot_probs.append(mean_boot_probs_cancer)
        
    
    
    return mean_point_probs_cancer, np.array(boot_probs)#np.quantile(boot_probs, 0.025), np.quantile(boot_probs, 0.975) 
    
    
    # TODO: cubic spline not fully convex... its like basically there but a little noisy
    plt.plot(np.arange(0, 20.5, 0.5), [spline_loaded(F(i, full_data["Align_Score"]), 1) for i in np.arange(0, 20.5, 0.5)] , color='red', label="Estimate", marker='x')
    # plt.plot(np.arange(0, 15.5, 0.5), [backdoor(i, data, treatment=treatment, outcome=outcome) for i in np.arange(0, 15.5, 0.5)], color='green', label="Backdoor", marker='+')
    
    plt.xlabel("Alignment Score (Å)")
    plt.ylabel("Probability of Cancer")
    plt.title("Dose Response Curve")
        
    plt.grid(True)
    plt.legend()
    plt.savefig("dose_response.png", format="png", dpi=300)

def calc_point_and_conf(intervention_data_filename, full_data_filename, spline_filename, main_treatment):
    with_point, with_boot = calc_effect(intervention_data_filename, full_data_filename, spline_filename=spline_filename,treatment="Align_Score_do_mut", main_treatment=main_treatment) 
    wo_point, wo_boot= calc_effect(intervention_data_filename, full_data_filename, spline_filename=spline_filename, treatment="Align_Score_do_no_mut", main_treatment=main_treatment)
    
    boots_stuff = with_boot - wo_boot
    
    return with_point - wo_point, np.quantile(boots_stuff, 0.025), np.quantile(boots_stuff, 0.975)
    

# Define the function to process each file
def process_files(project_name, treatment, log_fn):
    
    regex = f"outputs/{project_name}/intervention_data/*.csv"
    data = []
    # Find all CSV files matching the pattern
    csv_files = glob.glob(regex) # "post_treatment*.csv"

    for filename in csv_files:
        # Assuming `mut` and `spline_filename` are defined somewhere
        
        spl = filename.split("/")
        
        significance = spl[-1][-5]
        lab = 'Benign' if significance == 'b' else 'Pathogenic' if significance == 'p' else 'Unknown'
        # mutation = "_".join(spl[2:-1])
        mutation = spl[-1][:-6]
        point, bottom, top = calc_point_and_conf(filename, f'outputs/{project_name}/data.csv', f'outputs/{project_name}/pickles/spline.pkl', treatment)
        point_conf = f"Causal Effect: {point:.4f} ({bottom:.4f}, {top:.4f})"
        
        data.append((mutation, point, bottom, top, lab))
        log_fn(f"{mutation} [{lab}]: {point_conf}")

    return data
# Call the function
# process_files()
   
     

def main():
    pass
    

    
if __name__ == "__main__":
    main()