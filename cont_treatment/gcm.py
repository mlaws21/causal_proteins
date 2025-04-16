import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle
from scipy.optimize import linprog
import cvxpy as cp
from concurrent.futures import ThreadPoolExecutor, as_completed
from scipy.stats import norm
from scipy.interpolate import CubicSpline, PchipInterpolator
import sys

from mono import compute_ground

def F(a, data, treatment):
    # emperical distribution function
    # this models the cumulative distribution function by convention
    
    series = data[treatment]
    n = len(series)
    edf = np.sum(series <= a) / n
    return edf

def verify_spline_convex(spline):
    
    data = pd.read_csv("prion_data.csv")
    treatment = "Align_Score"
    
    myrange = np.linspace(data[treatment].min(), data[treatment].max(), 100)
    snd_derivs = np.array([spline(F(i, data, treatment=treatment), 2) for i in myrange])
    plt.plot(myrange, snd_derivs, color='green', label='Spline')
    plt.plot(myrange, [0 for i in myrange], color='red', label='Spline')
    
    
    plt.savefig("snd_derivs.png", format="png", dpi=300)
    
    
    plt.clf()
    
def calc_gcm(points_filename):
    with open(points_filename, "rb") as f:
        x, y = pickle.load(f)
    
    x_gcm, y_gcm = compute_gcm(x, y)

    # Create a dense set of x values for plotting the piecewise linear GCM.
    x_dense = np.linspace(x.min(), x.max(), 500)
    y_dense = piecewise_linear_interp(x_gcm, y_gcm, x_dense)

    # Plot the original data and the computed GCM.
    plt.figure(figsize=(8, 5))
    plt.plot(x, y, 'o', label='Data points')
    plt.plot(x_dense, y_dense, 'r-', label='GCM (piecewise linear)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Greatest Convex Minorant of Data')
    plt.legend()
    # plt.show()
    plt.savefig("snd_derivs.png", format="png", dpi=300)
    
    
    # verify_spline_convex(spline)
    
def get_left_deriv(point, gcm_x, gcm_y, data, treatment):
    # gcms are in sorted order
    #TODO
    if point == 0:
        return 0

    pointf = F(point, data, treatment)
    
    startpt = -1
    endpt = -1
    for i, ele in enumerate(gcm_x):
        if ele >= pointf:
            endpt = i
            startpt = i-1
            break
    
    # make sure we find something
    assert endpt != -1 and startpt != -1
    
    print(point, pointf)
    print(startpt, endpt)
    rise = gcm_y[endpt] - gcm_y[startpt]
    run = gcm_x[endpt] - gcm_x[startpt]
    slope = rise / run
    
    return slope
    
            
            
# def verify_points_convex(gcm_x, gcm_y, data, treatment):
#     myrange = np.linspace(data[treatment].min(), data[treatment].max(), 100)
    
    
    
    


def compute_gcm(x, y, data, treatment):
    """
    Compute the Greatest Convex Minorant (GCM) for points (x, y).
    Assumes that x is sorted in ascending order.
    Returns the x and y coordinates of the GCM knots.
    """
    # Ensure numpy arrays and sort by x
    x = np.asarray(x)
    y = np.asarray(y)
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    
    # Start with all indices as potential knots of the GCM.
    gcm_idx = []
    for i in range(len(x)):
        gcm_idx.append(i)
        # Check the convexity condition on the last three points.
        while len(gcm_idx) >= 3:
            k = gcm_idx[-1]
            j = gcm_idx[-2]
            i0 = gcm_idx[-3]
            # Calculate slopes between consecutive segments.
            slope1 = (y[j] - y[i0]) / (x[j] - x[i0])
            slope2 = (y[k] - y[j]) / (x[k] - x[j])
            # For convexity, slopes should be nondecreasing.
            if slope1 > slope2:
                # If the condition is violated, remove the middle point.
                gcm_idx.pop(-2)
            else:
                break
    return x[gcm_idx], y[gcm_idx]

def piecewise_linear_interp(x_knots, y_knots, x_vals):
    """
    Evaluate the piecewise linear interpolation defined by (x_knots, y_knots)
    at the points x_vals.
    """
    return np.interp(x_vals, x_knots, y_knots)


def main():
    
    with open("prion_spline.pkl", "rb") as f:
        x, y = pickle.load(f)
        
    data = pd.read_csv("prion_data.csv")
    treatment = "Align_Score"
    outcome = "Cancer"
    
    
    x_gcm, y_gcm = compute_gcm(x, y, data, treatment)
    
    # get_left_deriv(1, x_gcm, y_gcm)
    
    cubic_spline = PchipInterpolator(x, y)
    myrange = np.linspace(data[treatment].min(), data[treatment].max(), 50)
    
    # TODO: cubic spline not fully convex... its like basically there but a little noisy
    plt.plot(myrange, [cubic_spline(F(i, data, treatment=treatment), 1) for i in myrange] , color='Green', label="Spline", marker='+')
    plt.plot(myrange, [get_left_deriv(i, x_gcm, y_gcm, data, treatment) for i in myrange] , color='red', label="Estimate", marker='x')
    plt.plot(myrange, [compute_ground(data, i, treatment=treatment, outcome=outcome) for i in myrange], color='blue', label="Ground", marker='o')
    # plt.plot(np.arange(0, 15.5, 0.5), [backdoor(i, data, treatment=treatment, outcome=outcome) for i in np.arange(0, 15.5, 0.5)], color='green', label="Backdoor", marker='+')
    
    plt.xlabel("Alignment Score (Å)")
    plt.ylabel("Probability of Cancer")
    plt.title("Dose Response Curve")
        
    plt.grid(True)
    plt.legend()
    plt.savefig(f"prion_dose_response_custom.png", format="png", dpi=300)
    
if __name__ == "__main__":
    main()