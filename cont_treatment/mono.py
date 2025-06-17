
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle
from scipy.optimize import linprog, curve_fit
import cvxpy as cp
from concurrent.futures import ThreadPoolExecutor, as_completed
from scipy.stats import norm
from scipy.interpolate import CubicSpline, PchipInterpolator
from pygam import LinearGAM, s
import sys
import os
from bisect import bisect_right


def fa(data, treatment, log_fn=print):
    
    formula = f'{treatment} ~ 1'
    log_fn(f"Distribution of A Formula: {formula}")
    
    model = smf.glm(formula=formula, data=data, family=sm.families.Gaussian())
    result = model.fit()
    return result 
    

def mu(data, treatment, outcome, log_fn=print):
    # Ep[Y | A = a, W = w]
    
    confounders = [col for col in data.columns if col not in [treatment, outcome]]
    confounder_form = "+".join(confounders)
    formula = f'{outcome} ~ 1 + {treatment} + {confounder_form}'
    log_fn(f"Mu Formula: {formula}")
    
    # print(formula)
    model = smf.glm(formula=formula, data=data, family=sm.families.Binomial())
    result = model.fit()
    
    return result

def pi(data, treatment, outcome, log_fn=print):
    #pi(a|w)
    
    confounders = [col for col in data.columns if col not in [treatment, outcome]]
    confounder_form = "+".join(confounders)
    formula = f'{treatment} ~ 1 + {confounder_form}'
    log_fn(f"Pi Formula: {formula}")

    # gaussian bc continuous
    model = smf.glm(formula=formula, data=data, family=sm.families.Gaussian())
    result = model.fit()
    return result 

def gamma(a, data, ij_data, mu_est, pi_est, fa_est, treatment, outcome):

    
    # get nuissance vars
    mu_pred = mu_est.predict(data)
    pi_pred = pi_est.predict(data)
    fa_pred = fa_est.predict(data)
    
    Y = data[outcome]
    
    # get marginal distribution fn(a)
    std = np.std(data[treatment] - fa_pred)
    prob_a = norm.pdf(data[treatment], loc=fa_pred, scale=std)
    
    ind = data[treatment] <= a
    front = np.mean(ind * (Y- mu_pred)/(pi_pred/prob_a))
    
    # back part (double sum)
    indij = ij_data[treatment] <= a
    mu_predij = mu_est.predict(ij_data)
    back = np.mean(indij * mu_predij)
    
    return front + back
    
def F(a, data, treatment):
    # emperical distribution function
    # this models the cumulative distribution function by convention
    
    series = data[treatment]
    # print(series)
    n = len(series)
    edf = np.sum(series <= a) / n
    return edf

def theta(a, gcm_x, gcm_y, data, treatment="T"):
    Fa = F(a, data, treatment=treatment)
    
    # print(a)
    # print(Fa)
    ind = np.where(gcm_x == Fa)[0][0]
    
    slope_left = (gcm_y[ind] - gcm_y[ind - 1]) / (gcm_x[ind] - gcm_x[ind - 1])
    return slope_left
    # print(Fa)
    
def get_points(data, treatment, outcome, log_fn=print):
    # gets points for GCM in sorted order
    
    # start by adding (0, 0)
    xs = [0]
    ys = [0]


    # build estimators
    fa_est = fa(data, treatment, log_fn)
    pi_est = pi(data, treatment, outcome, log_fn)
    mu_est = mu(data, treatment, outcome, log_fn)
    
    # build dataframe that pairs all treatments with all confounders -- this allows us to vectorize later
    # use same cols as original data
    ij_df = pd.DataFrame(columns=data.columns)
    ij_df = ij_df.astype(data.dtypes)
    
    log_fn("Generating Giant Dataframe for Vectorization")
    # go through all rows
    for i, rowi in data.iterrows():
        if i % (len(data) // 10) == 0:
            log_fn(f"{(i / len(data)) * 100:.0f}% complete")
        
        # go through all rows
        for _, rowj in data.iterrows():
            
            cp = rowj.copy()
            cp[treatment] = rowi[treatment]
            # we make lots of new rows where 
            ij_df = pd.concat([ij_df, pd.DataFrame([cp])], ignore_index=True)

    # NOTE: can enforce uniqueness if needed -- but shouldnt need to bc continuous
    log_fn("Calculating Points")
    seen = set()
    for i, row in data.iterrows():
        if i % (len(data) // 10) == 0:
            log_fn(f"{(i / len(data)) * 100:.0f}% complete")
        if row[treatment] not in seen:
            xs.append(F(row[treatment], data, treatment))
            ys.append(gamma(row[treatment], data, ij_df, mu_est, pi_est, fa_est, treatment, outcome))
            seen.add(row[treatment])
    
    
    sorted_idx = np.argsort(xs)
    x_sorted = np.array(xs)[sorted_idx]
    y_sorted = np.array(ys)[sorted_idx]
    log_fn(f"Calculated {len(x_sorted)} Unique Points")
    return x_sorted, y_sorted
    
def non_smooth_gcm(x_sorted, y_sorted):
    
    
    # Ensure data is sorted by x
    n = len(x_sorted)

    # Define variables: f_i = f(x_i)
    f = cp.Variable(n)

    # Constraints:
    constraints = []

    # 1) f_i <= y_i
    for i in range(n):
        constraints.append(f[i] <= y_sorted[i])

    # 2) Convexity constraints: second difference >= 0
    #    i.e., f_{i-1} - 2 f_i + f_{i+1} >= 0
    for i in range(1, n-1):
        constraints.append(f[i-1] - 2*f[i] + f[i+1] >= 0)

    # Objective: maximize sum(f)
    objective = cp.Maximize(cp.sum(f))

    # Solve problem
    prob = cp.Problem(objective, constraints)
    result = prob.solve(solver=cp.SCS, verbose=False)

    if prob.status not in ["optimal", "optimal_inaccurate"]:
        raise ValueError(f"Convex minorant optimization failed with status: {prob.status}")

    f_values = f.value  # Optimal values for f_i

    return x_sorted, f_values

def is_convex(x, f):

    x = np.asarray(x, dtype=float)
    f = np.asarray(f, dtype=float)

    if len(x) != len(f):
        raise ValueError("x and f must have the same length.")
    if len(x) < 3:
        # With fewer than 3 points, you cannot have a "violation" of convexity
        # in the discrete second-difference sense.
        return True

    # Check slopes: slope(i-1 -> i) <= slope(i -> i+1)
    for i in range(1, len(x) - 1):
        slope_left = 0 if (x[i]   - x[i-1]) == 0 else (f[i]   - f[i-1]) / (x[i]   - x[i-1])
        slope_right = 0 if (x[i+1]   - x[i]) == 0 else (f[i+1] - f[i])   / (x[i+1] - x[i])
        if round(slope_left, 5) > round(slope_left, 5):
            return False

    return True


def compute_ground(data, treatment_value, treatment="T", outcome="Y"):
    # df1 = pd.read_csv("synthetic_data.csv")
    # df2 = pd.read_csv("synthetic_data_addition.csv")
    
    # data = df1
    
    # data = pd.concat([df1, df2], ignore_index=True)
    
    data[treatment] = treatment_value
    
    # b0, b1, b2, b3, b4 = -1, 0.1, -0.5, 0.3, 0.8  # Coefficients
    
    # logit = (
    #     b0
    #     + b1 * data["X1"]
    #     + b2 * data["X2"]
    #     + b3 * data["X3"]
    #     + b4 * data[treatment]
    # )
    
    #BRCA COEFFS
    # intercept, coef_old, coef_white, coef_unhealthy, coef_align = -3, 0.5, 0.1, 0.3, 0.3
    
    #PRION COEFFS
    intercept, coef_old, coef_white, coef_unhealthy, coef_align = -3, 0.5, 0.1, 0.3, 0.8  # Coefficients
    
    logit = (
        intercept
        + coef_old * data["Old"]
        + coef_white * data["White"]
        + coef_unhealthy * data["Unhealthy"]
        + coef_align * data["Align_Score"]
    )
    
    prob_Y = 1 / (1 + np.exp(-logit))  # Sigmoid function
    
    data[outcome] = np.random.binomial(1, prob_Y)
    
    return data[outcome].mean()

def compute_ground_truth(data, treatment_value, treatment, outcome, coeffs):

    data[treatment] = treatment_value
    
    
    
    columns = data.columns.tolist()
    columns.remove(outcome)
    logit = coeffs["Intercept"] + sum([coeffs[x] * data[x] for x in columns])
    
    prob_Y = 1 / (1 + np.exp(-logit))  # Sigmoid function
    
    data[outcome] = np.random.binomial(1, prob_Y)
    
    return data[outcome].mean()

def build_cubic_spline(x_data, f_data):

    x_data = np.asarray(x_data, dtype=float)
    f_data = np.asarray(f_data, dtype=float)

    # Build the cubic spline object
    cs = PchipInterpolator(x_data, f_data)
    return cs

def backdoor(a, data, treatment="T", outcome="Y"):
    
    
    confounders = [col for col in data.columns if col not in [outcome]]
    
    formula = f"{outcome} ~ 1 + {' + '.join(confounders)}"
    
    print("backdoor formula", formula)
    model = smf.glm(formula=formula, family=sm.families.Binomial(), data=data).fit()
    data_a = data.copy()
    data_a[treatment] = a
    data_0 = data.copy()
    data_0[treatment] = 0
    return np.mean(model.predict(data_a)) - np.mean(model.predict(data_0))
    
def verify(xs, ys, gcm_x, gcm_y):
    # verifies the GCM -- true if verifies
    
    # this enforces the minorant constrant
    
    for y, gy in zip(ys, gcm_y):

        if np.round(y, 3) < np.round(gy, 3):
            print("Minorant Constraint Violated")
            print(np.round(y, 3), np.round(gy, 3))
            return False
    
    return is_convex(gcm_x, gcm_y)

def verify_spline_convex(spline, data, treatment):
    
    myrange = np.linspace(data[treatment].min(), data[treatment].max(), 100)
    snd_derivs = np.array([spline(F(i, data, treatment=treatment), 2) for i in myrange])
    
    # assert(np.all(snd_derivs >= 0))
    # print(snd_derivs)
    

    # plt.savefig(f"dose_response{str(count)}.png", format="png", dpi=300)
    
    plt.clf()


def exp_func(x, a, b):
    return a * np.exp(b * x)

def exp_deriv(x, a, b):
    return a * b * np.exp(b * x)
    
def fit_exponential(xlist, ylist):
    """
    Fits an exponential curve y = a * exp(b * x) to the data.
    
    Parameters:
        xlist (list or np.ndarray): The x values.
        ylist (list or np.ndarray): The y values.
    
    Returns:
        a, b (float): Parameters for the best fit y = a * exp(b * x)
    """

    # Define the exponential model

    
    # Convert input to numpy arrays
    x = np.array(xlist)
    y = np.array(ylist)

    # Fit the curve
    params, _ = curve_fit(exp_func, x, y, p0=(1.0, 0.1))  # initial guess a=1, b=0.1
    return params[0], params[1]

#TODO Make sure the spline is convex -- but pretty close

def make_strict_convex_underestimator(points, delta=0):
    """
    Given points = [(x1,y1),…,(xn,yn)], returns a convex, piecewise-linear
    function f(x) with f(xi) < yi for all i.

    If delta is None, we pick a tiny shift based on the data range.
    """

    # 1) sort by x
    pts = sorted(points, key=lambda p: (p[0], p[1]))
    xs_all, ys_all = zip(*pts)

    # 2) build the *lower* convex hull of the points
    hull = []
    for x,y in pts:
        # while the last turn is not "convex downwards", pop
        while len(hull) >= 2:
            x1,y1 = hull[-2]
            x2,y2 = hull[-1]
            # cross‐product test for right‐turn or collinear
            if (x2 - x1)*(y - y1) - (y2 - y1)*(x - x1) <= 0:
                hull.pop()
            else:
                break
        hull.append((x,y))

    # 3) pick a small delta if not provided
    if delta is None:
        y_min, y_max = min(ys_all), max(ys_all)
        # use 1e-5 of the overall y‐range, or a tiny constant if flat
        delta = (y_max - y_min)*1e-5 if y_max > y_min else 1e-6

    # unpack hull into two lists for fast lookup
    xs_h, ys_h = zip(*hull)

    # 4) return the piecewise‐linear convex function f(x) = hull(x) - delta
    def f(x):
        # extrapolate with the first/last segment if outside data range
        if x <= xs_h[0]:
            i0, i1 = 0, 1
        elif x >= xs_h[-1]:
            i0, i1 = -2, -1
        else:
            # find the segment [xs_h[i], xs_h[i+1]] that contains x
            i0 = bisect_right(xs_h, x) - 1
            i1 = i0 + 1

        x0, y0 = xs_h[i0], ys_h[i0]
        x1, y1 = xs_h[i1], ys_h[i1]
        t = (x - x0) / (x1 - x0)
        # linear interp minus delta
        return (1-t)*y0 + t*y1 - delta
    
    def df(x):
        # same segment lookup
        if x <= xs_h[0]:

            i0, i1 = 0, 1
        elif x >= xs_h[-1]:

            
            i0, i1 = -2, -1
        else:
            i0 = bisect_right(xs_h, x) - 1
            i1 = i0 + 1

        x0, y0 = xs_h[i0], ys_h[i0]
        x1, y1 = xs_h[i1], ys_h[i1]
        # slope of the line
        return (y1 - y0) / (x1 - x0)

    return f, df

def generate_dr_curve(patient_numerical_data, project_name, coeffs, log_fn, treatment="Align_Score", outcome="Disease"):
    
    dr_pkl = f"pickles/{project_name}_spline.pkl"
    
    # TODO FIXME
    if os.path.exists(dr_pkl):
    # if False:
        log_fn(f"Reading in Dose Response Curve from {dr_pkl}")
        
        with open(dr_pkl, "rb") as f:
            xlist, ylist = pickle.load(f)
            
    else:
        # patient_numerical_data = patient_numerical_data.head(100)
        log_fn(f"Generating Dose Response Curve for {len(patient_numerical_data)} patients")
        
        xlist, ylist = get_points(patient_numerical_data, treatment, outcome, log_fn)   

        gcm_x, gcm_y = non_smooth_gcm(xlist, ylist)
        
        assert verify(xlist, ylist, gcm_x, gcm_y)
        
        
        with open(dr_pkl, "wb") as f:
            pickle.dump((xlist, ylist), f)
            
    log_fn("Building Spline")
    
    # a, b = fit_exponential(xlist, ylist)
    # print(f"Fitted: y = {a:.4f} * exp({b:.4f} * x)")
    plt.clf()
    pts = list(zip(xlist, ylist))
    mygcm, mygcm_deriv = make_strict_convex_underestimator(pts)
    cubic_spline = build_cubic_spline(xlist, ylist)
    # exp_ys = [exp_deriv(F(i, patient_numerical_data, treatment=treatment), a, b) for i in np.arange(0, 10.5, 0.5)]
    # print(exp_ys)
    
    if not os.path.exists(f"outputs/{project_name}/results"):
        os.makedirs(f"outputs/{project_name}/results")
    
    dr_graph_file = f"outputs/{project_name}/results/dr_curve.png"
    dr_eps_file =f"outputs/{project_name}/results/dr_curve.eps"
    
    x_vals = np.arange(0, 10.5, 0.5)
    
    # TODO: cubic spline not fully convex... its like basically there but a little noisy
    plt.plot(x_vals, [compute_ground_truth(patient_numerical_data.copy(), i, treatment, outcome, coeffs) for i in x_vals], color='blue', label="Ground", marker='o')
    plt.plot(x_vals, [cubic_spline(F(i, patient_numerical_data, treatment=treatment), 1) for i in x_vals] , color='red', label="Estimate", marker='x')
    # plt.plot(x_vals, [mygcm_deriv(F(i, patient_numerical_data, treatment=treatment)) for i in x_vals], color='black', label="gcm", marker='+')
    # plt.plot(x_vals, [F(i, patient_numerical_data, treatment=treatment) for i in x_vals], color='green', label="Exp", marker='+')
    # plt.plot(x_vals, exp_ys, color='green', label="Exp", marker='+')
    
    # plt.plot(np.arange(0, 15.5, 0.5), [backdoor(i, data, treatment=treatment, outcome=outcome) for i in np.arange(0, 15.5, 0.5)], color='green', label="Backdoor", marker='+')
    
    plt.xlabel("Misalignment Score (Å)")
    plt.ylabel("Probability of Disease")
    plt.title("Dose Response Curve")
        
    plt.grid(True)
    plt.legend()
    plt.savefig(dr_graph_file, format="png", dpi=300)
    plt.savefig(dr_eps_file, format="eps", dpi=300)
    
    
    log_fn(f"Saving Dose Response Graph to {dr_graph_file}")
    
     
    plt.clf()

    gcm_graph_file = f"outputs/{project_name}/results/gcm.png"
    
    plt.plot(np.arange(0, 1.0001, 0.001), [cubic_spline(x) for x in np.arange(0, 1.00001, 0.001)], color='green', label='Spline')
    plt.scatter(xlist, ylist, color='blue', label='Points')
    # plt.plot(np.arange(0, 1.0001, 0.001), [exp_func(x, a, b) for x in np.arange(0, 1.0001, 0.001)], color='red', label="Exp", marker='+')
    

    # Add labels and title
    plt.xlabel("F(Align Score)") # where F is the continous emperical distribution
    plt.ylabel("Gamma") # that crazy formula
    plt.title("Global Convex Minorant (GCM)")

    # Add grid and legend
    plt.grid(True)
    plt.legend()
    plt.savefig(gcm_graph_file, format="png", dpi=300)
    # plt.savefig(gcm_graph_file, format="eps", dpi=300)
    
    log_fn(f"Saving GCM Graph to {gcm_graph_file}")

    return cubic_spline

def main():
    pass
    
#     if len(sys.argv) < 2:
#         print("Usage: python mono.py [count]")
#         exit(-1)
#     count = int(sys.argv[1])
    
#     print("COUNT:", count)
#     # data = pd.read_csv("synthetic_data.csv").head(100)
#     data = pd.read_csv("prion_data.csv").head(count)
#     # data = pd.read_csv("cleaned_data.csv").head(200)
    
    
#     treatment = "Align_Score"
#     outcome = "Cancer"
    
#     xlist, ylist = get_points(data, treatment=treatment, outcome=outcome)   

#     gcm_x, gcm_y = non_smooth_gcm(xlist, ylist)
    
#     # cubic_spline = build_cubic_spline(gcm_x, gcm_y)
#     cubic_spline = build_cubic_spline(xlist, ylist)
    
#     with open("prion_spline.pkl", "wb") as f:
#         pickle.dump((xlist, ylist), f)
    
#     verify_spline_convex(cubic_spline, data, treatment)
    
#     # NOTE: this is an idea
#     # gam = LinearGAM(s(0, constraints="convex")).fit(gcm_x, gcm_y)
    
    
#     assert verify(xlist, ylist, gcm_x, gcm_y)
    
#     print(verify(xlist, ylist, xlist, ylist))
    
    
    
#     # NOTE: not using spline doesn't work because weird edge behaivor
#     # plt.plot(np.arange(0, 15.5, 0.5), [theta(i, gcm_x, gcm_y, data, treatment="Align_Score")for i in np.arange(0, 15.5, 0.5)] , color='red', label="Estimate", marker='x')
#     # plt.plot(np.arange(0, 15.5, 0.5), [theta(i, xlist, ylist, data, treatment="Align_Score")for i in np.arange(0, 15.5, 0.5)] , color='red', label="Estimate", marker='x')
    
    
#     # TODO: cubic spline not fully convex... its like basically there but a little noisy
#     plt.plot(np.arange(0, 15.5, 0.5), [cubic_spline(F(i, data, treatment=treatment), 1) for i in np.arange(0, 15.5, 0.5)] , color='red', label="Estimate", marker='x')
#     plt.plot(np.arange(0, 15.5, 0.5), [compute_ground(data, i, treatment=treatment, outcome=outcome) for i in np.arange(0, 15.5, 0.5)], color='blue', label="Ground", marker='o')
#     # plt.plot(np.arange(0, 15.5, 0.5), [backdoor(i, data, treatment=treatment, outcome=outcome) for i in np.arange(0, 15.5, 0.5)], color='green', label="Backdoor", marker='+')
    
#     plt.xlabel("Alignment Score (Å)")
#     plt.ylabel("Probability of Cancer")
#     plt.title("Dose Response Curve")
        
#     plt.grid(True)
#     plt.legend()
#     plt.savefig(f"prion_dose_response{str(count)}.png", format="png", dpi=300)
    
     
#     plt.clf()
    
#     # plt.plot(X_grid, deriv_1, color='green', label='Spline')
    
#     # plt.plot(gcm_x, gcm_y, color='red', label='GCM', marker='x')
#     plt.plot(np.arange(0, 1.0001, 0.001), [cubic_spline(x) for x in np.arange(0, 1.00001, 0.001)], color='green', label='Spline')
#     plt.scatter(xlist, ylist, color='blue', label='Points')

#     # Add labels and title
#     plt.xlabel("F(Align Score)") # where F is the continous emperical distribution
#     plt.ylabel("Gamma") # that crazy formula
#     plt.title("Global Convex Minorant (GCM)")

#     # Add grid and legend
#     plt.grid(True)
#     plt.legend()
#     plt.savefig(f"prion_gcm{str(count)}.png", format="png", dpi=300)


if __name__ == "__main__":
    main()