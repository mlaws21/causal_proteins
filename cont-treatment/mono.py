
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


def mu(data, treatment="T", outcome="Y"):
    # Ep[Y | A = a, W = w]
    
    confounders = [col for col in data.columns if col not in [treatment, outcome]]
    confounder_form = "+".join(confounders)
    formula = f'{outcome} ~ {treatment} + {confounder_form}'
    
    # print(formula)
    model = smf.glm(formula=formula, data=data, family=sm.families.Binomial())
    result = model.fit()
    
    return result

def g(data, treatment="T", outcome="Y"):
    #pi(a|w)/fp(a)
    # according to rohit i think we can ignore fp(a) -- at laest for now
    
    confounders = [col for col in data.columns if col not in [treatment, outcome]]
    confounder_form = "+".join(confounders)
    formula = f'{treatment} ~ {confounder_form}'
    
    # print(formula)
    # gaussian bc continuous
    model = smf.glm(formula=formula, data=data, family=sm.families.Gaussian())
    result = model.fit()
    return result 
    
#def fa(data, treatment="T")

def gamma(a, data, mu_est, g_est, treatment="T", outcome="Y"):
    
    n = len(data)
    
    confounders = [col for col in data.columns if col not in [treatment, outcome]]
    confounders_treatment = [treatment] + confounders 
    
    # g_est = g(data, treatment, outcome)
    # mu_est = mu(data, treatment, outcome)

    
    # front = sum([(row[outcome] - mu_est.predict(row[confounders_treatment]))/g_est.predict(row[confounders]) for _, row in data.iterrows() if row[treatment] <= a])/n
    
    # for _, row in data.iterrows():
    #     if row[treatment] <= a:
            
    #         # print(row[confounders_treatment][0])
    #         # print(row[confounders_treatment])
    #         mu_pred = mu_est.predict(row[confounders_treatment])
    #         #TODO figure out what is fp(a)
    #         g_pred = g_est.predict(row[confounders]) / F(a, data, treatment=treatment)
            
    #         # print(mu_pred.values[0], g_pred.values[0])
            
    #         front += (row[outcome] - mu_pred.values[0]) / g_pred.values[0]
            
    #         # print(mu_pred + g_pred)
            
    
    # print(front)
    front = 0
    back = 0
    for _, rowi in data.iterrows():
        
        if rowi[treatment] <= a:
            
            
            mu_pred = mu_est.predict(rowi[confounders_treatment])
            #TODO figure out what is fp(a)
            g_pred = g_est.predict(rowi[confounders]) #/ F(a, data, treatment=treatment)
    
            
            front += (rowi[outcome] - mu_pred.values[0]) / g_pred.values[0]
            
            for _, rowj in data.iterrows():
                
                new_data = rowj[confounders_treatment]
                new_data[treatment] = rowi[treatment]

                # print(mu_est.predict(new_data))
                predicted_values = mu_est.predict(new_data).values[0]
                back += predicted_values
                
    # print(back)
    # back /= (n*n)
    # print(back)
    
    total = (front / n) + (back / (n*n))
    # print(front + back)

    
    # print(type(front))
    
    return total
    
def F(a, data, treatment="T"):
    # emperical distribution function
    # this models the cumulative distribution function by convention
    
    series = data[treatment]
    # print(series)
    n = len(series)
    # print(a)
    edf = np.sum(series <= a) / n
    return edf
  

# def process_row(row, data, treatment, outcome):
  
#     val = row[treatment]
#     x = F(val, data, treatment)
#     y = gamma(val, data, treatment, outcome)
#     return x, y
  
def get_points(data, treatment="T", outcome="Y"):
    xs = [0]
    ys = [0]
    
    
    # with ThreadPoolExecutor(max_workers=64) as executor:
    #     # Submit tasks for each row
    #     futures = [
    #         executor.submit(process_row, row, data, treatment, outcome)
    #         for _, row in data.iterrows()
    #     ]

    #     # As each future completes, unpack the results
    #     for future in as_completed(futures):
    #         x, y = future.result()
    #         xs.append(x)
    #         ys.append(y)

    # return xs, ys
    

    g_est = g(data, treatment, outcome)
    mu_est = mu(data, treatment, outcome)

    
    for _, row in data.iterrows():
        xs.append(F(row[treatment], data, treatment))
        ys.append(gamma(row[treatment], data, mu_est, g_est, treatment, outcome))
        
    
    # xs = np.array(xs)
    
    
    
    return xs, ys
    
    
def is_below(point1, point2, query_point):
    """
    Determine if the query point is below or on the line defined by two points.

    Parameters:
    - point1: Tuple (x1, y1) representing the first point on the line.
    - point2: Tuple (x2, y2) representing the second point on the line.
    - query_point: Tuple (xq, yq) representing the query point.

    Returns:
    - True if the query point is below or on the line.
    - False if the query point is above the line.
    """
    # Extract coordinates
    x1, y1 = point1
    x2, y2 = point2
    xq, yq = query_point
    
    # Calculate the y-value on the line at xq
    y_on_line = y1 + (y2 - y1) * (xq - x1) / (x2 - x1)

    # Return True if query point is below or on the line, otherwise False
    return yq < y_on_line

    
def non_smooth_gcm(x, y):
    
    
        # Ensure data is sorted by x
    sorted_idx = np.argsort(x)
    x_sorted = np.array(x)[sorted_idx]
    y_sorted = np.array(y)[sorted_idx]
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
    # """
    # Computes the greatest convex minorant of the points (x, y).
    # """
    
    # n = len(x)
    
    # sorted_indices = np.argsort(x)
    # x = np.array(x)[sorted_indices]
    # y = np.array(y)[sorted_indices]


    # gcm_x = x[:2]
    # gcm_y = y[:2]
    
    # for i in range(2, n):
        
    #     print(i)
    #     print(gcm_x)

    #     if is_below((gcm_x[-2], gcm_y[-2]), (x[i], y[i]), (gcm_x[-1], gcm_y[-1])):
    #         print("below")
    #         gcm_x = np.append(gcm_x, gcm_x[-1])
    #         gcm_y = np.append(gcm_y, gcm_y[-1])
    #     else: 
    #         gcm_x = np.append(gcm_x, x[i])
    #         gcm_y = np.append(gcm_y, y[i])
            
    # return gcm_x, gcm_y
            

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
        
def theta(a, gcm_x, gcm_y, data, treatment="T"):
    Fa = F(a, data, treatment=treatment)
    
    # print(a)
    # print(Fa)
    ind = np.where(gcm_x == Fa)[0][0]
    
    slope_left = (gcm_y[ind] - gcm_y[ind - 1]) / (gcm_x[ind] - gcm_x[ind - 1])
    return slope_left
    # print(Fa)


def compute_mean_outcome(treatment_value, treatment="T", outcome="Y"):
    df1 = pd.read_csv("synthetic_data.csv")
    df2 = pd.read_csv("synthetic_data_addition.csv")
    # data = df1
    
    data = pd.concat([df1, df2], ignore_index=True)
    
    data[treatment] = treatment_value
    
    b0, b1, b2, b3, b4 = -1, 0.1, -0.5, 0.3, 0.8  # Coefficients
    
    logit = (
        b0
        + b1 * data["X1"]
        + b2 * data["X2"]
        + b3 * data["X3"]
        + b4 * data[treatment]
    )
    
    prob_Y = 1 / (1 + np.exp(-logit))  # Sigmoid function
    
    data[outcome] = np.random.binomial(1, prob_Y)
    
    return data[outcome].mean()

from scipy.interpolate import CubicSpline
def build_cubic_spline(x_data, f_data):
    """
    Build a cubic spline interpolant for the given data points (x_data, f_data).

    Parameters
    ----------
    x_data : array_like
        1D array of x-coordinates, assumed to be strictly increasing.
    f_data : array_like
        1D array of function values at each x_data point.

    Returns
    -------
    spline_function : callable
        A function that, given a scalar x, returns the interpolated f(x).
    spline_derivative : callable
        A function that, given a scalar x, returns the derivative f'(x).
    """
    # Convert inputs to numpy arrays (good practice, not always required)
    x_data = np.asarray(x_data, dtype=float)
    f_data = np.asarray(f_data, dtype=float)

    # Build the cubic spline object
    cs = CubicSpline(x_data, f_data)

    # Define a function to evaluate the spline at a given x
    def spline_function(x):
        """
        Evaluate the cubic spline interpolation at x.
        """
        return cs(x)

    # Define a function to evaluate the derivative of the spline at a given x
    def spline_derivative(x):
        """
        Evaluate the first derivative of the cubic spline interpolation at x.
        """
        # Option 1: Use cs.derivative()(x)
        return cs.derivative()(x)
        
        # Option 2: Alternatively, CubicSpline supports calling cs(x, nu=1)
        # which returns the first derivative:
        # return cs(x, nu=1)

    return spline_function, spline_derivative
  
def main():
    data = pd.read_csv("synthetic_data1.csv").head(20)

    xlist, ylist = get_points(data)
    

    # with open("xs.pkl", "wb") as file:
    #     pickle.dump(xlist, file)
        
    # with open("ys.pkl", "wb") as file:
    #     pickle.dump(ylist, file)
        
    # with open("xs.pkl", "rb") as file:
    #     xlist = pickle.load(file)
        
    # with open("ys.pkl", "rb") as file:
    #     ylist = pickle.load(file)
        
    sorted_idx = np.argsort(xlist)
    x_sorted = np.array(xlist)[sorted_idx]
    y_sorted = np.array(ylist)[sorted_idx]
    
    gcm_x, gcm_y = non_smooth_gcm(xlist, ylist)
    
    
    for i in range(len(x_sorted)):

        # print(x_sorted[i], ":", y_sorted[i], gcm_y[i])
        assert round(y_sorted[i], 3) >= round(gcm_y[i], 3)
        
    assert is_convex(gcm_x, gcm_y)
    
    # this number is the averate 
    
    # for i in np.arange(0, 1, 0.1):
    fx, f1_x = build_cubic_spline(gcm_x, gcm_y)
    
    
    
    for i in range(-20, 20):
    
        causal_effect = theta(i, gcm_x, gcm_y, data)
        
        
        Fi = F(i, data, treatment="T")
        print(f"Causal Effect (spline) {i:.2f}: {f1_x(Fi):.3f}")
    
        # print(f"Causal Effect {i:.2f}: {causal_effect:.3f}")
        print(f"Causal Effect {i:.2f}: {compute_mean_outcome(i):.3f}")
        

    
    
    plt.plot(range(-20, 20), [f1_x(F(i, data, treatment="T")) for i in range(-20, 20)] , color='red', label="Estimate", marker='x')
    plt.plot(range(-20, 20), [compute_mean_outcome(i) for i in range(-20, 20)], color='red', label="Ground", marker='o')
        
        
        
    
    # plt.plot(gcm_x, gcm_y, color='red', label='GCM', marker='x')
    # plt.plot(np.arange(0, 1.0001, 0.001), [fx(x) for x in np.arange(0, 1.00001, 0.001)], color='green', label='Spline')
    
    # plt.scatter(xlist, ylist, color='blue', label='Points')

    # Add labels and title
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.title("Global Convex Minorant (GCM)")

    # Add grid and legend
    plt.grid(True)
    plt.legend()
    plt.savefig("scatter_plot1.png", format="png", dpi=300)

    # Show the plot

    # gamma(0.4, data)
    


if __name__ == "__main__":
    main()