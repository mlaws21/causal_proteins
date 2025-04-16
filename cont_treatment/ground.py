import numpy as np
import pandas as pd

def compute_mean_outcome(treatment_value, treatment="T", outcome="Y"):
    df1 = pd.read_csv("synthetic_data.csv")
    # df2 = pd.read_csv("synthetic_data_addition.csv")
    data = df1
    
    # data = pd.concat([df1, df2], ignore_index=True)
    
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


# print(compute_mean_outcome(0.4))

for i in range(-20, 20):

    print(f"Causal Effect {i:.2f}: {compute_mean_outcome(i):.3f}")