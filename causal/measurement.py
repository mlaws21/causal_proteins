
# what do i need:
# i need data: some percentage chance of disease and 
# look at different types of noise
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from numpy.linalg import inv
from collections import defaultdict
"""
DAG:

A --> Y
|
|
V
A*

Need: 
P(A*| A)
P(x,y,z) = sum_w{I(z, w)P(x,y,w)}

"""
def sigmoid(x):
    return 1 / (1 + np.exp(-x))


def calc_joint_dist(df, param="Astar"):

    return {
        (0, 0): df[(df[param] == 0) & (df['Y'] == 0)].shape[0],
        (0, 1): df[(df[param] == 0) & (df['Y'] == 1)].shape[0],
        (1, 0): df[(df[param] == 1) & (df['Y'] == 0)].shape[0],
        (1, 1): df[(df[param] == 1) & (df['Y'] == 1)].shape[0]
      
    }

# epsilon = P(A* = 0 | A = 1)
# delta = P(A* = 1 | A = 0)
def matrix_adjustment(data, relational_data):
    
    epsilon = relational_data[(relational_data["Astar"] == 0) & (relational_data['A'] == 1)].shape[0] / relational_data[(relational_data['A'] == 1)].shape[0]
    delta = relational_data[(relational_data["Astar"] == 1) & (relational_data['A'] == 0)].shape[0] / relational_data[(relational_data['A'] == 0)].shape[0]
    
    print(f"epsilon: {epsilon:.2f}; delta: {delta:.2f}")

    M = np.array([[1 - delta, epsilon], [delta, 1- epsilon]])
    I = inv(M)
    joint_dist = calc_joint_dist(data)
    
    out_joint = defaultdict(int)
    for (astar, y), prob in joint_dist.items():
        for a in [0, 1]:
            out_joint[(a, y)] += I[a, astar] * prob
    
    return out_joint

def backdoor_adjustment(data: pd.DataFrame, a_name: str, y_name: str, c_names: list[str]) -> float:
    """
    Perform backdoor adjustment for a given treatment A and outcome Y using
    the covariates in Z
    """

    # make a regression formula
    c_names = ["1"] + c_names
    c_formula = " + ".join(c_names)
    regression_formula = f"{y_name} ~ {c_formula} + {a_name}"

    # fit a regression depending on whether Y is binary or not
    if set(data[y_name]) == {0, 1} or set(data[y_name]) == {True, False}:
        model = smf.glm(formula=regression_formula, family=sm.families.Binomial(), data=data).fit()
    else:
        model = smf.glm(formula=regression_formula, family=sm.families.Gaussian(), data=data).fit()

    data_a1 = data.copy() # make a copy for the interventional datasets
    data_a1[a_name] = 1
    data_a0 = data.copy()
    data_a0[a_name] = 0

    return round(np.mean(model.predict(data_a1) - model.predict(data_a0)), 3)

def data_from_joint(joint):
    
    data = []
    for k,v in joint.items():

        a, y = k
        for _ in range(int(v)):
            data.append({"A": a, "Y": y})
            
    df = pd.DataFrame(data)
    
    return df
    # for c in c_names:

def generate_data_confounders(n_samples=1000, relation_only=False):
    # Step 1: Generate confounders C1 and C2
    # C1 and C2 are normally distributed independent variables
    C1 = np.random.normal(0, 1, n_samples)  # Normally distributed with mean=0 and std=1
    C2 = np.random.normal(0, 1, n_samples)  # Normally distributed with mean=0 and std=1

    # Step 2: Generate A based on C1 and C2
    # Using a logistic model for A influenced by C1 and C2
    prob_A = 1 / (1 + np.exp(- (0.5 * C1 + 0.5 * C2)))  # Sigmoid function for probability
    A = np.random.binomial(1, prob_A, n_samples)  # A is binary, 0 or 1

    # Step 3: Generate A* as a noisy version of A
    # The probability of A* being equal to A depends on C1 and C2
    prob_Astar_equals_A = 0.9 + 0.05 * C1 - 0.05 * C2  # Adding small influence of C1 and C2 on noise
    prob_Astar_equals_A = np.clip(prob_Astar_equals_A, 0, 1)  # Clip to [0, 1] range
    Astar = np.array([a if np.random.uniform(0, 1) < p else 1 - a for a, p in zip(A, prob_Astar_equals_A)])

    # Step 4: Generate Y based on A, C1, and C2
    # Logistic model for Y, influenced by A, C1, and C2

    if relation_only:
        df = pd.DataFrame({

        "A": A,
        "Astar": Astar,
        
        })
        
        return df
    else:
        prob_Y = 1 / (1 + np.exp(- (1.0 * A + 0.3 * C1 + 0.3 * C2)))
        Y = np.array([np.random.choice([0, 1], p=[0.2, 0.8]) if a == 1 else np.random.choice([0, 1], p=[0.7, 0.3]) for a in A])
    
        
        df = pd.DataFrame({

        "A": A,
        "Astar": Astar,
        "C1": C1,
        "C2": C2,
        "Y": Y
        })
    
        return df

def generate_data(n_samples=1000, relation_only=False):
    # Step 1: Generate A (discrete variable)
    # Let's assume A has values 0 or 1 with equal probability
    
    A = np.random.choice([0, 1], size=n_samples, p=[0.5, 0.5])
    
    Astar = np.array([a if np.random.uniform(0, 1) < 0.9 else 1 - a for a in A])
    
    
    
    # Z = np.random.uniform(-5, 5, n_samples)  # A is a standard normal variable
    
    # # binarize 
    # A = (2 * Z + np.random.normal(0, 1, n_samples) > 0.5).astype(int)
    
    # # Step 3: Generate Y based on B
    
    # # error rate of 10%: right now A just depends of
    # # Astar = (A if np.random.uniform(0,1) < 0.9 else 1 - A) + Z/10 + np.random.normal(0, 1, n_samples) > 0
    # Astar = (([a if np.random.uniform(0, 1) < 0.9 else 1 - a for a in A]) + np.random.normal(0, 0.5, n_samples) > 0.5).astype(int)
    

    
    # Step 2: Generate Y based on A
    # Define conditional probabilities for Y given A
    # P(Y=1 | A=1) = 0.8, P(Y=1 | A=0) = 0.3
    # causal effect = 0.5
    if relation_only:
        df = pd.DataFrame({

        "A": A,
        "Astar": Astar,
        })
        
        return df
    else:
        Y = np.array([np.random.choice([0, 1], p=[0.2, 0.8]) if a == 1 else np.random.choice([0, 1], p=[0.7, 0.3]) for a in A])
    
        
        df = pd.DataFrame({

        "A": A,
        "Astar": Astar,
        "Y": Y
        })
    
        return df

# def generate_data(n_samples=1000):
#     # Step 1: Generate A
#     Z = np.random.uniform(-5, 5, n_samples)  # A is a standard normal variable
    
#     # binarize 
#     A = (2 * Z + np.random.normal(0, 1, n_samples) > 0.5).astype(int)
    
#     # Step 3: Generate Y based on B
    
#     # error rate of 10%: right now A just depends of
#     # Astar = (A if np.random.uniform(0,1) < 0.9 else 1 - A) + Z/10 + np.random.normal(0, 1, n_samples) > 0
#     Astar = ((A if np.random.uniform(0,1) < 0.9 else 1 - A) + np.random.normal(0, 0.5, n_samples) > 0.5).astype(int)
    

    
#     # Step 4: Generate B* based on both A and B
#     Y = (3 * A + np.random.normal(0, 1, n_samples)) > 1.5

    
#     return Z, A, Astar, Y


df = generate_data()
relation_df = generate_data(n_samples=10000, relation_only=True)



# print(df.head(15))

# but this doesnt work
print(backdoor_adjustment(df, "Astar", "Y", []))

# so we want this
print(backdoor_adjustment(df, "A", "Y", []))

# print((Z))

# print(np.random.normal(0, 1, 1000) )

joint = matrix_adjustment(df, relation_df)

print(joint)

synth_df = data_from_joint(joint)
print(synth_df)


print(backdoor_adjustment(synth_df, "A", "Y", []))
