import numpy as np
import pandas as pd

# Set random seed for reproducibility
np.random.seed(42)

# Number of samples
n = 9000

# Generate three binary covariates (X1, X2, X3)
X1 = np.random.binomial(1, 0.5, n)  # Binary with P(X1=1) = 0.5
X2 = np.random.binomial(1, 0.3, n)  # Binary with P(X2=1) = 0.3
X3 = np.random.binomial(1, 0.7, n)  # Binary with P(X3=1) = 0.7

# Generate a continuous treatment (T) as a function of the covariates (optional)
T = 5 * X1 + 3 * X2 - 2 * X3 + np.random.normal(0, 1, n)

# Generate the binary outcome (Y) using a logistic model
# Logistic regression formula: P(Y=1) = sigmoid(b0 + b1*X1 + b2*X2 + b3*X3 + b4*T)
b0, b1, b2, b3, b4 = -1, 0.1, -0.5, 0.3, 0.8  # Coefficients
logit = b0 + b1 * X1 + b2 * X2 + b3 * X3 + b4 * T
prob_Y = 1 / (1 + np.exp(-logit))  # Sigmoid functionxw
Y = np.random.binomial(1, prob_Y)  # Generate binary outcome based on probability

# Create a DataFrame for better visualization and analysis
data = pd.DataFrame({
    'X1': X1,
    'X2': X2,
    'X3': X3,
    'T': T,
    'Y': Y
})

# Display the first few rows
print(data.head())

print(sum(data["Y"]))
# Save to CSV if needed
data.to_csv("synthetic_data_addition.csv", index=False)
