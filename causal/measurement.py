
# what do i need:
# i need data: some percentage chance of disease and 
# look at different types of noise
import numpy as np

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def generate_data(n_samples=1000):
    # Step 1: Generate A
    A = np.random.normal(0, 1, n_samples)  # A is a standard normal variable
    
    # Step 2: Generate B based on A
    B = 2 * A + np.random.normal(0, 1, n_samples)  # B depends on A with some noise
    
    # Step 3: Generate Y based on B
    Y = sigmoid(3 * B + np.random.normal(0, 1, n_samples))  # Y depends on B with some noise
    
    # Step 4: Generate B* based on both A and B
    B_star = 1.5 * A + 0.5 * B + np.random.normal(0, 1, n_samples)  # B* depends on both A and B
    
    return A, B, B_star, Y

# Example usage:
A, B, B_star, Y = generate_data()
