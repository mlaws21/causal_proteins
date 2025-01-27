import pandas as pd
from sklearn.linear_model import LinearRegression, LogisticRegression
import numpy as np

def cont_treatment(data, treatment="T", outcome="Y"):
    a = data[treatment]
    y = data[outcome]
    l = data[[col for col in data.columns if col not in [treatment, outcome]]]
    
    a_min, a_max = min(a), max(a)
    a_vals = np.linspace(a_min, a_max, 1000)
    l_repeated = np.repeat(l, len(a_vals), axis=0)
    print(l_repeated.shape)
    a_repeated = np.tile(a_vals, 1000)
    la_new = np.column_stack((l_repeated, a_repeated))
    l_new = la_new[:, :-1]  # Remove the last column (a)
    
    print(l_new.shape)
    # Fit GLM (Linear Regression) for `a`
    glm_a = LinearRegression()
    glm_a.fit(l, a)  # Train on `l` to predict `a`
    pimod_vals = glm_a.predict(l_new)  # Predict `a` for new data
    sq_res = (a - glm_a.predict(l)) ** 2  # Squared residuals

    # Fit GLM (Linear Regression) for squared residuals
    glm_res = LinearRegression()
    glm_res.fit(l, sq_res)  # Train on `l` to predict squared residuals
    pi2mod_vals = glm_res.predict(l_new)  # Predict squared residuals for new data

    # Fit GLM (Logistic Regression) for binary outcome `y`
    la_train = np.column_stack((l, a))  # Combine `l` and `a` for training
    glm_y = LogisticRegression(solver='lbfgs', max_iter=1000, random_state=42)
    glm_y.fit(la_train, y)  # Train logistic regression to predict `y`
    muhat_vals = glm_y.predict_proba(la_new)[:, 1]  # Predict probabilities for new data

        
def main():
    data = pd.read_csv("synthetic_data.csv")
    
    cont_treatment(data)
    
    


if __name__ == "__main__":
    main()