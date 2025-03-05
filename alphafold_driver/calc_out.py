import pandas as pd
import numpy as np
# ID,Old,White,Unhealthy,Align Score,Cancer,Sequence,is_pathogenic
def compute_out():
    data = pd.read_csv("align.csv")
    
    intercept, coef_old, coef_white, coef_unhealthy, coef_align = -3, 0.5, 0.1, 0.3, 0.8  # Coefficients
    
    logit = (
        intercept
        + coef_old * data["Old"]
        + coef_white * data["White"]
        + coef_unhealthy * data["Unhealthy"]
        + coef_align * data["Align_Score"]
    )
    
    prob_Y = 1 / (1 + np.exp(-logit))  # Sigmoid function
    
    data["Cancer"] = np.random.binomial(1, prob_Y)
    
    print(np.mean(data["Cancer"]))
    

    data.to_csv("prion_data.csv", index=False)
    
compute_out()