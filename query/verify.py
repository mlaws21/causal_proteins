import pandas as pd
import numpy as np


data = pd.read_csv("pre_treatment.csv")

print(np.sum(data["is_pathogenic"] == True) / len(data))


