from __init__ import OUTPUT

import pandas as pd
import numpy as np
import os

ROOT = os.path.dirname(os.path.abspath(__file__))

print("MAIP")

df = pd.read_csv(os.path.join(OUTPUT, "data_11.csv"))
print(df.shape)
print(df)

maip = pd.read_csv(os.path.join(ROOT, "..", "maip", "maip_predictions.csv"))
print(maip.shape)
print(maip)

df["Maip"] = maip["model_score"]

cut = np.percentile(df["Maip"], 10)

df = df[df["Maip"] > cut]
print(df.shape)

df.to_csv(os.path.join(OUTPUT, "data_12.csv"), index = False)
