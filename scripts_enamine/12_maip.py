from __init__ import OUTPUT

import pandas as pd
import numpy as np
import os

ROOT = os.path.dirname(os.path.abspath(__file__))

print("MAIP")

df = pd.read_csv(os.path.join(OUTPUT, "data_11.csv"))
print(df.shape)

df.to_csv(os.path.join(OUTPUT, "data_12.csv"), index = False)
