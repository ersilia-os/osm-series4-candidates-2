from __init__ import OUTPUT

print("run _chemprop.sh in chemprop conda environment")
print("run _grover.sh in grover conda environment")

import pandas as pd
import numpy as np
import os

ROOT = os.path.dirname(os.path.abspath(__file__))

print("DEEP LEARNING PREDICTIONS")

df = pd.read_csv(os.path.join(OUTPUT, "data_12.csv"))
print(df.shape)

cp = pd.read_csv(os.path.join(ROOT, "_pred_chemprop.csv"))
gr = pd.read_csv(os.path.join(ROOT, "_pred_grover.csv"))

df = pd.concat([df, cp, gr], axis=1)

print(df.shape)

df.to_csv(os.path.join(OUTPUT, "data_13.csv"), index=False)
