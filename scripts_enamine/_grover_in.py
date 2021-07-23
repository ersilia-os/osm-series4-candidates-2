from __init__ import OUTPUT

import pandas as pd
import os


df = pd.read_csv(os.path.join(OUTPUT, "data_12.csv"))
smiles = list(df["Smiles"])

with open("_grover.csv", "w") as f:
    f.write("Smiles,Activity\n")
    for smi in smiles:
        f.write("{0},0\n".format(smi))
