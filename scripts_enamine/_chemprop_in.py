from __init__ import OUTPUT

import pandas as pd
import os


df = pd.read_csv(os.path.join(OUTPUT, "data_12.csv"))
smiles = list(df["Smiles"])

with open("_chemprop.csv", "w") as f:
    f.write("smiles\n")
    for smi in smiles:
        f.write("{0}\n".format(smi))
