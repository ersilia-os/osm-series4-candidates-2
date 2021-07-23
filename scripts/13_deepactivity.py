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

print("Filter based on chemprop")


cp_ic50_p50 = np.percentile(df["IC50Pred"], 50)
print(cp_ic50_p50)

ALPHA = 25

cp_clf_p25 = np.percentile(df["ActivityClfGraph"], ALPHA)
print(cp_clf_p25)

cp_reg_p25 = np.percentile(df["ActivityRegGraph"], ALPHA)
print(cp_reg_p25)

df = df[
    (df["IC50Pred"] < cp_ic50_p50) |
    (
        (df["ActivityClfGraph"] > cp_clf_p25) & (df["ActivityRegGraph"] > cp_reg_p25)
    )
        ]
print(df.shape)

ALPHA = 5

gr_clf_p25 = np.percentile(df["ActivityClfVocab"], ALPHA)
print(gr_clf_p25)

gr_reg_p25 = np.percentile(df["ActivityRegVocab"], ALPHA)
print(gr_reg_p25)

gr_ic_p75 = np.percentile(df["ActivityIcVocab"], 100-ALPHA)
print(gr_ic_p75)

df = df[
    (df["IC50Pred"] < cp_ic50_p50) |
    (
        (df["ActivityClfVocab"] > gr_clf_p25) & (df["ActivityRegVocab"] > gr_reg_p25) & (df["ActivityIcVocab"] < gr_ic_p75)
    )
        ]

print(df.shape)

df.to_csv(os.path.join(OUTPUT, "data_13.csv"), index=False)
