from __init__ import OUTPUT

import numpy as np
import pandas as pd
import os

df = pd.read_csv(os.path.join(OUTPUT, "data_10.csv"))

print("COMBOS")

ALPHA = 10

print("Accessibility")

ra_cut = np.percentile(df["RAScore"], ALPHA)
print("RAScore", ra_cut)

sa_cut = np.percentile(df["SAScore"], 100-ALPHA)
print("SAScore", sa_cut)

sy_cut = np.percentile(df["SybaScore"], ALPHA)
print("SybaScore", sy_cut)

df = df[(df["RAScore"] > ra_cut) & (df["SAScore"] < sa_cut) & (df["SybaScore"] > sy_cut)]

print("Predictions")

p1_cut = np.percentile(df["IC50Pred"], 100-ALPHA)
print("IC50Pred", p1_cut)

p2_cut = np.percentile(df["HighClassifier"], ALPHA)
print("HighClassifier", p2_cut)

p3_cut = np.percentile(df["WhalesDist3Act"], 100-ALPHA)
print("WhalesDist3Act", p3_cut)

df = df[(df["IC50Pred"] < p1_cut) & (df["HighClassifier"] > p2_cut) & (df["WhalesDist3Act"] < p3_cut)]
print(df.shape)

print("Families")

columns = ['TriazoloHeteroaryl',
 'TriazoloPhenyl',
 'TriazoloHeteroarylPara',
 'TriazoloHeteroarylMeta',
 'TriazoloHeteroarylOrto',
 'TriazoloNaphthalene',
 'PyrazineEther',
 'PyrazineAmide']

df_ = df[columns]

fam = []
for r in df_.values:
    fam += [np.sum(r)]

df = df[np.array(fam) > 0]
print(df.shape)

print("Some physicochemistry")

df = df[(df["SLogP"] < 4.5) | (df["NumRings"] <= 5)]

print(df.shape)

df.to_csv(os.path.join(OUTPUT, "data_11.csv"), index = False)
