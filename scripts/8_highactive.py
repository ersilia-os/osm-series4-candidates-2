from __init__ import OUTPUT

from tqdm import tqdm
import pandas as pd
from rdkit import Chem
import os, sys
import numpy as np
import joblib

print("CLASSIFIER HIGH ACTIVE")

df=pd.read_csv(os.path.join(OUTPUT, "data_7.csv"))

smiles = list(df["Smiles"])

ml_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "highpredictor", "results")

from mordred import Calculator, descriptors

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "predictor"))

from descriptors.utils import impute, normalize


IGNORE_3D = False

class Mordred(object):

    def __init__(self):
        pass

    def calc(self, mols):
        calc = Calculator(descriptors, ignore_3D=IGNORE_3D)
        df = calc.pandas(mols)
        return np.array(df, dtype=np.float)

desc = Mordred()

print("Loading transformation data")

with open(os.path.join(ml_folder, "medians.npy"), "rb") as f:
    medians = np.load(f)

with open(os.path.join(ml_folder, "mus.npy"), "rb") as f:
    mus = np.load(f)

with open(os.path.join(ml_folder, "sigmas.npy"), "rb") as f:
    sigmas = np.load(f)

with open(os.path.join(ml_folder, "sel1.npy"), "rb") as f:
    sel_idxs = np.load(f)

def chunker(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

print("Loading classifier model")
sel0 = joblib.load(os.path.join(ml_folder, "sel0.pkl"))
mdl = joblib.load(os.path.join(ml_folder, "classifier.pkl"))

yp = []
done = 0
print("Total", len(smiles))
for chunk in chunker(smiles, 1000):
    print("Calculating descriptors")
    mols = [Chem.MolFromSmiles(smi) for smi in tqdm(chunk)]
    X = desc.calc(mols)
    print(X.shape)
    X = impute(X, medians)
    X = normalize(X, mus, sigmas)
    print("Predicting")
    X = sel0.transform(X)
    print(X.shape)
    X = X[:,sel_idxs]
    print(X.shape)
    yp += list(mdl.predict_proba(X)[:,1])
    done += len(chunk)
    print("Done", done)

yp = np.array(yp)

print(df.shape)
df["HighClassifier"] = yp

print("Not very trustable predictor, only removing lowest quartile")
p25 = np.percentile(yp, 25)
print(p25)
df = df[df["HighClassifier"] > p25]

df.to_csv(os.path.join(OUTPUT, "data_8.csv"), index=False)
