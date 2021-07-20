from __init__ import OUTPUT

from tqdm import tqdm
import pandas as pd
from rdkit import Chem
import os, sys
import numpy as np
import joblib

df=pd.read_csv(os.path.join(OUTPUT, "data_3.csv"), nrows=3000)

smiles = list(df["Smiles"])

predictor_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "predictor")

sys.path.append(predictor_folder)

from descriptors.eosdescriptors.rdkit2d import Rdkit2d
from descriptors.utils import impute, normalize

desc = Rdkit2d()

print("Loading transformation data")
model_folder = os.path.join(predictor_folder, "models", "rdkit2d")

with open(os.path.join(model_folder, "medians.npy"), "rb") as f:
    medians = np.load(f)

with open(os.path.join(model_folder, "mus.npy"), "rb") as f:
    mus = np.load(f)

with open(os.path.join(model_folder, "sigmas.npy"), "rb") as f:
    sigmas = np.load(f)

def chunker(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

print("Loading classifier model")
sel0 = joblib.load(os.path.join(model_folder, "sel0.pkl"))
sel1 = joblib.load(os.path.join(model_folder, "sel1.pkl"))
sel2 = joblib.load(os.path.join(model_folder, "sel2.pkl"))
mdl = joblib.load(os.path.join(model_folder, "classifier.pkl"))

yp = []
done = 0
for chunk in chunker(smiles, 1000):
    print("Calculating descriptors")
    mols = [Chem.MolFromSmiles(smi) for smi in tqdm(chunk)]
    X = desc.calc(mols)
    X = impute(X, medians)
    X = normalize(X, mus, sigmas)
    print("Predicting")
    X = sel0.transform(X)
    X = sel1.transform(X)
    X = sel2.transform(X)
    yp += list(mdl.predict(X))
    done += len(chunk)
    print("Done", done)

df["Rdkit2dClassifier"]=yp

df=df[df["Rdkit2dClassifier"] == 1]

df.to_csv(os.path.join(OUTPUT, "data_4.csv"), index=False)


#df.to_csv(os.path.join(OUTPUT, "data_4.csv"), index=False)
