import sys
import os
import joblib
import json
import numpy as np
import pandas as pd
from tqdm import tqdm

from rdkit import Chem

from __init__ import OUTPUT

print("PREDICTOR")

ROOT = os.path.dirname(os.path.abspath(__file__))

predictorapp_folder = os.path.join(ROOT, "..", "predictorapp")
sys.path.append(predictorapp_folder)

from src.eosdescriptors.chembl import Chembl
from src.eosdescriptors.ecfp import Ecfp
from src.eosdescriptors.rdkit2d import Rdkit2d
from src.eosdescriptors.rdkitfpbits import RdkitFpBits

MODELS_DIR = os.path.join(predictorapp_folder, "model")

def get_models_files():
    models_files = {}
    with open(os.path.join(MODELS_DIR, "models.json"), "r") as f:
        all_models = json.load(f)
    for dn, tn in all_models:
        models_files[(dn, tn)] = os.path.join(MODELS_DIR, "{0}_{1}.pkl".format(dn, tn))
    return all_models, models_files

all_models, models_files = get_models_files()

def get_sessions():
    sessions = {}
    for k, v in models_files.items():
        sess = joblib.load(v)
        sessions[k] = sess
    return sessions

sessions = get_sessions()

necessary_descriptors = sorted(set([k[0] for k in all_models]))
necessary_tasks = sorted(set([k[1] for k in all_models]))

def get_descriptor_calculators():
    descriptor_calculators = {
        "chembl": Chembl(),
        "ecfp": Ecfp(),
        "rdkit2d": Rdkit2d(),
        "rdkitfpbits": RdkitFpBits()
    }
    return descriptor_calculators

descriptor_calculators = get_descriptor_calculators()

headers = []
for dn in necessary_descriptors:
    for tn in necessary_tasks:
        headers += [dn+"_"+tn]

def chunker(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def prediction(mols):
    preds = {}
    for dn in necessary_descriptors:
        fp = descriptor_calculators[dn]
        X = fp.calc(mols)
        for tn in necessary_tasks:
            sess = sessions[(dn, tn)]
            if tn == "classification":
                pred = sess.predict_proba(X)[:,1]
            else:
                pred = sess.predict(X)[:]
            preds[(dn, tn)] = pred
    Y = np.zeros((len(mols), len(preds)))
    i = 0
    for dn in necessary_descriptors:
        for tn in necessary_tasks:
            v = preds[(dn, tn)]
            Y[:,i] = v
            i += 1
    return Y


# Fit on series4 IC50 mean
print("Fitting metapredictor")

df = pd.read_csv(os.path.join(ROOT, "..", "predictor", "data", "data.csv"))
smiles = np.array(df["Smiles"])
y = np.array(df["Activity"])
mask = y < 10
y = y[mask]
smiles = smiles[mask]

mols = [Chem.MolFromSmiles(smi) for smi in smiles]
X = prediction(mols).astype(np.float32)

from sklearn.linear_model import LinearRegression
mdl = LinearRegression()

from sklearn.ensemble import GradientBoostingRegressor
# Set lower and upper quantile
LOWER_ALPHA = 0.1
UPPER_ALPHA = 0.9
# Each model has to be separate
lower_model = GradientBoostingRegressor(loss="quantile",
                                        alpha=LOWER_ALPHA)
# The mid model will use the default loss
mid_model = GradientBoostingRegressor(loss="ls")
upper_model = GradientBoostingRegressor(loss="quantile",
                                        alpha=UPPER_ALPHA)


mdl = mid_model

from sklearn.model_selection import ShuffleSplit
from sklearn.metrics import r2_score

splits = ShuffleSplit(n_splits=5, test_size=0.2)
scores = []

import matplotlib.pyplot as plt
for train_idx, test_idx in splits.split(X):
    mdl.fit(X[train_idx], y[train_idx])
    yp = mdl.predict(X[test_idx])
    yt = y[test_idx]
    scores += [r2_score(yp, yt)]
score = np.mean(scores)
print("R2 score", score)

upper_model.fit(X, y)
mid_model.fit(X, y)
lower_model.fit(X, y)

# Get smiles

df=pd.read_csv(os.path.join(OUTPUT, "data_5.csv")) # TODO CHANGE

smiles = list(df["Smiles"])

Y = None
for chunk in tqdm(chunker(smiles, 10000)):
    mols = [Chem.MolFromSmiles(smi) for smi in chunk]
    Y_ = prediction(mols).astype(np.float32)
    if Y is None:
        Y = Y_
    else:
        Y = np.vstack([Y, Y_])

y_u = upper_model.predict(Y)
y = mid_model.predict(Y)
y_l = lower_model.predict(Y)

df["IC50Pred"] = y
df["IC50PredUB"] = y_u
df["IC50PredLB"] = y_l

cutoff = np.mean(list(df["IC50Pred"]))

df = df[df["IC50Pred"] < cutoff]

df.to_csv(os.path.join(OUTPUT, "data_6.csv"), index=False)
