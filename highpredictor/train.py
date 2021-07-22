import os, sys
import joblib
from tqdm import tqdm
import shutil

ROOT = os.path.dirname(os.path.abspath(__file__))

ml_folder = os.path.join(ROOT, "results")

from mordred import Calculator, descriptors

IGNORE_3D = False

class Mordred(object):

    def __init__(self):
        pass

    def calc(self, mols):
        calc = Calculator(descriptors, ignore_3D=IGNORE_3D)
        df = calc.pandas(mols)
        return np.array(df, dtype=np.float)

desc = Mordred()

import pandas as pd
import numpy as np

from rdkit import Chem

df = pd.read_csv(os.path.join(ROOT, "..", "predictor", "data", "data.csv"))

smiles = np.array(df["Smiles"])

y = np.array(df["Activity"])

mask = y < 10
y = y[mask]
smiles = smiles[mask]

print(len(y))

mols = [Chem.MolFromSmiles(smi) for smi in smiles]

y_ = np.zeros(len(y))
y_[y < 1] = 1
y = np.array(y_)

print("Calculate descriptor")

X_file = os.path.join(ml_folder, "X.npy")
if os.path.exists(X_file):
    with open(X_file, "rb") as f:
        X = np.load(f)
else:
    X = desc.calc(mols)
    with open(X_file, "wb") as f:
        np.save(f, X)

print("Train model")

MAX_FEAT_INIT = 100

import os
from sklearn import svm
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import VarianceThreshold

sys.path.append(os.path.join(ROOT,".."))
from predictor.descriptors.utils import impute, normalize

print("Get medians")
medians = []
for i in range(X.shape[1]):
    vals = X[:,i]
    vals = vals[~np.isnan(vals)]
    if len(vals) == 0:
        medians += [0.]
    else:
        medians += [np.median(vals)]
medians = np.array(medians)
print(medians)
with open(os.path.join(ml_folder, "medians.npy"), "wb") as f:
    np.save(f, medians, allow_pickle=False)
X = impute(X, medians)

print("Normalizing descriptor")
mus = []
sigmas = []
for i in range(X.shape[1]):
    vals = X[:,i]
    mus += [np.mean(vals)]
    sigmas += [np.std(vals)]
mus = np.array(mus)
sigmas = np.array(sigmas)
with open(os.path.join(ml_folder, "mus.npy"), "wb") as f:
    np.save(f, mus, allow_pickle=False)
with open(os.path.join(ml_folder, "sigmas.npy"), "wb") as f:
    np.save(f, sigmas, allow_pickle=False)
X = normalize(X, mus, sigmas)

sel = VarianceThreshold(threshold=0.0)
sel.fit(X)
joblib.dump(sel, os.path.join(ml_folder, "sel0.pkl"))
X = sel.transform(X)

from featurewiz import featurewiz
print(X.shape)
feat_names = ["feat{0}".format(i) for i in range(X.shape[1])]
df_ = pd.DataFrame(X, columns=feat_names)
df_["y"] = y

sel1_file = os.path.join(ml_folder, "sel1.npy")
if os.path.exists(sel1_file):
    with open(sel1_file, "rb") as f:
        sel_idxs = np.load(f)
else:
    out1, out2 = featurewiz(df_, corr_limit=0.9, target="y", verbose=0)
    sel_idxs = []
    fw = set(out1)
    for i, fn in enumerate(feat_names):
        if fn in fw:
            sel_idxs += [i]
    sel_idxs = np.array(sel_idxs)
    with open(os.path.join(ml_folder, "sel1.npy"), "wb") as f:
        np.save(f, sel_idxs)
X = X[:,sel_idxs]
print(X.shape)

from flaml import AutoML
automl_settings = {
    "time_budget": 1000,  #  in seconds
    "task": "classification",
    "log_file_name": "automl.log",
}
automl = AutoML()
automl.fit(X_train=X, y_train=y, **automl_settings)
meta = {
    "estimator": automl.best_estimator,
    "loss": automl.best_loss,
}
print(meta)

mdl = automl.model
mdl.fit(X, y)
# remove catboost info folder if generated
cwd = os.getcwd()
catboost_info = os.path.join(cwd, "catboost_info")
if os.path.exists(catboost_info):
    shutil.rmtree(catboost_info)
if os.path.exists("automl.log"):
    os.remove("automl.log")

estimator = automl.model

from sklearn.metrics import auc
from sklearn.metrics import roc_curve

print("Cross-validation")
cv = StratifiedKFold(n_splits=10)
aucs = []
for i, (train, test) in tqdm(enumerate(cv.split(X, y))):
    estimator.fit(X[train], y[train])
    fpr, tpr, _ = roc_curve(y[test], estimator.predict_proba(X[test])[:,1])
    aucs += [auc(fpr, tpr)]
print(np.mean(aucs))

print("Final fit")
estimator.fit(X, y)
joblib.dump(estimator, os.path.join(ml_folder, "classifier.pkl"))
