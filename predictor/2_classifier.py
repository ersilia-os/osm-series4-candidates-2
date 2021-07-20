import pandas as pd
import numpy as np
import shutil
import joblib
import json
import os
from tqdm import tqdm
from descriptors.utils import impute, normalize

from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest

from sklearn import svm
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold

from sklearn.metrics import auc
from sklearn.metrics import roc_curve

ROOT = os.path.dirname(os.path.abspath(__file__))

MAX_FEAT_INIT = 100
RFE_STEPS = 5
CUTOFF = 2.5


use_descriptors = [
    "rdkit2d",
]

TO_NORMALIZE = {
    "rdkit2d",
    "mordred"
}

# Read activity data and binarize

y_ = np.array(pd.read_csv(os.path.join(ROOT, "data", "data.csv"))["Activity"])
y = np.zeros(len(y_))
y[y_ <= 2.5] = 1
print(np.sum(y), "actives")

def classifier_pipeline(descriptor_name):
    ml_folder = os.path.join(ROOT, "models", descriptor_name)
    print("Reading descriptor", descriptor_name)
    with open(os.path.join(ml_folder, "X.npy"), "rb") as f:
        X = np.load(f)
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
    with open(os.path.join(ml_folder, "medians.npy"), "wb") as f:
        np.save(f, medians, allow_pickle=False)
    print("Impute if necessary")
    X = impute(X, medians)
    if descriptor_name in TO_NORMALIZE:
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
    print("Filter zero variance")
    from sklearn.feature_selection import VarianceThreshold
    sel = VarianceThreshold(threshold=0.0)
    sel.fit(X)
    joblib.dump(sel, os.path.join(ml_folder, "sel0.pkl"))
    X = sel.transform(X)
    if X.shape[1] > MAX_FEAT_INIT*5:
        print("PCA")
        n_components = int(np.min([MAX_FEAT_INIT, X.shape[1]]))
        sel = PCA(n_components=n_components)
        sel.fit(X)
        joblib.dump(sel, os.path.join(ml_folder, "sel1.pkl"))
        X = sel.transform(X)
    else:
        print("KBest")
        n_features = int(np.min([MAX_FEAT_INIT, X.shape[1]]))
        sel = SelectKBest(k=n_features)
        sel.fit(X, y)
        joblib.dump(sel, os.path.join(ml_folder, "sel1.pkl"))
        X = sel.transform(X)
    print("Recursive feature elimination")
    estimator = svm.SVC(kernel='linear', probability=True,
                        random_state=42, class_weight="balanced")
    sel = RFECV(estimator, step=2, min_features_to_select=10, cv=5, verbose=3)
    sel.fit(X, y)
    joblib.dump(sel, os.path.join(ml_folder, "sel2.pkl"))
    X = sel.transform(X)
    print("Cross-validation")
    cv = StratifiedKFold(n_splits=10)
    estimator = svm.SVC(kernel='linear', probability=True,
                        random_state=42, class_weight="balanced")
    aucs = []
    for i, (train, test) in tqdm(enumerate(cv.split(X, y))):
        estimator.fit(X[train], y[train])
        fpr, tpr, _ = roc_curve(y[test], estimator.predict_proba(X[test])[:,1])
        aucs += [auc(fpr, tpr)]
    print(np.mean(aucs))

    print("Final fit")
    estimator.fit(X, y)
    joblib.dump(estimator, os.path.join(ml_folder, "classifier.pkl"))


for desc_name in use_descriptors:
    classifier_pipeline(desc_name)
