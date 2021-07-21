from __init__ import OUTPUT

from tqdm import tqdm
import pandas as pd
from rdkit import Chem
import os, sys

print("ACCESSIBILITY SCORES")

# SA Score

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
import utils.SA_Score.sascorer as sascorer

MAX_SASCORE = 6
df = pd.read_csv(os.path.join(OUTPUT, "data_2.csv"))

smiles = df["Smiles"].tolist()

mols = [Chem.MolFromSmiles(smi) for smi in tqdm(smiles)]

sa = [sascorer.calculateScore(mol) for mol in tqdm(mols)]

df["SAScore"] = sa

# RA Score

MIN_RASCORE = 0.9

import numpy as np
from rdkit.Chem import AllChem
import onnxruntime as rt


def ra_fingerprint(mol):
    """
    Converts SMILES into a counted ECFP6 vector with features.
    :param smiles: SMILES representation of the moelcule of interest
    :type smiles: str
    :return: ECFP6 counted vector with features
    :rtype: np.array
    """
    fp = AllChem.GetMorganFingerprint(mol, 3, useCounts=True, useFeatures=False)
    size = 2048
    arr = np.zeros((size,), np.int32)
    for idx, v in fp.GetNonzeroElements().items():
        nidx = idx % size
        arr[nidx] += int(v)
    return arr

def get_rascore_inference_session():
    sess = rt.InferenceSession("../utils/ra_model.onnx")
    return sess

RAScore = []

sess = get_rascore_inference_session()

for mol in tqdm(mols):
    ra_fp=ra_fingerprint(mol)
    ra_fp = np.array([ra_fp], dtype=np.float32)
    input_name = sess.get_inputs()[0].name
    label_name = sess.get_outputs()[1].name
    ra = sess.run([label_name], {input_name: ra_fp})[0][0][1]
    RAScore += [ra]

df["RAScore"]=RAScore

df = df[df["SAScore"] <= MAX_SASCORE]
df = df[df["RAScore"] >= MIN_RASCORE]


df.to_csv(os.path.join(OUTPUT, "data_3.csv"), index=False)
