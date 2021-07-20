import os
import pandas as pd
import random
import numpy as np
from standardiser import standardise
from rdkit import Chem

ROOT = os.path.dirname(os.path.abspath(__file__))

# Read molecules with activity

df = pd.read_csv(os.path.join(ROOT, "../data", "raw", "series4_processed.csv"))

y = np.array(df["activity"])
smiles = list(df["smiles"])
print(y.shape)
print(len(smiles))


class StandardizeSmiles(object):

    def __init__(self):
        pass

    def _smiles(self, mol):
        return Chem.MolToSmiles(mol)

    def _inchikey(self, mol):
        return Chem.inchi.MolToInchiKey(mol)

    def standardize(self, smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        try:
            mol = standardise.run(mol)
        except:
            return None
        if mol is None:
            return None
        ik = self._inchikey(mol)
        smi = self._smiles(mol)
        return ik, smi


std = StandardizeSmiles()

ik2smi = {}
ik2act = {}
for smi, y_ in zip(smiles, y):
    mol = std.standardize(smi)
    if mol is None:
        print(mol)
        continue
    ik, smi = mol[0], mol[1]
    print(ik, y_)
    if ik in ik2smi:
        continue
    else:
        ik2smi[ik] = smi
    if ik in ik2act:
        if ik2act[ik] > y_:
            continue
        else:
            ik2act[ik] = y_
    else:
        ik2act[ik] = y_


R = []
for k, v in ik2smi.items():
    R += [[k, v, ik2act[k]]]

df = pd.DataFrame(R, columns=["InchiKey", "Smiles", "Activity"])

df.to_csv(os.path.join(ROOT, "data", "data.csv"), index=False)
