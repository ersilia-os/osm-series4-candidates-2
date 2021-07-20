from __init__ import OUTPUT

from tqdm import tqdm
import pandas as pd
from rdkit import Chem
import os

from rdkit.Chem.Descriptors import qed
from rdkit.Chem.Descriptors import MolLogP
from rdkit.Chem.rdMolDescriptors import CalcNumRings
from rdkit.Chem.Descriptors import MolWt

MAX_MW = 550
MAX_SLOGP = 5
MAX_NUMRINGS = 7
MIN_QED = 0.35

df = pd.read_csv(os.path.join(OUTPUT, "data_1.csv"))

mols = [Chem.MolFromSmiles(smi) for smi in tqdm(list(df["Smiles"]))]


df["MolWt"] = [MolWt(mol) for mol in tqdm(mols)]
df["SLogP"] = [MolLogP(mol) for mol in tqdm(mols)]
df["NumRings"] = [CalcNumRings(mol) for mol in tqdm(mols)]
df["Qed"] = [qed(mol) for mol in tqdm(mols)]


df = df[df["MolWt"] <= MAX_MW]
df = df[df["SLogP"] <= MAX_SLOGP]
df = df[df["NumRings"] <= MAX_NUMRINGS]
df = df[df["Qed"] >= MIN_QED]

df.to_csv(os.path.join(OUTPUT, "data_2.csv"), index=False)
