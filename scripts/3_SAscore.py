from __init__ import OUTPUT

from tqdm import tqdm
import pandas as pd
from rdkit import Chem
import os, sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
import utils.SA_Score.sascorer as sascorer

MAX_SASCORE = 6
df = pd.read_csv(os.path.join(OUTPUT, "data_2.csv"))

smiles=df["Smiles"].tolist()

mols = [Chem.MolFromSmiles(smi) for smi in tqdm(smiles)]

sa = [sascorer.calculateScore(mol) for mol in tqdm(mols)]

df["SAScore"] = sa

df = df[df["SAScore"] <= MAX_SASCORE]

df.to_csv(os.path.join(OUTPUT, "data_3.csv"), index=False)
