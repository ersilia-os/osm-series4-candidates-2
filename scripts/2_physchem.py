from __init__ import OUTPUT

from tqdm import tqdm
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Descriptors import qed
from rdkit.Chem.Descriptors import MolLogP
from rdkit.Chem.rdMolDescriptors import CalcNumRings
from rdkit.Chem.Descriptors import MolWt

MAX_MW = 600
MAX_SLOGP = 5


df = pd.read_csv(os.path.join(OUTPUT, "data_1.csv"))

mols = [Chem.MolFromSmiles(smi) for smi in list(df["Smiles"])]

df["MolWt"] = [MolWt(mol) for mol in mols]
df["SLogP"] = [MolLogP(mol) for mol in mols]
df["NumRings"] = [CalcNumRings(mol) for mol in mols]
df["Qed"] = [qed(mol) for mol in mols]

df = df[df["MolWt"] <= MAX_MW]
df = df[df["SLogP"] <= MAX_SLOGP]

df.to_csv(os.path.join(OUTPUT, "data_2.csv"), index=False)
