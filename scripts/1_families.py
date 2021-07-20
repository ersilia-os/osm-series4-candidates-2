from __init__ import OUTPUT

from tqdm import tqdm
import pandas as pd
from rdkit import Chem

# Matching patterns

triazolo_patt = Chem.MolFromSmarts("*-c1nnc2cncc(-*)n12")
def is_triazolo(mol):
    if mol.HasSubstrMatch(triazolo_patt):
        return 1
    else:
        return 0

def is_heteroaryl(mol):
    if


def is_phenyl(mol):



df = pd.read_csv(os.path.join(OUTPUT, "data_0.csv"))
smiles = list(df["Smiles"])

mols = [Chem.MolFromSmiles(smi) for smi in tqdm(smiles)]

df["IsTriazolo"] = [is_triazolo(mol) for mol in mols]
df[""]

df = df[df["IsTriazolo"] == 1]
df.to_csv(os.path.join(OUTPUT, "data_1.csv"), index=False)
