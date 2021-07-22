from __init__ import OUTPUT

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from tqdm import tqdm

import collections


print("SCAFFOLDS")

df = pd.read_csv(os.path.join(OUTPUT, "data_6.csv"))
values = list(df["IC50Pred"])
smiles = list(df["Smiles"])

mols = [Chem.MolFromSmiles(smi) for smi in tqdm(smiles)]

print("Getting scaffolds")
scaffolds = [MurckoScaffold.GetScaffoldForMol(mol) for mol in tqdm(mols)]
print("... as inchikey")
scaffolds = [Chem.inchi.MolToInchiKey(sc) for sc in scaffolds]

sc2val = collections.defaultdict(list)
for i, sc in enumerate(scaffolds):
    sc2val[sc] += [(i, values[i])]

print("Sorting by value (lower first)")
sc2idx = {}
for k,v in sc2val.items():
    v = sorted(v, key=lambda x: x[1])
    idx = v[0][0]
    sc2idx[k] = idx

all_idxs = [i for i in range(df.shape[0])]
sel_idxs = set([v for _,v in sc2idx.items()])

print(len(sc2val))
print(df.shape)

df["Index"] = all_idxs
df = df[df["Index"].isin(sel_idxs)]
df.drop(columns=["Index"], inplace=True)

print("Saving")
df.to_csv(os.path.join(OUTPUT, "data_7.csv"), index = False)
