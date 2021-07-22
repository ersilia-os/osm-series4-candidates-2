from __init__ import OUTPUT

from tqdm import tqdm
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import BulkTanimotoSimilarity
import os, sys
import networkx as nx
import collections

print("NON REDUNDANT SET")

CUTOFF = 0.75

def mols_to_fingerprints(molecules, radius=3, useCounts=False, useFeatures=True):
    fingerprints = [AllChem.GetMorganFingerprint(
        mol,
        radius,
        useCounts=useCounts,
        useFeatures=useFeatures
    ) for mol in tqdm(molecules)]
    return fingerprints

df = pd.read_csv(os.path.join(OUTPUT, "data_9.csv"))
values = list(df["IC50Pred"])

smiles=df["Smiles"].tolist()
mols = [Chem.MolFromSmiles(smi) for smi in tqdm(smiles)]
fps=mols_to_fingerprints(mols)

G = nx.Graph()

for i in range(len(fps)):
    G.add_node(i)

for i in tqdm(range(len(fps)-1)):
    fp = fps[i]
    sims = BulkTanimotoSimilarity(fp, fps[(i+1):])
    for k, j in enumerate(range(i+1, len(fps))):
        if sims[k] >= CUTOFF:
            G.add_edge(i,j)

cliques = collections.defaultdict(list)
for i, clique in enumerate(nx.find_cliques(G)):
    cliques[i] += clique

cl2val = collections.defaultdict(list)
for k,v in cliques.items():
    cl2val[k] = [(i, values[i]) for i in v]

print("Sorting by value (lower first)")
cl2idx = {}
for k,v in cl2val.items():
    v = sorted(v, key=lambda x: x[1])
    idx = v[0][0]
    cl2idx[k] = idx

all_idxs = [i for i in range(df.shape[0])]
sel_idxs = set([v for _,v in cl2idx.items()])

print(len(cl2val))
print(df.shape)

df["Index"] = all_idxs
df = df[df["Index"].isin(sel_idxs)]
df.drop(columns=["Index"], inplace=True)

print(df.shape)

df.to_csv(os.path.join(OUTPUT, "data_10.csv"), index = False)
