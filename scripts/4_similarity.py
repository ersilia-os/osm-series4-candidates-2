from __init__ import OUTPUT

from tqdm import tqdm
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import BulkTanimotoSimilarity
import os, sys

print("SIMILARITY SCORES")

def mols_to_fingerprints(molecules, radius=3, useCounts=False, useFeatures=True):
    fingerprints = [AllChem.GetMorganFingerprint(
        mol,
        radius,
        useCounts=useCounts,
        useFeatures=useFeatures
    ) for mol in tqdm(molecules)]
    return fingerprints

raw_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "raw")
sys.path.append(raw_folder)

#get series4 molecules for tanimoto similarity
s4 = pd.read_csv(os.path.join(raw_folder, "series4_processed.csv"))
s4_smiles = s4["smiles"].tolist()
s4_mols = [Chem.MolFromSmiles(smi) for smi in s4_smiles]
ref_fps=mols_to_fingerprints(s4_mols)


df = pd.read_csv(os.path.join(OUTPUT, "data_3.csv"))
smiles=df["Smiles"].tolist()
mols = [Chem.MolFromSmiles(smi) for smi in tqdm(smiles)]
fps=mols_to_fingerprints(mols)
sims = []
for fp in tqdm(fps):
    sim=BulkTanimotoSimilarity(fp, ref_fps)
    maxsim = np.max(sim)
    sims += [maxsim]

df["Similarity"]=sims
df=df[df["Similarity"] <= 0.70]

df.to_csv(os.path.join(OUTPUT, "data_4.csv"), index = False)
