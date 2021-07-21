from __init__ import OUTPUT

import os
import shutil
import csv
import pandas as pd
from rdkit import Chem
from standardiser import standardise
from tqdm import tqdm

print("AGGREGATE")

# Output directory

if os.path.exists(OUTPUT):
    shutil.rmtree(OUTPUT)
os.mkdir(OUTPUT)

reinvent_dir="../data/reinvent"
mollib_dir="../data/mollib"

all_smiles=[]

for root, dirs, files in os.walk(reinvent_dir):
    for filename in files:
        path=os.path.join(root,filename)
        df=pd.read_csv(path)
        smiles=df["SMILES"].tolist()
        all_smiles += smiles
for root, dirs, files in os.walk(mollib_dir):
    for filename in files:
        if filename[-3:] == "txt":
            path=os.path.join(root,filename)
            with open(path, "r") as f:
                for line in f:
                    line=line.rstrip("\n")
                    all_smiles += [line]

all_smiles_file = os.path.join(OUTPUT, "all_smiles.txt")
textfile = open(all_smiles_file, "w")
for element in all_smiles:
    textfile.write(element + "\n")
textfile.close()

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



# Iterate over smiles

std = StandardizeSmiles()

ik2smiles = {}
for smi in tqdm(all_smiles[:10000]):
    mol = std.standardize(smi)
    if mol is None:
        continue
    ik2smiles[mol[0]] = mol[1]

# Known s4 molecules

s4smiles=pd.read_csv("../data/raw/series4_allsmiles.csv")
s4smiles_list=s4smiles["canonical"].tolist()
s4inchikey=[]
for smi in (s4smiles_list):
    s4mol=std.standardize(smi)
    if s4mol is None:
        continue
    s4inchikey.append(s4mol[0])
s4inchikey = set(s4inchikey)

# Save as dataframe

R = []
for k,v in ik2smiles.items():
    if k not in s4inchikey:
        R += [[k, v]]

df = pd.DataFrame(R, columns = ["InchiKey", "Smiles"])

df.to_csv(os.path.join(OUTPUT, "data_0.csv"), index=False)
