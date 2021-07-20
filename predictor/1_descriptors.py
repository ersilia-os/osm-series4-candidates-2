import os
import numpy as np
import pandas as pd
import shutil

from rdkit import Chem

# Calculate descriptors

ROOT = os.path.dirname(os.path.abspath(__file__))

df = pd.read_csv(os.path.join(ROOT, "data", "data.csv"))
mols = [Chem.MolFromSmiles(smi) for smi in list(df["Smiles"])]

from descriptors.eosdescriptors.rdkit2d import Rdkit2d
from descriptors.eosdescriptors.ecfp import Ecfp
from descriptors.eosdescriptors.chembl import Chembl
from descriptors.eosdescriptors.signaturizer import Signaturizer


descs = [
    #Avalon(),
    #Cddd(),
    Chembl(),
    Ecfp(),
    #Mordred(),
    Rdkit2d(),
    #RdkitFpBits(),
    Signaturizer()
]

def calculate_descriptor(mols, desc):
    name = desc.name
    print("Working on", name)
    X = desc.calc(mols)
    print("Done")
    output_dir = os.path.join(ROOT, "models", name)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    output_file = os.path.join(output_dir, "X.npy")
    with open(output_file, "wb") as f:
        np.save(f, X)


for desc in descs:
    calculate_descriptor(mols, desc)
