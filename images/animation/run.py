from rdkit import Chem
from rdkit.Chem import TemplateAlign
import numpy as np

import pandas as pd

df = pd.read_csv("../osm-series4-candidates-2/scripts/results/data_13.csv")

vals = np.array(df["IC50Pred"])
idxs = np.argsort(vals)

smiles = list(df["Smiles"])
smiles = [smiles[i] for i in idxs]

mols = [Chem.MolFromSmiles(smi) for smi in smiles]
values = np.array([i for i in range(len(mols))])/len(mols)

from mol2svg import Mol2svg

draw = Mol2svg("out")
draw.mols2png(mols, values)
