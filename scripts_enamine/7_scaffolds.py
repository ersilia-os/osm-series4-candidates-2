from __init__ import OUTPUT

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from tqdm import tqdm

import collections


print("SCAFFOLDS")

df = pd.read_csv(os.path.join(OUTPUT, "data_6.csv"))

print("Saving")
df.to_csv(os.path.join(OUTPUT, "data_7.csv"), index = False)
