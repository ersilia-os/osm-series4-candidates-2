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

df = pd.read_csv(os.path.join(OUTPUT, "data_9.csv"))

print(df.shape)

df.to_csv(os.path.join(OUTPUT, "data_10.csv"), index = False)
