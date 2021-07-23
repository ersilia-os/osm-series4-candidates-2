from __init__ import OUTPUT

from tqdm import tqdm
import pandas as pd
from rdkit import Chem
import os, sys
import numpy as np
import joblib

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(ROOT, "..", "whales"))

from whales.ChemTools import prepare_mol_from_sdf
from whales import do_whales
import whales.ChemTools as tools
from rdkit.Chem import PandasTools


print("HOPPING APPLICABILITY DOMAIN")

df=pd.read_csv(os.path.join(OUTPUT, "data_8.csv"))
print(df.shape)

pp = pd.DataFrame(df[["Smiles"]])

SDF_FILE = "tmp.sdf"

PandasTools.AddMoleculeColumnToFrame(pp,"Smiles",'molecule') # pp = doesn't work for me
PandasTools.WriteSDF(pp, SDF_FILE, molColName='molecule', properties=list(pp.columns))

print("Reading SDF")
mols = prepare_mol_from_sdf(SDF_FILE) # computes 3D geometry from a specified sdf file
os.remove(SDF_FILE)

print(len(mols))

print("Calculating WHALES")

whales_library = []
for mol in tqdm(mols): # runs over the library and updates WHALES
    whales_temp, lab = do_whales.whales_from_mol(mol)
    whales_library.append(whales_temp)
# convert the arrays into a pandas dataframe
df_whales_library = pd.DataFrame(whales_library,columns=lab)
df_whales_library.head() # library preview

X = np.array(df_whales_library)

scaler = joblib.load(os.path.join(ROOT, "..", "whales", "scaler.pkl"))

X = scaler.transform(X)

print("Nearest neighbors")

nn = joblib.load(os.path.join(ROOT, "..", "whales", "nn.pkl"))

distances, _ = nn.kneighbors(X)

df["WhalesDist3Act"] = distances[:,2]

print(df.shape)

WHALES_CUTOFF = 5

#df = df[df["WhalesDist3Act"] < WHALES_CUTOFF]
print(df.shape)

df.to_csv(os.path.join(OUTPUT, "data_9.csv"), index=False)
