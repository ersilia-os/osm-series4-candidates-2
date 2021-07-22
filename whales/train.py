from whales.ChemTools import prepare_mol_from_sdf
from whales import do_whales
import whales.ChemTools as tools

import os, sys
import pandas as pd
from rdkit.Chem import PandasTools
from tqdm import tqdm
import joblib
import numpy as np

from sklearn.preprocessing import RobustScaler

ROOT = os.path.dirname(os.path.abspath(__file__))

if os.path.exists("X.npy"):
    with open("X.npy", "rb") as f:
        X = np.load(f)
else:

    data_file = os.path.join(ROOT, "..", "predictor", "data", "data.csv")

    df = pd.read_csv(data_file)
    df = df[df["Activity"] < 1]

    print(df)
    print(df.shape)

    SDF_FILE = "act.sdf"

    print("Writing as SDF")
    # Smiles as SDF
    PandasTools.AddMoleculeColumnToFrame(df,"Smiles",'molecule') # pp = doesn't work for me
    PandasTools.WriteSDF(df, SDF_FILE, molColName='molecule', properties=list(df.columns))

    print("Reading SDF")
    mols = prepare_mol_from_sdf(SDF_FILE) # computes 3D geometry from a specified sdf file

    print("Calculating WHALES")

    whales_library = []
    for mol in tqdm(mols): # runs over the library and updates WHALES
        whales_temp, lab = do_whales.whales_from_mol(mol)
        whales_library.append(whales_temp)
    # convert the arrays into a pandas dataframe
    df_whales_library = pd.DataFrame(whales_library,columns=lab)
    df_whales_library.head() # library preview

    X = np.array(df_whales_library)

    scaler = RobustScaler()
    scaler.fit(X)
    X = scaler.transform(X)

    with open("X.npy", "wb") as f:
        np.save(f, X)

    joblib.dump(scaler, "scaler.pkl")


from sklearn.neighbors import NearestNeighbors

nn = NearestNeighbors(n_neighbors=4).fit(X)
distances, _ = nn.kneighbors(X)

joblib.dump(nn, "nn.pkl")

import matplotlib.pyplot as plt
plt.hist(distances[:,1], 100)
plt.hist(distances[:,2], 100)
plt.hist(distances[:,3], 100)

plt.savefig("dists.png", dpi=300)
