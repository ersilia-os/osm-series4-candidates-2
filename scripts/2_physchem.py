from __init__ import OUTPUT

from tqdm import tqdm
import pandas as pd
from rdkit import Chem
import os

from rdkit.Chem.Descriptors import qed
from rdkit.Chem.Descriptors import MolLogP
from rdkit.Chem.rdMolDescriptors import CalcNumRings
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.Lipinski import FractionCSP3
from rdkit.Chem.Lipinski import HeavyAtomCount
from rdkit.Chem.Lipinski import NumRotatableBonds
from rdkit.Chem.Lipinski import NumHeteroatoms
from rdkit.Chem.Fragments import fr_halogen
from rdkit.Chem.Fragments import fr_alkyl_halide

PRINT("PHYSCHEM")

MIN_MW = 300
MAX_MW = 550
MAX_SLOGP = 5
MAX_NUMRINGS = 7
MIN_QED = 0.35
MIN_CSP3 = 0.1
MAX_FRHALOGEN = 5
MIN_HEAVYATOM = 20
MAX_HEAVYATOM = 40
MIN_ROTATABLE = 4
MIN_HETEROATOMS = 6
MAX_HETEROATOMS = 14
MAX_FRALKYLHALIDE = 5

df = pd.read_csv(os.path.join(OUTPUT, "data_1.csv"))

mols = [Chem.MolFromSmiles(smi) for smi in tqdm(list(df["Smiles"]))]

df["MolWt"] = [MolWt(mol) for mol in tqdm(mols)]
df["SLogP"] = [MolLogP(mol) for mol in tqdm(mols)]
df["NumRings"] = [CalcNumRings(mol) for mol in tqdm(mols)]
df["Qed"] = [qed(mol) for mol in tqdm(mols)]
df["FractionCSP3"] = [CalcFractionCSP3(mol) for mol in tqdm(mols)]
df["FrHalogen"] = [fr_halogen(mol) for mol in tqdm(mols)]
df["HeavyAtom"] = [HeavyAtomCount(mol) for mol in tqdm(mols)]
df["Rotatable"] = [NumRotatableBonds(mol) for mol in tqdm(mols)]
df["Heteroatoms"] = [NumHeteroatoms(mol) for mol in tqdm(mols)]
df["FrAlkylHalide"] = [fr_alkyl_halide(mol) for mol in tqdm(mols)]

df = df[(df["MolWt"] >= MIN_MW)&(df["MolWt"] <= MAX_MW)]
df = df[df["SLogP"] <= MAX_SLOGP]
df = df[df["NumRings"] <= MAX_NUMRINGS]
df = df[df["Qed"] >= MIN_QED]
df = df[df["FractionCSP3"] >= MIN_CSP3]
df = df[df["FrHalogen"] <= MAX_FRHALOGEN]
df = df[(df["HeavyAtom"] >= MIN_HEAVYATOM)&(df["HeavyAtom"] <= MAX_HEAVYATOM)]
df = df[df["Rotatable"] >= MIN_ROTATABLE]
df = df[(df["Heteroatoms"] >= MIN_HETEROATOMS)&(df["Heteroatoms"] <= MAX_HETEROATOMS)]
df = df[df["FrAlkylHalide"] <= MAX_FRALKYLHALIDE]

df.to_csv(os.path.join(OUTPUT, "data_2.csv"), index=False)
