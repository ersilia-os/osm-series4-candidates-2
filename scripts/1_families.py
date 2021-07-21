from __init__ import OUTPUT

from tqdm import tqdm
import pandas as pd
from rdkit import Chem
import os

#substructures of interest:
#Triazolo core
pattern_triazolopyrazine = Chem.MolFromSmarts("*-c1nnc2cncc(-*)n12")

#Right Hand Side substituent with phenyl ring containing heteroatoms (or not)
pattern_triazoloheteroaryl = Chem.MolFromSmarts("*-c1cncc2nnc(n12)-[*]:1:[*]:[*]:[*]:[*]:[*]:1")

#RHS substituent with phenyl ring (no heteroatoms)
pattern_triazolophenyl = Chem.MolFromSmarts("[*]-c1cncc2nnc(-c3ccccc3)n12")

#RHS substituent with phenyl ring (can include heteroatoms) and para substituent on the RHS phenyl
pattern_triazoloheteroaryl_para = Chem.MolFromSmarts("*-c1cncc2nnc(n12)-[*]:1:[*]:c:[*](-[*]):[*]:[*]:1")

#RHS substituent with phenyl ring (can include heteroatoms) and para substituent on the RHS phenyl
pattern_triazoloheteroaryl_meta = Chem.MolFromSmarts("*-c1cncc2nnc(n12)-[*]:1:[*]:[*]:[*]:[*](-[*]):[*]:1")

#RHS substituent with phenyl ring (can include heteroatoms) and para substituent on the RHS phenyl
pattern_triazoloheteroaryl_orto = Chem.MolFromSmarts("*-c1cncc2nnc(n12)-[*]:1:[*]:[*]:[*]:[*]:[*]:1-[*]")

#RHS substituent with a napthalene
pattern_triazolonaphthalene = Chem.MolFromSmarts("[*]-c1cncc2nnc(-c3ccc4ccccc4c3)n12")

#LHS substituent has an ether
pattern_pyrazineether = Chem.MolFromSmarts("[*]-[#8]-c1cncc2nnc(-[*])n12")

#LHS substituent has an amide
pattern_pyrazineamide = Chem.MolFromSmarts("[*]-[#7]-[#6](=O)-c1cncc2nnc(-[*])n12")

def is_triazolopyrazine(mol):
    if mol.HasSubstructMatch(pattern_triazolopyrazine):
        return 1
    else:
        return 0

def is_triazoloheteroaryl(mol):
    if mol.HasSubstructMatch(pattern_triazoloheteroaryl):
        return 1
    else:
        return 0

def is_triazolophenyl(mol):
    if mol.HasSubstructMatch(pattern_triazolophenyl):
        return 1
    else:
        return 0

def is_triazoloheteroaryl_para(mol):
    if mol.HasSubstructMatch(pattern_triazoloheteroaryl_para):
        return 1
    else:
        return 0

def is_triazoloheteroaryl_meta(mol):
    if mol.HasSubstructMatch(pattern_triazoloheteroaryl_meta):
        return 1
    else:
        return 0

def is_triazoloheteroaryl_orto(mol):
    if mol.HasSubstructMatch(pattern_triazoloheteroaryl_orto):
        return 1
    else:
        return 0

def is_triazolonaphthalene(mol):
    if mol.HasSubstructMatch(pattern_triazolonaphthalene):
        return 1
    else:
        return 0

def is_pyrazineether(mol):
    if mol.HasSubstructMatch(pattern_pyrazineether):
        return 1
    else:
        return 0

def is_pyrazineamide(mol):
    if mol.HasSubstructMatch(pattern_pyrazineamide):
        return 1
    else:
        return 0


df=pd.read_csv(os.path.join(OUTPUT, "data_0.csv"))
smiles = list(df["Smiles"])

mols = [Chem.MolFromSmiles(smi) for smi in tqdm(smiles)]

df["IsTriazoloPyrazine"] = [is_triazolopyrazine(mol) for mol in tqdm(mols)]
df["TriazoloHeteroaryl"] = [is_triazoloheteroaryl(mol) for mol in tqdm(mols)]
df["TriazoloPhenyl"] = [is_triazolophenyl(mol) for mol in tqdm(mols)]
df["TriazoloHeteroarylPara"] = [is_triazoloheteroaryl_para(mol) for mol in tqdm(mols)]
df["TriazoloHeteroarylMeta"] = [is_triazoloheteroaryl_meta(mol) for mol in tqdm(mols)]
df["TriazoloHeteroarylOrto"] = [is_triazoloheteroaryl_orto(mol) for mol in tqdm(mols)]
df["TriazoloNaphthalene"] = [is_triazolonaphthalene(mol) for mol in tqdm(mols)]
df["PyrazineEther"] = [is_pyrazineether(mol) for mol in tqdm(mols)]
df["PyrazineAmide"] = [is_pyrazineamide(mol) for mol in tqdm(mols)]

df = df[df["IsTriazoloPyrazine"] == 1]
df.to_csv(os.path.join(OUTPUT, "data_1.csv"), index=False)
