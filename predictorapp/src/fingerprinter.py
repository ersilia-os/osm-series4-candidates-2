import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

radius = 3
useCounts = True
useFeatures = True


class Ecfp(object):
    def __init__(self):
        self.name = "ecfp"
        self.radius = radius
        self.useCounts = useCounts
        self.useFeatures = useFeatures

    def calc(self, mols):
        fps = [
            AllChem.GetMorganFingerprint(
                mol, self.radius, useCounts=self.useCounts, useFeatures=self.useFeatures
            )
            for mol in mols
        ]
        size = 2048
        nfp = np.zeros((len(fps), size), np.int32)
        for i, fp in enumerate(fps):
            for idx, v in fp.GetNonzeroElements().items():
                nidx = idx % size
                nfp[i, nidx] += int(v)
        return np.array(nfp, dtype=np.int32)


def get_fingerprints(smiles):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    fingerprinter = Ecfp()
    return fingerprinter.calc(mols)


def mols_to_fingerprints(molecules, radius=3, useCounts=True, useFeatures=True):
    fingerprints = [AllChem.GetMorganFingerprint(
        mol,
        radius,
        useCounts=useCounts,
        useFeatures=useFeatures
    ) for mol in molecules]
    return fingerprints

def ra_fingerprint(mol):
    """
    Converts SMILES into a counted ECFP6 vector with features.
    :param smiles: SMILES representation of the moelcule of interest
    :type smiles: str
    :return: ECFP6 counted vector with features
    :rtype: np.array
    """
    fp = AllChem.GetMorganFingerprint(mol, 3, useCounts=True, useFeatures=False)
    size = 2048
    arr = np.zeros((size,), np.int32)
    for idx, v in fp.GetNonzeroElements().items():
        nidx = idx % size
        arr[nidx] += int(v)
    return arr
    