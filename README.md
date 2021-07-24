## OSM Series 4 Candidates with Deep Generative Models - Round 2

A new round of series 4 candidates for the [Open Source Malaria Project](https://github.com/opensourcemalaria), including molecules generated from low-data generative models (adapted from the [ETH Modlab](https://github.com/ETHmodlab/virtual_libraries)) and molecules generated in a second round using the [Reinvent 2.0](https://github.com/MolecularAI/Reinvent) generative model with improved activity predictors. See Open Source Malaria discussion [#34](https://github.com/OpenSourceMalaria/Series4_PredictiveModel/issues/34).

We have tried to complement OSM issue [#29](https://github.com/OpenSourceMalaria/Series4_PredictiveModel/issues/34) opened by Evariste Technologies. In that issue, a specific region of the chemical space is **exploited** to identify highly active compounds. Here, we provide a rougher **exploration** of the chemical space with the hope to identify alternatie lead compounds:

![exploration_vs_exploitation](https://user-images.githubusercontent.com/19725330/126868711-f834d617-6a0d-44c5-927f-11abd36541b7.png)


## Data

- All 405,766 molecules generated (with duplicates eliminated) can be found here: [data_0.csv](https://github.com/ersilia-os/osm-series4-candidates-2/blob/main/scripts/results/data_0.csv)
- A selection of the best 557 candidates according to the pipeline below, rendered the following molecules: [data_13.csv](https://github.com/ersilia-os/osm-series4-candidates-2/blob/main/scripts/results/data_13.csv)
- A final list of the **best** 90 candidates based on activity can be found here: [eosi_s4_candidates_90.csv](https://github.com/ersilia-os/osm-series4-candidates-2/blob/main/scripts/results/eosi_s4_candidates_90.csv)
- Explore the 90 candidates in [this app](https://share.streamlit.io/ersilia-os/osm-series4-candidates-2/main/app.py)!

## Results columns

Candidate molecules are listed along with the following columns:

### Identifiers

- EosID: ID number from Ersilia Open Source Initiative
- InchiKey
- Smiles

### Activity predictions

- IC50Pred: Activity prediction based on multiple descriptors and classical ML models. The lower the better. It is probably biased towards high values, so hopefully it is a conservative estimate. Includes confidence interval (Upper Bound (UB) and Lower Bound (LB))
- DeepActivity: Activity prediction based on deep learning models (Grover and ChemProp). The higher the better. It is a composite z-score between several deep learning scores (chemprop, grover; trained on classification and regression tasks). Includes confidence interval (Upper Bound (UB) and Lower Bound (LB))
- Maip: Blood-stage antimalarial activity prediction using the [MAIP tool](https://www.ebi.ac.uk/chembl/maip/)

### Applicability domain

- WhalesDist3Act: WHALES descriptors distance to the top-3 actives in the training set. These
descriptors are used for scaffold hopping
- Similarity: Tanimoto similarity to known series 4 compounds

### Accessibility:

- SAScore: Synthetic accessibility. The lower the better.
- RAScore: Retrosynthetic accessibility as predicted by the [Reymond lab](https://github.com/reymond-group/RAscore). The higher the better.
- SybaScore: Fragment-based accessibility score by [Voršilák et al](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00439-2). The higher the better.

### Physicochemical properties

- SLogP: Solubility. Should be < 5.
- QED: Drug-likeness. The higher the better.
- NumRings: Number of rings in the molecule.
- FractionCSP3: Number of tertiary carbons.
- FrHalogen: Number of halogen groups.
- HeavyAtom: Heavy atom count.
- Rotatable: Number of rotatable bonds.
- Heteroatoms: Number of Heteroatoms.
- FrAlkylHalide: Fragments containing Halides.

### Chemotype

- TriazoloHeteroaryl: Contains an heteroaryl ring in the RHS.
- TriazoloPhenyl: Contains a phenyl (no heteroatoms) in the RHS.
- TriazoloHeteroaryl - Para / - Meta / - Orto: Contain substituents in para, meta or orto positions

## Molecule generation steps

1. A first batch of molecules were generated in May2021 using [Reinvent 2.0](https://github.com/MolecularAI/Reinvent). A detailed explanation as well as results analysis of this first round can be found in our GitHub repo [ersilia-os/osm-series4-candidates](https://github.com/ersilia-os/osm-series4-candidates). We generated 116,728 new series 4 candidates. [Code](https://github.com/ersilia-os/osm-series4-candidates).
2. A second batch of molecules (209,310) has been generated for the purposes of this analysis using Reinvent 2.0 in exploration mode and optimizing for activity based on a simple QSAR model build with RDKIT descriptors. [Code](https://drive.google.com/drive/folders/1YDfnBz8EEKw5bB6q5htM0cVbOf8A8N_n?usp=sharing).
3. A third batch of molecules (150,365) has been generated using a [low-data generative model](https://github.com/ETHmodlab/virtual_libraries) taking as pre-training populations the ChEMBL and a large fragments library. [Code](https://drive.google.com/drive/folders/1YCj5l2jzpXlyB6aPK5JjgTyfzttNbC3k?usp=sharing).

## Selection of best candidates

All unique final molecules (405,766) have undergone a recursive selection process based on physicochemical properties, synthetic accessibility and predicted activity as follows:

![](images/selection01.png)

Scripts for the filterings applied can be found in the `scripts` folder in this repository.

## Run pipeline

For transparency and reproducibility, we provide code to run the full pipeline for candidate selection. Please download and uncompress the following folders and files:

* [chemprop](https://drive.google.com/file/d/1WDN3NRTC4T98f-6St9YT8wDXO8foZOg5/view?usp=sharing)
* [grover](https://drive.google.com/file/d/11_zSh1635KcP6GGgiVTozmE96A1N-z-U/view?usp=sharing)
* [predictorapp](https://drive.google.com/file/d/1skShCUFMrpkLFJvYqvxbQU5DpPsC86Ii/view?usp=sharing)
* [syba.pkl](https://drive.google.com/file/d/1tPA1vprB7gEwxMy_25Cz_PqDIEtzOBXK/view?usp=sharing) (save it in `utils/syba.pkl`)
* [ra_model.onnx](https://drive.google.com/file/d/1x_Y5oOZOnxkb1hHjs8B9a8wlbF8izlGf/view?usp=sharing) (save it in `utils/ra_model.onnx`)

The notebook with the process to select the best 90 candidates can be found [here](https://deepnote.com/project/Open-Source-Malaria-Series-4-Round-2-Zq8tjyh_Q4qjsK0NKdSk0A/%2Feosi_s4_candidates_90.csv).
