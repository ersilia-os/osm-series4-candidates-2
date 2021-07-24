import streamlit as st
import os
import pandas as pd
import mols2grid
import streamlit.components.v1 as components

ROOT = os.path.dirname(os.path.abspath(__file__))
data = os.path.join(ROOT, "scripts", "results", "eosi_s4_candidates_90.csv")


st.set_page_config(layout="wide")# force wide display
st.title("OSM Series 4 candidates - Round 2")
st.markdown("Open Source Malaria Issue [#34](https://github.com/OpenSourceMalaria/Series4_PredictiveModel/issues/34) | Learn about [Ersilia Open Source Initiative](https://ersilia.io) | Pipeline [code](https://github.com/ersilia-os/osm-series4-candidates-2)")
st.markdown("Full dataset: [eosi_s4_candidates_90.csv](https://github.com/ersilia-os/osm-series4-candidates-2/blob/main/scripts/results/eosi_s4_candidates_90.csv)")

@st.cache
def load_data():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv(data)
    return df

df=load_data()

columns = ["EosId", "IC50Pred", "DeepActivity", "Maip", "WhalesDist3Act", "Similarity", "SAScore", "RAScore", "SybaScore", "SLogP", "Qed"]

selection = st.selectbox("View top 30 molecules by", columns)
more_cols = st.checkbox("More columns")

ascending = {
   "EosId",
   "IC50Pred",
   "WhalesDist3Act",
   "Similarity",
   "SAScore",
   "SLogP"
}

for col in columns:
    if col in ascending:
        asc = True
    else:
        asc = False
    if col == selection:
        df_ = df.sort_values(col, ascending = asc).head(30)
        break

descriptions = {
   "EosId": "ID number from Ersilia Open Source Initiative",
   "IC50Pred": "Activity prediction based on multiple descriptors and classical ML models. The lower the better. It is probably biased towards high values, so hopefully it is a conservative estimate. Includes confidence interval (Upper Bound (UB) and Lower Bound (LB))",
   "DeepActivity": "Activity prediction based on deep learning models (Grover and ChemProp). The higher the better. It is a composite z-score between several deep learning scores (chemprop, grover; trained on classification and regression tasks). Includes confidence interval (Upper Bound (UB) and Lower Bound (LB))",
   "Maip": "Blood-stage antimalarial activity prediction using the MAIP tool",
   "WhalesDist3Act": "WHALES descriptors distance to the top-3 actives in the training set. These descriptors are used for scaffold hopping" ,
   "Similarity": "Tanimoto similarity to known series 4 compounds",
   "SAScore": "Synthetic accessibility. The lower the better",
   "RAScore": "Retrosynthetic accessibility as predicted by the Reymond lab. The higher the better",
   "SybaScore": "Fragment-based accessibility score. The higher the better.",
   "SLogP": "RDKit Crippen's calculation",
   "Qed": "Drug-likeness. The higher the better",
}

st.text(descriptions[selection])
if asc:
    st.text("Sorted by {0} in ASCENDING order".format(selection))
else:
    st.text("Sorted by {0} in DESCENDING order".format(selection))

_cols = [c for c in list(df_.columns) if c not in {"InchiKey", "Smiles"}]
if not more_cols:
    _cols = [c for c in _cols if c in columns]

df_display = df_[_cols]

formats = dict((c, "{:.2f}") for c in _cols[1:])

try:
    st.dataframe(df_display.style.format(formats))
except:
    st.dataframe(df_display)

df_ = df_.round(2)

df_["Selection"] = ["{0}: {1}".format(selection, v) for v in list(df_[selection])]

n_cols = 6
n_rows = 3
raw_html = mols2grid.display(df_, smiles_col = "Smiles",
subset=["EosId", "img", "Selection"], tooltip=["Smiles", "EosId"],
n_cols = n_cols,
n_rows = n_rows,
selection=False)._repr_html_()
components.html(raw_html, width=None, height=600, scrolling=False)
