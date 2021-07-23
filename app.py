import streamlit as st
import os
import pandas as pd
import mols2grid
import streamlit.components.v1 as components

ROOT = os.path.dirname(os.path.abspath(__file__))
data = os.path.join(ROOT, "scripts", "results", "eosi_s4_candidates_90.csv")


st.set_page_config(layout="wide")# force wide display
st.title("Molecules for OSM Series 4 by Ersilia Open Source Initiative")
st.markdown("Open Source Malaria Issue [#34](https://github.com/OpenSourceMalaria/Series4_PredictiveModel/issues/34) | Learn about [Ersilia Open Source Initiative](https://ersilia.io)")

@st.cache
def load_data():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv(data)
    return df

df=load_data()

#reorder columns to show activity values first
df=df.round(3) #round activity to 4 decimals

columns = ["EosId", "IC50Pred", "DeepActivity", "Maip", "WhalesDist3Act", "Similarity", "SAScore", "RAScore", "SybaScore", "SLogP", "Qed"]

selection = st.selectbox("View Top 20 molecules by", columns)

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
        df_ = df.sort_values(col, ascending = asc).head(20)

if asc:
    st.text("Low is better")
else:
    st.text("High is better")

_cols = columns[:1] + columns[3:]
    
st.write(df_[_cols])

n_cols = 10
raw_html = mols2grid.display(df_, smiles_col = "Smiles",
subset=["img"], tooltip=["Smiles", "EosId"],
selection=False, n_cols=n_cols)._repr_html_()
components.html(raw_html, width=None, height=600, scrolling=True)
