import streamlit as st
import os
import pandas as pd
import mols2grid
import streamlit.components.v1 as components

ROOT = os.path.dirname(os.path.abspath(__file__))
data = os.path.join(ROOT, "scripts", "results", "processed.csv")


st.set_page_config(layout="wide")# force wide display
st.title("Molecules for OSM Series 4 by Ersilia Open Source Initiative")
st.text("Learn about [Ersilia Open Source Initiative](https://ersilia.io)")

@st.cache
def load_data():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv(data)
    return df

df=load_data()

#reorder columns to show activity values first
df=df[["EosId", "InchiKey", "Smiles", "Activity", "Accessibility", "Similarity", "MolWt", "SLogP", "SAScore", "RAScore", "SybaScore", "IC50Pred", "IC50PredUB", "IC50PredLB", "HighClassifier", "WhalesDist3Act", "Maip", "Qed", "NumRings", "FractionCSP3", "FrHalogen", "HeavyAtom", "Rotatable", "Heteroatoms", "FrAlkylHalide", "TriazoloHeteroaryl", "TriazoloPhenyl", "TriazoloHeteroarylPara", "TriazoloHeteroarylMeta", "TriazoloHeteroarylOrto", "TriazoloNaphthalene", "PyrazineEther", "PyrazineAmide"]]
df=df.round(3) #round activity to 4 decimals

selection = st.sidebar.selectbox("View Top 20 molecules by", ["Activity", "Accessibility", "Similarity", "None"])

if selection == "Activity":
    df = df.sort_values('Activity', ascending = False).head(20)
elif selection == "Accessibility":
    df = df.sort_values('Accessibility', ascending = False).head(20)
elif selection == "Similarity":
    df = df.sort_values('Similarity', ascending = False).head(20)
else:
    df = df

st.write(df)

n_cols = 10
raw_html = mols2grid.display(df.head(100), smiles_col = "Smiles",
subset=["img"], tooltip=["Smiles", "EosId"],
selection=False, n_cols=n_cols)._repr_html_()
components.html(raw_html, width=None, height=600, scrolling=True)


filt_step = st.sidebar.selectbox("Filtering step", ["1", "2", "3", "4", "5", "6", "7", "8"])
data_file = os.path.join(ROOT, "scripts", "results", "data_{0}.csv".format(filt_step))
df2 = pd.read_csv(data_file, nrows=10000)
df2 = df2.sample(1000)

st.write(df2)

n_cols = 10
raw_html = mols2grid.display(df2.head(100), smiles_col = "Smiles",
subset=["img"], tooltip=["Smiles"],
selection=False, n_cols=n_cols)._repr_html_()
components.html(raw_html, width=None, height=600, scrolling=True)
