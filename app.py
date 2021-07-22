import streamlit as st
import os
import pandas as pd
import mols2grid
import streamlit.components.v1 as components

ROOT = os.path.dirname(os.path.abspath(__file__))

filt_step = st.sidebar.selectbox("Filtering step", ["1", "2", "3", "4", "5", "6", "7", "8"])
data_file = os.path.join(ROOT, "scripts", "results", "data_{0}.csv".format(filt_step))
df = pd.read_csv(data_file, nrows=10000)
df = df.sample(1000)

st.write(df)

n_cols = 10
raw_html = mols2grid.display(df.head(100), smiles_col = "Smiles",
subset=["img"], tooltip=["Smiles"],
selection=False, n_cols=n_cols)._repr_html_()
components.html(raw_html, width=None, height=600, scrolling=True)
