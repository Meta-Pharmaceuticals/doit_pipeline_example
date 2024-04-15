import streamlit as st
import pandas as pd
import json
from pathlib import Path
import sys

def get_file(ps):
    return " ".join([Path(p).name for p in ps])

def load_history():
    data = json.load(open('./history'))
    col = ['status', 'script', 'input', 'target', 'params', 'message']
    df = pd.DataFrame.from_dict(data)
    df.set_index('id', inplace=True)
    df.sort_index(inplace=True, ascending = False)
    df = df[col]

    df['script'] = df['script'].map(lambda p: Path(p).name)
    df['input'] = df['input'].map(get_file)
    df['target'] = df['target'].map(get_file)
    df['params'] = df['params'].map(get_file)
    return df

st.subheader('Current Task Dependencies')
plotfile = Path('./tasks.png')
if plotfile.exists():
    st.image(str(plotfile.absolute()))

st.subheader('Current Task Definition')
code = open('./dodo.py').read()
st.code(code, language='python')

st.subheader('Current Task History')
history = load_history()
history.astype(str)
st.dataframe(history)





