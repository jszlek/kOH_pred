import streamlit as st
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw


def load_view():
    st.title('Home Page')
    smiles = 'COC(=O)c1c[nH]c2cc(OC(C)C)c(OC(C)C)cc2c1=O'
    mol = Chem.MolFromSmiles(smiles)
    im = Draw.MolToImage(mol)

    st.write('Hello World!')
    st.image(im)