import h2o
import streamlit as st
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from import_model import my_model, my_dummy_data


# page and sidebar style
# st.set_page_config(
#     page_title="Reaction OH constant prediction",
#     page_icon=":microscope:",
# )

st.markdown(
    """
    <style>
    [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {
        width: 250px;
        height: 48rem;
        margin-bottom: 24rem;
    }
    [data-testid="stSidebar"][aria-expanded="false"] > div:first-child {
        width: 250px;
        margin-left: -250px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


# defined functions
def get_smiles():
    st.write(st.session_state.smiles_key)


# content
st.title('View and predict single molecule')

# input SMILES
smiles_input = st.text_input(label='Please enter SMILES',
                             value='COC(=O)c1c[nH]c2cc(OC(C)C)c(OC(C)C)cc2c1=O',
                             on_change=get_smiles,
                             key='smiles_key')
# print out current SMILES
st.write('The current SMILES is \n', smiles_input)

# process smiles into mol image
try:
    mol = Chem.MolFromSmiles(smiles_input)
    im = Draw.MolToImage(mol)
    st.image(im)
except ValueError:
    st.markdown('__Cannot process SMILES into MOL.__')
    st.markdown('__Wrong SMILES format.__')


# calculate dummy data of provided by user

#my_model = main_page.my_model
#my_dummy_data = main_page.my_dummy_data

res = my_model.predict(h2o.H2OFrame(my_dummy_data))

st.write(float(res.as_data_frame().iloc[0]))
