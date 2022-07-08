import h2o
import pandas as pd
import streamlit as st
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from import_model import my_model, my_dummy_data
from decimal import Decimal
from mordred import Calculator, PathCount, Autocorrelation, MoeType

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


def calculate_descriptors(smiles_format):
    mol_format = Chem.MolFromSmiles(smiles_format)
    # create descriptor calculator with all descriptors
    calc1 = Calculator()
    calc1.register(PathCount.PathCount(10, False, False, False))
    calc1.register(Autocorrelation.GATS(order=3, prop='c'))
    calc1.register(Autocorrelation.GATS(order=3, prop='dv'))
    calc1.register(Autocorrelation.AATS(order=4, prop='dv'))
    calc1.register(Autocorrelation.AATS(order=4, prop='d'))
    calc1.register(Autocorrelation.AATS(order=7, prop='pe'))
    calc1.register(MoeType.EState_VSA(2))
    desc = calc1(mol_format)
    return desc


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

    my_df = calculate_descriptors(smiles_input)
    # st.write(my_df.keys())
    # st.write(my_df.values())
    # st.write(my_df.asdict())

except ValueError:
    st.markdown('__Cannot process SMILES into MOL.__')
    st.markdown('__Wrong SMILES format.__')

# st.markdown(my_df)
my_df = pd.DataFrame([my_df.asdict()])
my_df['Reaction_temp_K'] = 298

st.dataframe(my_df)

# calculate dummy data of provided by user
res = my_model.predict(h2o.H2OFrame(my_df))
my_number = float(res.as_data_frame().iloc[0]*(10**9))
my_decimal = '%.4E' % Decimal(my_number)
st.markdown('## prediction: ' + str(my_decimal))


# # if data changes recalculate descriptors, form new dataframe and make predictions
# if st.session_state.smiles_key:
#     mol = Chem.MolFromSmiles(smiles_input)
#     im = Draw.MolToImage(mol)
#     st.image(im)
