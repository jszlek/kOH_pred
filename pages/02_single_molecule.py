import h2o
import pandas as pd
import streamlit as st
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from import_model import my_model
from decimal import Decimal
from mordred import Calculator, PathCount, Autocorrelation, MoeType


# st.markdown(f'''
#     <style>
#         section[data-testid="stSidebar"] .css-ng1t4o {{width: 10rem;}}
#         section[data-testid="stSidebar"] .css-1d391kg {{width: 10rem;}}
#     </style>
# ''',unsafe_allow_html=True)

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
    col1.write(st.session_state.smiles_key)


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
# create column layout
col1, col2 = st.columns(2)


# input SMILES
smiles_input = col1.text_input(label='Please enter SMILES',
                             value='CCCC(CCO)CCCCO',
                             on_change=get_smiles,
                             key='smiles_key')
# input temperature in K
temp_input = col1.slider(label='Please select temperature [K]',
                       min_value=278, max_value=328, value=298,
                       step=None, format=None, key=None, help=None,
                       on_change=None, args=None, kwargs=None, disabled=False)

# print out current SMILES
with col1:
    st.write('The current SMILES is \n', smiles_input)

# process smiles into mol image
try:
    mol = Chem.MolFromSmiles(smiles_input)
    im = Draw.MolToImage(mol, size=(300,200))
    with col1:
        st.image(im)
    my_df = calculate_descriptors(smiles_input)

except ValueError:
    with col1:
        st.markdown('__Cannot process SMILES into MOL.__')
        st.markdown('__Wrong SMILES format.__')

# st.markdown(my_df)
my_df = pd.DataFrame([my_df.asdict()])
my_df['Reaction_temp_K'] = temp_input

with col2:
    st.markdown('### Input variables: ')
    st.dataframe(my_df)

# calculate dummy data of provided by user
res = my_model.predict(h2o.H2OFrame(my_df))
my_number = float(res.as_data_frame().iloc[0]*(10**9))
my_decimal = '%.4E' % Decimal(my_number)

with col2:
    st.markdown('### Prediction: ' + str(my_decimal))


# # if data changes recalculate descriptors, form new dataframe and make predictions
# if st.session_state.smiles_key:
#     mol = Chem.MolFromSmiles(smiles_input)
#     im = Draw.MolToImage(mol)
#     st.image(im)
