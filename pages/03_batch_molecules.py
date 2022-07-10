import h2o
import base64
import pandas as pd
import streamlit as st
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from import_model import my_model
from io import BytesIO
from decimal import Decimal
from mordred import Calculator, PathCount, Autocorrelation, MoeType

# page and sidebar style
# st.set_page_config(
#     page_title="Reaction OH constant prediction",
#     page_icon=":microscope:",
# )


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


def make_predictions(smiles_format, temperature):

    # try:
    mol = Chem.MolFromSmiles(smiles_format)
    im = Draw.MolToImage(mol, size=(100, 50))
    my_df = calculate_descriptors(smiles_format)

    # except ValueError:
    #     print('__Cannot process SMILES into MOL.__')
    #     print('__Wrong SMILES format.__')

    my_df = pd.DataFrame([my_df.asdict()])
    my_df['Reaction_temp_K'] = temperature
    res = my_model.predict(h2o.H2OFrame(my_df))
    f_prediction = float(res.as_data_frame().iloc[0])
    return {'Image': im, 'SMILES': smiles_format, 'Prediction': f_prediction}


def image_base64(im):
    with BytesIO() as buffer:
        im.save(buffer, 'png')
        return base64.b64encode(buffer.getvalue()).decode()


def image_formatter(im):
    return f'<img src="data:image/png;base64,{image_base64(im)}">'

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


st.markdown("# Batch molecules")
st.markdown('Please upload data table with two columns: `SMILES` and `Reaction temperature`.')
st.markdown('Please see an example below:')

uploaded_file = st.file_uploader("Choose a file")

if uploaded_file is not None:
    # Load data table:
    dataframe = pd.read_csv(uploaded_file)
    # process dataframe rows into mol images and predictions
    # Iterating over two columns with use of `zip`

    prediction = [make_predictions(x, y) for x, y in zip(dataframe.iloc[:, 0], dataframe.iloc[:, 1])]
    # st.write(prediction)
    # result_dict = {'Image': image, 'SMILES': smiles, 'Prediction': prediction}
    result = pd.DataFrame.from_records(prediction)
    predict_dataframe = pd.DataFrame(result)
    st.write(predict_dataframe.to_html(formatters={'Image': image_formatter}, escape=False), unsafe_allow_html=True)
    st.download_button(
        label="Download data as CSV",
        data=predict_dataframe.to_csv().encode('utf-8'),
        file_name='result.csv',
        mime='text/csv',
    )
    # predict_dataframe.to_excel('tmp.xlsx')
