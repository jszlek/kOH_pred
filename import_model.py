import streamlit as st
import h2o
import pandas as pd
from PATHS import MODEL

# run H2O server
@st.cache
def run_server_and_load_model():
    h2o.init()
    my_path = MODEL
    my_model = h2o.load_model(my_path)
    return my_model


@st.cache
def load_test_data():
    my_data = pd.DataFrame(
        [[0, 1.01833781668656, 6.60688196451292, 0, 1.10679611650485, 1, 0, 278]],
        columns=['MPC10', 'GATS3c', 'EState_VSA2', 'AATS4dv', 'GATS3dv', 'AATS4d', 'ATSC7pe', 'Reaction_temp_K']
    )
    return my_data


# read standard data
my_model = run_server_and_load_model()
my_dummy_data = load_test_data()
