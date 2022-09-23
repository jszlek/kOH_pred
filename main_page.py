# kOH_pred.py
# kOH_pred, prediction of aqueous OH kinetics of saturated alcohols.
# Copyright (C) 2022 Jakub Szlęk
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version. This program is
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details. You should have received a copy of
# the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

import streamlit as st


st.set_page_config(
    page_title="Reaction OH constant prediction",
    page_icon=":microscope:",
)

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

# page content
st.markdown("### Main page")
st.markdown('Application kOH_pred is used to predict an aqueous OH kinetics of saturated alcohols.')
st.markdown('This application is one of the results presented in publication: ')
st.markdown('_Temperature-dependent aqueous OH kinetics of saturated alcohols; new structure-activity relationship models and atmospheric lifetimes_')
st.markdown('by: _Bartłomiej Witkowski, Priyanka Jain, Jakub Szlęk, Beata Wileńska, and Tomasz Gierczak_')
st.markdown('')
st.markdown('Please proceed to other pages (sidebar)')
st.markdown('### Instructions')
HtmlFile = open("pages/Instructions.html", 'r', encoding='utf-8')
source_code = HtmlFile.read()
st.markdown(source_code,unsafe_allow_html=True)
st.write('### Summary')
st.markdown('kOH_pred, prediction of aqueous OH kinetics of saturated alcohols.')
st.markdown('Copyright (C) 2022 Jakub Szlęk')
st.markdown('This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.')
st.markdown('The engine of the application uses a model developed by [H2O AutoML](https://docs.h2o.ai/h2o/latest-stable/h2o-docs/automl.html) software.')
st.markdown('')


st.session_state.update(st.session_state)
