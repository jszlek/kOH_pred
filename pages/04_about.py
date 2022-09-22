# 04_about.py
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

# page and sidebar style
# st.set_page_config(
#     page_title="Reaction OH constant prediction",
#     page_icon=":microscope:",
# )
# st.write(st.state.SessionState.filtered_state.fget)


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


st.markdown("# About")
st.markdown('Software developer and maintainer: Jakub Szlęk, <j.szlek@uj.edu.pl>')
HtmlFile = open("pages/GNU_GPL_License_v3.html", 'r', encoding='utf-8')
source_code = HtmlFile.read()
st.markdown(source_code,unsafe_allow_html=True)
