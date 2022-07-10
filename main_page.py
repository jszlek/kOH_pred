# kOH_pred.py
# License
# etc.

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
st.markdown("# Main page 🎈")
st.markdown('I am here!')
st.markdown('Main page content')
st.markdown('Please proceed to other pages (sidebar)')
st.write("### Summary of model used")



# import streamlit as st

# from SessionState import _get_state

# st.write(st.session_state['page_config'])


st.session_state.update(st.session_state)


# if 'page_config' not in st.session_state:
#     st.session_state.page_config.page_title = 'Reaction OH constant prediction'
#     st.session_state.page_config.page_icon = ":microscope:"

# state = st.session_state._get_state()
#
#
# state.page_config = st.set_page_config(
#     page_title="Reaction OH constant prediction",
#     page_icon=":microscope:",
# )
#
# # rest of the code
#
# state.sync()