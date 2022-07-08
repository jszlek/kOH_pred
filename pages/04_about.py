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
st.sidebar.markdown("# About")