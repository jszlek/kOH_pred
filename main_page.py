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
st.write("# Welcome to Streamlit! 👋")
st.markdown(
    """
    Streamlit is an open-source app framework built specifically for
    Machine Learning and Data Science projects.
    **👈 Select a demo from the sidebar** to see some examples
    of what Streamlit can do!
    ### Want to learn more?
    - Check out [streamlit.io](https://streamlit.io)
    - Jump into our [documentation](https://docs.streamlit.io)
    - Ask a question in our [community
      forums](https://discuss.streamlit.io)
    ### See more complex demos
    - Use a neural net to [analyze the Udacity Self-driving Car Image
      Dataset](https://github.com/streamlit/demo-self-driving)
    - Explore a [New York City rideshare dataset](https://github.com/streamlit/demo-uber-nyc-pickups)
"""
)
#    )


# if __name__ == "__main__":
#     run()




# import streamlit as st
# from rdkit.Chem import AllChem as Chem
# from rdkit.Chem import Draw
# import utils as utl
# from views import home, about, analysis, options, configuration
#
# st.set_page_config(layout="wide", page_title='Reaction constant prediction')
# st.set_option('deprecation.showPyplotGlobalUse', False)
# utl.inject_custom_css()
# utl.navbar_component()
#
#
# def navigation():
#     route = utl.get_current_route()
#     if route == "home":
#         home.load_view()
#         sidebar_std()
#     elif route == "about":
#         about.load_view()
#         sidebar_std()
#     elif route == "analysis":
#         analysis.load_view()
#     elif route == "options":
#         options.load_view()
#     elif route == "configuration":
#         configuration.load_view()
#     elif route is None:
#         home.load_view()
#         sidebar_std()
#
#
# def sidebar_std():
#     # Using object notation
#     file_picker = st.sidebar.file_uploader("Choose a file")
#
#     if file_picker is not None:
#         # To read file as bytes:
#         bytes_data = file_picker.getvalue()
#         st.write(bytes_data)
#
#
# navigation()
#
#
