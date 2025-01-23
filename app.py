import streamlit as st
import pandas as pd
import math
from pathlib import Path

# Set the title and favicon that appear in the Browser's tab bar.
st.set_page_config(
    page_title='Protein Interaction Predictor',
    page_icon=':earth_americas:', # This is an emoji shortcode. Could be a URL too.
)

# -----------------------------------------------------------------------------
# Declare some useful functions.

@st.cache_data
def get_protein_data():
    """Grab GDP data from a CSV file.

    This uses caching to avoid having to read the file every time. If we were
    reading from an HTTP endpoint instead of a file, it's a good idea to set
    a maximum age to the cache with the TTL argument: @st.cache_data(ttl='1d')
    """

    # Instead of a CSV on disk, you could read from an HTTP endpoint here too.
    #DATA_FILENAME = Path(__file__).'github/gisthuband/protein_interaction_predictor/condensed_feature_engineered.csv'
    protein_df = pd.read_csv('https://github.com/gisthuband/Protein_Interaction_Predictor/blob/main/condensed_feature_engineered.csv')

   

    # The data above has columns like:
    ###protein 1 & 2:
    #sequence
    #len
    #phobic count
    #philic count
    #acidic count
    #basic count
    #aromatic count
    #sulfur count
    #and whether or not the proteins interact 
    #
    # So let's remove the sequence and extrac only numeric information
    

    # Convert years from string to integers
    

    return protein_df

usable_df = get_protein_data()

# -----------------------------------------------------------------------------
# Draw the actual page

# Set the title that appears at the top of the page.
'''
# :earth_americas: Protein Interaction Predictor (PIP)

This is my protein interaction predictor.  Just input protein 1's sequence and protein 2's sequence, and click 'submit'.

'''

# Add some spacing
''
''
usable_df.head()