import streamlit as st
import pandas as pd
import math
from pathlib import Path
! pip install scikit-learn
from sklearn.model_selection import GridSearchCV
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn.metrics import classification_report, confusion_matrix, roc_curve, auc

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
    DATA_FILENAME = Path(__file__).parent/'data/condensed_feature_engineered.csv'
    protein_df = pd.read_csv(DATA_FILENAME)


    

    return protein_df

usable_df = get_protein_data()

st.write(usable_df.head())
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


user_text1 = st.text_area('input protein sequence 1')
st.write(user_text1)

user_text2 = st.text_area('input protein sequence 2')
st.write(user_text2)

seqs = [user_text1, user_text2]

''
''

def seq_features(seqs):
    
    #structure = {'protein_1_seq':'p', 'protein_1_len':'p', '1_phobic_count':'p', '1_philic_count':'p', '1_basic_count': 'p', '1_acidic_count': 'p', '1_aromatic_count':'p', '1_sulfur_count': 'p', 'protein_2_seq':'p', 'protein_2_len':'p', '2_phobic_count':'p', '2_philic_count':'p', '2_basic_count': 'p', '2_acidic_count': 'p', '2_aromatic_count':'p', '2_sulfur_count': 'p'}
    #new_df = pd.DataFrame(data=structure, index=[0])
    
        
    p1_list = []
    p2_list = []
        
    for i in range(2):
            
        amino_acids = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q', 'R','S','T','V','W','Y']
    
        ind_counts = {}
        for a in amino_acids:
            ind_counts[a] = 0

        phobic = {'A':0, 'F':0, 'G':0, 'I':0, 'L':0, 'M':0, 'P':0, 'V':0, 'W':0, 'Y':0, 'phobic_total':0}
        philic = {'C':0, 'D':0, 'E':0, 'H':0, 'K':0, 'N':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'philic_total':0}
        basic = {'H':0, 'K':0, 'R':0, 'basic_total':0}
        acidic = {'D':0, 'E':0, 'acidic_total':0}
        aromatic = {'F':0, 'H':0, 'W':0, 'Y':0, 'aromatic_total':0}
        sulfur = {'C':0 , 'M':0, 'sulfur_total':0}
            
            
        seq = seqs[i]
        seq_split = [y for y in seq]
        seq_len = len(seq_split)
            
        for z in seq_split:
                
            if z in ind_counts:
                ind_counts[z] += 1
                    
            if z in phobic:
                phobic[z] += 1
                phobic['phobic_total'] += 1
                    
            if z in philic:
                philic[z] += 1
                philic['philic_total'] += 1
                    
            if z in basic:
                basic[z] += 1
                basic['basic_total'] += 1
                    
            if z in acidic:
                acidic[z] += 1
                acidic['acidic_total'] += 1
                    
            if z in aromatic:
                aromatic[z] += 1
                aromatic['aromatic_total'] += 1
                    
            if z in sulfur:
                sulfur[z] += 1
                sulfur['sulfur_total'] += 1
                    
                    
        if i == 0:
                
            p1_list = [seq, seq_len, phobic['phobic_total'], philic['philic_total'], basic['basic_total'], acidic['acidic_total'], aromatic['aromatic_total'], sulfur['sulfur_total']]
                
        elif i != 0:
                
            p2_list = [seq, seq_len, phobic['phobic_total'], philic['philic_total'], basic['basic_total'], acidic['acidic_total'], aromatic['aromatic_total'], sulfur['sulfur_total']]
                
    p3_list = p1_list + p2_list
        
        
    return p3_list


merged = seq_features(seqs)

''
''
st.write(merged)