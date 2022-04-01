def mols_to_NHA(x): # input mmpa dataframe, grab smarts/smiles output number of heavy atoms
    return Chem.MolFromSmarts(x).GetNumHeavyAtoms()

def clean_mmpa_pairs_len(mmpa_df):
    temp=pd.DataFrame() # temp dataframe
    if 'LHS' not in mmpa_df.columns: # add LHS and RHS if not present
        mmpa_df = master_functions.split_transition(mmpa_df, 'smirks')     # produce LHS and RHS
    else:
        temp['common_core_HA'] = mmpa_df['common_core'].apply(mols_to_NHA) # produce number of heavy atoms
        temp['LHS_HA'] = mmpa_df['LHS'].apply(mols_to_NHA)
        temp['RHS_HA'] = mmpa_df['LHS'].apply(mols_to_NHA)
        
        temp['len_check'] = np.where((temp['LHS_HA'] >= temp['common_core_HA']) & (temp['RHS_HA'] >= temp['common_core_HA'])
                     , 'fail', 'pass') # compare lengths of heavy atoms
        
        mmpa_df = mmpa_df.drop(temp[temp['len_check']=='fail'].index) # drop index that failed length check
        # return temp to debug
    return mmpa_df # retun mmpa dataframe without len_check == fail
