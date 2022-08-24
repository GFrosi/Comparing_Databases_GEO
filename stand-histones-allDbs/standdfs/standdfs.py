import pandas as pd
import numpy as np
from functools import reduce
import sys



def read_csv(path):
    """Receives a path to
    a csv file and returns
    a dataframe"""

    df = pd.read_csv(path)
    
    return df


def rename_col_drop_unnamed(df, col_old, col_new):
    '''Receives a df, old col name and new col name.
    Returns a df with no Unnamed column, and renamed
    col (old by new)'''
    
    df1 = df.copy()
    if 'Unnamed: 0' in df1.columns:
        df1.drop(columns=['Unnamed: 0'], inplace=True)
    
    else:
        print('No unnamed column found.')
        
    df1.rename(columns={col_old: col_new}, inplace = True)
    
    return df1


def check_GSM(df, col):
    '''Receives a df and
    column name (GSM).
    Returns a df with just
    samples with a GSM ID.
    '''

    df1 = df.copy()
    df1 = df1[df1[col].str.contains('GSM', na=False)]
    
    return df1


def list_dict_map(df, col_target):
    '''Receives a dataframe and a column
    target name. Returns a dictionary with
    keys as the typos of histones and inputs,
    and values with the standardized term.'''
    
    #maybe here pass a dict of regex as for GEO standardization!
    list_map = ['h3k4me3','h3k4me1','h3k9me3','h3k36me3','h3k27ac','h3k27me3','h3k9ac','h3k9/14ac','igg','wce','input','input control','control'] #input terms
    df_map = df[df[col_target].str.contains('|'.join(list_map), case=False, na=False)] #df to get the typos
    inp = list(set(df_map[col_target].tolist())) #list of hist to create the dict keys
    dict_map = dict.fromkeys(inp) #dict with each histone typo as keys
    dict_map = {k.lower().strip(): v for k,v in dict_map.items()} #standardazing to lower case
    
    for ele in list_map:
        for k,v in dict_map.items():
            if ele.lower() in k and not ele.startswith('h3'): #inputs
                dict_map[k] = 'input'

            if ele.lower() in k and ele.startswith('h3'): #histones
                dict_map[k] = ele

    return dict_map


def map_dict(df, col_hist):
    '''Receives a dataframe 
    and a column target name.
    Returns a dataframe with a
    new column standardized target
    column.'''

    dict_hist = list_dict_map(df, col_hist) #creating dict to map
    df_map = df.copy() #copying the original df
    new_col = col_hist+'_stand' #new target_stand mapped column
    df_map[new_col] = df_map[col_hist].str.lower().map(dict_hist).fillna(df_map[col_hist]) #creating new column
    # df_map.to_csv('test_map_ngs.csv')
    
    return df_map


def merge_df_outer(list_df):
    '''Receives a list of df
    (df_maps). Returns a merged
    df including all columns from
    metadata databases.'''

    df_merged = reduce(lambda left,right: pd.merge(left,right, on=['GSM'], how='outer').fillna('----'), list_df) #all Histones/inputs from all dbs
    df_merged.drop(columns=['SRX_y', 'GSE_y'], inplace=True)
    
    
    df_merged.rename(columns={'GSE_x':'GSE_GEO', 'Cell_type_x':'Cell_type_GEO', 'SRX_x':'SRX_GEO',
                          'Organism_x':'Organism_GEO', 
                         'Cell_type_y':'Cell_type_ChIP_Atlas', 
                         'target':'Target_CistromeDB',
                         'cell_type':'Cell_type_CistromeDB',
                         'Organism_y':'Organism_NGS-QC',
                         }, inplace=True)

    return df_merged


def filter_histones_inputs(df, col_target, col_gse): #generate them for all dbs, then merge outer
        '''Receives a metadata df (e.g GEO), target
        and GSE column name. Returns a df containig
        samples identified as Histones of interest and
        their respective inputs (input samples from
        the same series).'''

        #sorting values in GSE columns recover all inputs/histones for superseries (i.g NGS GSE88957, GSE88955 / GSE88955, GSE88957)
        df[col_gse] = df[col_gse].str.strip().sort_values().apply(lambda x: ",".join(sorted(str(x).split(","))))

        #Testing problematic example ngs-qc
        # a = df[df['GSM'].str.contains('GSM2356346')] #why?
        # c = df[df['GSM'].str.contains('GSM2356341')]
        
        # print('c')
        # frosi = c[col_gse].str.split(',')
        # print(frosi)

        # print('a')
        # problem = a[col_gse].str.split(',')
        # print(problem)


        #generating hist + input dfs
        list_hist = ['h3k4me3','h3k4me1','h3k9me3','h3k36me3','h3k27ac','h3k27me3','h3k9ac']
        df_hist = df[df[col_target].str.contains('|'.join(list_hist), case=False, na=False)] #maybe change to str.match
        list_gse_hist = list(set(df_hist[col_gse].tolist()))

        df_gse = df[df[col_gse].isin(list_gse_hist)]
        df_input = df_gse[df_gse[col_target].str.contains('input', case=False, na=False)] #change just for input. Map function including list of inputs

        # b = df_input[df_input[col_gse].str.contains('GSE88957')]
        # print(b)

        df_hist_input = pd.concat([df_hist, df_input]) #concatenating hist+inputs


        return df_hist_input


def histone_df(df1, df2): #adjust here to get all dbs
    '''Receives a df1 (df_histone_input)
    and df2 (merged_df). Returns a merged_df
    filtered by samples in df1.'''

    df_result = df1[['GSM']].merge(df2, how='left', on='GSM')
    df_result.drop_duplicates(subset=['GSM'], inplace=True)
 
    return df_result
