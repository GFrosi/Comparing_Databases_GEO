import pandas as pd
import numpy as np
from functools import reduce
import sys
import re



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
    list_map = ['h3k4me3','h3k4me1','h3k9me3','h3k36me3','h3k27ac','h3k27me3','igg','wce','input','input control','control'] #input terms
    # list_map = ['h3k4me3','h3k4me1','h3k9me3','h3k36me3','h3k27ac','h3k27me3','h3k9ac','h3k9/14ac','igg','wce','input','input control','control'] #input terms
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

    return df_map


def map_dict_catalog(df, col_hist_cat):

    #Adjusted running (e.g ab8580 (Lot GR - GSM2576925)
    #consider lower case when map
    #consider if k in column to replace it (check if it is possible using map)

    dict_antb = {
    '12-371':'input','ab46540':'input','G4018':'input','i2511':'input','sc-2027':'input',
    '07-360':'h3k27ac','ab4729':'h3k27ac','8173':'h3k27ac','39135':'h3k27ac','39133':'h3k27ac',
    'C15410196':'h3k27ac','C15410184':'h3k27ac',
    'ab8895':'h3k4me1','8895':'h3k4me1','39297':'h3k4me1','39635':'h3k4me1',
    '07-436':'h3k4me1','C15410194':'h3k4me1','cs-037-100':'h3k4me1','pAb-037-050':'h3k4me1',
    '07-473':'h3k4me3','ab8580':'h3k4me3','39159':'h3k4me3','39-159':'h3k4me3','305-34819':'h3k4me3',
    '9751S':'h3k4me3','9751':'h3k4me3','17-614':'h3k4me3','307-34813':'h3k4me3','C15200152':'h3k4me3','C15410003':'h3k4me3',
    'ab1012':'h3k4me3','pAb-003-050':'h3k4me3','pAb-MEHAHS-024':'h3k4me3','05-745R':'h3k4me3','04-745':'h3k4me3',
    '07-449':'h3k27me3','ab6002':'h3k27me3','CS200603':'h3k27me3','CS20060':'h3k27me3','C36B11':'h3k27me3','C15410069':'h3k27me3',
    'ABE44':'h3k27me3','39155':'h3k27me3','192985':'h3k27me3','17-622':'h3k27me3','9733':'h3k27me3','39536':'h3k27me3','61017':'h3k27me3',
    'C15410195':'h3k27me3','61017':'h3k27me3','301-95253':'h3k27me3',
    '300-95289':'h3k36me3','CS-058-100':'h3k36me3','61101':'h3k36me3','61021':'h3k36me3','302-95283':'h3k36me3','ab9050':'h3k36me3','C15410192':'h3k36me3',
    'ab8898':'h3k9me3','39765':'h3k9me3','07-442':'h3k9me3','17-10242':'h3k9me3','17-625':'h3k9me3','9754S':'h3k9me3','ab1186':'h3k9me3',
    '39161':'h3k9me3','C15410193':'h3k9me3','pAb-056-050':'h3k9me3'
    }

    #keys to lower case to avoid mismatches while looking at col values
    dict_antb_lower = {k.lower():v for k,v in dict_antb.items()} 

    #Creating regex patterns for keys (several times keys are subsrings of col values, i.e ab8580 (Lot GR - GSM2576925 )

    pat = '|'.join(r"\b{}\b".format(x) for x in dict_antb_lower.keys())
    df[col_hist_cat] = df[col_hist_cat].str.lower().str.extract('('+ pat + ')', expand=False).map(dict_antb_lower)


    return df



def merge_df_outer(list_df):
    '''Receives a list of df
    (df_maps). Returns a merged
    df including all columns from
    metadata databases.'''

    df_merged = reduce(lambda left,right: pd.merge(left,right, on=['GSM'], how='outer').fillna('----'), list_df) #all Histones/inputs from all dbs
    # cols = df_merged.columns.values.tolist()
    # print('df_merged cols:', cols)
    # sys.exit()

    #organizing columns    - CHECK COLUMNS TOMORROW!
    df_merged.drop(columns=['SRX_y', 'GSE_y','Unnamed: 0'], inplace=True)
    
    df_merged.rename(columns={'Organism_x':'Organism_GEO',
                        'GSE_x':'GSE-GEO',
                        'SRX_x': 'SRX_GEO',
                        'Organism_y':'Organism-NGS-QC',
                        'cell_line_x':'cell-line-CA',
                        'tissue_type_x':'tissue-type-CA',
                        'cell_line_y':'cell-line-CistromeDB',
                        'tissue_type_y':'tissue-type-CistromeDB',
                        'target':'Target-CistromeDB',
                        'target_stand':'Target-CistromeDB_stand'}, inplace=True)


    #standardizing col names
    # df_merged.drop(columns=['Unnamed: 0'], inplace=True)
    df_merged.columns = df_merged.columns.str.replace(' ', '-')
    df_merged.columns = df_merged.columns.str.replace('_', '-')
    # df_merged.columns = df_merged.columns.str.replace('_', '-')
    df_merged.columns = df_merged.columns.str.capitalize() #fisrt letter


    cols = df_merged.columns.values.tolist()
    print('df_merged cols:', cols)
    # sys.exit()

    return df_merged


def filter_histones_inputs(df, col_target):
    '''Receives a metadata df (e.g GEO) and
    target column name. Returns a df containig
    samples identified as Histones of interest
    and inputs'''

    list_targets = ['H3K4me3','H3K4me1','H3K9me3','H3K36me3','H3K27ac','H3K27me3','Input'] #stand columns
    
    return df[df[col_target].str.contains('|'.join(list_targets), case=False, na=False)]    


# def filter_histones_inputs(df, col_target, col_gse): #generate them for all dbs, then merge outer
#         '''Receives a metadata df (e.g GEO), target
#         and GSE column name. Returns a df containig
#         samples identified as Histones of interest and
#         their respective inputs (input samples from
#         the same series).'''

#         #sorting values in GSE columns recover all inputs/histones for superseries (i.g NGS GSE88957, GSE88955 / GSE88955, GSE88957)
#         df[col_gse] = df[col_gse].str.strip().sort_values().apply(lambda x: ",".join(sorted(str(x).split(","))))

#         #generating hist + input dfs
#         list_hist = ['h3k4me3','h3k4me1','h3k9me3','h3k36me3','h3k27ac','h3k27me3']
#         # list_hist = ['h3k4me3','h3k4me1','h3k9me3','h3k36me3','h3k27ac','h3k27me3','h3k9ac']
#         df_hist = df[df[col_target].str.contains('|'.join(list_hist), case=False, na=False)] #maybe change to str.match
#         list_gse_hist = list(set(df_hist[col_gse].tolist()))

#         df_gse = df[df[col_gse].isin(list_gse_hist)]
#         df_input = df_gse[df_gse[col_target].str.contains('input', case=False, na=False)] #change just for input. Map function including list of inputs

#         df_hist_input = pd.concat([df_hist, df_input]) #concatenating hist+inputs
        
#         return df_hist_input




def histone_df(df1, df2): #adjust here to get all dbs
    '''Receives a df1 (df_histone_input)
    and df2 (merged_df). Returns a merged_df
    filtered by samples in df1.'''

    df_result = df1[['GSM']].merge(df2, how='left', left_on='GSM', right_on='Gsm')

    # cols = df_result.columns.values.tolist()

    # for i in cols:
    #     print(i)

    df_result.drop_duplicates(subset=['GSM'], inplace=True)
 
    return df_result
