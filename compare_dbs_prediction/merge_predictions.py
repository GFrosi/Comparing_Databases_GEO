import pandas as pd 
import sys
import os
import glob
from functools import reduce



def merge_pred(list_dfs):

    df = reduce(lambda df1,df2: pd.merge(df1,df2,on='md5sum'), list_dfs)
    
    df.to_csv('C-A_merged_predictions_54k_total_2023.tsv', sep="\t",index=False)



def clean_dfs(path_csv):

    csv_files = glob.glob(os.path.join(path_csv, "*.csv"))
    list_dfs = [os.path.basename(file).split('.csv')[0] for file in csv_files]

    cols_del = ['2nd pred class','1rst/2nd prob diff','1rst/2nd prob ratio']

    clean_dfs = []

    for file,df in zip(csv_files,list_dfs):

        df = pd.read_csv(file)
        df.drop(cols_del, axis=1, inplace=True)
        clean_dfs.append(df)

    return clean_dfs

def main():

    path_csv = sys.argv[1]
    clean_dfs_list = clean_dfs(path_csv)
    merge_pred(clean_dfs_list)

    # df_sex = pd.read_csv(sys.argv[1])
    # df_assay11 = pd.read_csv(sys.argv[2])
    # df_assay13 = pd.read_csv(sys.argv[3])
    # df_bio = pd.read_csv(sys.argv[4])
    # df_life = pd.read_csv(sys.argv[5])
    # df_disease = pd.read_csv(sys.argv[6])
    # df_paired = pd.read_csv(sys.argv[7])


if __name__ == "__main__":

    main()
