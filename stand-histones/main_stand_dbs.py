import sys
import argparse
import pandas as pd
from functools import reduce
from standdfs import standdfs as sd 
from datetime import datetime



def main():

    print('Starting...')

    df_geo = sd.read_csv(args.geo)
    df_ngs = sd.read_csv(args.ngs)
    df_ca = sd.read_csv(args.ca)
    df_cistrome = sd.read_csv(args.cistrome)
    df_ngs = sd.rename_col_drop_unnamed(df_ngs, 'Public ID', 'GSM') 

    #filtering samples with GSM ID
    df_geo_gsm = sd.check_GSM(df_geo, 'GSM')
    df_ngs_gsm = sd.check_GSM(df_ngs, 'GSM')
    df_ca_gsm = sd.check_GSM(df_ca, 'GSM')
    df_cistrome_gsm = sd.check_GSM(df_cistrome, 'GSM')

    #creating stand target column - histones
    df_map_geo = sd.map_dict(df_geo_gsm, 'Target-GEO')
    df_map_ngs = sd.map_dict(df_ngs_gsm, 'Target molecule')
    df_map_ca = sd.map_dict(df_ca_gsm, 'Antigen')
    df_map_cistrome = sd.map_dict(df_cistrome_gsm, 'target')
    # sys.exit()   

    #list df to merge all dfs
    list_df = [df_map_geo, df_map_ngs, df_map_ca, df_map_cistrome]
    merge_df_left = sd.merge_df_left(list_df)
    print('Your merged df has ' + str(len(merge_df_left)))
    # print(merge_df_left['Target molecule_stand']) - correct
    #saving df_merged 
    # merge_df_left.to_csv('Compiled_dbs_2022_'+ str(len(merge_df_left)) + '.csv', index=False) # here is ok, the stand data

    #Histones + Inputs df (the problem is here...)
    df_hist_inp_geo = sd.filter_histones_inputs(df_map_geo, 'Target-GEO_stand', 'GSE') #ok
    
    df_result_histones_dbs = sd.histone_df(df_hist_inp_geo, merge_df_left)
    print('Your Histone_merged df has ' + str(len(df_result_histones_dbs)))
    #saving df_merged_histones based on GEO classification
    date = datetime.now().strftime("%Y_%m_%d") 
    df_result_histones_dbs.to_csv('Histones_basedGEO_alldbs_dbs_'+ str(len(df_result_histones_dbs)) + '_' + date +'.tsv', index=False, sep='\t')

    print('Files successfully saved!')





if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="A script to merge several dataframes from different databases, standardize and filter the Histones and Input samples"
    )

    parser.add_argument(
        '-g', '--geo', action="store",
        help="GEO metadata csv file generated by GEO-Metadata script",
        required=True
    )

    parser.add_argument(
        '-n', '--ngs',
         help="NGS-QC metadata csv file generated by NGS-QC-extraction script",
         required=True

    )
    parser.add_argument(
        '-c', '--ca',
         help="ChIP-Atlas metadata csv file generated by ChIP-Atla-extraction script",
         required=True

    )


    parser.add_argument(
        '-C', '--cistrome',
         help="Cistrome metadata csv file generated by Cistrome-extraction script",
         required=True
    )

    args = parser.parse_args()
    main()