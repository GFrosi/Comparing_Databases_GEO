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
  
    date = datetime.now().strftime("%Y_%m_%d") 


    #list df to merge all dfs
    list_df = [df_map_geo, df_map_ngs, df_map_ca, df_map_cistrome] #samples with GSM
    df_merge_outer = sd.merge_df_outer(list_df)
    df_merge_outer_unique = df_merge_outer.drop_duplicates('GSM', keep='first')
    print('Your merged df has ' + str(len(df_merge_outer_unique)))
    # print(merge_df_left['Target molecule_stand']) - correct
    #saving df_merged 
    # print(df_merge_outer.nunique())
    # df_merge_outer_unique.to_csv('Compiled_dbs_2022_outer_'+ str(len(df_merge_outer_unique)) + '_' + date +'.csv', index=False) # here is ok, the stand data

    #Histones + Inputs df 
    # df_hist_inp_geo = sd.filter_histones_inputs(df_map_geo, 'Target-GEO_stand', 'GSE') #ok
    df_hist_inp_ngs = sd.filter_histones_inputs(df_map_ngs, 'Target molecule_stand', 'Study ID') #ok
    # df_hist_inp_ca = sd.filter_histones_inputs(df_map_ca, 'Antigen_stand', 'GSE') #ok testing all ca samples
    # df_hist_inp_cistrome = sd.filter_histones_inputs(df_map_cistrome, 'target_stand', 'GSE') #ok

    print(len(df_hist_inp_ngs))

    sys.exit()


    df_hist_inp_geo.to_csv('df_hist_inp_geo.csv', index=False)
    df_hist_inp_ngs.to_csv('df_hist_inp_ngs.csv', index=False)
    df_hist_inp_ca.to_csv('df_hist_inp_ca.csv', index=False)
    df_hist_inp_cistrome.to_csv('df_hist_inp_cistrome.csv', index=False)

    #concatenate all dfs based on histones per db
    df_to_concate_geo = sd.histone_df(df_hist_inp_geo, df_merge_outer_unique)
    df_to_concate_ngs = sd.histone_df(df_hist_inp_ngs, df_merge_outer_unique)
    df_to_concate_ca = sd.histone_df(df_hist_inp_ca, df_merge_outer_unique)
    df_to_concate_cistrome = sd.histone_df(df_hist_inp_cistrome, df_merge_outer_unique)

    df_histone_complete = pd.concat([df_to_concate_geo,df_to_concate_ngs,df_to_concate_ca,df_to_concate_cistrome]).drop_duplicates('GSM', keep='first').reset_index(drop=True)
    print('Your concatenated histone file has ' + str(len(df_histone_complete)))
    df_histone_complete.to_csv('Histones_basedDBs_'+ str(len(df_histone_complete)) + '_' + date +'.tsv', index=False, sep='\t')

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