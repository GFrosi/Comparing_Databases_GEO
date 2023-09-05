import pandas as pd
import sys





def main():

    df_total = pd.read_csv(sys.argv[1], sep="\t")

    list_col = df_total.columns.values.tolist()
    # print(len(list_col))
    # print(list_col)

    to_rm = ['Tmrs', 'Qc-stamps', 'Qc-reports', 'Local-qc-indicators', 
             '2.5-pcdenqc', '2.5-pcqcsim', '5-pcdenqc', '5-pcqcsim', '10-pcdenqc', '10-pcqcsim',
             'Is-correcting', 'Paper-reference', 'Paper-pmid', 'Paper-journal-name', 'Name',
             'Paper-lab', 'Sign', 'Judge-map', 'Judge-peaks', 'Judge-fastqc', 'Judge-frip', 'Judge-pbc',
            'Judge-motif', 'Judge-dhs', 'Map-perc', 'Peaks', 'Control-number', 'Fastqc', 'Frip', 'Sample',
            'Meta', 'Map-number', 'Pbc-qc', 'Dhs', 'Raw',
            'h3k27ac', 'h3k27me3', 'h3k36me3', 'h3k4me1', 'h3k4me3', 'h3k9me3', 'input', 'mrna_seq', 
            'rna_seq', 'wgbs-pbat', 'wgbs-standard', 'biomaterial_type_y',
            'cell line', 'primary cell', 'primary cell culture', 'primary tissue', 'Antigen_class_paired',
            'FALSE', 'TRUE', 'sex_sex', 'female', 'male', 'disease_disease', 
            'Cancer', 'Disease_y', 'Healthy/None'
             ]



    df_clean = df_total.drop(to_rm, axis=1)

    print('length df:',len(df_clean.columns.values.tolist()))
    print(df_clean.columns.values.tolist())
    

    print('df_clean saved with length of cols: ', len(df_clean.columns.values.tolist()))
    df_clean.to_csv('Hist_DBs_clean_CA_pred_gdrive.tsv', sep="\t")


if __name__ == "__main__":


    main()
