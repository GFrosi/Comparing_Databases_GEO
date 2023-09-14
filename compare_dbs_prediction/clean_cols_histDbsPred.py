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
            'Cancer', 'Disease_y', 'Healthy/None', 'Gsm', 'Address','Pmid', 'Processing-logs','Species','Adress', 'Srx'
             ]


    df_clean = df_total.drop(to_rm, axis=1)

    print('length dfafter clean columns:',len(df_clean.columns.values.tolist()))
    # print(df_clean.columns.values.tolist())
    
    cols_reorder = ['GSM', 'Release-date', 'Library-strategy', 'Organism-geo', 'Gpl', 'Gpl-title', 'Gse-geo', 'Gse-title', 'Gsm-title', 'Cell', 
                    'Disease_x', 'Sex-geo', 'Source', 'Chip-antibody-catalog', 'Target', 'Target-interest', 'Target-target', 'Target-catalog', 
                    'Target-geo', 'Cl-target', 'Srx-geo', 'Srr', 'Srr-count', 'Target-geo-stand', 'Target-target-stand', 'Target-catalog-stand','ENCODE_GSM', 'ENCODE_GSE',
                    'Experimental-id', 'Genome-assembly', 'Antigen-class', 'Antigen', 'Cell-type-class', 'Cell-type', 'Cell-type-description', 
                    'Title', 'Metadata', 'Age', 'Biomaterial-type', 'Cell-line-ca', 'Disease.1', 'Disease-ontology-uri', 'Donor-health-status',
                    'Donor-health-status-ontology-uri', 'Donor-id', 'Origin-sample-ontology-uri', 'Sample-ontology-uri', 'Sex', 'Sra', 'Tissue-type-ca','Antigen-stand',
                    'Gse', 'Cell-line-cistromedb', 'Target-cistromedb', 'Disease-state-name', 'Cell-type.1', 'Tissue-type-cistromedb', 'Target-cistromedb-stand',
                    'Study-id', 'Organism-ngs-qc', 'Data-type', 'Target-molecule', 'Seq-plataform', 'Submission-date', 'Target-molecule-stand',
                    'md5sum', 'True class_assay', 'Predicted class_assay', 'Same?_assay', 'Max pred_assay', 'True class_biomat', 'Predicted class_biomat', 'Same?_biomat',
                    'Max pred_biomat', 'True class_paired', 'Predicted class_paired', 'Same?_paired', 'Max pred_paired', 'True class_sex', 'Predicted class_sex', 
                    'Same?_sex', 'Max pred_sex', 'True class_disease', 'Predicted class_disease', 'Same?_disease', 'Max pred_disease',
                    'number_dbs_targets','number_dbs_agreement', 'dbs_agreement', 'Identical_4', 'Identical_3', 'Consensus_dbs_4', 'Consensus_dbs_3',
                    'dbs_target_match_pred', 'number_dbs_target_match_pred', 'number_dbs_target_no_match_pred',
                    'GEO', 'CA', 'CistromeDB', 'NGSQC', 'ENCODE', 'Merged-valid-sample', 'Merged-valid-sample_ENC',
                    'Geo-venn', 'CA-venn', 'Cistrome-venn', 'Ngs-venn', 'ENCODE-venn', 'Merged-venn', 'Merged-venn_ENC',
                    'groups', 'groups_ENCODE', 'groupsDbsEnc'
    ]

    df_clean_reorder = df_clean[cols_reorder]

    df_clean_reorder.to_csv('Hist_DBs_clean_CA_pred_reorder_cols_2023-09-07.tsv', sep="\t", index=False)


if __name__ == "__main__":


    main()
