# Comparing_Databases_GEO
Compiling and comparing databases (GEO_based), and checking the metadata validation


# Requiremets

```python > 3.0``` and ```create a env with pandas```

# stand-histones

## main_stand_dbs.py

### Usage
```
A script to merge several dataframes from different databases, standardize and
filter the Histones and Input samples

optional arguments:
  -h, --help            show this help message and exit
  -g GEO, --geo GEO     GEO metadata csv file generated by GEO-Metadata script
  -n NGS, --ngs NGS     NGS-QC metadata csv file generated by NGS-QC-
                        extraction script
  -c CA, --ca CA        ChIP-Atlas metadata csv file generated by ChIP-Atla-
                        extraction script
  -C CISTROME, --cistrome CISTROME
                        Cistrome metadata csv file generated by Cistrome-
                        extraction script
```

## merge_prediction.py

### Usage
```
Script to merge the table generated by main_stand_dbs.py and the EpiLaP prediction table.

python merge_prediction.py Histones_basedGEO_alldbs_dbs.tsv ChIP_Atlas_pred_EpiLaP_22k.csv
```


# compare_dbs_prediction/

## main.py

```
Script to generate a table including columns containing the information of how many databases agree/disagree compared to EpiLaP prediction, and how many in each situation above agree/disagree among them.  

python main.py Histones_basedGEO_alldbs_dbs.tsv Histones_allDBs_CA_pred_comparisonDBs.tsv
```