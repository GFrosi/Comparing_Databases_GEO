#!/bin/bash
#SBATCH --time=24:00:00  
#SBATCH --account=
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=comparing-dbs-stand

module load python/3.6
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip

pip install --no-index -r requirements.txt

echo "Starting"

python Comparing_Databases_GEO/stand-histones-allDbs/main_stand_dbs.py -g GEO_metadata_2023_91930_stand.csv -n NGS-QC/NGS_HS_ChipSeq_nodup_33233.csv -c ChIP-Atlas/CA_hg38_Hs_GSM_GSE_2022_antigenclass_2023_01_24.csv -C Cistrome/Cistrome_filter_human_noENC.csv

echo "Done"


