###
srun -n 4 --mem=120G --time=4:01:00 --pty bash
module load Python/3.8.6-GCCcore-10.2.0
module load libxml2/2.9.8-GCCcore-6.4.0
module load OpenSSL/1.1
module load R/4.1.3-foss-2021b
module load GMP/6.2.1-GCCcore-11.2.0
module load BCFtools/1.14-GCC-11.2.0
source /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/env/bin/activate
###

import pandas as pd
import numpy as np
import tempfile
import json,os
import argparse


pd.set_option('display.float_format', '{:.2E}'.format)
new_temp_dir = '/edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/'
tempfile.tempdir = new_temp_dir
print(tempfile.gettempdir())

##File location /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Disease_Phenotype_GWAS/  ; daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.trios ; Proff Todd shared
path=os.getcwd()
ID="PGC3_SCZ_NoUKB"

filename="daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.trios" #location /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Disease_Phenotype_GWAS/
out_prefix=filename
ID="PGC_SCZ_noukbiobank"


fdf=pd.read_csv(filename,sep="\t")
fdf3=fdf[['CHR', 'SNP', 'BP', 'A1', 'A2','FRQ_U_93456', 'INFO', 'OR', 'SE', 'P', 'Nca', 'Nco']]
fdf3[['BP','Nca', 'Nco']]=fdf3[['BP','Nca', 'Nco']].astype("int")
fdf3["BETA"]=np.log(fdf3["OR"])
fdf3[['FRQ_U_93456', 'INFO', 'BETA', 'SE', 'P']]=fdf3[['FRQ_U_93456', 'INFO', 'BETA', 'SE', 'P']].astype("float")
fdf3=fdf3[['CHR','BP','SNP','A1','A2',"BETA",'SE','Nco','Nca','P','FRQ_U_93456','INFO' ]]

paramsdict={"chr_col": 0,
    "pos_col": 1,
    "snp_col": 2,
    "ea_col": 3,
    "oa_col": 4,
    "beta_col": 5,
    "se_col": 6,
    "ncontrol_col": 7,
    "ncase_col":8,
    "pval_col": 9,
    "eaf_col": 10,
    "imp_info_col":11,
    "delimiter": "\t",
    "header": "true",
    "build": "GRCh37"}


with open(f'{out_prefix}.txt', 'w') as f:
  json.dump(paramsdict, f)

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{out_prefix}.vcf \
    --data {path}/{out_prefix}.tsv \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    --json {out_prefix}.txt \
    --dbsnp /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/dbsnp_nochr.v153.hg37_sortd.vcf.gz \
    --id {ID} > {out_prefix}.error 2>&1 ''' )
