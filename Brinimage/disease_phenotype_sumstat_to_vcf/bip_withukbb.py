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

##File location /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Disease_Phenotype_GWAS/  ; https://ipsych.dk/en/research/downloads
os.system("wget https://figshare.com/ndownloader/files/26603681")
os.system("mv 26603681 pgc-bip2021-all.vcf.tsv.gz")
os.system(f'''zgrep -v "##" pgc-bip2021-all.vcf.tsv.gz >pgc-bip2021-all.tsv''')


filename="pgc-bip2021-all.tsv" #location /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Disease_Phenotype_GWAS/

out_prefix=filename[:-4]
ID="BIP_PGC3"

fdf=pd.read_csv(f"{out_prefix}.tsv",sep="\t")
fdf3=fdf[['#CHROM', 'ID', 'POS', 'A1', 'A2','FCON', 'IMPINFO','BETA', 'SE', 'PVAL','NCAS', 'NCON']]
fdf3[['POS','NCAS', 'NCON']]=fdf3[['POS','NCAS', 'NCON']].astype("int")
fdf3[['FCON', 'IMPINFO', 'BETA', 'SE', 'PVAL']]=fdf3[['FCON', 'IMPINFO', 'BETA', 'SE', 'PVAL']].astype("float")

fdf3=fdf3[['#CHROM','POS','ID','A1','A2',"BETA",'SE','NCON','NCAS','PVAL','FCON','IMPINFO' ]]
fdf3.to_csv(f'{out_prefix}.tsv',sep="\t",index=None)

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
    "build": "GRCh37"
    }

with open(f'{out_prefix}.txt', 'w') as f:
  json.dump(paramsdict, f)

##Convert to vcf format 
path=os.getcwd()

os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
    --out {path}/{out_prefix}.vcf \
    --data {path}/{out_prefix}.tsv \
    --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    --json {out_prefix}.txt \
    --dbsnp /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/dbsnp_nochr.v153.hg37_sortd.vcf.gz \
    --id {ID} > {out_prefix}.error 2>&1 ''' )

