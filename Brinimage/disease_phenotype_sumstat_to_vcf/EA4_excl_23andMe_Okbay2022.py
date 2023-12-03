
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

##File location /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Disease_Phenotype_GWAS/  ; https://ipsych.dk/en/research/downloads
## Downloaded from https://thessgac.com/papers/
os.system(f'''zgrep -v "##" EA4_additive_excl_23andMe.txt.gz >EA4_additive_excl_23andMe.tsv''')


filename="EA4_additive_excl_23andMe.tsv" #location /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Disease_Phenotype_GWAS/

out_prefix=filename[:-4]
ID="EA4_excl_23andMe_Okbay2022"


fdf=pd.read_csv(f"{out_prefix}.tsv",sep="\t")
fdf["N"]=765283
fdf2=fdf[['rsID', 'Chr', 'BP', 'Effect_allele', 'Other_allele', 'EAF_HRC', 'Beta','SE', 'P', 'N']]
fdf3=fdf2[~fdf2["rsID"].isna()]
fdf3[['BP','N']]=fdf3[['BP','N']].astype("int")
fdf3[['EAF_HRC', 'Beta', 'SE', 'P']]=fdf3[['EAF_HRC', 'Beta', 'SE', 'P']].astype("float")
fdf3=fdf3[['Chr','BP','rsID','Effect_allele','Other_allele',"Beta",'SE','N','P','EAF_HRC' ]]
fdf3.to_csv(f'{out_prefix}.tsv',sep="\t",index=None)

paramsdict={"chr_col": 0,
    "pos_col": 1,
    "snp_col": 2,
    "ea_col": 3,
    "oa_col": 4,
    "beta_col": 5,
    "se_col": 6,
    "ncontrol_col": 7,
    "pval_col": 8,
    "eaf_col": 9,
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
