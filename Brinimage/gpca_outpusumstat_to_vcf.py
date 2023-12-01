
import glob
import pandas as pd
import numpy as np
import tempfile
import json,os
import argparse


pd.set_option('display.float_format', '{:.2E}'.format)
new_temp_dir = '/edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/'
tempfile.tempdir = new_temp_dir
print(tempfile.gettempdir())

filenames=glob.glob("/edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/*results.txt.gz")

for filename in filenames:
    out_preefix=f'{filename[:-7]}'
    fdf=pd.read_csv(filename,sep="\t")
    fdf=fdf[['CHR','BP','SNPID','EA','OA','BETA','SE','N_eff','PVAL','EAF']]
    fdf3=fdf[~fdf["CHR"].isna()]
    fdf3["CHR"]=fdf3["CHR"].astype("int").astype("str")
    fdf3[['BP','N_eff']]=fdf3[['BP','N_eff']].astype("int")
    fdf3[['EAF', 'BETA','SE','PVAL']]=fdf3[['EAF', 'BETA','SE','PVAL']].astype("float")
    fdf3.to_csv(f'{out_preefix}.tsv',sep="\t",index=None)
    
    ## CrossMap.py vcf GRCh38_to_GRCh37.chain.gz dbsnp_nochr.v153.hg38.vcf.gz Homo_sapiens.GRCh37.dna.primary_assembly.fa dbsnp_nochr.v153.hg37.vcf
    ## bgzip dbsnp_nochr.v153.hg37.vcf
    ## tabix -p vcf dbsnp_nochr.v153.hg37.vcf.gz
    
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
    
    with open(f'{out_preefix}_params.txt', 'w') as f:
        json.dump(paramsdict, f)
    
    ##Convert to vcf format 
    path=os.getcwd()
    ID="Cognition_Meta"
    
    os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
        --out {path}/{out_preefix}.vcf \
        --data {path}/{out_preefix}.tsv \
        --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
        --dbsnp /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/dbsnp_nochr.v153.hg37_sortd.vcf.gz \
        --json {out_preefix}_params.txt \
        --id {ID} > {out_preefix}.error 2>&1 ''' )
