
###  do it from unix 
srun -n 4 --mem=120G --time=24:01:00 --pty bash
module load Python/3.8.6-GCCcore-10.2.0
module load libxml2/2.9.8-GCCcore-6.4.0
module load OpenSSL/1.1
module load R/4.1.3-foss-2021b
module load GMP/6.2.1-GCCcore-11.2.0
module load BCFtools/1.14-GCC-11.2.0
source /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/env/bin/activate 

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

filenames=glob.glob("*txt.gz")
filenames=[x for x in filenames if "variants" not in x]

info=pd.read_csv("variants.txt.gz",sep="\s")
info["chr1"]=info["chr"].str.replace("^0","",regex=True)
info.drop("chr",inplace=True,axis=1)
info=info.rename(columns={"chr1":"chr"})
info.columns=['rsid','BP','EA', 'OA', 'EAF','info', 'CHR']
info["BP"]=info["BP"].astype("int")
info["CHR"]=info["CHR"].astype("str")
info=info[["CHR",'rsid',"BP","EA","OA","info",'EAF']]



for filename in filenames:
    print(f"{filename} started running")
    out_prefix=f'IDP_{filename[:-7]}'
    ID=out_prefix.replace(".txt","")
    fdf=pd.read_csv(filename,sep=" ")
    fdf["chr1"]=fdf["chr"].str.replace("^0X","X",regex=True)
    fdf['chr']=np.where(fdf['chr1'].isna(),fdf['chr'],fdf['chr1'])
    fdf['N_eff']=33000
    fdf=fdf.drop('chr1',axis=1)
    fdf=fdf.rename(columns={"chr":"CHR","pos":"BP","a1":"EA","a2":"OA","beta":"BETA","se":"SE",'pval(-log10)':"P"})
    fdf["CHR"]=fdf["CHR"].astype("str")
    fdf[['BP','N_eff']]=fdf[['BP','N_eff']].astype("int")
    fdf3=pd.merge(fdf,info,on=['CHR','BP','rsid','EA','OA'])
    fdf3=fdf3.rename(columns={"rsid":"SNPID","P":"PVAL"})
    fdf3[['EAF', 'BETA','SE','PVAL']]=fdf3[['EAF', 'BETA','SE','PVAL']].astype("float")
    fdf3=fdf3[['CHR','BP','SNPID','EA','OA','BETA','SE','N_eff','PVAL','EAF','info']]
    
    ## Conveert the -log10 p value to P value
    if fdf[fdf["P"]>1].shape[0]>1:
        fdf["P"]=np.power(10,-fdf["P"])
    
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
        "imp_info_col":10,
        "delimiter": "\t",
        "header": "true",
        "build": "GRCh37"
        }
    with open(f'{out_prefix}_params.txt', 'w') as f:
        json.dump(paramsdict, f)
    
    ##Convert to vcf format 
    path=os.getcwd()
    print(f"{filename} started vcf cration")
    os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
        --out {path}/{out_prefix}.vcf \
        --data {path}/{out_prefix}.tsv \
        --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
        --dbsnp /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/dbsnp_nochr.v153.hg37_sortd.vcf.gz \
        --json {out_prefix}_params.txt \
        --id {ID} > {out_prefix}.error 2>&1 ''' )
