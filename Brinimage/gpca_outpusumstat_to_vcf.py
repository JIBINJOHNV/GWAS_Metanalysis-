
###  do it from unix 
srun -n 4 --mem=120G --time=24:01:00 --pty bash
module load Python/3.8.6-GCCcore-10.2.0
module load libxml2/2.9.8-GCCcore-6.4.0
module load OpenSSL/1.1
module load R/4.1.3-foss-2021b
module load GMP/6.2.1-GCCcore-11.2.0
module load BCFtools/1.14-GCC-11.2.0
source /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/env/bin/activate 
###


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

#filenames=glob.glob("*Correlation.N_weighted_GWAMA.results.txt.gz")
filenames=glob.glob("*covariates.N_weighted_GWAMA.results.txt.gz")

info=pd.read_csv("variants.txt.gz",sep="\s")
info["chr1"]=info["chr"].str.replace("^0","",regex=True)
info.drop("chr",inplace=True,axis=1)
info=info.rename(columns={"chr1":"chr"})
info['SNPID']=info["chr"]+"_"+info["pos"].astype("str")+"_"+info["a1"]+"_"+info["a2"]
info.columns=['rsid','BP','EA', 'OA', 'EAF','info', 'CHR','SNPID']
info["BP"]=info["BP"].astype("int")
info["CHR"]=info["CHR"].astype("str")
info=info[["CHR","BP","EA","OA","SNPID","info"]]



for filename in filenames:
    print(f"{filename} started running")
    out_preefix=f'{filename[:-7]}'
    ID=out_preefix.replace(".N_weighted_GWAMA.results","")
    ID=ID.replace("dbscan_clust","Cluster").replace("GenomicPCA_Correlation","GPCA_CORR")
    ID=ID.replace("dbscan_clust","Cluster").replace("GenomicPCA_covariates","GPCA_COV")
    fdf=pd.read_csv(filename,sep="\t")
    fdf=fdf[['CHR','BP','SNPID','EA','OA','BETA','SE','N_eff','PVAL','EAF']]
    fdf3=fdf[~fdf["CHR"].isna()]
    fdf3["CHR"]=fdf3["CHR"].astype("int").astype("str")
    fdf3[['BP','N_eff']]=fdf3[['BP','N_eff']].astype("int")
    fdf3[['EAF', 'BETA','SE','PVAL']]=fdf3[['EAF', 'BETA','SE','PVAL']].astype("float")
    fdf3["CHR"]=fdf3["CHR"].str.replace("23","X")
    fdf3['SNPID']=fdf3["CHR"]+"_"+fdf3["BP"].astype("str")+"_"+fdf3["EA"]+"_"+fdf3["OA"]
    fdf3=pd.merge(fdf3,info,on=['CHR', 'BP', 'EA', 'OA', 'SNPID'])
    fdf3.to_csv(f'{out_preefix}.tsv',sep="\t",index=None)
    
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
    
    with open(f'{out_preefix}_params.txt', 'w') as f:
        json.dump(paramsdict, f)
    
    ##Convert to vcf format 
    path=os.getcwd()
    print(f"{filename} started vcf cration")
    os.system(f'''python /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/main.py \
        --out {path}/{out_preefix}.vcf \
        --data {path}/{out_preefix}.tsv \
        --ref /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
        --dbsnp /edgehpc/dept/human_genetics/users/jjohn1/gwas_vcf_reffiles/dbsnp_nochr.v153.hg37_sortd.vcf.gz \
        --json {out_preefix}_params.txt \
        --id {ID} > {out_preefix}.error 2>&1 ''' )
    

