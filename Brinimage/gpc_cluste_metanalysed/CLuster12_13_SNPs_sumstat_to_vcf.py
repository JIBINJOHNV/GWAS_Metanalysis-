import glob
import pandas as pd
import numpy as np
import tempfile
import json,os
import argparse

parser=argparse.ArgumentParser(description="for converting Brain image cluter12 & 13 metaanalysed GenomicPC sumstat to vcf convertion")
parser.add_argument('-inputsumstat','--inputsumstat', help="GenomicPCA output sumstat", required=True)
args=parser.parse_args()


filename=args.inputsumstat
sample_size=12181081

info=pd.read_csv("variants.txt.gz",sep="\s")
info["chr1"]=info["chr"].str.replace("^0","",regex=True)
info.drop("chr",inplace=True,axis=1)
info=info.rename(columns={"chr1":"chr"})
info.columns=['rsid','BP','EA','OA','EAF','info', 'CHR']
info["BP"]=info["BP"].astype("int")
info["CHR"]=info["CHR"].astype("str")
info=info[["CHR","BP","EA","OA","info",'EAF']]
info['SNPID']=info["CHR"]+"_"+info["BP"].astype("str")+"_"+info["EA"]+"_"+info["OA"]


print(f"{filename} started running")
out_prefix=filename.replace("_METAANALYSIS_1.tbl","")
ID=out_prefix.replace("IDP_","")
out_prefix=out_prefix.replace("IDP_","")

df=pd.read_csv(filename,sep="\t")
df=df.rename(columns={"MarkerName":"SNPID","Allele1":"EA","Allele2":"OA","Effect":"BETA",
                    "StdErr":"SE","P-value":'PVAL'})[['SNPID', 'EA', 'OA', 'BETA', 'SE', 'PVAL']]
df["EA"]=df["EA"].str.upper()
df["OA"]=df["OA"].str.upper()
df2=pd.merge(df,info,on=['SNPID'])
df2["BETA2"]=np.where(df2["EA_x"]!=df2["EA_y"],-df2["BETA"],df2["BETA"])
df2=df2[~df2["EA_x"].isna()]
df2=df2.drop(["EA_x","OA_x","BETA"],axis=1).rename(columns={"EA_y":"EA","OA_y":"OA","BETA2":"BETA"})
df2["N_eff"]=sample_size 
df2=df2[['CHR', 'BP','SNPID','EA', 'OA','BETA','SE',"N_eff",'PVAL','EAF','info']]
df2[['BP','N_eff']]=df2[['BP','N_eff']].astype("int")

df2['EAF']=df2['EAF'].astype("float")
df2['BETA']=df2['BETA'].astype("float")
df2['SE']=df2['SE'].astype("float")
df2['PVAL']=df2['PVAL'].astype("float")
df2["CHR"]=df2["CHR"].str.replace("23","X")
df2=df2.sort_values(by=["CHR","BP"])
df2.to_csv(f'{out_prefix}.tsv',sep="\t",index=None)

paramsdict={"chr_col": 0,
        "pos_col": 1,
        "snp_col": 2,
        "oa_col": 3,
        "ea_col": 4,
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
