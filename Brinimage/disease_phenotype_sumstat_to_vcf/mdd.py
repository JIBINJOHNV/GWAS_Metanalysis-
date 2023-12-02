###
srun -n 4 --mem=120G --time=24:01:00 --pty bash
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
print(tempfile.gettempdir())

##File location /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Disease_Phenotype_GWAS/  ; https://ipsych.dk/en/research/downloads
os.system("wget https://ipsych.dk/fileadmin/ipsych.dk/Downloads/daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01.gz")
os.system("wget https://ipsych.dk/fileadmin/ipsych.dk/Downloads/Readme_excl_UKB.pdf")
os.system(f'''zgrep -v "##" daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01.gz >daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01.tsv''')


filename="daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01.tsv" #location /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Disease_Phenotype_GWAS/

out_prefix=filename[:-4]
ID="Depression_iPSYCH_2023"
fdf=pd.read_csv(f"{out_prefix}.tsv",sep=" ")


fdf3=fdf[['CHR', 'SNP', 'BP', 'A1', 'A2','FRQ_U_507679', 'INFO', 'OR', 'SE', 'P', 'Nca', 'Nco']]
fdf3=fdf3[~fdf3["CHR"].isna()]
fdf3["CHR"]=fdf3["CHR"].astype("int").astype("str")

fdf3[['BP','Nca', 'Nco']]=fdf3[['BP','Nca', 'Nco']].astype("int")
fdf3[['FRQ_U_186843', 'INFO', 'OR', 'SE', 'P']]=fdf3[['FRQ_U_507679', 'INFO', 'OR', 'SE', 'P']].astype("float")
fdf3=fdf3[['CHR','BP','SNP','A1','A2',"OR",'SE','Nco','Nca','P','FRQ_U_507679','INFO' ]]
fdf3[fdf3["OR"]<0]

fdf3["OR"]=np.log(fdf3["OR"])
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