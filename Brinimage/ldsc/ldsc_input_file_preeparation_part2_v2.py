import pandas as pd
import os,subprocess
import numpy as np

vcf_location="/edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/"
input_file=f"dbscan_clust_1_15_GenomicPCA_Covar_Corr_MDDwoBP_BIP_pgc3_SCZ_PGC3_ClutMeta_SNPOnly_NoMAV.vcf.gz"
output_prefix=f"dbscan_clust_1_15_GenomicPCA_Covar_Corr_MDDwoBP_BIP_pgc3_SCZ_PGC3_ClutMeta_SNPOnly_NoMAV"

os.system(f'''zgrep -v ".:.:.:." {vcf_location}{input_file}  | bgzip -c > {output_prefix}_Nomissing.vcf.gz''')
os.system(f'''tabix -f -p vcf {output_prefix}_Nomissing.vcf.gz''')


command = ['bcftools', 'query', '-l', f'{output_prefix}_Nomissing.vcf.gz']
sumstats_files = subprocess.run(command, capture_output=True, text=True)
sumstats_files=sumstats_files.stdout.split("\n")
sumstats_files=[x for x in sumstats_files if x!=""]

for sample in sumstats_files:
    os.system(f'''bcftools view -s {sample} {output_prefix}_Nomissing.vcf.gz | bgzip -c > {sample}.vcf.gz''')
    os.system(f'''tabix -p vcf {sample}.vcf.gz''')
    os.system(f'''gatk VariantsToTable  \
                -V {sample}.vcf.gz \
                -F CHROM -F POS -F ID -F REF -F ALT -GF ES -GF SE -GF LP -GF SS -GF NC -GF AF -GF SI \
                -O {sample}.tsv''') # -GF AF, -GF SI


binary_sumstat_types=['PGC3_SCZ','PGC3_SCZ_NoUKB','BIP_PGC3', 'BIP_noukbiobank',
                      'Depression_iPSYCH_2023','ADHD2022_iPSYCH_deCODE']


for sample in sumstats_files:
    df=pd.read_csv(f'{sample}.tsv',sep="\t")
    df.columns=[x.replace(f'{sample}.',"") for x in df.columns ]
    df["P"]=10**(-df["LP"])
    df['Z']=df['ES']/df['SE']
    
    ## Adding dummy values to for Allle frequency and INFo score
    df['AF']=np.where(df['AF'].isna(),0.1,df['AF'])
    df['SI']=np.where(df['SI'].isna(),0.9,df['SI'])
    
    if sample in binary_sumstat_types:
        v=df['NC']/(df['NC']+df['SS'])
        df['N']=4*v*(1-v)*(df['NC']+df['SS'])
        df['N']=df['N'].round().astype("int")
        df=df.rename(columns={
            'ID':"SNP",'CHROM':"CHR",'POS':"POS",'ALT':"A1",'REF':"A2",'AF':"eaf_A1",
            'ES':"effect",'SE':"SE",'SI':"INFO",'P':"P"})
        df=df[['CHR', 'POS', 'SNP', 'A2', 'A1', 'eaf_A1', 'effect', 'SE','Z','INFO','N','P']]
    
    elif sample not in binary_sumstat_types::
        df=df.rename(columns={
            'ID':"SNP",'CHROM':"CHR",'POS':"POS",'ALT':"A1",'REF':"A2",'AF':"eaf_A1",
            'ES':"effect",'SE':"SE",'SI':"INFO","SS":"N",'P':"P"})
        df=df[['CHR', 'POS', 'SNP', 'A2', 'A1', 'eaf_A1', 'effect', 'SE','Z','INFO','N','P']]
    
    df.to_csv(f"{sample}_munge_input.tsv",index=None,sep="\t")
