import pandas as pd
import os,subprocess

gpcatype="CORR" # COV
gpcatyp2="covariates" # Correlation
vcf_location="/edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/"
input_file=f"dbscan_clust_1_15_GenomicPCA_{gpcatyp2}_Disease_Phenotype_SNPOnly.vcf.gz"
output_prefix=f"dbscan_clust_1_15_GenomicPCA_{gpcatyp2}_Disease_Phenotype_SNPOnly"



os.system(f'''zgrep -v ".:.:.:." {vcf_location}{input_file}  | bgzip -c > {output_prefix}_Nomissing.vcf.gz''')
os.system(f'''tabix -f -p vcf {output_prefix}_Nomissing.vcf.gz''')


command = ['bcftools', 'query', '-l', f'{output_prefix}_Nomissing.vcf.gz']
# Run the command and capture its output
sumstats_files = subprocess.run(command, capture_output=True, text=True)
sumstats_files=sumstats_files.stdout.split("\n")
sumstats_files=[x for x in sumstats_files if x!=""]

for sample in sumstats_files:
    os.system(f'''bcftools view -s {sample} {output_prefix}_Nomissing.vcf.gz | bgzip -c > {sample}.vcf.gz''')
    os.system(f'''tabix -p vcf {sample}.vcf.gz''')
    os.system(f'''gatk VariantsToTable  \
                -V {sample}.vcf.gz \
                -F CHROM -F POS -F ID -F REF -F ALT -GF AF -GF ES -GF SE -GF LP -GF SS -GF SI -GF NC \
                -O {sample}.tsv''')


sumstat_types={f'Cluster_1_GPCA_{gpcatype}':"quant", f'Cluster_2_GPCA_{gpcatype}':"quant", 
 f'Cluster_3_GPCA_{gpcatype}':"quant", f'Cluster_4_GPCA_{gpcatype}':"quant", 
 f'Cluster_5_GPCA_{gpcatype}':"quant", f'Cluster_6_GPCA_{gpcatype}':"quant", 
 f'Cluster_7_GPCA_{gpcatype}':"quant", f'Cluster_8_GPCA_{gpcatype}':"quant", 
 f'Cluster_9_GPCA_{gpcatype}':"quant", f'Cluster_10_GPCA_{gpcatype}':"quant", 
 f'Cluster_11_GPCA_{gpcatype}':"quant", f'Cluster_12_GPCA_{gpcatype}':"quant", 
 f'Cluster_13_GPCA_{gpcatype}':"quant", f'Cluster_14_GPCA_{gpcatype}':"quant", 
 f'Cluster_15_GPCA_{gpcatype}':"quant", 
 'Lam_2021_CTP':"quant",'EA4_excl_23andMe_Okbay2022':"quant",
 'PGC3_SCZ':"binary", 'PGC3_SCZ_NoUKB':"binary", 
 'BIP_PGC3':"binary", 'BIP_noukbiobank':"binary", 
 'Depression_iPSYCH_2023':"binary", 
 'ADHD2022_iPSYCH_deCODE':"binary"}


for sample in sumstat_types.keys():
    df=pd.read_csv(f'{sample}.tsv',sep="\t")
    df.columns=[x.replace(f'{sample}.',"") for x in df.columns ]
    df["P"]=10**(-df["LP"])
    df['Z']=df['ES']/df['SE']
    
    if sumstat_types[sample]=="binary":
        v=df['NC']/(df['NC']+df['SS'])
        df['N']=4*v*(1-v)*(df['NC']+df['SS'])
        df['N']=df['N'].round().astype("int")
        df=df.rename(columns={
            'ID':"SNP",'CHROM':"CHR",'POS':"POS",'ALT':"A1",'REF':"A2",'AF':"eaf_A1",
            'ES':"effect",'SE':"SE",'SI':"INFO",'P':"P"})
        df=df[['CHR', 'POS', 'SNP', 'A2', 'A1', 'eaf_A1', 'effect', 'SE','Z','INFO','N','P']]
    
    elif sumstat_types[sample]=="quant":
        df=df.rename(columns={
            'ID':"SNP",'CHROM':"CHR",'POS':"POS",'ALT':"A1",'REF':"A2",'AF':"eaf_A1",
            'ES':"effect",'SE':"SE",'SI':"INFO","SS":"N",'P':"P"})
        df=df[['CHR', 'POS', 'SNP', 'A2', 'A1', 'eaf_A1', 'effect', 'SE','Z','INFO','N','P']]
    
    df.to_csv(f"{sample}_munge_input.tsv",index=None,sep="\t")
