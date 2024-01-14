import pandas as pd
import numpy as np
import pandas as pd
from scipy.stats import norm

con_sp=pd.read_csv("PLEIO_Output_with_Concordant_Specific_witheffectsize.tsv",sep="\t") 

const_cols=['CHROM', 'POS','SNP', 'REF', 'ALT', 'PLEIO_META_PLEIO_P', 'FRQ2']
file_names=['SCZ_PGC3_Primary_22_EUR', 'CTP_LAM_21_EUR', 'EDU_Ex23ME_LEE_18_EUR']

for file_name in file_names:
    df=con_sp[const_cols+[file_name]]
    df['SE']=abs(df[file_name]/ norm.ppf(df['PLEIO_META_PLEIO_P'] / 2))
    df['N']=100000
    df=df[['CHROM', 'POS', 'SNP','ALT', 'REF',file_name,'SE','N','PLEIO_META_PLEIO_P', 'FRQ2']]
    df.columns=['CHROM', 'POS', 'SNP','ALT', 'REF','Beta','SE','N','P', 'FRQ']
    df.to_csv(f'PLEIO_Output_with_Concordant_Specific_{file_name}.tsv',sep="\t",index=None)



##Disconcordant
con_sp=pd.read_csv("PLEIO_Output_with_Disconcordant_witheffectsize.tsv",sep="\t") 

const_cols=['CHROM', 'POS','SNP', 'REF', 'ALT', 'PLEIO_META_PLEIO_P', 'FRQ2']
file_names=['SCZ_PGC3_Primary_22_EUR', 'CTP_LAM_21_EUR', 'EDU_Ex23ME_LEE_18_EUR']

for file_name in file_names:
    df=con_sp[const_cols+[file_name]]
    df['SE']=abs(df[file_name]/ norm.ppf(df['PLEIO_META_PLEIO_P'] / 2))
    df['N']=100000
    df=df[['CHROM', 'POS', 'SNP','ALT', 'REF',file_name,'SE','N','PLEIO_META_PLEIO_P', 'FRQ2']]
    df.columns=['CHROM', 'POS', 'SNP','ALT', 'REF','Beta','SE','N','P', 'FRQ']
    df.to_csv(f'PLEIO_Output_with_Disconcordant_Specific_{file_name}.tsv',sep="\t",index=None)


##Both
con_sp=pd.read_csv("PLEIO_Output_with_Both_ConcordantDiscordant_Specific_withPLeioBeta_frq.tsv",sep="\t") 
const_cols=['CHROM', 'POS','SNP', 'REF', 'ALT', 'PLEIO_META_PLEIO_P', 'FRQ']
file_names=['SCZ_PGC3_Primary_22_EUR', 'CTP_LAM_21_EUR', 'EDU_Ex23ME_LEE_18_EUR']

for file_name in file_names:
    df=con_sp[const_cols+[file_name]]
    df['SE']=abs(df[file_name]/ norm.ppf(df['PLEIO_META_PLEIO_P'] / 2))
    df['N']=100000
    df=df[['CHROM', 'POS', 'SNP','ALT', 'REF',file_name,'SE','N','PLEIO_META_PLEIO_P', 'FRQ']]
    df.columns=['CHROM', 'POS', 'SNP','ALT', 'REF','Beta','SE','N','P', 'FRQ']
    df.to_csv(f'PLEIO_Output_with_Both_ConcordantDiscordant_Specific_{file_name}.tsv',sep="\t",index=None)

