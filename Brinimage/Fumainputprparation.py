import pandas as pd 

info=pd.read_csv("variants.txt.gz",sep="\s")
info["chr1"]=info["chr"].str.replace("^0","",regex=True)
info.drop("chr",inplace=True,axis=1)
info=info.rename(columns={"chr1":"chr"})
info['SNPID']=info["chr"]+"_"+info["pos"].astype("str")+"_"+info["a1"]+"_"+info["a2"]
info.columns=['rsid','BP','EA', 'OA', 'EAF','info', 'CHR','SNPID']
info["BP"]=info["BP"].astype("int")
info["CHR"]=info["CHR"].astype("str")


file="dbscan_clust_1_GenomicPCA_Correlation.N_weighted_GWAMA.results.txt.gz"

df=pd.read_csv(file,sep="\t")
df["BP"]=df["BP"].astype("int")
df["CHR"]=df["CHR"].astype("str")
df["CHR"]=df["CHR"].str.replace("23","X")


df['SNPID']=df["CHR"]+"_"+df["BP"].astype("str")+"_"+df["EA"]+"_"+df["OA"]


merged=pd.merge(df,info,on=['SNPID', 'CHR','BP','EA', 'OA', 'EAF'])

mergd2=merged[['CHR', 'BP','SNPID','EA', 'OA','EAF','info','N_eff','BETA', 'SE','PVAL']]

indels=mergd2[ (mergd2["EA"].str.len()>1) | (mergd2["OA"].str.len()>1) ]
snps=mergd2[ (mergd2["EA"].str.len()==1) & (mergd2["OA"].str.len()==1) ]

snps.to_csv(f"{file[:-7]}SNP.txt",sep="\t",index=None)
mergd2.to_csv(f"{file[:-7]}_Sumstat.txt",sep="\t",index=None)

