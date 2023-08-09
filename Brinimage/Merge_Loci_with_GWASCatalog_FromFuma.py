import glob
import pandas as pd
from pybedtools import BedTool
import numpy as np
import argparse



parser = argparse.ArgumentParser(description="This script is for Merginng Fuma identified GWAS catlog from multiple Fuma output")
parser.add_argument('-fumafilesufix',
                    '--fumafilesufix',help="sufix of fiuma output files eg: GenomicRiskLoci.txt",required=True)

fumafile=args.fumafilesufix
fumafile="*"+fumafile
files=glob.glob(fumafile)



connccat_df=pd.DataFrame()
for file in files:
    df1=pd.read_csv(file,sep="\t")
    df1=df1[[ 'chr', 'bp','Trait','PMID']].drop_duplicates()
    df1['start']=df1["bp"]-1
    df1['end']=df1["bp"]+1
    df1=df1.drop("bp",axis=1)
    df1["GWAS_CAtalog_Info"]=df1["PMID"].astype("str")+":"+df1["Trait"]
    df1=df1.drop(["Trait","PMID"],axis=1)
    df1=df1.groupby(["chr","start","end"]).agg({x:",".join for x in df1.iloc[:,3:].columns}).reset_index()
    connccat_df=pd.concat([connccat_df,df1])

connccat_df=connccat_df.drop_duplicates()
connccat_df=connccat_df.groupby(["chr","start","end"]).agg({x:",".join for x in df1.iloc[:,3:].columns}).reset_index()
connccat_df_bed=BedTool.from_dataframe(connccat_df)

fuma_merged_loci="/home/jjohn41/DBSCAN/Genomic_risk_loci/Fuma_Genomic_Loci_Comparison.csv"
fuma_merged_loci_df=pd.read_csv(fuma_merged_loci)
fuma_merged_loci_df=fuma_merged_loci_df.rename(columns={"chrom":"chr"})
master_loci=fuma_merged_loci_df[["chr","start","end"]].drop_duplicates()

master_locibed = BedTool.from_dataframe(master_loci)


master_loci_concat_merged=master_locibed.intersect(connccat_df_bed, wa=True, wb=True,loj=True)

master_loci_concat_merged_df=master_loci_concat_merged.to_dataframe()
master_loci_concat_merged_df=master_loci_concat_merged_df.drop(["name","score","strand"],axis=1).drop_duplicates()

master_loci_concat_merged_df=master_loci_concat_merged_df.groupby(["chrom","start","end"]).agg({x:",".join for x in master_loci_concat_merged_df.iloc[:,3:].columns}).reset_index()
master_loci_concat_merged_df.columns
master_loci_concat_merged_df.rename(columns={"thickStart":"GWAS_Catalog","chrom":"chr"},inplace=True)

Final_results=pd.merge(fuma_merged_loci_df,master_loci_concat_merged_df,on=["chr","start","end"])
Final_results["MosTest"]=np.where(Final_results["GWAS_Catalog"].str.contains("MOSTest"),"MOSTest","No_MOSTest")

Final_results.to_csv("Fuma_Genomic_Loci_Comparison_With_GWAS_catlogInfo.csv",index=None)
