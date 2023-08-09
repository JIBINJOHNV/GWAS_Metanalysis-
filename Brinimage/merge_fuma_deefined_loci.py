import glob
import pandas as pd
from pybedtools import BedTool
import numpy as np
import argparse


parser = argparse.ArgumentParser(description="This script is for preparing inputs for the meta analysis")
parser.add_argument('-fumafilesufix',
                    '--fumafilesufix',help="sufix of fiuma output files eg: GenomicRiskLoci.txt",required=True)

fumafile=args.fumafilesufix
fumafile="*"+fumafile
files=glob.glob(fumafile)


fumafile="GenomicRiskLoci.txt"

MainDf_list=list()
MainDf=pd.DataFrame()

for file in files:
    df1=pd.read_csv(file,sep="\t")
    df1a=df1[['chr','start', 'end']]
    MainDf_list.append(df1a)

MainDf=pd.concat(MainDf_list)
MainDf=MainDf.sort_values(by=["chr","start","end"]).drop_duplicates()
MainDfbed = BedTool.from_dataframe(MainDf)
MainDf_merged=MainDfbed.merge()

FinalResulktsdf=MainDf_merged.to_dataframe().copy()


for file in files:
    df1=pd.read_csv(file,sep="\t")
    colprefix=file.split("/")[1].replace("_GenomicRiskLoci.txt","")+":"
    df1a=df1[['chr','start', 'end','nIndSigSNPs','nLeadSNPs','rsID','p']]
    df1abed = BedTool.from_dataframe(df1a)
    df1bed_merged=MainDf_merged.intersect(df1abed, wa=True, wb=True,loj=True)
    resultdf=df1bed_merged.to_dataframe()
    resultdf.columns=["chrom","start","end"]+[colprefix+x for x in df1a.columns  ]
    FinalResulktsdf=pd.merge(FinalResulktsdf,resultdf,on=['chrom', 'start', 'end'],how="outer")
    FinalResulktsdf=FinalResulktsdf.drop_duplicates()


df2=FinalResulktsdf.copy()
df2=df2.astype("str")
df3=df2.groupby(["chrom","start","end"]).agg({x:",".join for x in df2.iloc[:,3:].columns}).reset_index()

for Col in list(df3.iloc[:,3:].columns):
    df3[Col]=df3[Col].apply(lambda x: ";".join(set(x.split(","))))

df3=df3.sort_values(by=["chrom","start","end"]).drop_duplicates()
df3.to_csv("Fuma_Genomic_Loci_Comparison.csv",index=None)


#For_UpSetR_Analysis
df3=pd.read_csv("Fuma_Genomic_Loci_Comparison.csv")
df5=df3[['chrom', 'start', 'end']+[x for x in df3.columns if x.endswith(":p")]]


for col in [x for x in df3.columns if x.endswith(":p")]:
    df5[col]=np.where(df5[col]==".",0,1)

df5[[ x for x in df5.columns if x.endswith(":p")]].to_csv("For_UpSetR_Analysis.csv",index=None)



#For_Vendiagram_Analysis
df3=pd.read_csv("Fuma_Genomic_Loci_Comparison.csv")
df5=df3[['chrom', 'start', 'end']+[x for x in df3.columns if x.endswith(":p")]]

for col in [x for x in df3.columns if x.endswith(":p")]:
    df5[col]=np.where(df5[col]==".",np.nan,df5['chrom'].astype('str')+"_"+df5['start'].astype('str')+"_"+df5['end'].astype("str"))

df5[[ x for x in df5.columns if x.endswith(":p")]].to_csv("For_Vendiagram_Analysis.csv",index=None)

