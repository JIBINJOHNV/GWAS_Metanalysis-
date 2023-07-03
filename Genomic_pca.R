library(data.table)
library(devtools)
library(GenomicSEM)
library(knitr)
library(psych)
library(glue)


source("/home/jjohn41/Softwares/Mycustomscript/N_weighted_GWAMA.function.1_2_6.R")
my_GWAMA<-multivariate_GWAMA

##Reference
#1)https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html
#2)https://rpubs.com/AnnaFurtjes/genomicPCA
#3)https://rpubs.com/AnnaFurtjes/genomicPCAevaluation

datafile="GenomicPCA_file.tsv"

# as reference SNP list use 1000G reference genome (HapMap 3)
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
hm3<-"/home/jjohn41/Softwares/Resourses/GenomicSEM/w_hm3.snplist"
ld<- "/home/jjohn41/Softwares/Resourses/GenomicSEM/eur_w_ld_chr/"
wld<- "/home/jjohn41/Softwares/Resourses/GenomicSEM/eur_w_ld_chr/"


datadf=fread("GenomicPCA_file.tsv")
gwas_files<-datadf$filename
trait.names<-datadf$traitname
sample.prev<-datadf$sampleprevalence
population.prev<-datadf$populationprevalence
N<-404291


#creating inputs for mungesumstat & associatio testing
for (index in 1:nrow(datadf)){

    f_name<-as.character( datadf[index,"filename"])
    o_name<-glue(as.character( datadf[index,"traitname"]),"_munge_inputs.txt")

    df<-fread(f_name)
    df<-df[,c("MarkerName","CHR","POS","A1","A2","EAF_A1","Beta","SE","Pval")]
    names(df)<-c("SNP","CHR","POS","A1","A2","eaf_A1","beta","se","p")
    df$N<-N
    fwrite(df,file=o_name, quote=FALSE,col.names=TRUE,row.names=F,sep=" ")

    ##Fil for association testing
    o_name2<-glue(as.character( datadf[index,"traitname"]),"_GenomicPCA_inputs.tsv")
    df$Z<-df$beta/df$se
    df<-df[,c("SNP","CHR","POS","A1","A2","eaf_A1","N","Z","p")]
    names(df)<-c("SNPID","CHR","BP","EA","OA","EAF","N","Z","P")
    fwrite(df,file=o_name2, quote=FALSE,col.names=TRUE,row.names=F,sep=" ") 
    }


##Munge function calling using the files saved inthe above step
munge_inputfile<-paste(datadf$traitname, "_munge_inputs.txt", sep = "")
munge_traitname<-datadf$traitname
munge(files=munge_inputfile,hm3=hm3,trait.names = munge_traitname)



########## Perform linkage disequilibrium score regression (LDSC)
traits<-paste(datadf$traitname, ".sumstats.gz", sep = "")

LDSCoutput_risk_taking<-ldsc(traits=traits,
                             sample.prev=datadf$sampleprevalence,population.prev=datadf$populationprevalence,
                             ld=ld,wld=wld,
                             trait.names = trait.names,stand = T)

save(LDSCoutput_risk_taking, file="genomicPCA_LDSC.RData")
load("genomicPCA_LDSC.RData")




# Make an empty list with length equal to the number of input files
dat<-vector("list")

gwama_files<-paste(datadf$traitname, "_GenomicPCA_inputs.tsv", sep = "")
all_files<-as.data.frame(cbind(gwama_files, file_names = datadf$traitname))
  
# Read in the data 
n<-0
for(file in 1:nrow(all_files)){
  n<-n+1
  print(n)
  dat[[all_files$file_names[n]]]<-fread(all_files$gwama_files[n],data.table = F)
}



# select the matrix of LDSC intercepts 
# this matrix should be in the same order as data read into dat but I have added some code to specifically make sure they take the same order
order<-names(dat)
dimnames(LDSCoutput_risk_taking$I)[[1]]<-dimnames(LDSCoutput_risk_taking$S)[[2]]
dimnames(LDSCoutput_risk_taking$I)[[2]]<-dimnames(LDSCoutput_risk_taking$S)[[2]]
CTI<-as.matrix(LDSCoutput_risk_taking$I[order,order])

# load correlation matrix
dimnames(LDSCoutput_risk_taking$S_Stand)[[1]]<-dimnames(LDSCoutput_risk_taking$S)[[2]]
dimnames(LDSCoutput_risk_taking$S_Stand)[[2]]<-dimnames(LDSCoutput_risk_taking$S)[[2]]
cormatrix<-LDSCoutput_risk_taking$S_Stand[order,order]

# eigen decomposition
eigenvectors<-eigen(cormatrix)$vectors
eigenvalues<-eigen(cormatrix)$values

# eigen decomposition
eigenvectors<-eigen(cormatrix)$vectors
eigenvalues<-eigen(cormatrix)$values

# calculate standardised loadings
loadings<-as.vector(eigenvectors%*%sqrt(diag(eigenvalues))[,1])


# run modified GWAMA function
my_GWAMA(x=dat,cov_Z=CTI,h2=loadings,
          out=".",name="GenomicPCA_Correlation",
          output_gz=F,check_columns=F)
