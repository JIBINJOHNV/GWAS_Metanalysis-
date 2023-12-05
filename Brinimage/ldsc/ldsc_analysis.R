#srun -n 4 --mem=120G --time=24:01:00 --pty bash
#module load R/4.3.0-foss-2021b

.libPaths("/home/jjohn1/modulefiles/R4.3.0-foss-2021b")

library(GenomicSEM)
library(data.table)
library(devtools)
library(GenomicSEM)
library(knitr)
library(psych)
library(glue)
library(magrittr)
library(dplyr)
library(doSNOW)


# as reference SNP list use 1000G reference genome (HapMap 3)
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
hm3<-"/edgehpc/dept/human_genetics/users/jjohn1/Software/ldsc_reference_Files/GenomicSEMFiles/eur_w_ld_chr/w_hm3.snplist"
ld<- "/edgehpc/dept/human_genetics/users/jjohn1/Software/ldsc_reference_Files/GenomicSEMFiles/eur_w_ld_chr/"
wld<- "/edgehpc/dept/human_genetics/users/jjohn1/Software/ldsc_reference_Files/GenomicSEMFiles/eur_w_ld_chr/"


datafile="GenomicPCA_CORR_COV_file.tsv"



datadf=fread(datafile)
gwas_files<-datadf$filename
trait.names<-datadf$traitname
sample.prev<-datadf$sampleprevalence
population.prev<-datadf$populationprevalence


##Munge function calling using the files saved inthe above step
munge_inputfile<-datadf$filename
munge_traitname<-datadf$traitname
munge(files=munge_inputfile,hm3=hm3,trait.names = munge_traitname,info.filter = 0.3,maf.filter=0.01)

########## Perform linkage disequilibrium score regression (LDSC)
traits<-paste(datadf$traitname, ".sumstats.gz", sep = "")

LDSCoutput_risk_taking<-ldsc(traits=traits,
                             sample.prev=datadf$sampleprevalence,population.prev=datadf$populationprevalence,
                             ld=ld,wld=wld,
                             trait.names = trait.names,stand = T)

save(LDSCoutput_risk_taking, file="genomicPCA_Disease_LDSC.RData")
load("genomicPCA_Disease_LDSC.RData")


#Genetic covariance ; S	 ; header prsent ; 37
dimnames(LDSCoutput_risk_taking$S)[[1]]<-dimnames(LDSCoutput_risk_taking$S)[[2]]
covariance_df<-LDSCoutput_risk_taking$S
write.csv(covariance_df,"genomicPCA_Disease_LDSC_covariance.csv")

#Genetic covariance ; S	 ; header prsent ; 37 ; S_Stand
dimnames(LDSCoutput_risk_taking$S_Stand)[[1]]<-dimnames(LDSCoutput_risk_taking$S_Stand)[[2]]
covariance_stand_df<-LDSCoutput_risk_taking$S_Stand
write.csv(covariance_stand_df,"genomicPCA_Disease_LDSC_covariance_Standerdised.csv")


#intercepts I;  header prsent ; 37
dimnames(LDSCoutput_risk_taking$I)[[1]]<-dimnames(LDSCoutput_risk_taking$I)[[2]]
intervept_df<-LDSCoutput_risk_taking$I
write.csv(intervept_df,"genomicPCA_Disease_LDSC_Intercept.csv")


#variance covariance matrix of the parameter estimates in S ; 703 ; no header
LDSCoutput_risk_taking$V


#variance covariance matrix of the parameter estimates in S ; 703 ; no header
LDSCoutput_risk_taking$V_Stand





hdl.covstruct <- hdl(traits,
                     sample.prev = datadf$sampleprevalence,
                     population.prev = datadf$populationprevalence ,
                     trait.names=trait.names,
                     LD.path="/edgehpc/dept/human_genetics/users/jjohn1/Software/UKB_imputed_hapmap2_SVD_eigen99_extraction/", method = "jackknife")


save(hdl.covstruct, file="hdl_covstruct.RData")

hdl.covstruct.2 <- hdl(traits,
                     sample.prev = datadf$sampleprevalence,
                     population.prev = datadf$populationprevalence ,
                     trait.names=trait.names,
                     LD.path="/edgehpc/dept/human_genetics/users/jjohn1/Software/UKB_imputed_hapmap2_SVD_eigen99_extraction/", method = "piecewise")

save(hdl.covstruct.2, file="hdl_covstruct2.RData")

