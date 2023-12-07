## Correlation
bcftools merge \
     dbscan_clust_1_15_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz \
     /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Orignal_IDP/IDP_Brainimag.vcf.gz \
     /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Disease_Phenotype_GWAS/Disease_Pheenotype_Sumstat.vcf.gz | \
     bgzip -c > dbscan_clust_1_15_GenomicPCA_Correlation_Disease_Phenotype_IDP.vcf.gz

tabix -f -p vcf dbscan_clust_1_15_GenomicPCA_Correlation_Disease_Phenotype_IDP.vcf.gz

bcftools view --types snps dbscan_clust_1_15_GenomicPCA_Correlation_Disease_Phenotype_IDP.vcf.gz | \
       bgzip -c >dbscan_clust_1_15_GenomicPCA_Correlation_Disease_Phenotype_IDP_SNPOnly.vcf.gz
tabix -f -p vcf dbscan_clust_1_15_GenomicPCA_Correlation_Disease_Phenotype_IDP_SNPOnly.vcf.gz


#bcftools norm --multiallelics -any dbscan_clust_1_15_GenomicPCA_Correlation_Disease_Phenotype_SNPOnly.vcf.gz |\
#|awk -F"\t" '++seen[$1,$2] == 1'


## covariates
bcftools merge \
     dbscan_clust_2_15_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz \
     /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Orignal_IDP/IDP_Brainimag.vcf.gz \
     /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Disease_Phenotype_GWAS/Disease_Pheenotype_Sumstat.vcf.gz | \
     bgzip -c > dbscan_clust_2_15_GenomicPCA_covariates_Disease_Phenotype_IDP.vcf.gz

tabix -f -p vcf dbscan_clust_2_15_GenomicPCA_covariates_Disease_Phenotype_IDP.vcf.gz 

bcftools view --types snps dbscan_clust_2_15_GenomicPCA_covariates_Disease_Phenotype_IDP.vcf.gz | \
       bgzip -c >dbscan_clust_2_15_GenomicPCA_covariates_Disease_Phenotype_IDP_SNPOnly.vcf.gz

tabix -f -p vcf dbscan_clust_2_15_GenomicPCA_covariates_Disease_Phenotype_IDP_SNPOnly.vcf.gz




