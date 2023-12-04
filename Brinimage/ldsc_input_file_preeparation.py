

bcftools merge \
     dbscan_clust_1_15_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz \
     /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Disease_Phenotype_GWAS/Disease_Pheenotype_Sumstat.vcf.gz | \
     bgzip -c > dbscan_clust_1_15_GenomicPCA_Correlation_Disease_Phenotype.vcf.gz

tabix -p vcf dbscan_clust_1_15_GenomicPCA_Correlation_Disease_Phenotype.vcf.gz


bcftools merge \
     dbscan_clust_2_15_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz \
     /edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Disease_Phenotype_GWAS/Disease_Pheenotype_Sumstat.vcf.gz | \
     bgzip -c > dbscan_clust_1_15_GenomicPCA_covariates_Disease_Phenotype.vcf.gz

tabix -p vcf dbscan_clust_1_15_GenomicPCA_covariates_Disease_Phenotype.vcf.gz



gatk VariantsToTable  \
    -V PGC3_SCZ_wave3.european.autosome.public.v3.vcf.gz \
    -F CHROM -F POS -F ID -F REF -F ALT -GF AF -GF ES -GF SE -GF LP -GF SS -GF SI -GF NC \
    -O output.table

sed -i 's/PGC3_SCZ_EUR.//g' output.table
