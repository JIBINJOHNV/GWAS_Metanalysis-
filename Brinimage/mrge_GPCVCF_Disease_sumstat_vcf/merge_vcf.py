
#### GenomicPCA_covariates
bcftools merge \
        dbscan_clust_2_15_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz \
        Disease_Phenotype_GWAS/MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz | bgzip -c > dbscan_clust_2_15_GenomicPCA_covariates_MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz 

tabix -p vcf dbscan_clust_2_15_GenomicPCA_covariates_MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz

bcftools view \
      --types snps dbscan_clust_2_15_GenomicPCA_covariates_MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz | bgzip -c >dbscan_clust_2_15_GenomicPCA_covariates_MDDwoBP_BIP_pgc3_SCZ_PGC3_SNPOnly.vcf.gz

tabix -p vcf dbscan_clust_2_15_GenomicPCA_covariates_MDDwoBP_BIP_pgc3_SCZ_PGC3_SNPOnly.vcf.gz





####GenomicPCA_Correlation

bcftools merge \
      dbscan_clust_1_15_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz \
      Disease_Phenotype_GWAS/MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz | bgzip -c > dbscan_clust_1_15_GenomicPCA_Correlation_MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz 

tabix -p vcf dbscan_clust_1_15_GenomicPCA_Correlation_MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz


bcftools view \
      --types snps dbscan_clust_1_15_GenomicPCA_Correlation_MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz | bgzip -c >dbscan_clust_1_15_GenomicPCA_Correlation_MDDwoBP_BIP_pgc3_SCZ_PGC3_SNPOnly.vcf.gz

tabix -p vcf dbscan_clust_1_15_GenomicPCA_Correlation_MDDwoBP_BIP_pgc3_SCZ_PGC3_SNPOnly.vcf.gz
