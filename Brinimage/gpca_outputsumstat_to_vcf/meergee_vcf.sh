bcftools merge \
      dbscan_clust_1_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz dbscan_clust_2_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz \
      dbscan_clust_3_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz dbscan_clust_4_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz \
      dbscan_clust_5_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz dbscan_clust_6_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz \
      dbscan_clust_7_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz dbscan_clust_8_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz \
      dbscan_clust_9_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz dbscan_clust_10_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz \
      dbscan_clust_11_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz dbscan_clust_12_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz \
      dbscan_clust_13_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz dbscan_clust_14_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz \
      dbscan_clust_15_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz | bgzip -c > dbscan_clust_1_15_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz

tabix -p vcf dbscan_clust_1_15_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz

bcftools merge \
    dbscan_clust_2_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz dbscan_clust_3_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz \
    dbscan_clust_4_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz dbscan_clust_5_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz \
    dbscan_clust_6_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz dbscan_clust_7_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz \
    dbscan_clust_8_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz dbscan_clust_9_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz \
    dbscan_clust_10_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz dbscan_clust_11_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz \
    dbscan_clust_12_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz dbscan_clust_13_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz \
    dbscan_clust_14_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz \
    dbscan_clust_15_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz | bgzip -c > dbscan_clust_2_15_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz

tabix -p vcf dbscan_clust_2_15_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz
