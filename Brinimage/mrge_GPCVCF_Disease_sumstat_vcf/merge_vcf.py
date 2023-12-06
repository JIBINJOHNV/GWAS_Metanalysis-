
#### GenomicPCA_covariates
bcftools merge \
        dbscan_clust_2_15_GenomicPCA_covariates.N_weighted_GWAMA.results.vcf.gz \
        Disease_Phenotype_GWAS/Disease_Pheenotype_Sumstat.vcf.gz \
        /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Orignal_IDP/IDP_Brainimag.vcf.gz | \
        bgzip -c > dbscan_clust_2_15_GenomicPCA_covariates_MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz 

tabix -p vcf dbscan_clust_2_15_GenomicPCA_covariates_Disease_Phenotype.vcf.gz



##Select SNPS only from the vcf file
bcftools view \
      --types snps dbscan_clust_2_15_GenomicPCA_covariates_MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz | bgzip -c >dbscan_clust_2_15_GenomicPCA_covariates_MDDwoBP_BIP_pgc3_SCZ_PGC3_SNPOnly.vcf.gz

tabix -p vcf dbscan_clust_2_15_GenomicPCA_covariates_MDDwoBP_BIP_pgc3_SCZ_PGC3_SNPOnly.vcf.gz


## Remove multi alleleic variants
bcftools norm --multiallelics +any \
        dbscan_clust_2_15_GenomicPCA_covariates_MDDwoBP_BIP_pgc3_SCZ_PGC3_SNPOnly.vcf.gz |bcftools view --max-alleles 2 \
        bgzip -c > dbscan_clust_2_15_GenomicPCA_covariates_MDDwoBP_BIP_pgc3_SCZ_PGC3_SNPOnly_NoMAV.vcf.gz
        



####GenomicPCA_Correlation
bcftools merge \
      dbscan_clust_1_15_GenomicPCA_Correlation.N_weighted_GWAMA.results.vcf.gz \
      Disease_Phenotype_GWAS/MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz | bgzip -c > dbscan_clust_1_15_GenomicPCA_Correlation_MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz 

tabix -p vcf dbscan_clust_1_15_GenomicPCA_Correlation_MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz


##Select SNPS only from the vcf file
bcftools view \
      --types snps dbscan_clust_1_15_GenomicPCA_Correlation_MDDwoBP_BIP_pgc3_SCZ_PGC3.vcf.gz | bgzip -c >dbscan_clust_1_15_GenomicPCA_Correlation_MDDwoBP_BIP_pgc3_SCZ_PGC3_SNPOnly.vcf.gz

tabix -p vcf dbscan_clust_1_15_GenomicPCA_Correlation_MDDwoBP_BIP_pgc3_SCZ_PGC3_SNPOnly.vcf.gz




