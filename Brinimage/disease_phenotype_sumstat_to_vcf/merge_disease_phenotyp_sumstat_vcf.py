

bcftools merge \
       PGC3_SCZ_wave3.european.autosome.public.v3.vcf.gz daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.trios.vcf.gz \
       pgc-bip2021-all.vcf.gz daner_bip_pgc3_nm_noukbiobank.vcf.gz \
       daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01.vcf.gz \
       Lam_et_al_2021_CognitiveTaskPerformance.vcf.gz \
       ADHD2022_iPSYCH_deCODE_PGC.meta.vcf.gz | bgzip -c > Disease_Pheenotype_Sumstat.vcf.gz

tabix -p vcf Disease_Pheenotype_Sumstat.vcf.gz
