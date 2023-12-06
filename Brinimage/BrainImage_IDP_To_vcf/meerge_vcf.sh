
bcftools merge \
            IDP_0356.vcf.gz IDP_0512.vcf.gz \
            IDP_0553.vcf.gz IDP_0628.vcf.gz  \
            IDP_1580.vcf.gz IDP_2089.vcf.gz | bgzip -c > IDP_Brainimag.vcf.gz
