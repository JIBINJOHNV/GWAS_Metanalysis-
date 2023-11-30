

gatk VariantsToTable  \
    -V PGC3_SCZ_wave3.european.autosome.public.v3.vcf.gz \
    -F CHROM -F POS -F ID -F REF -F ALT -GF AF -GF ES -GF SE -GF LP -GF SS -GF SI -GF NC \
    -O output.table

sed -i 's/PGC3_SCZ_EUR.//g' output.table
