

for file in $(cat GenomicPCA_CORR_COV_file.tsv | sed '1d'| cut -f1|tr "\n" " "); do

files1=$(grep -v $file GenomicPCA_CORR_COV_file.tsv| sed '1d' | cut -f1 | sed 's/GPCA_//g' | sed 's/_munge_input.tsv/.sumstats.gz/g' | tr "\n" ",") 
fil0=$(echo $file | sed 's/GPCA_//g' | sed 's/_munge_input.tsv/.sumstats.gz/g')
files="${fil0},${files1}"
files="${files::-1}"

s_pre0=$(grep $file GenomicPCA_CORR_COV_file.tsv| cut -f3 | sed 's/NA/nan/g' | tr "\n" ",") 
s_pre1=$(grep -v $file GenomicPCA_CORR_COV_file.tsv| sed '1d' | cut -f3 | sed 's/NA/nan/g' | tr "\n" ",") 
s_pre="${s_pre0}${s_pre1}"
s_pre="${s_pre::-1}"

p_pre0=$(grep $file GenomicPCA_CORR_COV_file.tsv| cut -f4 | sed 's/NA/nan/g' | tr "\n" ",") 
p_pre1=$(grep -v $file GenomicPCA_CORR_COV_file.tsv| sed '1d' | cut -f4 | sed 's/NA/nan/g' | tr "\n" ",")
prev="${p_pre0}${p_pre1}"
p_prev="${prev::-1}"


ldsc.py --rg ${files} \
--ref-ld-chr /edgehpc/dept/human_genetics/users/jjohn1/Software/ldsc_reference_Files/GenomicSEMFiles/eur_w_ld_chr/ \
--w-ld-chr /edgehpc/dept/human_genetics/users/jjohn1/Software/ldsc_reference_Files/GenomicSEMFiles/eur_w_ld_chr/ \
--out ${file} \
--samp-prev ${s_pre} \
--pop-prev  ${p_prev}


done 

