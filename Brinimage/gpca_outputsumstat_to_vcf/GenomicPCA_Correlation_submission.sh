#!/bin/bash

#SBATCH --array=1-15
#SBATCH --ntasks=4
#SBATCH --mem=120GB
#SBATCH --partition=cpu
#SBATCH --time=10:00:00

#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

module load Python/3.8.6-GCCcore-10.2.0
module load libxml2/2.9.8-GCCcore-6.4.0
module load OpenSSL/1.1
module load R/4.1.3-foss-2021b
module load GMP/6.2.1-GCCcore-11.2.0
module load BCFtools/1.14-GCC-11.2.0


source /edgehpc/dept/human_genetics/users/jjohn1/Software/gwas2vcf/env/bin/activate


python brain_image_gpca_sumstat_to_vcf.py -inputsumstat=dbscan_clust_${SLURM_ARRAY_TASK_ID}_GenomicPCA_Correlation.N_weighted_GWAMA.results.txt.gz
