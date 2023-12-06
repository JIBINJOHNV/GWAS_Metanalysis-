#!/bin/bash

#SBATCH --array=1-6
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
# Set the input files
input_files=("0356.txt.gz" "0512.txt.gz" "0553.txt.gz" "0628.txt.gz" "1580.txt.gz" "2089.txt.gz")

# Set the input file for the current array index
input_file=${input_files[$SLURM_ARRAY_TASK_ID-1]}


python BrainImage_IDP_To_vcf.py -inputsumstat="$input_file"


##Location /mnt/depts/dept04/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/Orignal_IDP
