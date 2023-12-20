#!/bin/bash

#SBATCH --array=1-2
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


# Define your Python script and input file directory 
PYTHON_SCRIPT="CLuster12_13_SNPs_sumstat_to_vcf.py" 
INPUT_DIR="/edgehpc/dept/human_genetics/users/jjohn1/Outcome_GWAS/All_GPCA_MetaSumstat/cluster_16_17_Meta" # Define the input files 
INPUT_FILES=("CLuster12_13_SNPs_Correlation_METAANALYSIS_1.tbl" "CLuster12_13_SNPs_Covariates_METAANALYSIS_1.tbl") # Get the input file for the current array task 
INPUT_FILE=${INPUT_FILES[${SLURM_ARRAY_TASK_ID}-1]} # Run your Python script with the selected input file 

python ${PYTHON_SCRIPT} --inputsumstat ${INPUT_DIR}/${INPUT_FILE}


