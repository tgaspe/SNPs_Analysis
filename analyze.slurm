#!/bin/bash -l

#SBATCH --account=lp_edu_large_omics
#SBATCH --cluster=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=chromosome1_analysis
#SBATCH --time=15:00   # Walltime in HH:MM:SS
#SBATCH --mem=15G     # Memory in gigabytes
#SBATCH --cpus-per-task=1   # Number of CPU cores
#SBATCH --output=snps_analysis.out   # Standard output file
#SBATCH --error=snps_analysis.err   # Standard error file

# Make the analysis script executable
#chmod u+x analyze_GFF_features.sh

export PATH=/lustre1/project/stg_00079/teaching/I0U19a_conda_2024/bin/:$PATH

# Run the SNPs analysis script
snakemake -c1 -s "./Snakefile_SNP_Call"

