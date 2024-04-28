#!/bin/bash -l

#SBATCH --account=lp_edu_large_omics
#SBATCH --cluster=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=snps_analysis
#SBATCH --time=60:00   # Walltime in HH:MM:SS
#SBATCH --mem=15G     # Memory in gigabytes
#SBATCH --cpus-per-task=8   # Number of CPU cores
#SBATCH --output=slurm_job.out   # Standard output file
#SBATCH --error=slurm_job.err   # Standard error file
#SBATCH --mail-type=END,FAIL # Specify the type of events to receive notifications for (ALL, BEGIN, END, FAIL, REQUEUE, etc.)
#SBATCH --mail-user=theodoro.gasperinterracamargo@student.kuleuven.be

export PATH=/lustre1/project/stg_00079/teaching/I0U19a_conda_2024/bin/:$PATH

# Run the SNPs analysis script
snakemake -c8 -s "./Snakefile_SNP_Call"
snakemake -c8 -s "./Snakefile_SNP_Call" --report report.html
zip submission.zip Snakefile_SNP_Call plot_graphs.py analysis_snps.slurm snp_call_snakemake/genes.vcf snp_call_snakemake/report.html