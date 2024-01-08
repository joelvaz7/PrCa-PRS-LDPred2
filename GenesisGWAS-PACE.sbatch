#!/bin/bash
#SBATCH -J GENESIS-GWAS-GWide
#SBATCH -A gts-jlachance6-biocluster
#SBATCH -N2 --ntasks-per-node=24
#SBATCH --mem-per-cpu=20G
#SBATCH -t2880
#SBATCH -q inferno
#SBATCH -o GENESIS-GWAS-GWide.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jvaz7@gatech.edu

# change directory to the directory you submitted the script from
cd $SLURM_SUBMIT_DIR
# load anaconda and other software you will need
module load anaconda3/2022.05
conda activate r-environment
module load r/4.2.1-bio

Rscript ./GenesisGWAS-PACE.R
