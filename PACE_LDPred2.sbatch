#!/bin/bash
#SBATCH -J LDPRed2_all_PrCa
#SBATCH -A gts-jlachance6-biocluster
#SBATCH -C cpu-amd
#SBATCH --ntasks=1 --cpus-per-task=32
#SBATCH --mem=0
#SBATCH -t7200
#SBATCH -q inferno
#SBATCH -o LDPRed2_all_PrCa.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jvaz7@gatech.edu

# change directory to the directory you submitted the script from
cd $SLURM_SUBMIT_DIR
# load anaconda and other software you will need
module load anaconda3/2022.05
module load r/4.2.1-bio
#conda create --name r-environment r-essentials r-base r-remotes r-magrittr r-data.table r-devtools r-fmsb -y
#conda create --name r-environment -y
conda activate r-environment
#conda install -c r-remotes -y

Rscript ./LDPred2_GWide_PACE.R
