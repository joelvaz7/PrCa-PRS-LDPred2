# PRS for prostate cancer - LDPred2

## Files
* Ethnic-grouping.R: Convert Level 1 ethnicity to top-level in UKBB.
* GenesisGWAS-PACE.sbatch: Running GenesisGWAS.R on PACE.
* GenesisGWAS.R: Run Genesis to generate summary statistics - genome-wide.
* GenesisGWAS_loop.R: Run Genesis to generate summary statistics - chromosome-by-chromosome.
* LDPred2_GWide_PACE.R: Genome-wide LDPred2 to generate PRS for the target population (UKBB PrCa here), using summary statistics of the MADCaP population.
* LDPred2_PACE.R: LDPred2 on PACE.
* PACE_LDPred2.sbatch: sbatch file for LDPred2 on PACE.
* PRScs-mods.R: Converting Genesis summ-stats to LDPred2-readable format.

## Data
**All files can be found on the SMB drive**
* UKBB PrCa data: //biolachance.biology.gatech.edu/projects/UKBB/PrCa_all
* MaDCaP data: //biolachance.biology.gatech.edu/projects/MADCaP_combined
