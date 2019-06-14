The files described in here were used for analysis of the resulting imputation pipeline experiments.
Descriptions for each file found below.

/MitoImputePrep/scripts/R/ANALYSIS/calculate_95CI.R
This script contains the formula/function to calculate the 95% CI intervals.
While we used it on certain metrics stated in the paper, this can be used on any numeric values.

/MitoImputePrep/scripts/R/ANALYSIS/HiMC/check_haplogroup_concordance.R
This script iteratively reads through all the experiments (KHAP, MAF, MCMC).
First the haplogroups are assiged for each sample in the truth set, then haplogroups are assigned using HiMC.
Then for each experiment and each chip, haplogroups are assigned.
The same happens for macro-haplogroups, where for non-L haplogroups the haplogroup assigned is truncated to the major haplogroup (ie H), then compared to the truth set.
For L haplogroups, the haplogroup is truncated to the first subclade (ie L3).

/MitoImputePrep/scripts/R/ANALYSIS/HiMC/Cleanup_concordance_tables_HiMC.R
Following on from the check_haplogroup_concordance.R script, this script calculates haplogroup and macrohaplogroup concordance.
Summary statistics are then calculated.
Each individual experiment is saved, then all are combined into one big data table.
Plots are able to be produced from this script, however the HiMC_1kGP_plots described later produce the most up to date plots.

/MitoImputePrep/scripts/R/ANALYSIS/HiMC/concordance_tables_ADNI.R
Similar to check_haplogroup_concordance.R, this script assigns haplogroups for the ADNI imputed dataset.
The samples from the Whole Genome re-Sequenced ADNI3 dataset have their haplogroups assigned and used as the truthset.
The samples from the genotyped-only ADNI1 dataset have their haplogroups assigned.
Haplogroup (and macrohaplogroup) assignment is then compared and the concordance calculated.

/MitoImputePrep/scripts/R/ANALYSIS/MCC/Cleanup_concordance_tables_MCC.R
This script is similar in function to Cleanup_concordance_tables_HiMC.R, however it performs Matthew's correlation coefficient calculations instead of 

/MitoImputePrep/scripts/R/ANALYSIS/MCC/MCC_emmeans.R

/MitoImputePrep/scripts/R/Plotting/HiMC_1kGP_plots.R
/MitoImputePrep/scripts/R/Plotting/MCC_concordance_tables.R