Files used in the imputation pipeline for 1000 Genomes Project data (in silico microarrays)
This script is used to iteratively submit one job per in silico microarray to the RAIJIN cluster at NCI (National Computational Infrastructure) in Canberra, Australia.
This script takes one argument, which is the reference panel being used from which missing variants will be imputed.

BASH file: /MitoImputePrep/scripts/BASH/DATA_PROCESSING/ThousandGenomes_Imputation/MassDeploy_ThousandGenomes_imputeFrom_RefPan.sh
This file is the parent script.

USAGE:
$ sh /MitoImputePrep/scripts/BASH/DATA_PROCESSING/ThousandGenomes_Imputation/MassDeploy_ThousandGenomes_imputeFrom_RefPan.sh <REFERENCE_PANEL>
Where <REFERENCE_PANEL> is a variable. You need to note the reference panel that will be used (exactly as it is in the reference panel file in XXX).
You may modify the MCMC, BURNIN, and KHAP settings as required for your particular analyses. However, these are set to these defaults for the 'recommended settings':
MCMC = 1
BURNIN = 0
KHAP = 500

SCRIPTS USED WITHIN:
BASH:
/MitoImputePrep/scripts/BASH/ThousandGenomes_imputeFrom_RefPan.sh

PYTHON:
/MitoImputePrep/scripts/PYTHON/pickFirstAlt.py
/MitoImputePrep/scripts/PYTHON/vcf2fasta_rCRS.py
/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py

R:
/MitoImputePrep/scripts/R/assign_sex_label.R 
/MitoImputePrep/scripts/R/FixSamplesFile_raijin.R
/MitoImputePrep/scripts/R/plink_sites_map.R

LISTS USED WITHIN:
/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt

MODULES AND APPLICATIONS CALLED UPON BY DAUGHTER SCRIPTS:
python v2.7.11
R v3.4.3
bcftools v1.8
plink v1.9
impute2 v2.3.2
vt
java/jdk v1.8.0_60
HaploGrep v2.1.19


DETAILS:
/MitoImputePrep/scripts/BASH/ThousandGenomes_imputeFrom_RefPan.sh
This script is used to submit jobs to the RAIJIN cluster.
Using the qsub command, it takes 6 arguments:
REFpanel – the reference panel from which missing variants will be imputed.
MtPlatforms – the name of the in silico microarray whose missing variants will be imputed.
mcmc – length of the Markov chain Monte Carlo
burn – burn-in length of the Markov chain Monte Carlo
khap – The number of reference haplotypes to be used
ne – the effective population size.
The REFpanel option is set by the argument passed to the parent BASH script 
The MtPlatforms option is what is being iteratively changed by the parent script.
All other options remain constant and are set in the parent BASH script.
