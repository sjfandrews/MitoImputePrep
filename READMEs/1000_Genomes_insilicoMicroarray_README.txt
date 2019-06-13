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
/MitoImputePrep/scripts/R/MCC_Genotypes.R

LISTS USED WITHIN:
/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt

MODULES AND APPLICATIONS CALLED UPON BY DAUGHTER SCRIPTS:
python v2.7.11
R v3.4.3
bcftools v1.8
plink v1.9
impute2 v2.3.2
vt v0.57721
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

/MitoImputePrep/scripts/PYTHON/pickFirstAlt.py
[ASK BRIAN OR RUSSELL]

/MitoImputePrep/scripts/PYTHON/vcf2fasta_rCRS.py
This script takes in a VCF file for mtDNA and converts it to a FASTA formatted multiple sequence alignment file.
Sites are numbered according to the revised Cambridge Reference Sequence.
There are four options: --vcf_file, --out_file, --include-rCRS, and --verbose
--vcf_file is the input VCF file. VCF files can be haploid or diploid, but if they are diploid this programme will select the genotype on the left of / or |. So if any heteroplasmy if observed, this is not the script for you.
--out_file is the output FASTA file. If no file is specified it will be output to the same location as the VCF file, albeit with .vcf / .vcf.gz replaced with .fasta
--include-rCRS includes the revised Cambridge Reference sequence in the FASTA file.
--verbose is verbose mode.
What this script does in essence is duplicate the rCRS n times (n being the number of sequences in the VCF file).
Then, for each sequence in the VCF file, for each site in the VCF file, if the genotype that sequence has at that site matches the reference allele, pass. If is the alternative allele, replace the reference nucleotide with the alternative nucleotide.

/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py
This script takes FASTA formatted multiple sequence alignments of mtDNA and converts them to a VCF format.
Reference alleles are set to the revised Cambridge Reference Sequence.
There are 9 options.
--infile is the input FASTA file. Preferably curated to the rCRS numbering system (s=16569)
--outfile OUTFILE is the output VCF file. Includes all sites, even invariants. If no file is specified it will be output to the same location as the VCF file, albeit with .fasta replaced with .vcf
--gap2missing turn gaps (-) to missing (N) character states.
--diploid creates diploid VCF file instead of haploid (ie 0|0 instead of 0).
--ID tags the ID column as MT<POS>.
--quality tags for QUAL column (default: 999).
--filt tags for FILTER column (default: PASS).
--verbose turns on verbose mode.
--add_alt forces the VCF to always have an alternative allele.

/MitoImputePrep/scripts/R/assign_sex_label.R 
This script takes a .SAMPLES file and appends a column with "M" to denote male sex label.
There are two user inputs required:
ARG1 = input .samples file.
ARG2 = output .samples file.

/MitoImputePrep/scripts/R/FixSamplesFile_raijin.R
This script files the .samples file by changing the 4th element of the 1st row to a "D". This makes it compatible with IMPUTE2.
ARG1 = inpute .samples file.

/MitoImputePrep/scripts/R/plink_sites_map.R
This script fixes the site map .map file generated by plink. It creates a new first column with "MT" at all sites.
ARG1 = input .map file

/MitoImputePrep/scripts/R/MCC_Genotypes.R
This script calculates the Matthew's correlation coefficient for imputed + genotyped VCF file and the genotyped only file.
It takes 5 arguments.
ARG1 = The VCF file for the truthset (whole-molecule resequencing of the 1000 Genomes Phase 3 mtDNA data set.
ARG2 = The VCF file for the genotyped only dataset.
ARG3 = The VCF file for the genotyped + imputed dataset
ARG4 = The INFO file that resulted from the IMPUTE2 process.
ARG5 = The output file.


#
