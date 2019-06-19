Files used in the imputation pipeline for Alzheimer's Disease Neuroimaging Initiative datasets.
This script takes one argument, which is the reference panel being used from which missing variants will be imputed.


BASH files:
	- /MitoImputePrep/scripts/BASH/DATA_PROCESSING/ADNI_imputation/Impute_ADNI_redo.sh # For the 258 genotyped samples that appear in both ADNI1 and ADNI3
	- /MitoImputePrep/scripts/BASH/DATA_PROCESSING/ADNI_imputation/Impute_ADNI_redo_noReSeq.sh # For the 499 genotyped samples that appear only in ADNI1
	- /MitoImputePrep/scripts/BASH/DATA_PROCESSING/ADNI_imputation/Impute_ADNI_12GO.sh # For the 1199 samples that appear in ADNI GO
These file are the parent scripts. They all do the same thing, but for different data sets (see below).

USAGE:
$ /MitoImputePrep/scripts/BASH/DATA_PROCESSING/create_RefPans/Impute_ADNI_redo.sh <REFERENCE_PANEL_VERSION>
Where <REFERENCE_PANEL_VERSION> is a variable.
<REFERENCE_PANEL_VERSION> is reference panel version and the minor allele frequency (ie ReferencePanel_v1_0.01 where ReferencePanel_v1 is the reference panel is 0.01 is the minor allele frequency).

SCRIPTS USED WITHIN:
BASH:
N/A

PYTHON:
/MitoImputePrep/scripts/PYTHON/fix_vcf_names.py

R:
/MitoImputePrep/scripts/R/DATA_PROCESSING/assign_sex_label.R
/MitoImputePrep/scripts/R/FixSamplesFile_raijin.R
/MitoImputePrep/scripts/R/HiMC_haplogroup_assignment.R
/MitoImputePrep/scripts/R/MCC_Genotypes.R

LISTS USED WITHIN:
N/A

MODULES AND APPLICATIONS CALLED UPON BY DAUGHTER SCRIPTS:
python v2.7.11
R v3.4.3
bcftools v1.4.1
bcftools v1.9
plink v1.9
vt v0.57721


DETAILS:
<BASH SCRIPTS>
<WHAT THEY ARE CALLED, WHAT THEY DO>
<ARGUEMENTS>

<PYTHON SCRIPTS>

/MitoImputePrep/scripts/PYTHON/fix_vcf_names.py
This file takes in a VCF file and fixes the names of the samples in the header.
This is done so that the same ADNI1 and ADNI3 samples can be directly and easily compared later.
This script takes 4 arguments: --vcf_file, --out_file, --csv_file, --verbose.
--vcf_file is the input VCF file.
--out_file is the output VCF file with the sample names fixed.
--csv_file is the csv file containing file names and comparisons.
--verbose turns on verbose mode.

<R SCRIPTS>

/MitoImputePrep/scripts/R/DATA_PROCESSING/assign_sex_label.R
This script takes in a .sample file and adds a new column with an "M" sex label.
This script takes one argument:
ARG1 = The input .sample file.

/MitoImputePrep/scripts/R/FixSamplesFile_raijin.R
This script files the .samples file by changing the 4th element of the 1st row to a "D". This makes it compatible with IMPUTE2.
ARG1 = inpute .samples file.

/MitoImputePrep/scripts/R/HiMC_haplogroup_assignment.R
This script takes in .ped and .map files from the imputed data set and assigns haplogroups according PhyloTree 17 via the HiMC package.
This script takes 3 arguments:
ARG1 = Input .ped file.
ARG2 = Input .map file.
ARG3 = Output .csv file.

/MitoImputePrep/scripts/R/MCC_Genotypes.R
This script calculates the Matthew's correlation coefficient for imputed + genotyped VCF file and the genotyped only file.
It takes 5 arguments.
ARG1 = The VCF file for the truthset (whole-molecule resequencing of the 1000 Genomes Phase 3 mtDNA data set.
ARG2 = The VCF file for the genotyped only dataset.
ARG3 = The VCF file for the genotyped + imputed dataset
ARG4 = The INFO file that resulted from the IMPUTE2 process.
ARG5 = The output file.


#
