Files used to create different versions of the reference panel for the MitoImpute pipeline.
This script takes the reference panel fasta multiple sequence alignment file as input and produces a VCF file of the reference panel, as well as PLINK and OXFORD formats.
This script:
	- Converts the ambiguous characters to missing character states.
	- Converts the fasta file to a VCF file.
	- Performs QC by removing sequences we define to be low quality (missing character states > 7, including gaps as they would have been converted to missing).
	- Filters to a Minor Allele Frequency specified in the script.
	- Extracts sample IDs and adds a male sex label to them.
	- Converts the reference panel VCF to PLINK, OXFORD, GEN, SAMPLE formats.
	- Creates a recombination map (all sites r = 0, as no recombination as assumed).
	- Copies the relevant files to the Git repository.


BASH file: /MitoImputePrep/scripts/BASH/DATA_PROCESSING/create_RefPans/make_RefPan.sh
This file is the parent script.

USAGE:
$ /MitoImputePrep/scripts/BASH/DATA_PROCESSING/create_RefPans/make_RefPan.sh <PATH_TO_REFERENCE_PANEL_FASTA> <MINOR_ALLELE_FREQUENCY>
Where <PATH_TO_REFERENCE_PANEL_FASTA> <Minor_Allele_Frequency> are variables. <PATH_TO_REFERENCE_PANEL_FASTA> is the path to the reference panel fasta files, and <MINOR_ALLELE_FREQUENCY> is the minor allele frequency it will be filtered to.
If you include a path to a reference panel that doesn't exist (or contains typos, etc), it will default to the current reference panel.
Minor allele frequency needs to be set between >0.0 and <1.0.

SCRIPTS USED WITHIN:
BASH:
N/A

PYTHON:
/MitoImputePrep/scripts/PYTHON/ambiguous2missing.py
/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py

R:
/MitoImputePrep/scripts/R/DATA_PROCESSING/removeLowQuality_cmdline.R
/MitoImputePrep/scripts/R/DATA_PROCESSING/assign_sex_label.R
/MitoImputePrep/scripts/R/DATA_PROCESSING/mt_recombination_map.R

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
<WHAT THEY ARE CALLED, WHAT THEY DO>
<ARGUEMENTS>

/MitoImputePrep/scripts/PYTHON/ambiguous2missing.py
/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py

<R SCRIPTS>
<WHAT THEY ARE CALLED, WHAT THEY DO>
<ARGUEMENTS>

/MitoImputePrep/scripts/R/DATA_PROCESSING/removeLowQuality_cmdline.R
/MitoImputePrep/scripts/R/DATA_PROCESSING/assign_sex_label.R
/MitoImputePrep/scripts/R/DATA_PROCESSING/mt_recombination_map.R

#
