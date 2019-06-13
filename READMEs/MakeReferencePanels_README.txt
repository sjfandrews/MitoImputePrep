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
N/A

<PYTHON SCRIPTS>

/MitoImputePrep/scripts/PYTHON/ambiguous2missing.py
This script takes in a FASTA formatted multiple sequence alignment file and converts ambiguous character states to missing character states ("N").
One option allows for gap character states ("-") also be converted to missing character states.
There are four options: --infile, --outfile, --gap2missing, --verbose.
--infile is the input FASTA file.
--outfile is the output FASTA file.
--gap2missing turns on the mode whereby gap character states ("-") will be converted to missing character states.
--verbose is verbose mode.

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

<R SCRIPTS>

/MitoImputePrep/scripts/R/DATA_PROCESSING/removeLowQuality_cmdline.R
This script takes in a FASTA formatted multiple sequence alignment file and removes sequences according to user specified quality control criteria.
This QC is on the basis of the number of allowed missing character states and the number of ambiguous character states.
This script takes four arguments:
ARG1 = The input FASTA file.
ARG2 = The output text file containing only high quality sequences.
ARG3 = The maximum number of missing character states allowed.
ARG4 = The maximum number of gap character states allowed.

/MitoImputePrep/scripts/R/DATA_PROCESSING/assign_sex_label.R
This script takes in a .sample file and adds a new column with an "M" sex label.
This script takes one argument:
ARG1 = The input .sample file.

/MitoImputePrep/scripts/R/DATA_PROCESSING/mt_recombination_map.R
This file takes in a VCF file and produces a recombination map for mitochondrial DNA.
Because we assume no recombination, the recombination rate for all sites is r = 0.
This file takes three arguments:
ARG1 = The input VCF file.
ARG2 = The output recombination map map file.
ARG3 = The output strand file.

#
