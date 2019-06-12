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
	- 
	-

Files used in the imputation pipeline for 1000 Genomes Project data (in silico microarrays)
This script is used to iteratively submit one job per in silico microarray to the RAIJIN cluster at NCI (National Computational Infrastructure) in Canberra, Australia.
This script takes one argument, which is the reference panel being used from which missing variants will be imputed.

BASH file: /MitoImputePrep/scripts/BASH/DATA_PROCESSING/create_RefPans/make_RefPan_v${NUM}.sh (Where ${NUM} is the version)
This file is the parent script.

USAGE:
<USEAGE>
<DESCRIPTION>
<OPTIONS>

SCRIPTS USED WITHIN:
BASH:
<LIST OF BASH SCRIPTS>

PYTHON:
<LIST OF PYTHON SCRIPTS>

R:
<LIST OF R SCRIPTS>

LISTS USED WITHIN:
<LISTS ITERATED OVER>

MODULES AND APPLICATIONS CALLED UPON BY DAUGHTER SCRIPTS:
<MODULE | VERSION>


DETAILS:
<BASH SCRIPTS>
<WHAT THEY ARE CALLED, WHAT THEY DO>
<ARGUEMENTS>

<PYTHON SCRIPTS>
<WHAT THEY ARE CALLED, WHAT THEY DO>
<ARGUEMENTS>

<R SCRIPTS>
<WHAT THEY ARE CALLED, WHAT THEY DO>
<ARGUEMENTS>

#
