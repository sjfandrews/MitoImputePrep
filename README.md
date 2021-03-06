# MitoImputePrep
A snakemake pipeline for preparing the datafiles required for [MitoImpute](https://github.com/sjfandrews/MitoImpute) and [MitoImputeShiny](https://github.com/sjfandrews/MitoImputeShiny).

## Reference Panel
For mitochondrial imputaiton, a custom reference panel was constructed. 44,299 mitochondrial sequences were downloaded from NCBI using the search term:

```(016500[SLEN]:016600[SLEN]) AND Homo[Organism] AND mitochondrion[FILT] AND complete genome NOT (Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus OR ancient human remains OR shotgun)```

The 44,299 sequences were aligned to Easteal, Jermiin, and Ott mitochondrial master alignment (~8,000 manually aligned sequneces...) using Geneious 10.2.6.

The **ReferencePanel.smk** files provids a pipeline for processing a FASTA alignment of mitochondrial sequnces for use as a reference panel in IMPUTE2. The steps in the pipeline are:
1. Convert [ambigous nucleotide positions](https://en.wikipedia.org/wiki/Nucleic_acid_notation) to missing (-)
2. Convert FASTA files to VCF
3. Identify/remove samples with more then eight gaps in sequence
4. Remove SNPs
   - decompose multiallelic variants
   - fill info field on vcf
   - remove SNPs where the alternate allele is missing (-)
   - remove SNPs with a MAF < 0.01
5. Assign Male sex to Reference samples
   - Extract a list of sample names from vcf
   - Add additional column assigning male sex
6. Convert VCF files to:
   - gen/sample format used by IMPUTE2
   - hap/legend/sample format used by IMPUTE2
7. Construct a fine-scale recombination map for the mitochondria (i.e 0)

## Usage for testing and validation prior to SnakeMake.
Before the SnakeMake pipeline was constructed, the pipeline was written in bash scripting and tested on a cluster computer.
This was also done to test how the pipeline performs in terms of imputation accuracy.
This section provides details on that initial pipeline. Henceforth, this shall be referred to as the test pipeline.

### Creation of the Reference Panel
Files used to create different versions of the reference panel for the MitoImpute pipeline.
This script takes the reference panel fasta multiple sequence alignment file as input and produces a VCF file of the reference panel, as well as PLINK and OXFORD formats.
This script:
*	Converts the ambiguous characters to missing character states.
*	Converts the fasta file to a VCF file.
*	Performs QC by removing sequences we define to be low quality (missing character states > 7, including gaps as they would have been converted to missing).
*	Filters to a Minor Allele Frequency specified in the script.
*	Extracts sample IDs and adds a male sex label to them.
*	Converts the reference panel VCF to PLINK, OXFORD, GEN, SAMPLE formats.
*	Creates a recombination map (all sites r = 0, as no recombination is assumed).
*	Copies the relevant files to the Git repository.


BASH file: [make_RefPan.sh](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/BASH/DATA_PROCESSING/create_RefPans/make_RefPan.sh)
This file is the parent script.

USAGE:
```
$ make_RefPan.sh <PATH_TO_REFERENCE_PANEL_FASTA> <MINOR_ALLELE_FREQUENCY>
```
Where `<PATH_TO_REFERENCE_PANEL_FASTA>` and `<Minor_Allele_Frequency>` are variables.
`<PATH_TO_REFERENCE_PANEL_FASTA>` is the path to the reference panel fasta files, and `<MINOR_ALLELE_FREQUENCY>` is the minor allele frequency it will be filtered to.
If you include a path to a reference panel that doesn't exist (or contains typos, etc), it will default to the current reference panel.
Minor allele frequency needs to be set between >0.0 and <1.0.

#### SCRIPTS USED WITHIN:
##### BASH:
*	N/A

##### PYTHON:
* 	[ambiguous2missing.py](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/PYTHON/ambiguous2missing.py)
* 	[fasta2vcf_mtDNA.py](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/PYTHON/fasta2vcf_mtDNA.py)

##### R:
*	[removeLowQuality_cmdline.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/DATA_PROCESSING/removeLowQuality_cmdline.R)
*	[assign_sex_label.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/DATA_PROCESSING/assign_sex_label.R)
*	[mt_recombination_map.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/DATA_PROCESSING/mt_recombination_map.R)

##### LISTS USED WITHIN:
*	N/A

##### MODULES AND APPLICATIONS CALLED UPON BY DAUGHTER SCRIPTS:
*	python v2.7.11
*	R v3.4.3
*	bcftools v1.4.1
*	bcftools v1.9
*	plink v1.9
*	vt v0.57721


#### DETAILS:
##### BASH SCRIPTS
*	N/A

##### PYTHON SCRIPTS

[ambiguous2missing.py](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/PYTHON/ambiguous2missing.py)
This script takes in a FASTA formatted multiple sequence alignment file and converts ambiguous character states to missing character states ("N").
One option allows for gap character states ("-") also be converted to missing character states.
There are four options:
*	`--infile` is the input FASTA file.
*	`--outfile` is the output FASTA file.
*	`--gap2missing` turns on the mode whereby gap character states ("-") will be converted to missing character states.
*	`--verbose` is verbose mode.

[fasta2vcf_mtDNA.py](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/PYTHON/fasta2vcf_mtDNA.py)
This script takes FASTA formatted multiple sequence alignments of mtDNA and converts them to a VCF format.
Reference alleles are set to the revised Cambridge Reference Sequence.
There are 9 options.
*	`--infile` is the input FASTA file. Preferably curated to the rCRS numbering system (s=16569)
*	`--outfile` OUTFILE is the output VCF file. Includes all sites, even invariants. If no file is specified it will be output to the same location as the VCF file, albeit with .fasta replaced with .vcf
*	`--gap2missing` turn gaps (-) to missing (N) character states.
*	`--diploid` creates diploid VCF file instead of haploid (ie 0|0 instead of 0).
*	`--ID` tags the ID column as MT<POS>.
*	`--quality` tags for QUAL column (default: 999).
*	`--filt` tags for FILTER column (default: PASS).
*	`--verbose` turns on verbose mode.
*	`--add_alt` forces the VCF to always have an alternative allele.

##### R SCRIPTS

[removeLowQuality_cmdline.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/DATA_PROCESSING/removeLowQuality_cmdline.R)
This script takes in a FASTA formatted multiple sequence alignment file and removes sequences according to user specified quality control criteria.
This QC is on the basis of the number of allowed missing character states and the number of ambiguous character states.
This script takes four arguments:
*	`ARG1` = The input FASTA file.
*	`ARG2` = The output text file containing only high quality sequences.
*	`ARG3` = The maximum number of missing character states allowed.
*	`ARG4` = The maximum number of gap character states allowed.

[assign_sex_label.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/DATA_PROCESSING/assign_sex_label.R)
This script takes in a .sample file and adds a new column with an "M" sex label.
This script takes one argument:
*	`ARG1` = The input .sample file.

[mt_recombination_map.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/DATA_PROCESSING/mt_recombination_map.R)
This file takes in a VCF file and produces a recombination map for mitochondrial DNA.
Because we assume no recombination, the recombination rate for all sites is r = 0.
This file takes three arguments:
*	`ARG1` = The input VCF file.
*	`ARG2` = The output recombination map map file.
*	`ARG3` = The output strand file.

### Imputation of missing variants on in silico microarrays (1,000 Genomes Project phase 3 data)
Files used in the imputation pipeline for 1000 Genomes Project data (in silico microarrays)
This script is used to iteratively submit one job per in silico microarray to the RAIJIN cluster at NCI (National Computational Infrastructure) in Canberra, Australia.
This script takes one argument, which is the reference panel being used from which missing variants will be imputed.

BASH file: [MassDeploy_ThousandGenomes_imputeFrom_RefPan.sh](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/BASH/DATA_PROCESSING/ThousandGenomes_Imputation/MassDeploy_ThousandGenomes_imputeFrom_RefPan.sh)
This file is the parent script.

USAGE:
```
$ sh MassDeploy_ThousandGenomes_imputeFrom_RefPan.sh <REFERENCE_PANEL>
```
Where `REFERENCE_PANEL` is a variable. You need to note the reference panel that will be used (exactly as it is in the reference panel file in XXX).
You may modify the MCMC, BURNIN, and KHAP settings as required for your particular analyses. However, these are set to these defaults for the 'recommended settings':
*	`MCMC = 1`
*	`BURNIN = 0`
*	`KHAP = 500`

SCRIPTS USED WITHIN:
BASH:
*	[MassDeploy_ThousandGenomes_imputeFrom_RefPan.sh](MassDeploy_ThousandGenomes_imputeFrom_RefPan.sh)
*	[ThousandGenomes_imputeFrom_RefPan.sh](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/BASH/ThousandGenomes_imputeFrom_RefPan.sh)

PYTHON:
*	[pickFirstAlt.py](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/PYTHON/pickFirstAlt.py)
*	[vcf2fasta_rCRS.py](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/PYTHON/vcf2fasta_rCRS.py)
*	[fasta2vcf_mtDNA.py](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/PYTHON/fasta2vcf_mtDNA.py)


R:
*	[assign_sex_label.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/assign_sex_label.R)
*	[FixSamplesFile_raijin.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/FixSamplesFile_raijin.R)
*	[plink_sites_map.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/plink_sites_map.R)
*	[MCC_Genotypes.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/MCC_Genotypes.R)

LISTS USED WITHIN:
*	[b37_platforms.txt](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/INFORMATION_LISTS/b37_platforms.txt)

MODULES AND APPLICATIONS CALLED UPON BY DAUGHTER SCRIPTS:
*	python v2.7.11
*	R v3.4.3
*	bcftools v1.8
*	plink v1.9
*	impute2 v2.3.2
*	vt v0.57721
*	java/jdk v1.8.0_60
*	HaploGrep v2.1.19


DETAILS:
[ThousandGenomes_imputeFrom_RefPan.sh](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/BASH/ThousandGenomes_imputeFrom_RefPan.sh)
This script is used to submit jobs to the RAIJIN cluster.
Using the qsub command, it takes 6 arguments:
*	`REFpanel` – the reference panel from which missing variants will be imputed.
*	`MtPlatforms` – the name of the in silico microarray whose missing variants will be imputed.
*	`mcmc` – length of the Markov chain Monte Carlo
*	`burn` – burn-in length of the Markov chain Monte Carlo
*	`khap` – The number of reference haplotypes to be used
*	`ne` – the effective population size.
The REFpanel option is set by the argument passed to the parent BASH script.
The MtPlatforms option is what is being iteratively changed by the parent script.
All other options remain constant and are set in the parent BASH script.

[pickFirstAlt.py](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/PYTHON/pickFirstAlt.py)
[ASK BRIAN OR RUSSELL]

[vcf2fasta_rCRS.py](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/PYTHON/vcf2fasta_rCRS.py)
This script takes in a VCF file for mtDNA and converts it to a FASTA formatted multiple sequence alignment file.
Sites are numbered according to the revised Cambridge Reference Sequence.
There are four options:
*	`--vcf_file` is the input VCF file. VCF files can be haploid or diploid, but if they are diploid this programme will select the genotype on the left of / or |. So if any heteroplasmy if observed, this is not the script for you.
*	`--out_file` is the output FASTA file. If no file is specified it will be output to the same location as the VCF file, albeit with .vcf / .vcf.gz replaced with .fasta
*	`--include-rCRS` includes the revised Cambridge Reference sequence in the FASTA file.
*	`--verbose` is verbose mode.
What this script does in essence is duplicate the rCRS n times (n being the number of sequences in the VCF file).
Then, for each sequence in the VCF file, for each site in the VCF file, if the genotype that sequence has at that site matches the reference allele, pass. If is the alternative allele, replace the reference nucleotide with the alternative nucleotide.

[fasta2vcf_mtDNA.py](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/PYTHON/fasta2vcf_mtDNA.py)
This script takes FASTA formatted multiple sequence alignments of mtDNA and converts them to a VCF format.
Reference alleles are set to the revised Cambridge Reference Sequence.
There are 9 options.
*	`--infile` is the input FASTA file. Preferably curated to the rCRS numbering system (s=16569)
*	`--outfile` OUTFILE is the output VCF file. Includes all sites, even invariants. If no file is specified it will be output to the same location as the VCF file, albeit with .fasta replaced with .vcf
*	`--gap2missing` turn gaps (-) to missing (N) character states.
*	`--diploid` creates diploid VCF file instead of haploid (ie 0|0 instead of 0).
*	`--ID tags` the ID column as MT<POS>.
*	`--quality` tags for QUAL column (default: 999).
*	`--filt tags` for FILTER column (default: PASS).
*	`--verbose` turns on verbose mode.
*	`--add_alt` forces the VCF to always have an alternative allele.

[assign_sex_label.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/assign_sex_label.R)
This script takes a .SAMPLES file and appends a column with "M" to denote male sex label.
There are two user inputs required:
*	`ARG1` = input .samples file.
*	`ARG2` = output .samples file.

[FixSamplesFile_raijin.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/FixSamplesFile_raijin.R)
This script files the .samples file by changing the 4th element of the 1st row to a "D". This makes it compatible with IMPUTE2.
*	`ARG1` = inpute .samples file.

[plink_sites_map.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/plink_sites_map.R)
This script fixes the site map .map file generated by plink. It creates a new first column with "MT" at all sites.
*	`ARG1` = input .map file

[MCC_Genotypes.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/MCC_Genotypes.R)
This script calculates the Matthew's correlation coefficient for imputed + genotyped VCF file and the genotyped only file.
It takes 5 arguments.
*	`ARG1` = The VCF file for the truthset (whole-molecule resequencing of the 1000 Genomes Phase 3 mtDNA data set.
*	`ARG2` = The VCF file for the genotyped only dataset.
*	`ARG3` = The VCF file for the genotyped + imputed dataset
*	`ARG4` = The INFO file that resulted from the IMPUTE2 process.
*	`ARG5` = The output file.

### Imputation pipeline tested on the Alzheimer's Disease Neuroimaging Initiative (ADNI) 1 and 3 datasets.
Files used in the imputation pipeline for Alzheimer's Disease Neuroimaging Initiative datasets.
This script takes one argument, which is the reference panel being used from which missing variants will be imputed.


BASH files:
*	[Impute_ADNI_redo.sh](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/BASH/DATA_PROCESSING/ADNI_imputation/Impute_ADNI_redo.sh) # For the 258 genotyped samples that appear in both ADNI1 and ADNI3
*	[Impute_ADNI_redo_noReSeq.sh](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/BASH/DATA_PROCESSING/ADNI_imputation/Impute_ADNI_redo_noReSeq.sh) # For the 499 genotyped samples that appear only in ADNI1
*	[Impute_ADNI_12GO.sh](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/BASH/DATA_PROCESSING/ADNI_imputation/Impute_ADNI_12GO.sh) # For the 1199 samples that appear in ADNI GO
These file are the parent scripts. They all do the same thing, but for different data sets (see below).

USAGE:
```
$ Impute_ADNI_redo.sh <REFERENCE_PANEL_VERSION>
```
Where <REFERENCE_PANEL_VERSION> is a variable.
<REFERENCE_PANEL_VERSION> is reference panel version and the minor allele frequency (ie ReferencePanel_v1_0.01 where ReferencePanel_v1 is the reference panel is 0.01 is the minor allele frequency).

SCRIPTS USED WITHIN:
BASH:
N/A

PYTHON:
[fix_vcf_names.py](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/PYTHON/fix_vcf_names.py)

R:
[assign_sex_label.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/DATA_PROCESSING/assign_sex_label.R)
[FixSamplesFile_raijin.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/FixSamplesFile_raijin.R)
[HiMC_haplogroup_assignment.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/HiMC_haplogroup_assignment.R)
[MCC_Genotypes.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/MCC_Genotypes.R)

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

[fix_vcf_names.py](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/PYTHON/fix_vcf_names.py)
This file takes in a VCF file and fixes the names of the samples in the header.
This is done so that the same ADNI1 and ADNI3 samples can be directly and easily compared later.
This script takes 4 arguments:
*	`--vcf_file` is the input VCF file.
*	`--out_file` is the output VCF file with the sample names fixed.
*	`--csv_file` is the csv file containing file names and comparisons.
*	`--verbose` turns on verbose mode.

<R SCRIPTS>

[assign_sex_label.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/DATA_PROCESSING/assign_sex_label.R)
This script takes in a .sample file and adds a new column with an "M" sex label.
This script takes one argument:
*	`ARG1` = The input .sample file.

[FixSamplesFile_raijin.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/FixSamplesFile_raijin.R)
This script files the .samples file by changing the 4th element of the 1st row to a "D". This makes it compatible with IMPUTE2.
*	`ARG1` = inpute .samples file.

[HiMC_haplogroup_assignment.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/HiMC_haplogroup_assignment.R)
This script takes in .ped and .map files from the imputed data set and assigns haplogroups according PhyloTree 17 via the HiMC package.
This script takes 3 arguments:
*	`ARG1` = Input .ped file.
*	`ARG2` = Input .map file.
*	`ARG3` = Output .csv file.

[MCC_Genotypes.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/MCC_Genotypes.R)
This script calculates the Matthew's correlation coefficient for imputed + genotyped VCF file and the genotyped only file.
It takes 5 arguments.
*	`ARG1` = The VCF file for the truthset (whole-molecule resequencing of the 1000 Genomes Phase 3 mtDNA data set.
*	`ARG2` = The VCF file for the genotyped only dataset.
*	`ARG3` = The VCF file for the genotyped + imputed dataset
*	`ARG4` = The INFO file that resulted from the IMPUTE2 process.
*	`ARG5` = The output file.

### Analyses for the pre-SnakeMake pipeline
The files described in here were used for analysis of the resulting imputation pipeline experiments.
Descriptions for each file found below.

[calculate_95CI.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/ANALYSIS/calculate_95CI.R)
This script contains the formula/function to calculate the 95% CI intervals.
While we used it on certain metrics stated in the paper, this can be used on any numeric values.

[check_haplogroup_concordance.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/ANALYSIS/HiMC/check_haplogroup_concordance.R)
This script iteratively reads through all the experiments (KHAP, MAF, MCMC).
First the haplogroups are assiged for each sample in the truth set, then haplogroups are assigned using HiMC.
Then for each experiment and each chip, haplogroups are assigned.
The same happens for macro-haplogroups, where for non-L haplogroups the haplogroup assigned is truncated to the major haplogroup (ie H), then compared to the truth set.
For L haplogroups, the haplogroup is truncated to the first subclade (ie L3).

[Cleanup_concordance_tables_HiMC.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/ANALYSIS/HiMC/Cleanup_concordance_tables_HiMC.R)
Following on from the check_haplogroup_concordance.R script, this script calculates haplogroup and macrohaplogroup concordance.
Summary statistics are then calculated.
Each individual experiment is saved, then all are combined into one big data table.
Plots are able to be produced from this script, however the HiMC_1kGP_plots.R described later produce the most up to date plots.
This also script performs statistical tests on the MCC datasets.
It performs linear mixed model evaluations using ANOVA estimated marginal means.
The results from these evaluations are then output to CSV format.

[concordance_tables_ADNI.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/ANALYSIS/HiMC/concordance_tables_ADNI.R)
Similar to check_haplogroup_concordance.R, this script assigns haplogroups for the ADNI imputed dataset.
The samples from the Whole Genome re-Sequenced ADNI3 dataset have their haplogroups assigned and used as the truthset.
The samples from the genotyped-only ADNI1 dataset have their haplogroups assigned.
Haplogroup (and macrohaplogroup) assignment is then compared and the concordance calculated.

[MCC_emmeans.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/ANALYSIS/MCC/MCC_emmeans.R)
This script performs statistical tests on the MCC datasets.
It performs linear mixed model evaluations using ANOVA estimated marginal means.
The results from these evaluations are then output to CSV format.

[Cleanup_concordance_tables_MCC.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/ANALYSIS/MCC/Cleanup_concordance_tables_MCC.R)
This script is similar in function to Cleanup_concordance_tables_HiMC.R, however it performs Matthew's correlation coefficient calculations instead of haplogroup concordance.
This is conducted iteratively through all the experiments (KHAP, MAF, MCMC).
Summary statistics are calculated and saved.

[HiMC_1kGP_plots.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/Plotting/HiMC_1kGP_plots.R)
This script produces box and whisker plots of the results of the HiMC experiments.

[MCC_concordance_tables.R](https://github.com/sjfandrews/MitoImputePrep/blob/master/scripts/R/Plotting/MCC_concordance_tables.R)
his script produces box and whisker plots of the results of the MCC experiments.

# Snakemake Implementation
Snakemake pipeline requires further validation


## Thousand Genomes Validation
To validate the MitoImpute pipeline, we evaluate the imputation perforamnce using mitochondrial sequences from Thousand Genomes. Using [Strand Files](http://www.well.ox.ac.uk/~wrayner/strand/), which provide the strand orientation and position of variants of the most common genotyping platform on build 37, a list of mtSNPs for each platform is created. For each platform, the mtSNPs are then extracted from the Thousand Genome mitochondrial Sequences. Using these 'Typed SNPs' mitochondrial SNPs from the refernce panel are then imputed. Impuation quality is then evaluated using  [MitoImputeShiny](https://github.com/sjfandrews/MitoImputeShiny). This pipeline consistes for two Snakemake files.


**PlatformStrandFiles.smk:** Script to downloand the [Strand Files](http://www.well.ox.ac.uk/~wrayner/strand/). The steps in the pipeline include:
1. Download Strand files
2. Extract list of mitochondrial SNPs on each genotyping platforms
3. Write a summary file with the number of mtSNPs on each platform and a list of platforms with mtSNPs

**ThousandGenomes.smk:** Pipline to download and process mitochondrial sequences from Thousand Genomes and creat 'in silico' microarrys based on the downloaded strand files. The steps in the pipeline include:
1. Download Thousand Genomes mitochondrial sequences
2. Normalize the Thousand Genomes vcf
   - split multialleic sites
   - remove indels and mnps
   - join multialleic sites
   - fill info field on vcf
3. Decompose multiallelic variants - other alternate allels are set to missing
   - Pick first alternate allele in multiallelic sites (ie most frequent alt allele)
4. Convert Thousand Genomes VCF to plink (.map/.ped)
5. Assign Male sex to Thousand Genomes samples
   - Extract a list of sample names from vcf
   - Add additional column assigning male sex
6. Extract Platform specific mtSNPs from Thousand Genomes Sequences
7. Convert Pltaform files from VCF to:
   - gen/sample format used by IMPUTE2
   - Plink format
8. Runs the chromosome X IMPUTE2 imputation protocol.
9. Fixes chromosome label on the IMPUTE2 output
10. Converts the Imputed files to:
    - Plink format
    - vcf format
11. Generates a html rmarkdown report


## MitoImpute
**mtImpute.smk:** Pipeline for imputing untyped mitochondrial SNPs in a study dataset from a custom mitochondrial reference panel. Used by [MitoImpute](https://github.com/sjfandrews/MitoImpute).

A custom reference panel for imputation can be found in the ```MitoImpute/DerivedData/ReferencePanel/``` directory. The key files consist of:
1. -h: A file of known haplotypes ```(ReferencePanel.hap.gz)```.
2. -l: Legend file(s) with information about the SNPs in the -h file ```(ReferencePanel.legend.gz)```
3. -m: A fine-scale recombination map for the region to be analyzed ```(MtMap.txt)```

setting REFDATA in the ```mtImpute_config.yaml``` file to ```path/to/MitoImpute/DerivedData/ReferencePanel``` will automaticlay call these files.

## Getting Started
### Installation
Be sure to download and install the latest versions of the following software packages:
1. [Python 3](https://www.python.org/downloads/)
2. [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
3. [R](https://cran.r-project.org/)
4. [PLINK](https://www.cog-genomics.org/plink2)
5. [BCFtools](http://samtools.github.io/bcftools/howtos/install.html)
6. [vt](https://genome.sph.umich.edu/wiki/Vt#Installation)
7. [Impute2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#download)

The following R packages are also required:
1. [tidyverse](https://www.tidyverse.org/packages/)
2. [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html)
3. [Hi-MC](https://github.com/vserch/himc)
4. [ggforce](https://github.com/thomasp85/ggforce)

Note that the development versions of ggforce (required for plotting alluvial diagrams) and Hi-MC (required for mitochondrial haplogroup assignment) are required. These packages can be isntalled directly from github using devtools (see their respective pages).

The following Python modules are required:
1. [pysam](http://pysam.readthedocs.io/en/latest/api.html)

Once all the prerequiste software is isntalled, MitoImputePrep can be installed on a git-enabled machine by typeing:

```bash
git clone https://github.com/sjfandrews/MitoImputePrep
```

### Usage Overview
#### Reference Panel
To construct the reference panel run the following code:

```bash
snakemake -s ReferencePanel.smk
```

Options for the snakemake file are set in the corresponding config file ```ReferencePanel_config.yaml``` file. The avaliable options are:

```bash
FileName: 'name of fasta file'
DataIn: 'path/to/input/directory'
DataOut: 'path/to/output/directory'
```

#### Thousand Genomes Validation
To run the Thousand Genomes validation pipeline, run the following code:

```bash
snakemake -s PlatformStrandFiles.smk
snakemake -s ThousandGenomes.smk
```

By default, the reference panel is set to the example reference panel.

#### MitoImpute
To impute mitochondrial SNPs in a study dataset, run the following code:

```bash
snakemake -s mtImpute.smk
```

Options for the snakemake file are set in the corresponding config file ```mtImpute_config.yaml``` file. The avaliable options are:

```bash
SAMPLE: 'name of binary plink file'
DATAIN: 'path/to/plink/file'
DATAOUT: 'path/to/output/file'
REFDATA: 'path/to/reference/panel'
```

The default options are for the example dataset.
