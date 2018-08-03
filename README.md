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
