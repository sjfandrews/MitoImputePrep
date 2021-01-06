
========================
= REFERENCE ALIGNMENTS =
========================

The reference alignments included in this repository are:
*	'McInerney_Master_Alignment_July18_2018.fasta.gz'
*	'hsapiensCRS7k.fasta.gz'

'McInerney_Master_Alignment_July18_2018.fasta.gz' is the novel reference alignment constructed on in 2018 from the sequences downloaded on the 18th of July, 2018. It contains 44,299 aligned complete mitochondrial DNA sequences. These sequences are all 16,569 DNA nucleotide states long (to match the numbering conventions of the revised Cambridge Reference Sequence - Andrews et al., 1999). From this alignment the Reference Panels were filtered down to 36,960 sequences and filtered to thresholds detailed in the manuscript.

'hsapiensCRS7k.fasta.gz' is the previous reference alignment constructed in 2011 by Simon Easteal and Lars Jermiin. It contains 7,747 aligned complete mitochondrial DNA sequences. These sequences are all 16,569 DNA nucleotide states long (to match the numbering conventions of the revised Cambridge Reference Sequence - Andrews et al., 1999). This curated alignment was used to align the sequences downloaded on the 18th of July 2018. Novel sequences were aligned in batches of 2,500 sequences. Any gaps forced into 'hsapiensCRS7k.fasta.gz' were removed via the custom python script 'MitoImputePrep/scripts/PYTHON/curate_rCRS_FASTA.py'. The 'McInerney_Master_Alignment_July18_2018.fasta.gz' reference alignment can be recreated by downloading the complete mitochondrial DNA sequences available on GenBank on the 18th of July 2018 and aligning them to 'hsapiensCRS7k.fasta.gz' as detailed in the manuscript.

Both alignments are contained in the directory 'MitoImpute/resources/alignments/'.

Both alignments are compressed in the .gz format to allow upload to GitHub. Many programs designed to handle FASTA formatted multiple sequence alignments can parse .gz. However, these files can be decompressed in a BASH terminal using the following command:
	$ tar -zxvf McInerney_Master_Alignment_July18_2018.fasta.gz

====================
= REFERENCE PANELS =
====================

The reference panels included in this repository are:
*	ReferencePanel_v1_0.01
*	ReferencePanel_v1_0.005
*	ReferencePanel_v1_0.001

Each of these corresponds to a filtering of sites to a minor allele frequency of 1%, 0.5%, and 0.01%, respectively. All reference panels contain variant information for the 36,960 sequences from the 'McInerney_Master_Alignment_July18_2018.fasta.gz' reference alignment.

All references panels can be found in the directory: 'MitoImpute/resources/'. Each specific reference panel further has its own subdirectory. Within these subdirectories, additional files can be found, such as *.gen.gz, .*hap.gz, *.legend.gz, *.vcf.gz, *_sampleList_sex.txt, and *.ped. These files will be necessary for using a reference panel for genotype imputation. As these files were created using BCFtools and PLINK, please refer to https://samtools.github.io/bcftools/bcftools.html and https://www.cog-genomics.org/plink2, respectively for additional information.









