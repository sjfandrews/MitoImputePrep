'''Snakefile for Construction of Reference Panel 0.1'''
# snakemake -s ReferencePanel.smk
# snakemake -s ReferencePanel.smk --dag | dot -Tsvg > dag_ReferencePanel.svg

configfile: "ReferencePanel_config.yaml"
DATAIN = config['DataIn']
DATAOUT = config['DataOut']
FILENAME = config['FileName']

BPLINK = ["bed", "bim", "fam"]

rule all:
    input:
        expand(DATAOUT + "/ReferencePanel.{ext}", ext = ['hap.gz', 'legend.gz']),
        expand(DATAOUT + "/ReferencePanel.{ext}", ext = ['ped', 'map']),
        expand(DATAOUT + "/ReferencePanel.{ext}", ext = ['gen.gz', 'samples']),
        expand(DATAOUT + "/MtMap.txt"),
        expand(DATAOUT + "/MtStrand.txt"),
        expand(DATAOUT + "/ReferenceSNPs.txt")

## 1. Run the ambiguous2missing.py script to change ambiguous character states to missing data:
rule ambiguous2missing:
    input:
        in_script = "scripts/PYTHON/ambiguous2missing.py",
        in_fasta = DATAIN + "/{Reference}.fasta"
    output: temp(DATAOUT + "/{Reference}_ambig2missing.fasta")
    shell:
        'python {input.in_script} -i {input.in_fasta} -o {output.out} -v'

# 2. Identify samples with highQuality sequences
rule LowQualitySequences:
    input:
        in_script = "scripts/R/removeLowQuality_cmdline.R",
        in_fasta = DATAOUT + "/{Reference}_ambig2missing.fasta"
    output: DATAOUT + "/ReferencePanel_highQualitySequences.txt"
    shell:
        'Rscript {input.in_script} {input.in_fasta} {output}'

## 3a. Run the fasta2vcf_mtDNA.py script
rule fasta2vcf:
    input:
        in_script = "scripts/PYTHON/fasta2vcf_mtDNA.py",
        in_fasta = DATAOUT + "/{Reference}_ambig2missing.fasta"
    output:
        out_vcf = temp(DATAOUT + "/{Reference}_ambig2missing.vcf.gz")
    shell:
        'python {input.in_script} -i {input.in_fasta} -o {output.out_vcf} -v'

# 3b. Pass the resulting VCF through BCFTOOLS to make sure it conforms to all standards and index it
rule VcfCheck:
    input: expand(DATAOUT + "/{Reference}_ambig2missing.vcf.gz", Reference=FILENAME),
    output: temp(DATAOUT + "/Reference_panal.vcf.gz"),
    shell:
        'bcftools view -Oz -o {output} {input}'
        'bcftools index {output}'

## 4. Remove low quality sequences from VCF
rule RemoveLowQuality:
    input:
        qual = DATAOUT + "/ReferencePanel_highQualitySequences.txt",
        in_vcf = DATAOUT + "/Reference_panal.vcf.gz",
    output: temp(DATAOUT + "/ReferencePanel_highQual.vcf.gz"),
    shell:
        '''
bcftools view --force-samples -S {input.qual} -Oz -o {output} {input.in_vcf}
bcftools index {output}
'''

## 5. Apply filtration criteria
rule SiteFiltration:
    input: DATAOUT + "/ReferencePanel_highQual.vcf.gz",
    output: DATAOUT + "/ReferencePanel.vcf.gz",
    shell:
        '''
vt decompose {input} | bcftools +fill-tags | \
  bcftools view -i \'ALT!="-" \' | \
  bcftools view -q 0.01 -Q 0.99 -Oz -o {output}
'bcftools index {output}')
'''

## 6a. Extract sample names from Reference Panel
rule RefSampleNames:
    input:
        in_vcf = DATAOUT + "/ReferencePanel.vcf.gz",
    output:
        out_samples = DATAOUT + "/RefSampleList.txt",
    shell:
        'bcftools query -l {input.in_vcf} > {output.out_samples}'

## 6b. Assign M sex label to reference Samples
rule RefSampleSex:
    input: DATAOUT + "/RefSampleList.txt"
    output: DATAOUT + "/RefSampleList_sex.txt"
    shell:
        '''awk -F "\t" '$2 = "M"' {input} {output}'''

## 7. Convert to Oxford format
rule Oxford:
    input:
        in_vcf = DATAOUT + "/ReferencePanel.vcf.gz",
        in_sex = DATAOUT + "/RefSampleList_sex.txt"
    output:
        expand(DATAOUT + "/ReferencePanel.{ext}", ext = ['hap.gz', 'legend.gz'])
    params:
        out = DATAOUT + "/ReferencePanel"
    shell:
        'bcftools convert --haplegendsample {params.out} {input.in_vcf} --sex {input.in_sex}'

## 8. Generate .ped and .map files
rule Plink:
    input:
        in_vcf = DATAOUT + "/ReferencePanel.vcf.gz",
    output:
        expand(DATAOUT + "/ReferencePanel.{ext}", ext = ['ped', 'map'])
    params:
        out = DATAOUT + "/ReferencePanel"
    shell:
        'plink --vcf {input.in_vcf} --recode --double-id --keep-allele-order --out {params.out}'

## 9. Generate .gen and .sample files
rule GenSample:
    input:
        in_vcf = DATAOUT + "/ReferencePanel.vcf.gz",
        in_sex = DATAOUT + "/RefSampleList_sex.txt"
    output:
        expand(DATAOUT + "/ReferencePanel.{ext}", ext = ['gen.gz', 'samples'])
    params:
        out = DATAOUT + "/ReferencePanel"
    shell:
        'bcftools convert --gensample {params.out} {input.in_vcf} --sex {input.in_sex}'

## 10. Construct .map file for IMPUTE2
rule MakeMapFile:
    input:
        in_script = "scripts/R/mt_recombination_map.R",
        in_vcf = DATAOUT + "/ReferencePanel.vcf.gz"
    output:
        out_map = DATAOUT + "/MtMap.txt",
        out_strand = DATAOUT + "/MtStrand.txt"
    shell:
        'Rscript {input.in_script} {input.in_vcf} {output.out_map} {output.out_strand}'

## 11. Output list of ReferenceSNPs
rule RefSNPs:
    input:
        in_vcf = DATAOUT + "/ReferencePanel.vcf.gz",
    output:
        out =  DATAOUT + "/ReferenceSNPs.txt"
    shell:
        'bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n" {input.in_vcf} -o {output.out}'
