'''Snakefile for Cleaning Thousand Genomes'''
# snakemake -s ThousandGenomes.smk
# snakemake -s ThousandGenomes.smk --dag | dot -Tsvg > dag_ThousandGenomes.svg
# snakejob -s ThousandGenomes.smk -j 8000 --max-jobs-per-second 1 --keep-going

import os
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

shell.prefix('module load plink/1.90 R/3.4.3 curl/7.61.0 vt bcftools impute2; ')

#with open('data/platforms/Mt_platforms.txt', "r") as f:
#    MtPlatforms = [x.rstrip() for x in f]
MtPlatforms = ['GSA-24v1-0_A2-b37', 'Human610-Quadv1_B-b37', 'NeuroX_15036164_A-b37']

#REFDATA = "example/ReferencePanel"
REFDATA = "DerivedData/ReferencePanel"

## Parameters for IMPUTE2
configfile: 'ThousandGenomes_config.yaml'
MCMC = [str(x).split(".") for x in config['MCMC']]
MCMC = {'.'.join(x): {'iter': x[0], 'burnin': x[1]} for x in MCMC}
# 1KGP MCMC params {'1.0': {'burnin': '0', 'iter': '1'}, ...}
KHAP = config['KHAP']

FTP = FTPRemoteProvider()
RWD = os.getcwd()

TG_derived = 'DerivedData/ThousandGenomes' #Directory for derived 1kgp files
rule all:
    input:
        expand(TG_derived + '/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.{ext}',
               ext=['gen.gz', 'samples', 'ped', 'map'],
               MtPlatforms=MtPlatforms),
        expand(TG_derived + '/chrMT_1kg_norm_firstAlt.{ext}',
               ext = ['ped', 'map', '.vcf.gz']),
        expand(TG_derived + '/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/' +
               'chrMT_1kg_{MtPlatforms}_imputed{tail}',
               MtPlatforms=MtPlatforms, KHAP=KHAP, MCMC=MCMC,
               tail=['','_samples','.ped', '.map', '.vcf', '_QC.html'])

# 1. Pull down 1000 genomes mitochondrial vcf file from ftp
TG_releasedir = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
TG_vcf = 'ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz'

rule Get1kgMT_vcf:
    input:
        vcf = FTP.remote(TG_releasedir + TG_vcf, keep_local=True),
        tbi = FTP.remote(TG_releasedir + TG_vcf + '.tbi', keep_local=True)
    output:
        vcf = 'data/ThousandGenomes/' + TG_vcf,
        tbi = 'data/ThousandGenomes/' + TG_vcf + '.tbi'
    shell:
        "mv {input.vcf} {output.vcf}; "
        "mv {input.tbi} {output.tbi}"

rule NormaliseVcf:
    input:
        vcf = rules.Get1kgMT_vcf.output.vcf,
        fasta = "data/ReferencePanel/rCRS.fasta"
    output:
        vcf = "DerivedData/ThousandGenomes/chrMT_1kg_norm_firstAlt.vcf.gz"
    shell:
        '''
vt {input.vcf} | bcftools norm -f {input.fasta} | bcftools view -v snps | \
   bcftools norm -d all | bcftools +fill-tags -o {output.vcf} -Oz
'''

rule wgs_vcf2Plink:
    input:
        vcf = "DerivedData/ThousandGenomes/chrMT_1kg_norm_firstAlt.vcf.gz",
    output:
        expand("DerivedData/ThousandGenomes/chrMT_1kg_norm_firstAlt.{ext}", ext = ['ped', 'map'])
    params:
        out = "DerivedData/ThousandGenomes/chrMT_1kg_norm_firstAlt"
    shell:
        'plink --vcf {input.vcf} --recode --double-id --keep-allele-order --out {params.out}'

## 6a. Extract sample names from Reference Panel
rule SampleNames1kg:
    input:
        "DerivedData/ThousandGenomes/chrMT_1kg_norm_firstAlt.vcf.gz",
    output:
        "DerivedData/ThousandGenomes/SampleList1kg.txt",
    shell:
        'bcftools query -l {input} > {output}'

## 6b. Assign M sex label to  Samples
rule SampleSex1kg:
    input:
        in_samples = "DerivedData/ThousandGenomes/SampleList1kg.txt",
        in_script = "scripts/R/assign_sex_label.R"
    output:
        "DerivedData/ThousandGenomes/SampleList1kg_sex.txt",
    shell:
        'Rscript {input.in_script} {input.in_samples} {output}'

rule ExtractPlatformMTsnps:
    input:
        MTSnps = 'data/platforms/{MtPlatforms}/{MtPlatforms}_MT_snps.txt',
        vcf_1kg = "DerivedData/ThousandGenomes/chrMT_1kg_norm_firstAlt.vcf.gz"
    output:
        vcf = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.vcf.gz"
    shell:
        "bcftools view -R {input.MTSnps} {input.vcf_1kg} -Oz -o {output.vcf}; "
        "bcftools index {output.vcf}"

rule vcf2gensample:
    input:
        vcf = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.vcf.gz",
        sex = "DerivedData/ThousandGenomes/SampleList1kg_sex.txt",
        script = "scripts/R/FixSamplesFile.R"
    output:
        expand("DerivedData/ThousandGenomes/{{MtPlatforms}}/chrMT_1kg_{{MtPlatforms}}.{ext}", ext = ['gen.gz', 'samples']),
    params:
        out = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}"
    shell:
        'bcftools convert --gensample {params.out} {input.vcf} --sex {input.sex}; '
        'Rscript {input.script} {params.out}.samples'

rule vcf2Plink:
    input:
        vcf = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.vcf.gz",
    output:
        expand("DerivedData/ThousandGenomes/{{MtPlatforms}}/chrMT_1kg_{{MtPlatforms}}.{ext}", ext = ['ped', 'map'])
    params:
        out = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}"
    shell:
        'plink --vcf {input.vcf} --recode --double-id --keep-allele-order --out {params.out}'

rule Impute2:
    input:
        m = expand('{RefData}/MtMap.txt', RefData=REFDATA),
        h = expand('{RefData}/ReferencePanel.hap.gz', RefData=REFDATA),
        l = expand('{RefData}/ReferencePanel.legend.gz', RefData=REFDATA),
        g = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.gen.gz",
        sample = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.samples",
    output:
        'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed',
        'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed_info',
        'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed_samples'
    params:
        out = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed',
        iter = lambda wildcards: MCMC[wildcards['MCMC']]['iter'],
        burnin = lambda wildcards: MCMC[wildcards['MCMC']]['burnin'],
        khap = '{KHAP}'
    shell:
        'impute2 -chrX -m {input.m} -h {input.h} -l {input.l} -g {input.g} \
        -sample_g {input.sample} -int 1 16569 -Ne 20000 -o {params.out} \
        -iter {params.iter} -burnin {params.burnin} -k_hap {params.khap}'


rule FixChromName:
    input:
        InFile = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed'
    output:
        OutFile = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed_ChromFixed'
    shell:"""
        awk '{{$1 = "26"; print}}' {input.InFile} > {output.OutFile}
    """

rule oxford2ped:
    input:
        gen = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed_ChromFixed',
        sample = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed_samples'
    output:
        expand("DerivedData/ThousandGenomes/{{MtPlatforms}}/KHAP_{{KHAP}}/MCMC_{{MCMC}}/chrMT_1kg_{{MtPlatforms}}_imputed.{ext}", ext = ['ped', 'map'])
    params:
        out = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed'
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --output-chr 26 --recode --out {params.out}'

rule oxford2vcf:
    input:
        gen = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed_ChromFixed',
        sample = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed_samples'
    output:
        "DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed.vcf"
    params:
        out = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed'
    shell:
        'plink --gen {input.gen} --sample {input.sample} --hard-call-threshold 0.49 \
        --keep-allele-order --output-chr 26 --recode vcf --out {params.out}'

rule Imputation_QC_Report:
    input:
        script = 'scripts/R/MT_imputation_QC.Rmd',
        wgs_map = 'DerivedData/ThousandGenomes/chrMT_1kg_norm_firstAlt.map',
        wgs_ped = 'DerivedData/ThousandGenomes/chrMT_1kg_norm_firstAlt.ped',
        wgs_vcf = 'DerivedData/ThousandGenomes/chrMT_1kg_norm_firstAlt.vcf.gz',
        typ_map = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.map",
        typ_ped = "DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.ped",
        typ_vcf = 'DerivedData/ThousandGenomes/{MtPlatforms}/chrMT_1kg_{MtPlatforms}.vcf.gz',
        imp_map = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed.map',
        imp_ped = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed.ped',
        imp_vcf = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed.vcf',
        imp_info = 'DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_imputed_info',
    output:
        "DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/chrMT_1kg_{MtPlatforms}_mtImputed_QC.html"
    params:
        rwd = RWD,
        output_dir = "DerivedData/ThousandGenomes/{MtPlatforms}/KHAP_{KHAP}/MCMC_{MCMC}/",
        info_cut = '0'

    shell:
        "R -e 'rmarkdown::render("
        """"{input.script}", output_file = "{output}", output_dir = "{params.output_dir}", \
params = list(rwd = "{params.rwd}", info.cut = "{params.info_cut}", \
wgs.map = "{input.wgs_map}", wgs.ped = "{input.wgs_ped}", wgs.vcf = "{input.wgs_vcf}", \
typ.map = "{input.typ_map}", typ.ped = "{input.typ_ped}", typ.vcf = "{input.typ_vcf}", \
imp.map = "{input.imp_map}", imp.ped = "{input.imp_ped}", imp.vcf= "{input.imp_vcf}", \
imp.info = "{input.imp_info}"))' --slave
        """
