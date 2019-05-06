'''Snakefile for downloading microarray strand datasets'''
# snakemake -s PlatformStrandFiles.smk

import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import re
import requests as rq

shell.executable("/bin/bash")

wrayner = {'src': rq.get('https://www.well.ox.ac.uk/~wrayner/strand/').text,
           're': '(?<=href=").+?b37(?=-strand.zip)'}
wrayner['noAffy'] = re.split('Affymetrix data files', wrayner['src'])[0]
wrayner['platforms'] = re.findall(wrayner['re'], wrayner['noAffy'])

HTTP = HTTPRemoteProvider()

rule all:
    input:
        expand('data/platforms/{platform}/{platform}_MT_snps.txt',
               platform = wrayner['platforms']),
        'data/platforms/Nsnps_Mt_platforms.txt',
        'data/platforms/Mt_platforms.txt'


'''
Download, unzip and copy non-Affy strands from Will Rayner's website at Oxford.
Download the zips, extract and copy the file ending in ".strand" to
"platform.strand" in data/platforms/[platform name] from b37_strandfiles.txt
'''
rule StrandFiles:
    input:
        HTTP.remote("https://www.well.ox.ac.uk/~wrayner/strand/{platform}-strand.zip", allow_redirects=True)
    output:
        'data/platforms/{platform}/platform.strand'
    params:
        directory = 'data/platforms/{platform}'
    shell:
        "unzip {input} -d {params.directory} *.strand; "
        "mv {params.directory}/*.strand {params.directory}/platform.strand"

'''
Generate files with the format {CHR}\t{POS} and no headers containing the MtSNPs
for each platform.
'''
rule StrandFilesMT:
    input:
        script = 'scripts/R/StrandFiles_ExtractMTsnps.R',
        strands = 'data/platforms/{platform}/platform.strand'
    output:
        out = 'data/platforms/{platform}/{platform}_MT_snps.txt'
    shell:
        'Rscript {input.script} {input.strands} {output.out}'

'''
Summarize which platforms have MtSNPs and how many snps on each platform.
'''
rule strandSummary:
    input:
        StrandFiles = expand('data/platforms/{platform}/platform.strand',
                             platform = wrayner['platforms']),
        script = "scripts/R/StrandFiles_ExtractMTSummary.R",
    output:
        sum_nSNPS = 'data/platforms/Nsnps_Mt_platforms.txt',
        sum_MTplats = 'data/platforms/Mt_platforms.txt'
    params:
        directory = 'data/platforms'
    shell:
        'Rscript {input.script} {params.directory} {output.sum_nSNPS} {output.sum_MTplats}'
