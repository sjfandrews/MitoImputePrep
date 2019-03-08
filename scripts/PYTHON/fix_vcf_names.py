#!/usr/bin/env python

import time
import gzip
import argparse
import os
import sys
import math
import numpy
import csv
import fcntl, termios, struct

from tqdm import *

samples_csv = "/Users/u5015730/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH.csv"
WGS_vcf = "/Volumes/TimMcInerney/MitoImpute/data/ADNI_raijin/ADNI_reseq/adni_mito_genomes_180214_n258_subset.vcf.gz"
relabelled_vcf = "/Volumes/TimMcInerney/MitoImpute/data/ADNI_raijin/ADNI_reseq/adni_mito_genomes_180214_n258_subset_relabelled.vcf"

csv_lines = []

with open(samples_csv, 'r') as sc:
    for line in sc:
        if line.startswith("SAMPLE_NUMBER"):
            line = line.strip("\n")
            line = line.split(",")
            csv_header = line
        else:
            line = line.strip("\n")
            line = line.split(",")
            csv_lines.append(line)

seqs = {}
for i in range(len(csv_lines)):
    seqs[csv_lines[i][0]] = csv_lines[i][5]
    
meta_lines = []
var_lines = []

with gzip.open(WGS_vcf, 'r') as vcf:
    for line in vcf:
        if line.startswith("##"):
            meta_lines.append(line)
        elif line.startswith("#C"):
            line = line.strip("\n")
            line = line.split("\t")
            headers = line
        else:
            var_lines.append(line)

info_headers = headers[:9]
seq_headers = headers[9:]
            
for i in range(len(seq_headers)):
    seq_headers[i] = seqs[seq_headers[i]]
    
new_headers = info_headers + seq_headers
    
with open(relabelled_vcf, 'w') as o:
    for i in meta_lines:
        o.write(i)
    o.write("\t".join(new_headers) + "\n")
    for i in var_lines:
        o.write(i)