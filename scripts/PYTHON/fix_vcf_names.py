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

def main():
    
    start_time = time.time()
    start_display = 'Process started at ' + time.strftime("%d/%m/%Y") + " " + time.strftime("%H:%M:%S")
    
    #DEFINE INPUT ARGUMENTS
    parser=argparse.ArgumentParser(description='Fix the names in ADNI VCF file so they match',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
      
    parser.add_argument('-i', '--vcf_file', dest='vcf_file', type=str, required=True, help='input vcf file')
    parser.add_argument('-o', '--out_file', dest='out_file', type=str, required=False, help='output vcf file')
    parser.add_argument('-c', '--csv_file', dest='csv_file', type=str, required=True, help='sample names file')
    parser.add_argument('-v', '--verbose', dest='verbose', action="store_true", required=False, help='turn on verbose mode')
      
    args=parser.parse_args()    
    
    samples_csv = args.csv_file
    WGS_vcf = args.vcf_file
    relabelled_vcf = args.out_file  
    verbose = args.verbose
    
    if relabelled_vcf is None:
        if WGS_vcf.endswith(".vcf.gz"):
            relabelled_vcf = WGS_vcf[:-7] + "_relabelled.vcf"
        else:
            relabelled_vcf = WGS_vcf[:-4] + "_relabelled.vcf"
    
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
    #
    stop_time = time.time()
    stop_display = 'Process completed at ' + time.strftime("%d/%m/%Y") + " " + time.strftime("%H:%M:%S")
    mem_used = memory_usage_resource()
    
    print
    print '*\t' + str(start_display)
    print '*\t' + str(stop_display)
    time_adj = time.time() - start_time
    if time_adj < 60:
        print('*\tProcess completed in %s seconds' % (round(time_adj, 2)))
    if time_adj >= 60 and time_adj < 3600:
        print('*\tProcess completed in %s minutes' % (round((time_adj / 60), 2)))
    if time_adj >= 3600 and time_adj < 86400:
        print('*\tProcess completed in %s hours' % (round(((time_adj / 60) / 60), 2)))
    if time_adj >= 86400:
        print('*\tProcess completed in %s days' % (round((((time_adj / 60) / 60) / 24), 2)))
    if mem_used < 1024:
        print '*\tProcess used %s MB of memory' % ('%.2f' % (mem_used))
    if mem_used >= 1024:
        print '*\tProcess used %s GB of memory' % ('%.2f' % (mem_used / 1024))
    
    if verbose:
        print
        print(' FINISHED ADDING INFERRED ANCESTRAL SEQUENCE '.center(int(terminal_size()[0]), '='))
        print        
        
if __name__=="__main__":
    main()      

'''
samples_csv = "/Users/TimMcInerney/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH.csv"
#WGS_vcf = "/Volumes/TimMcInerney/MitoImpute/data/ADNI_raijin/ADNI_reseq/adni_mito_genomes_180214_n258_subset.vcf.gz"
#relabelled_vcf = "/Volumes/TimMcInerney/MitoImpute/data/ADNI_raijin/ADNI_reseq/adni_mito_genomes_180214_n258_subset_relabelled.vcf"
WGS_vcf="/Volumes/TimMcInerney/MitoImpute/data/ADNI_REDO/WGS/VCF/adni_mito_genomes_180214_fixed_n258_biallelic.vcf.gz"
relabelled_vcf="/Volumes/TimMcInerney/MitoImpute/data/ADNI_REDO/WGS/VCF/adni_mito_genomes_180214_fixed_n258_biallelic_relabelled.vcf"
'''