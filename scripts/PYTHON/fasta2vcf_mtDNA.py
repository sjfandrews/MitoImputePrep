#!/usr/bin/env python3

## IMPORT NECESSARY PACKAGES
import gzip
import argparse
import re
import os
import time
from tqdm import *
import pandas as pd
import numpy as np
from collections import OrderedDict

# Fill in the CHROM, POS, REF, ALT fields and calculate AC
def alt_alleles(sites_list, ref_fasta, chrom='MT', add_alt=False, verbose=False):
    def count_sort(ref, alts, add_alt):
        alts_set = set(alts) - {'?', 'X', 'N', ref}
        if alts_set:
            counts = {k: alts.count(k) for k in alts_set}
            counts = OrderedDict(sorted(counts.items(),
                                        key=lambda t: t[1],
                                        reverse=True))
            return {'alts': list(counts.keys), 'AC': list(counts.values)}
        else: # if the site is monomorphic
            if add_alt:
                alt = ['C'] if ref == 'A' else ['A']
            else:
                # If no alternative alleles found,
                #  set alternative allele to missing
                alt = ['.']
            return {'alts': alt, 'AC': ['.']}

    def update_df(df, pos, ref):
        counts = count_sort(ref, alts, add_alt)
        df.at[idx, 'POS'] = idx
        df.at[idx, 'REF'] = ref
        df.at[idx, 'ALT'] = counts.alts
        df.at[idx, 'n_alt'] = counts.AC
        return(df)

    ref_list = list(ref_fasta) #reference assembly as char list
    rl = len(ref_list) #length of reference assembly
    assert rl == len(sites_list), "Ref length does not match seq length"
    df = pd.DataFrame(columns=['#CHROM', 'POS', 'REF', 'ALT', 'n_alt'],
                      index=list(range(1, rl + 1)))
    idx = 0
    if verbose:
        pbar = tqdm(initial=0, total=rl - 1)
    for ref, alts in zip(ref_list, sites_list):
        idx += 1
        df = update_df(df, idx, ref, alts)
        if verbose:
            pbar.update(1)
    df['#CHROM'] = chrom
    return df

# Fill in the GT fields
def proc_snps(sites_list, df_site, labels, diploid=False, dip_symbol='/'):
    ref_and_alts = [[r] + a for r,a in zip(df_site['REF'], df_site['ALT'])]

    #preallocate numpy array allowing up to 2 digit diploid or 5 digit haploid
    GT = np.zeros((len(sites_list), len(labels)), dtype='<U5')

    idx = 0
    for site, ra in zip(ref_and_alts, sites_list):
        # Replace bases with index of alt allele (or "." if not an alt)
        GT_ = [ra.index(samp) if samp in ra else '.' for samp in site]
        if diploid: #double the haploid genotype if requested
            GT_ = [str(x)+dip_symbol+str(x) for x in GT_]
        GT[idx, :] = GT_ #save to array
        idx += 1

    #Turn into DF and return
    return pd.DataFrame(GT, index=range(1, len(labels) + 1), columns=labels)

def other_cols(args, df_site):
    INFO = [','.join([str(i) for i in AC]) for AC in df_site['n_alt']]
    df = df_site
    df['INFO'] = INFO
    df['ID'] = ['MT_{}'.format(x) for x in df['POS']] if args.ID else '.'
    df['ALT'] = [','.join(x) for x in df['ALT']
    df['QUAL'] = args.quality
    df['FILTER'] = args.filt
    df['FORMAT'] = 'GT'
    return df[['#CHROM', 'POS', 'ID', 'REF', 'ALT',
               'QUAL', 'FILTER', 'INFO', 'FORMAT']]

## DEFINE MAIN SCRIPT/FUNCTION
def main():
    #START THE TIMER
    start_time = time.time()
    #DEFINE INPUT ARGUMENTS
    parser = argparse.ArgumentParser(
       description='fasta2vcf',
       formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--infile', type=str, required=True, help='input fasta file (.fasta)')
    parser.add_argument('-o', '--outfile', type=str, required=False, help='output VCF file (.vcf or .vcf.gz)')
    #parser.add_argument('-g', '--gap2missing', action="store_true", required=False, help='turn gaps to missing (N)')
    parser.add_argument('-d', '--diploid', action="store_true", required=False, help='create diploid VCF file')
    parser.add_argument('-id', '--ID', action="store_true", required=False, help='tag the ID column as MT<POS>')
    parser.add_argument('-q', '--quality', type=str, required=False, default=".", help='tag for QUAL column')
    parser.add_argument('-f', '--filt', type=str, required=False, default="PASS", help='tag for FILTER column')
    #parser.add_argument('-c', '--chromosome', type=str, default='MT', required=False, help='specify the chromosome')
    parser.add_argument('-v', '--verbose', action="store_true", required=False, help='turn on verbose mode')
    parser.add_argument('-r', '--ref_fasta', type=argparse.FileType('r'), default='rCRS.fasta', help='reference fasta file')
    parser.add_argument('-a', '--add_alt', action="store_true", required=False, help='always have an alternative allele')

    args = parser.parse_args()

    verbose = args.verbose

    dip_symbol = "/"

    # READ IN THE FILES

    ## Input file
    if verbose:
        print('READ IN INFILE:\t\t{!s}'.format(args.infile))
    labels = [] # CREATE AN EMPTY LIST TO STORE SAMPLE LABELS
    seqs = [] # CREATE AN EMPTY LIST TO STORE SEQUENCES
    with open(args.infile, 'r') as f:
        for line in f:
            line = line.strip() # STRIP WHITE SPACE
            if line.startswith('>'):
                labels.append(line.strip('>')) # APPEND SAMPLE LABELS
            else:
                seqs.append(line.strip()) # APPEND SEQUENCE INFORMATION

    ## Output file
    outfile = args.outfile
    if outfile is None:
        print(True)
        outfile = os.path.splitext(args.infile)[0] + '.vcf.gz'
    if verbose:
        print('READ IN OUTFILE:\t\t{!s}'.format(outfile))
    if outfile.endswith('.gz'):
        outfile = gzip.open(outfile,'wt')
    else:
        outfile = open(outfile,'wt')

    ## Reference file
    # look for a ">MT" in the reference fasta and read the next line as the ref
    if verbose:
        print('READ IN REFERENCE:\t\t{!s}'.format(args.ref_fasta.name))
    lastline = ''
    while lastline != '>MT':
        lastline = args.ref_fasta.readline().strip()
    ref_fasta = args.ref_fasta.readline().strip()
    args.ref_fasta.close()

    # CREATE THE METADATA AT THE BEGINNING OF THE VCF
    meta = '''##fileformat=VCFv4.1
    ##contig=<ID=MT,length=16569>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=AC,Number=.,Type=Integer,Description="Alternate allele counts, comma delimited when multiple">'''

    # ASSIGN REFERENCE ALLELES FROM THE ref_fasta TO DICTIONARY
    if verbose:
        print('ASSIGNING REFERENCE ALLELES')

    # PULL OUT INFORMATION ABOUT ALTERNATIVE ALLELES

    seqs_mat = np.array([list(x) for x in seqs])
    sites_list = [list(i) for i in seqs_mat.T] #list of alleles for each site
    samps_list = [list(i) for i in seqs_mat] #list of alleles for each sample

    df_site = alt_alleles(sites_list, ref_fasta, chrom='MT',
                          args.add_alt, verbose)
    df_snps = proc_snps(sites_list, df_site, labels, args.diploid, dip_symbol)
    df_site = other_cols(args, df_site)
    df_out = pd.concat([df_site, df_snps], axis=1)

    outfile.write(meta) # WRITE METADATA LINES TO FILE
    if verbose:
        print('METADATA WRITTEN TO FILE')

    df_out.to_csv(outfile, sep='\t', quotechar='')
    if verbose:
        print('DATA WRITTEN TO FILE')

    outfile.close() # CLOSE OUTFILE
    if verbose:
        print('FILE CLOSED\n')

    if verbose:
        time_adj = time.time() - start_time
        if time_adj < 60:
            print('*\tProcess completed in %s seconds' % (round(time_adj, 2)))
        if time_adj >= 60 and time_adj < 3600:
            print('*\tProcess completed in %s minutes' % (round((time_adj / 60), 2)))
        if time_adj >= 3600 and time_adj < 86400:
            print('*\tProcess completed in %s hours' % (round(((time_adj / 60) / 60), 2)))
        if time_adj >= 86400:
            print('*\tProcess completed in %s days\n' % (round((((time_adj / 60) / 60) / 24), 2)))

if __name__=="__main__":
    main()
