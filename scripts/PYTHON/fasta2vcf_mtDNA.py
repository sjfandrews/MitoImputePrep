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

    #reference assembly as np array of ascii codes
    ref_list = np.array([ord(x) for x in ref_fasta], dtype=np.uint8)
    rl = len(ref_fasta) #length of reference assembly
    assert rl == sites_list.shape[0], "Ref length does not match seq length"
    removeme = {ord(x) for x in ['?', 'X', 'N']}
    alt_alleles = []
    allele_counts = []

    if verbose:
        pbar = tqdm(initial=1, total=rl)

    for idx in range(rl):
        alts = sites_list[idx, :]
        rmv = removeme | {ref_list[idx]}
        alts_set = list(set(np.unique(alts)) - rmv)
        if alts_set:
            counts = [(k, list(alts).count(k)) for k in alts_set]
            counts = sorted(counts, key=lambda t: t[1], reverse=True)
            alt_alleles.append([x[0] for x in counts])
            allele_counts.append([x[1] for x in counts])
            if verbose:
                pbar.update(1)
            continue
        elif add_alt:
            alt = [ord('C')] if ref_list[idx] == ord('A') else ord(['A'])
        else:
            # If no alternative alleles found,
            #  set alternative allele to missing
            alt = [ord('.')]
        alt_alleles.append(alt)
        allele_counts.append([ord('.')])
        if verbose:
            pbar.update(1)
    print(idx)
    return {'REF': list(ref_list),'ALT': alt_alleles, 'AC': allele_counts}

# Fill in the GT fields
def proc_snps(sites_list, df_site, labels, diploid=False, dip_symbol='/', verbose=True):
    ref_and_alts = [[r] + a for r,a in zip(df_site['REF'], df_site['ALT'])]

    #preallocate numpy array allowing up to 2 digit diploid or 5 digit haploid
    GT = np.zeros((len(sites_list), len(labels)), dtype='<U5')
    rl = sites_list.shape[0] #length of reference assembly
    if verbose:
        pbar = tqdm(initial=1, total=rl)
    for idx in range(rl):
        ra = ref_and_alts[idx]
        ra = {x: y for x,y in zip(ra, range(len(ra)))}
        # Replace bases with index of alt allele (or "." if not an alt)
        GT_ = [ra.get(samp, '.') for samp in sites_list[idx, :]]
        if diploid: #double the haploid genotype if requested
            GT_ = [str(x)+dip_symbol+str(x) for x in GT_]
        GT[idx, :] = GT_ #save to array
        if verbose:
            pbar.update(1)

    #Turn into DF and return
    return pd.DataFrame(GT, index=range(1, rl + 1),
                        columns=labels)

def other_cols(args, df_site):
    INFO = [','.join([str(i) for i in AC]) for AC in df_site['ALT']]
    df = df_site
    df['INFO'] = INFO
    df['#CHROM'] = 'MT'
    df['POS'] = list(range(1, df.shape[0] + 1))
    df['ID'] = ['MT_{}'.format(x) for x in df['POS']] if args.ID else '.'
    df['ALT'] = [','.join([chr(a) for a in x]) for x in df['ALT']]
    df['REF'] = [chr(x) for x in df['REF']]
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

    ## Reference file
    # look for a ">MT" in the reference fasta and read the next line as the ref
    if verbose:
        print('READ IN REFERENCE:\t\t{!s}'.format(args.ref_fasta.name))
    lastline = ''
    while lastline != '>MT':
        lastline = args.ref_fasta.readline().strip()
    ref_fasta = args.ref_fasta.readline().strip()
    args.ref_fasta.close()

    ## Input file
    if verbose:
        print('READ IN INFILE:\t\t{!s}'.format(args.infile))
    labels = [] # CREATE AN EMPTY LIST TO STORE SAMPLE LABELS
    seqs = bytearray(b'') # CREATE AN EMPTY BYTEARRAY TO STORE SEQUENCES
    with open(args.infile, 'r') as f:
        for line in f:
            line = line.strip() # STRIP WHITE SPACE
            if line.startswith('>'):
                labels.append(line.strip('>')) # APPEND SAMPLE LABELS
            else:
                seqlength = len(ref_fasta)
                assert len(line) == seqlength, 'Unaligned sequence'
                seqs += line.strip().encode() # APPEND SEQUENCE INFORMATION

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

    # CREATE THE METADATA AT THE BEGINNING OF THE VCF
    meta = '''##fileformat=VCFv4.1
    ##contig=<ID=MT,length=16569>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=AC,Number=.,Type=Integer,Description="Alternate allele counts, comma delimited when multiple">'''

    # PULL OUT INFORMATION ABOUT ALTERNATIVE ALLELES
    if verbose:
        print('TRANSPOSING GENOTYPES')
    #numpy array with genotype rows and sample columns
    seqs_mat = np.array(seqs).reshape((len(labels), seqlength)).T
    seqs = None #remove sequences bytearray

    if verbose:
        print('ASSIGNING ALTERNATE ALLELES')
    site_dict = alt_alleles(seqs_mat, ref_fasta, 'MT', args.add_alt, verbose)
    df_site = pd.DataFrame(site_dict, index=range(1, seqlength + 1))

    if verbose:
        print('TRANSLATING GENOTYPES')
    df_snps = proc_snps(seqs_mat, df_site, labels,
                        args.diploid, dip_symbol, verbose)

    if verbose:
        print('CONSTRUCTING OTHER SITE COLUMNS')
    df_allcols = other_cols(args, df_site)

    df_out = pd.concat([df_allcols, df_snps], axis=1)

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
    import pdb; pdb.set_trace()

if __name__=="__main__":
    main()
