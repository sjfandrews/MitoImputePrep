#!/usr/bin/env python

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

def alt_alleles(seqs_mat, ref_fasta, chrom='MT', add_alt=False, verbose=False):
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
        df.at[idx, 'pos'] = idx
        df.at[idx, 'ref'] = ref
        df.at[idx, 'alt'] = counts.alts
        df.at[idx, 'n_alt'] = counts.AC
        return(df)

    ref_list = list(ref_fasta) #reference assembly as char list
    alt_list_list = [list(i) for i in seqs_mat.T] #list of alts for each site
    rl = len(ref_list) #length of reference assembly
    assert rl == len(alt_list_list), "Ref length does not match seq length"
    df = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'n_alt'],
                      index=list(range(1, rl + 1)))
    idx = 0
    if verbose:
        pbar = tqdm(initial=0, total=rl - 1)
    for ref, alts in zip(ref_list, alt_list_list):
        idx += 1
        df = update_df(df, idx, ref, alts)
        if verbose:
            pbar.update(1)
    df['CHROM'] = chrom
    return df

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
    parser.add_argument('-g', '--gap2missing', action="store_true", required=False, help='turn gaps to missing (N)')
    parser.add_argument('-d', '--diploid', action="store_true", required=False, help='create diploid VCF file')

    parser.add_argument('-id', '--ID', action="store_true", required=False, help='tag the ID column as MT<POS>')
    parser.add_argument('-q', '--quality', type=str, required=False, default="999", help='tag for QUAL column')
    parser.add_argument('-f', '--filt', type=str, required=False, default="PASS", help='tag for FILTER column')
    #parser.add_argument('-c', '--chromosome', type=str, default='MT', required=False, help='specify the chromosome')
    parser.add_argument('-v', '--verbose', action="store_true", required=False, help='turn on verbose mode')
    parser.add_argument('-r', '--ref_fasta', type=argparse.FileType('r'), default='rCRS.fasta', help='reference fasta file')
    parser.add_argument('-a', '--add_alt', action="store_true", required=False, help='always have an alternative allele')

    args = parser.parse_args()

    infile = args.infile
    outfile = args.outfile
    gap2missing = args.gap2missing
    ID = args.ID
    quality = args.quality
    filt = args.filt
    diploid = args.diploid
    verbose = args.verbose
    add_alt = args.add_alt

    dip_symbol = "/"

    # The revised Cambridge Reference Sequence:
    # Andrews et al. (1999) Reanalysis and revision of the Cambridge reference sequence for human mitochondrial DNA. Nature Genetics 23(2): 147-147.
    # http://doi.org/10.1038/13779

    # READ IN THE FILES
    if verbose:
        print('READ IN INFILE:\t\t{!s}'.format(infile))
    if outfile is None:
        print(True)
        outfile = os.path.splitext(infile)[0] + '.vcf.gz'
    infile = open(infile, 'r')
    if verbose:
        print('READ IN OUTFILE:\t\t{!s}'.format(outfile))
    if outfile.endswith('.gz'):
        outfile = gzip.open(outfile,'wt')
    else:
        outfile = open(outfile,'wt')

    # look for a ">MT" in the reference fasta and read the next line as the ref
    if verbose:
        print('READ IN REFERENCE:\t\t{!s}'.format(args.ref_fasta.name))
    lastline = ''
    while lastline != '>MT':
        lastline = args.ref_fasta.readline().strip()
    ref_fasta = args.ref_fasta.readline().strip()
    args.ref_fasta.close()

    meta = '##fileformat=VCFv4.1\n##contig=<ID=MT,length=16569>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' # CREATE THE METADATA AT THE BEGINNING OF THE VCF
    header = ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'] # OPEN THE HEADER COLUMNS
    seqs = [] # CREATE AN EMPTY LIST TO STORE SEQUENCES
    labels = [] # CREATE AN EMPTY LIST TO STORE SAMPLE LABELS
    ref_sites = {} # CREATE AN EMPTY DICTIONARY TO STORE REFERENCE ALLELES
    alt_sites = {} # CREATE AN EMPTY DICTIONARY TO STORE ALTERNATIVE ALLELES
    count_d = {} # CREATE AN EMPTY DICTIONARY TO STORE ALTERNATIVE ALLELE COUNTS
    alt_order = {} # CREATE AN EMPTY DICTIONARY TO STORE ALTERNATIVE ALLELE ORDER

    # ASSIGN REFERENCE ALLELES FROM THE ref_fasta TO DICTIONARY
    if verbose:
        print('ASSIGNING REFERENCE ALLELES')
    ref_sites = {k + 1: v for k,v in zip(range(len(ref_fasta)), ref_fasta)}

    # ASSIGN INFORMATION FROM FASTA FILE TO LISTS
    if verbose:
        print('PULLING INFORMATION FROM FASTA FILE')

    for line in infile:
        line = line.strip() # STRIP WHITE SPACE
        if line.startswith('>'):
            header.append(line.strip('>')) # APPEND SAMPLE LABELS TO HEADERS LIST
            labels.append(line.strip('>')) # APPEND SAMPLE LABELS TO LABELS LIST
        else:
            seqs.append(line.strip()) # APPEND SEQUENCE INFORMATION TO SEQS LIST

    # PULL OUT INFORMATION ABOUT ALTERNATIVE ALLELES

    seqs_mat = np.array([list(x) for x in seqs])
    sites_list = [list(i) for i in seqs_mat.T]

    if verbose:
        print('ASSIGNING ALTERNATIVE ALLELES')
        sites_ = tqdm(range(len(ref_fasta)))
    else:
        sites_ = range(len(ref_fasta))

    for site in sites_:
        tmp_lst = [] # OPEN TEMPORARY LIST TO STORE ALLELES FOR GIVEN SITE
        alleleC = {} # OPEN TEMPORARY DICTIONARY TO STORE ALLELE COUNTS FOR GIVEN SITE
        tmp_lst = []
        for sample in range(len(labels)):
            tmp_lst.append(seqs[sample][site])
            # For each sample in the fasta file, append the allele at site to temp list
        tmp_set = set(tmp_lst) # CREATE SET OF ALL ALLELES AT GIVEN SITE
        if ref_sites[str(site + 1)] in tmp_set:
            tmp_set.remove(ref_sites[str(site + 1)]) # REMOVE THE REFERENCE ALLELE FROM SET
        else:
            pass
        if 'N' in tmp_set:
            tmp_set.remove('N') # REMOVE 'N' (MISSING) ALLELE FROM SET
        else:
            pass
        if 'X' in tmp_set:
            tmp_set.remove('X') # REMOVE 'X' (MISSING) ALLELE FROM SET
        else:
            pass
        if '?' in tmp_set:
            tmp_set.remove('?') # REMOVE '?' (MISSING) ALLELE FROM SET
        else:
            pass
        if len(tmp_set) == 0:
            if add_alt:
                if ref_sites[str(site + 1)] == "A":
                    alt_sites[str(site + 1)] = "C"
                    alt_order[str(site + 1)] = "C"
                else:
                    alt_sites[str(site + 1)] = "A"
                    alt_order[str(site + 1)] = "A"
            else:
                # If no alternative alleles found, set alternative allele to missing
                alt_sites[str(site + 1)] = '.'
                alt_order[str(site + 1)] = '.'
        else:
            # Otherwise
            alt_sites[str(site + 1)] = tmp_set # ADD SET OF ALLELES TO DICTIONARY
            for allele in tmp_set:
                if allele == 'N' or allele == 'X' or allele == '?':
                    pass # SKIP MISSING DATA ALLELES
                else:
                    alleleC[allele] = tmp_lst.count(allele) # ADD COUNT OF GIVEN ALLELE AT GIVEN SITE TO DICTIONARY
            count_d[str(site + 1)] = alleleC # ADD SITE-SPECIFIC ALLELE COUNT DICTIONARY TO OVERALL COUNT DICTIONARY
            tmp_order = sorted(count_d[str(site + 1)], key=count_d[str(site + 1)].get, reverse=True) # ORDER SAID DICTIONARY
            alt_order[str(site + 1)] = ','.join(tmp_order) # ADD ORDERED DICTIONARY TO OVERALL DICTIONARY

    header = '\t'.join(header) # JOIN ALL THE HEADER COLUMNS TOGETHER WITH A TAB-DELIMITER
    outfile.write(meta) # WRITE METADATA LINES TO FILE
    if verbose:
        print('METADATA WRITTEN TO FILE')
    outfile.write(header + '\n') # WRITE HEADER LINES TO FILE
    if verbose:
        print('COLUMNS WRITTEN TO FILE')

    # START CREATING THE SITE-WISE LINES OF THE OUTPUT VCF FILE
    if verbose:
        print('WRITING SITE-WISE INFORMATION TO FILE')
        for site in tqdm(range(16569)):
            # SITE INFORMATION COLUMNS FIRST
            tmp_list = ['MT'] # CHROM column
            tmp_list.append(str(site + 1)) # POS column
            if ID:
                tmp_list.append("MT_" + str(site + 1))
            else:
                tmp_list.append('.') # ID column
            tmp_list.append(ref_sites[str(site + 1)]) # REF column
            #tmp_list.append(','.join(alt_sites[str(site + 1)])) # ALT column
            tmp_list.append(alt_order[str(site + 1)]) # ALT column
            tmp_list.append(quality) # QUAL column
            #tmp_list.append('.') # QUAL column
            tmp_list.append(filt) # FILTER column
            #tmp_list.append('.') # FILTER column
            tmp_list.append('.') # INFO column
            tmp_list.append('GT') # FORMAT column

            # ALLELE INFORMATION COLUMNS FOR EACH SAMPLE NEXT
            for sample in range(len(labels)):
                # Has this many if statements to accomodate every IUPAC ambiguity code
                if seqs[sample][site] == 'N':
                    if diploid:
                        tmp_list.append("." + dip_symbol + ".")
                    else:
                        tmp_list.append(".")
                    #tmp_list.append(".")
                elif seqs[sample][site] == ref_sites[str(site + 1)]:
                    if diploid:
                        tmp_list.append("0" + dip_symbol + "0")
                    else:
                        tmp_list.append("0")
                    #tmp_list.append('0')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[0]:
                    if diploid:
                        tmp_list.append("1" + dip_symbol + "1")
                    else:
                        tmp_list.append("1")
                    #tmp_list.append('1')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[1]:
                    if diploid:
                        tmp_list.append("2" + dip_symbol + "2")
                    else:
                        tmp_list.append("2")
                    #tmp_list.append('2')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[2]:
                    if diploid:
                        tmp_list.append("3" + dip_symbol + "3")
                    else:
                        tmp_list.append("3")
                    #tmp_list.append('3')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[3]:
                    if diploid:
                        tmp_list.append("4" + dip_symbol + "4")
                    else:
                        tmp_list.append("4")
                    #tmp_list.append('4')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[4]:
                    if diploid:
                        tmp_list.append("5" + dip_symbol + "5")
                    else:
                        tmp_list.append("5")
                    #tmp_list.append('5')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[5]:
                    if diploid:
                        tmp_list.append("6" + dip_symbol + "6")
                    else:
                        tmp_list.append("6")
                    #tmp_list.append('6')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[6]:
                    if diploid:
                        tmp_list.append("7" + dip_symbol + "7")
                    else:
                        tmp_list.append("7")
                    #tmp_list.append('7')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[7]:
                    if diploid:
                        tmp_list.append("8" + dip_symbol + "8")
                    else:
                        tmp_list.append("8")
                    #tmp_list.append('8')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[8]:
                    if diploid:
                        tmp_list.append("9" + dip_symbol + "9")
                    else:
                        tmp_list.append("9")
                    #tmp_list.append('9')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[9]:
                    if diploid:
                        tmp_list.append("10" + dip_symbol + "10")
                    else:
                        tmp_list.append("10")
                    #tmp_list.append('10')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[10]:
                    if diploid:
                        tmp_list.append("11" + dip_symbol + "11")
                    else:
                        tmp_list.append("11")
                    #tmp_list.append('11')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[11]:
                    if diploid:
                        tmp_list.append("12" + dip_symbol + "12")
                    else:
                        tmp_list.append("12")
                    #tmp_list.append('12')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[12]:
                    if diploid:
                        tmp_list.append("13" + dip_symbol + "13")
                    else:
                        tmp_list.append("13")
                    #tmp_list.append('13')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[13]:
                    if diploid:
                        tmp_list.append("14" + dip_symbol + "14")
                    else:
                        tmp_list.append("14")
                    #tmp_list.append('14')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[14]:
                    if diploid:
                        tmp_list.append("15" + dip_symbol + "15")
                    else:
                        tmp_list.append("15")
                    #tmp_list.append('15')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[15]:
                    if diploid:
                        tmp_list.append("16" + dip_symbol + "16")
                    else:
                        tmp_list.append("16")
                    #tmp_list.append('16')
                else:
                    if diploid:
                        tmp_list.append("." + dip_symbol + ".")
                    else:
                        tmp_list.append(".")
                    #tmp_list.append('.')

            tmp_list = '\t'.join(tmp_list) # JOIN ALL COLUMNS TOGETHER WITH A TAB-DELIMITER
            outfile.write(tmp_list + '\n') # WRITE LINE TO FILE
    else:
        for site in range(16569):
            # SITE INFORMATION COLUMNS FIRST
            tmp_list = ['MT'] # CHROM column
            tmp_list.append(str(site + 1)) # POS column
            tmp_list.append('.') # ID column
            tmp_list.append(ref_sites[str(site + 1)]) # REF column
            #tmp_list.append(','.join(alt_sites[str(site + 1)])) # ALT column
            tmp_list.append(alt_order[str(site + 1)]) # ALT column
            tmp_list.append('.') # QUAL column
            tmp_list.append('.') # FILTER column
            tmp_list.append('.') # INFO column
            tmp_list.append('GT') # FORMAT column

            # ALLELE INFORMATION COLUMNS FOR EACH SAMPLE NEXT
            for sample in range(len(labels)):
                # Has this many if statements to accomodate every IUPAC ambiguity code
                if seqs[sample][site] == 'N':
                    tmp_list.append('.')
                elif seqs[sample][site] == ref_sites[str(site + 1)]:
                    tmp_list.append('0')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[0]:
                    tmp_list.append('1')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[1]:
                    tmp_list.append('2')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[2]:
                    tmp_list.append('3')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[3]:
                    tmp_list.append('4')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[4]:
                    tmp_list.append('5')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[5]:
                    tmp_list.append('6')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[6]:
                    tmp_list.append('7')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[7]:
                    tmp_list.append('8')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[8]:
                    tmp_list.append('9')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[9]:
                    tmp_list.append('10')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[10]:
                    tmp_list.append('11')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[11]:
                    tmp_list.append('12')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[12]:
                    tmp_list.append('13')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[13]:
                    tmp_list.append('14')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[14]:
                    tmp_list.append('15')
                elif seqs[sample][site] == alt_order[str(site + 1)].split(',')[15]:
                    tmp_list.append('14')
                else:
                    tmp_list.append('.')

            tmp_list = '\t'.join(tmp_list) # JOIN ALL COLUMNS TOGETHER WITH A TAB-DELIMITER
            outfile.write(tmp_list + '\n') # WRITE LINE TO FILE

    infile.close() # CLOSE INFILE
    outfile.close() # CLOSE OUTFILE
    if verbose:
        print('FILES CLOSED\n')


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
