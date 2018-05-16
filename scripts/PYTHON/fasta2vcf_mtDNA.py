#!/usr/bin/env python2

## IMPORT NECESSARY PACKAGES
import gzip
import argparse
import re
import os
import time
from tqdm import *
from collections import OrderedDict

## DEFINE MAIN SCRIPT/FUNCTION
def main():
    #START THE TIMER
    start_time = time.time()
    #DEFINE INPUT ARGUMENTS
    parser=argparse.ArgumentParser(description='ambiguous2missing',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile', type=str, required=True, help='input fasta file (.fasta)')
    parser.add_argument('-o', '--outfile', dest='outfile', type=str, required=False, help='output VCF file (.vcf or .vcf.gz)')
    #parser.add_argument('-c', '--chromosome', dest='chromosome', type=str, default='MT', required=False, help='specify the chromosome')
    parser.add_argument('-v', '--verbose', dest='verbose', action="store_true", required=False, help='turn on verbose mode')

    args = parser.parse_args()

    infile = args.infile
    outfile = args.outfile
    verbose = args.verbose

    # The revised Cambridge Reference Sequence:
    # Andrews et al. (1999) Reanalysis and revision of the Cambridge reference sequence for human mitochondrial DNA. Nature Genetics 23(2): 147-147.
    # http://doi.org/10.1038/13779
    rCRS = 'GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAAACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCACTTTTAACAGTCACCCCCCAACTAACACATTATTTTCCCCTCCCACTCCCATACTACTAATCTCATCAATACAACCCCCGCCCATCCTACCCAGCACACACACACCGCTGCTAACCCCATACCCCGAACCAACCAAACCCCAAAGACACCCCCCACAGTTTATGTAGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACAAATAGGTTTGGTCCTAGCCTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATCCCCGTTCCAGTGAGTTCACCCTCTAAATCACCACGATCAAAAGGAACAAGCATCAAGCACGCAGCAATGCAGCTCAAAACGCTTAGCCTAGCCACACCCCCACGGGAAACAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGATACCCCACTATGCTTAGCCCTAAACCTCAACAGTTAAATCAACAAAACTGCTCGCCAGAACACTACGAGCCACAGCTTAAAACTCAAAGGACCTGGCGGTGCTTCATATCCCTCTAGAGGAGCCTGTTCTGTAATCGATAAACCCCGATCAACCTCACCACCTCTTGCTCAGCCTATATACCGCCATCTTCAGCAAACCCTGATGAAGGCTACAAAGTAAGCGCAAGTACCCACGTAAAGACGTTAGGTCAAGGTGTAGCCCATGAGGTGGCAAGAAATGGGCTACATTTTCTACCCCAGAAAACTACGATAGCCCTTATGAAACTTAAGGGTCGAAGGTGGATTTAGCAGTAAACTAAGAGTAGAGTGCTTAGTTGAACAGGGCCCTGAAGCGCGTACACACCGCCCGTCACCCTCCTCAAGTATACTTCAAAGGACATTTAACTAAAACCCCTACGCATTTATATAGAGGAGACAAGTCGTAACATGGTAAGTGTACTGGAAAGTGCACTTGGACGAACCAGAGTGTAGCTTAACACAAAGCACCCAACTTACACTTAGGAGATTTCAACTTAACTTGACCGCTCTGAGCTAAACCTAGCCCCAAACCCACTCCACCTTACTACCAGACAACCTTAGCCAAACCATTTACCCAAATAAAGTATAGGCGATAGAAATTGAAACCTGGCGCAATAGATATAGTACCGCAAGGGAAAGATGAAAAATTATAACCAAGCATAATATAGCAAGGACTAACCCCTATACCTTCTGCATAATGAATTAACTAGAAATAACTTTGCAAGGAGAGCCAAAGCTAAGACCCCCGAAACCAGACGAGCTACCTAAGAACAGCTAAAAGAGCACACCCGTCTATGTAGCAAAATAGTGGGAAGATTTATAGGTAGAGGCGACAAACCTACCGAGCCTGGTGATAGCTGGTTGTCCAAGATAGAATCTTAGTTCAACTTTAAATTTGCCCACAGAACCCTCTAAATCCCCTTGTAAATTTAACTGTTAGTCCAAAGAGGAACAGCTCTTTGGACACTAGGAAAAAACCTTGTAGAGAGAGTAAAAAATTTAACACCCATAGTAGGCCTAAAAGCAGCCACCAATTAAGAAAGCGTTCAAGCTCAACACCCACTACCTAAAAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACCAATCTATCACCCTATAGAAGAACTAATGTTAGTATAAGTAACATGAAAACATTCTCCTCCGCATAAGCCTGCGTCAGATTAAAACACTGAACTGACAATTAACAGCCCAATATCTACAATCAACCAACAAGTCATTATTACCCTCACTGTCAACCCAACACAGGCATGCTCATAAGGAAAGGTTAAAAAAAGTAAAAGGAACTCGGCAAATCTTACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCATAATCACTTGTTCCTTAAATAGGGACCTGTATGAATGGCTCCACGAGGGTTCAGCTGTCTCTTACTTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCGGGCATAACACAGCAAGACGAGAAGACCCTATGGAGCTTTAATTTATTAATGCAAACAGTACCTAACAAACCCACAGGTCCTAAACTACCAAACCTGCATTAAAAATTTCGGTTGGGGCGACCTCGGAGCAGAACCCAACCTCCGAGCAGTACATGCTAAGACTTCACCAGTCAAAGCGAACTACTATACTCAATTGATCCAATAACTTGACCAACGGAACAAGTTACCCTAGGGATAACAGCGCAATCCTATTCTAGAGTCCATATCAACAATAGGGTTTACGACCTCGATGTTGGATCAGGACATCCCGATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAAAGTCCTACGTGATCTGAGTTCAGACCGGAGTAATCCAGGTCGGTTTCTATCTACNTTCAAATTCCTCCCTGTACGAAAGGACAAGAGAAATAAGGCCTACTTCACAAAGCGCCTTCCCCCGTAAATGATATCATCTCAACTTAGTATTATACCCACACCCACCCAAGAACAGGGTTTGTTAAGATGGCAGAGCCCGGTAATCGCATAAAACTTAAAACTTTACAGTCAGAGGTTCAATTCCTCTTCTTAACAACATACCCATGGCCAACCTCCTACTCCTCATTGTACCCATTCTAATCGCAATGGCATTCCTAATGCTTACCGAACGAAAAATTCTAGGCTATATACAACTACGCAAAGGCCCCAACGTTGTAGGCCCCTACGGGCTACTACAACCCTTCGCTGACGCCATAAAACTCTTCACCAAAGAGCCCCTAAAACCCGCCACATCTACCATCACCCTCTACATCACCGCCCCGACCTTAGCTCTCACCATCGCTCTTCTACTATGAACCCCCCTCCCCATACCCAACCCCCTGGTCAACCTCAACCTAGGCCTCCTATTTATTCTAGCCACCTCTAGCCTAGCCGTTTACTCAATCCTCTGATCAGGGTGAGCATCAAACTCAAACTACGCCCTGATCGGCGCACTGCGAGCAGTAGCCCAAACAATCTCATATGAAGTCACCCTAGCCATCATTCTACTATCAACATTACTAATAAGTGGCTCCTTTAACCTCTCCACCCTTATCACAACACAAGAACACCTCTGATTACTCCTGCCATCATGACCCTTGGCCATAATATGATTTATCTCCACACTAGCAGAGACCAACCGAACCCCCTTCGACCTTGCCGAAGGGGAGTCCGAACTAGTCTCAGGCTTCAACATCGAATACGCCGCAGGCCCCTTCGCCCTATTCTTCATAGCCGAATACACAAACATTATTATAATAAACACCCTCACCACTACAATCTTCCTAGGAACAACATATGACGCACTCTCCCCTGAACTCTACACAACATATTTTGTCACCAAGACCCTACTTCTAACCTCCCTGTTCTTATGAATTCGAACAGCATACCCCCGATTCCGCTACGACCAACTCATACACCTCCTATGAAAAAACTTCCTACCACTCACCCTAGCATTACTTATATGATATGTCTCCATACCCATTACAATCTCCAGCATTCCCCCTCAAACCTAAGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGCTTAAACCCCCTTATTTCTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAATTAATCCCCTGGCCCAACCCGTCATCTACTCTACCATCTTTGCAGGCACACTCATCACAGCGCTAAGCTCGCACTGATTTTTTACCTGAGTAGGCCTAGAAATAAACATGCTAGCTTTTATTCCAGTTCTAACCAAAAAAATAAACCCTCGTTCCACAGAAGCTGCCATCAAGTATTTCCTCACGCAAGCAACCGCATCCATAATCCTTCTAATAGCTATCCTCTTCAACAATATACTCTCCGGACAATGAACCATAACCAATACTACCAATCAATACTCATCATTAATAATCATAATAGCTATAGCAATAAAACTAGGAATAGCCCCCTTTCACTTCTGAGTCCCAGAGGTTACCCAAGGCACCCCTCTGACATCCGGCCTGCTTCTTCTCACATGACAAAAACTAGCCCCCATCTCAATCATATACCAAATCTCTCCCTCACTAAACGTAAGCCTTCTCCTCACTCTCTCAATCTTATCCATCATAGCAGGCAGTTGAGGTGGATTAAACCAAACCCAGCTACGCAAAATCTTAGCATACTCCTCAATTACCCACATAGGATGAATAATAGCAGTTCTACCGTACAACCCTAACATAACCATTCTTAATTTAACTATTTATATTATCCTAACTACTACCGCATTCCTACTACTCAACTTAAACTCCAGCACCACGACCCTACTACTATCTCGCACCTGAAACAAGCTAACATGACTAACACCCTTAATTCCATCCACCCTCCTCTCCCTAGGAGGCCTGCCCCCGCTAACCGGCTTTTTGCCCAAATGGGCCATTATCGAAGAATTCACAAAAAACAATAGCCTCATCATCCCCACCATCATAGCCACCATCACCCTCCTTAACCTCTACTTCTACCTACGCCTAATCTACTCCACCTCAATCACACTACTCCCCATATCTAACAACGTAAAAATAAAATGACAGTTTGAACATACAAAACCCACCCCATTCCTCCCCACACTCATCGCCCTTACCACGCTACTCCTACCTATCTCCCCTTTTATACTAATAATCTTATAGAAATTTAGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTGCAATACTTAATTTCTGTAACAGCTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTACTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCTAAGCACCCTAATCAACTGGCTTCAATCTACTTCTCCCGCCGCCGGGAAAAAAGGCGGGAGAAGCCCCGGCAGGTTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAATCACCTCGGAGCTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCACTCAGCCATTTTACCTCACCCCCACTGATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACCTATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATCTCCCATATTGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTATGGTCTGAGCTATGATATCAATTGGCTTCCTAGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCATATTTCACCTCCGCTACCATAATCATCGCTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACTCCACGGAAGCAATATGAAATGATCTGCTGCAGTGCTCTGAGCCCTAGGATTCATCTTTCTTTTCACCGTAGGTGGCCTGACTGGCATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACGTTGTAGCCCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATTTCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCCATTTCACTATCATATTCATCGGCGTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCCTATCCGGAATGCCCCGACGTTACTCGGACTACCCCGATGCATACACCACATGAAACATCCTATCATCTGTAGGCTCATTCATTTCTCTAACAGCAGTAATATTAATAATTTTCATGATTTGAGAAGCCTTCGCTTCGAAGCGAAAAGTCCTAATAGTAGAAGAACCCTCCATAAACCTGGAGTGACTATATGGATGCCCCCCACCCTACCACACATTCGAAGAACCCGTATACATAAAATCTAGACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAAAAAGGTATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAATTATAGGCTAAATCCTATATATCTTAATGGCACATGCAGCGCAAGTAGGTCTACAAGACGCTACTTCCCCTATCATAGAAGAGCTTATCACCTTTCATGATCACGCCCTCATAATCATTTTCCTTATCTGCTTCCTAGTCCTGTATGCCCTTTTCCTAACACTCACAACAAAACTAACTAATACTAACATCTCAGACGCTCAGGAAATAGAAACCGTCTGAACTATCCTGCCCGCCATCATCCTAGTCCTCATCGCCCTCCCATCCCTACGCATCCTTTACATAACAGACGAGGTCAACGATCCCTCCCTTACCATCAAATCAATTGGCCACCAATGGTACTGAACCTACGAGTACACCGACTACGGCGGACTAATCTTCAACTCCTACATACTTCCCCCATTATTCCTAGAACCAGGCGACCTGCGACTCCTTGACGTTGACAATCGAGTAGTACTCCCGATTGAAGCCCCCATTCGTATAATAATTACATCACAAGACGTCTTGCACTCATGAGCTGTCCCCACATTAGGCTTAAAAACAGATGCAATTCCCGGACGTCTAAACCAAACCACTTTCACCGCTACACGACCGGGGGTATACTACGGTCAATGCTCTGAAATCTGTGGAGCAAACCACAGTTTCATGCCCATCGTCCTAGAATTAATTCCCCTAAAAATCTTTGAAATAGGGCCCGTATTTACCCTATAGCACCCCCTCTACCCCCTCTAGAGCCCACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAGAGAACCAACACCTCTTTACAGTGAAATGCCCCAACTAAATACTACCGTATGGCCCACCATAATTACCCCCATACTCCTTACACTATTCCTCATCACCCAACTAAAAATATTAAACACAAACTACCACCTACCTCCCTCACCAAAGCCCATAAAAATAAAAAATTATAACAAACCCTGAGAACCAAAATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCCTAGGCCTACCCGCCGCAGTACTGATCATTCTATTTCCCCCTCTATTGATCCCCACCTCCAAATATCTCATCAACAACCGACTAATCACCACCCAACAATGACTAATCAAACTAACCTCAAAACAAATGATAACCATACACAACACTAAAGGACGAACCTGATCTCTTATACTAGTATCCTTAATCATTTTTATTGCCACAACTAACCTCCTCGGACTCCTGCCTCACTCATTTACACCAACCACCCAACTATCTATAAACCTAGCCATGGCCATCCCCTTATGAGCGGGCACAGTGATTATAGGCTTTCGCTCTAAGATTAAAAATGCCCTAGCCCACTTCTTACCACAAGGCACACCTACACCCCTTATCCCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTACGCCTAACCGCTAACATTACTGCAGGCCACCTACTCATGCACCTAATTGGAAGCGCCACCCTAGCAATATCAACCATTAACCTTCCCTCTACACTTATCATCTTCACAATTCTAATTCTACTGACTATCCTAGAAATCGCTGTCGCCTTAATCCAAGCCTACGTTTTCACACTTCTAGTAAGCCTCTACCTGCACGACAACACATAATGACCCACCAATCACATGCCTATCATATAGTAAAACCCAGCCCATGACCCCTAACAGGGGCCCTCTCAGCCCTCCTAATGACCTCCGGCCTAGCCATGTGATTTCACTTCCACTCCATAACGCTCCTCATACTAGGCCTACTAACCAACACACTAACCATATACCAATGATGGCGCGATGTAACACGAGAAAGCACATACCAAGGCCACCACACACCACCTGTCCAAAAAGGCCTTCGATACGGGATAATCCTATTTATTACCTCAGAAGTTTTTTTCTTCGCAGGATTTTTCTGAGCCTTTTACCACTCCAGCCTAGCCCCTACCCCCCAATTAGGAGGGCACTGGCCCCCAACAGGCATCACCCCGCTAAATCCCCTAGAAGTCCCACTCCTAAACACATCCGTATTACTCGCATCAGGAGTATCAATCACCTGAGCTCACCATAGTCTAATAGAAAACAACCGAAACCAAATAATTCAAGCACTGCTTATTACAATTTTACTGGGTCTCTATTTTACCCTCCTACAAGCCTCAGAGTACTTCGAGTCTCCCTTCACCATTTCCGACGGCATCTACGGCTCAACATTTTTTGTAGCCACAGGCTTCCACGGACTTCACGTCATTATTGGCTCAACTTTCCTCACTATCTGCTTCATCCGCCAACTAATATTTCACTTTACATCCAAACATCACTTTGGCTTCGAAGCCGCCGCCTGATACTGGCATTTTGTAGATGTGGTTTGACTATTTCTGTATGTCTCCATCTATTGATGAGGGTCTTACTCTTTTAGTATAAATAGTACCGTTAACTTCCAATTAACTAGTTTTGACAACATTCAAAAAAGAGTAATAAACTTCGCCTTAATTTTAATAATCAACACCCTCCTAGCCTTACTACTAATAATTATTACATTTTGACTACCACAACTCAACGGCTACATAGAAAAATCCACCCCTTACGAGTGCGGCTTCGACCCTATATCCCCCGCCCGCGTCCCTTTCTCCATAAAATTCTTCTTAGTAGCTATTACCTTCTTATTATTTGATCTAGAAATTGCCCTCCTTTTACCCCTACCATGAGCCCTACAAACAACTAACCTGCCACTAATAGTTATGTCATCCCTCTTATTAATCATCATCCTAGCCCTAAGTCTGGCCTATGAGTGACTACAAAAAGGATTAGACTGAACCGAATTGGTATATAGTTTAAACAAAACGAATGATTTCGACTCATTAAATTATGATAATCATATTTACCAAATGCCCCTCATTTACATAAATATTATACTAGCATTTACCATCTCACTTCTAGGAATACTAGTATATCGCTCACACCTCATATCCTCCCTACTATGCCTAGAAGGAATAATACTATCGCTGTTCATTATAGCTACTCTCATAACCCTCAACACCCACTCCCTCTTAGCCAATATTGTGCCTATTGCCATACTAGTCTTTGCCGCCTGCGAAGCAGCGGTGGGCCTAGCCCTACTAGTCTCAATCTCCAACACATATGGCCTAGACTACGTACATAACCTAAACCTACTCCAATGCTAAAACTAATCGTCCCAACAATTATATTACTACCACTGACATGACTTTCCAAAAAACACATAATTTGAATCAACACAACCACCCACAGCCTAATTATTAGCATCATCCCTCTACTATTTTTTAACCAAATCAACAACAACCTATTTAGCTGTTCCCCAACCTTTTCCTCCGACCCCCTAACAACCCCCCTCCTAATACTAACTACCTGACTCCTACCCCTCACAATCATGGCAAGCCAACGCCACTTATCCAGTGAACCACTATCACGAAAAAAACTCTACCTCTCTATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGAACTAATCATATTTTATATCTTCTTCGAAACCACACTTATCCCCACCTTGGCTATCATCACCCGATGAGGCAACCAGCCAGAACGCCTGAACGCAGGCACATACTTCCTATTCTACACCCTAGTAGGCTCCCTTCCCCTACTCATCGCACTAATTTACACTCACAACACCCTAGGCTCACTAAACATTCTACTACTCACTCTCACTGCCCAAGAACTATCAAACTCCTGAGCCAACAACTTAATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTAAAACTAGGCGGCTATGGTATAATACGCCTCACACTCATTCTCAACCCCCTGACAAAACACATAGCCTACCCCTTCCTTGTACTATCCCTATGAGGCATAATTATAACAAGCTCCATCTGCCTACGACAAACAGACCTAAAATCGCTCATTGCATACTCTTCAATCAGCCACATAGCCCTCGTAGTAACAGCCATTCTCATCCAAACCCCCTGAAGCTTCACCGGCGCAGTCATTCTCATAATCGCCCACGGGCTTACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTACGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACTAATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAACCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAACACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAACACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTACGACCCCTTATTTACCGAGAAAGCTCACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGCACACTACTATAACCACCCTAACCCTGACTTCCCTAATTCCCCCCATCCTTACCACCCTCGTTAACCCTAACAAAAAAAACTCATACCCCCATTATGTAAAATCCATTGTCGCATCCACCTTTATTATCAGTCTCTTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTTATTATCTCGAACTGACACTGAGCCACAACCCAAACAACCCAGCTCTCCCTAAGCTTCAAACTAGACTACTTCTCCATAATATTCATCCCTGTAGCATTGTTCGTTACATGGTCCATCATAGAATTCTCACTGTGATATATAAACTCAGACCCAAACATTAATCAGTTCTTCAAATATCTACTCATCTTCCTAATTACCATACTAATCTTAGTTACCGCTAACAACCTATTCCAACTGTTCATCGGCTGAGAGGGCGTAGGAATTATATCCTTCTTGCTCATCAGTTGATGATACGCCCGAGCAGATGCCAACACAGCAGCCATTCAAGCAATCCTATACAACCGTATCGGCGATATCGGTTTCATCCTCGCCTTAGCATGATTTATCCTACACTCCAACTCATGAGACCCACAACAAATAGCCCTTCTAAACGCTAATCCAAGCCTCACCCCACTACTAGGCCTCCTCCTAGCAGCAGCAGGCAAATCAGCCCAATTAGGTCTCCACCCCTGACTCCCCTCAGCCATAGAAGGCCCCACCCCAGTCTCAGCCCTACTCCACTCAAGCACTATAGTTGTAGCAGGAATCTTCTTACTCATCCGCTTCCACCCCCTAGCAGAAAATAGCCCACTAATCCAAACTCTAACACTATGCTTAGGCGCTATCACCACTCTGTTCGCAGCAGTCTGCGCCCTTACACAAAATGACATCAAAAAAATCGTAGCCTTCTCCACTTCAAGTCAACTAGGACTCATAATAGTTACAATCGGCATCAACCAACCACACCTAGCATTCCTGCACATCTGTACCCACGCCTTCTTCAAAGCCATACTATTTATGTGCTCCGGGTCCATCATCCACAACCTTAACAATGAACAAGATATTCGAAAAATAGGAGGACTACTCAAAACCATACCTCTCACTTCAACCTCCCTCACCATTGGCAGCCTAGCATTAGCAGGAATACCTTTCCTCACAGGTTTCTACTCCAAAGACCACATCATCGAAACCGCAAACATATCATACACAAACGCCTGAGCCCTATCTATTACTCTCATCGCTACCTCCCTGACAAGCGCCTATAGCACTCGAATAATTCTTCTCACCCTAACAGGTCAACCTCGCTTCCCCACCCTTACTAACATTAACGAAAATAACCCCACCCTACTAAACCCCATTAAACGCCTGGCAGCCGGAAGCCTATTCGCAGGATTTCTCATTACTAACAACATTTCCCCCGCATCCCCCTTCCAAACAACAATCCCCCTCTACCTAAAACTCACAGCCCTCGCTGTCACTTTCCTAGGACTTCTAACAGCCCTAGACCTCAACTACCTAACCAACAAACTTAAAATAAAATCCCCACTATGCACATTTTATTTCTCCAACATACTCGGATTCTACCCTAGCATCACACACCGCACAATCCCCTATCTAGGCCTTCTTACGAGCCAAAACCTGCCCCTACTCCTCCTAGACCTAACCTGACTAGAAAAGCTATTACCTAAAACAATTTCACAGCACCAAATCTCCACCTCCATCATCACCTCAACCCAAAAAGGCATAATTAAACTTTACTTCCTCTCTTTCTTCTTCCCACTCATCCTAACCCTACTCCTAATCACATAACCTATTCCCCCGAGCAATCTCAATTACAATATATACACCAACAAACAATGTTCAACCAGTAACTACTACTAATCAACGCCCATAATCATACAAAGCCCCCGCACCAATAGGATCCTCCCGAATCAACCCTGACCCCTCTCCTTCATAAATTATTCAGCTTCCTACACTATTAAAGTTTACCACAACCACCACCCCATCATACTCTTTCACCCACAGCACCAATCCTACCTCCATCGCTAACCCCACTAAAACACTCACCAAGACCTCAACCCCTGACCCCCATGCCTCAGGATACTCCTCAATAGCCATCGCTGTAGTATATCCAAAGACAACCATCATTCCCCCTAAATAAATTAAAAAAACTATTAAACCCATATAACCTCCCCCAAAATTCAGAATAATAACACACCCGACCACACCGCTAACAATCAATACTAAACCCCCATAAATAGGAGAAGGCTTAGAAGAAAACCCCACAAACCCCATTACTAAACCCACACTCAACAGAAACAAAGCATACATCATTATTCTCGCACGGACTACAACCACGACCAATGATATGAAAAACCATCGTTGTATTTCAACTACAAGAACACCAATGACCCCAATACGCAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTACTCACCAGACGCCTCAACCGCCTTTTCATCAATCGCCCACATCACTCGAGACGTAAATTATGGCTGAATCATCCGCTACCTTCACGCCAATGGCGCCTCAATATTCTTTATCTGCCTCTTCCTACACATCGGGCGAGGCCTATATTACGGATCATTTCTCTACTCAGAAACCTGAAACATCGGCATTATCCTCCTGCTTGCAACTATAGCAACAGCCTTCATAGGCTATGTCCTCCCGTGAGGCCAAATATCATTCTGAGGGGCCACAGTAATTACAAACTTACTATCCGCCATCCCATACATTGGGACAGACCTAGTTCAATGAATCTGAGGAGGCTACTCAGTAGACAGTCCCACCCTCACACGATTCTTTACCTTTCACTTCATCTTGCCCTTCATTATTGCAGCCCTAGCAACACTCCACCTCCTATTCTTGCACGAAACGGGATCAAACAACCCCCTAGGAATCACCTCCCATTCCGATAAAATCACCTTCCACCCTTACTACACAATCAAAGACGCCCTCGGCTTACTTCTCTTCCTTCTCTCCTTAATGACATTAACACTATTCTCACCAGACCTCCTAGGCGACCCAGACAATTATACCCTAGCCAACCCCTTAAACACCCCTCCCCACATCAAGCCCGAATGATATTTCCTATTCGCCTACACAATTCTCCGATCCGTCCCTAACAAACTAGGAGGCGTCCTTGCCCTATTACTATCCATCCTCATCCTAGCAATAATCCCCATCCTCCATATATCCAAACAACAAAGCATAATATTTCGCCCACTAAGCCAATCACTTTATTGACTCCTAGCCGCAGACCTCCTCATTCTAACCTGAATCGGAGGACAACCAGTAAGCTACCCTTTTACCATCATTGGACAAGTAGCATCCGTACTATACTTCACAACAATCCTAATCCTAATACCAACTATCTCCCTAATTGAAAACAAAATACTCAAATGGGCCTGTCCTTGTAGTATAAACTAATACACCAGTCTTGTAAACCGGAGATGAAAACCTTTTTCCAAGGACAAATCAGAGAAAAAGTCTTTAACTCCACCATTAGCACCCAAAGCTAAGATTCTAATTTAAACTATTCTCTGTTCTTTCATGGGGAAGCAGATTTGGGTACCACCCAAGTATTGACTCACCCATCAACAACCGCTATGTATTTCGTACATTACTGCCAGCCACCATGAATATTGTACGGTACCATAAATACTTGACCACCTGTAGTACATAAAAACCCAATCCACATCAAAACCCCCTCCCCATGCTTACAAGCAAGTACAGCAATCAACCCTCAACTATCACACATCAACTGCAACTCCAAAGCCACCCCTCACCCACTAGGATACCAACAAACCTACCCACCCTTAACAGTACATAGTACATAAAGCCATTTACCGTACATAGCACATTACAGTCAAATCCCTTCTCGTCCCCATGGATGACCCCCCTCAGATAGGGGTCCCTTGACCACCATCCTCCGTGAAATCAATATCCCGCACAAGAGTGCTACTCTCCTCGCTCCGGGCCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGTCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATG'

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

    meta = '##fileformat=VCFv4.1\n##contig=<ID=MT,length=16569>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' # CREATE THE METADATA AT THE BEGINNING OF THE VCF
    header = ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'] # OPEN THE HEADER COLUMNS
    seqs = [] # CREATE AN EMPTY LIST TO STORE SEQUENCES
    labels = [] # CREATE AN EMPTY LIST TO STORE SAMPLE LABELS
    ref_sites = {} # CREATE AN EMPTY DICTIONARY TO STORE REFERENCE ALLELES
    alt_sites = {} # CREATE AN EMPTY DICTIONARY TO STORE ALTERNATIVE ALLELES
    count_d = {} # CREATE AN EMPTY DICTIONARY TO STORE ALTERNATIVE ALLELE COUNTS
    alt_order = {} # CREATE AN EMPTY DICTIONARY TO STORE ALTERNATIVE ALLELE ORDER

    # ASSIGN REFERENCE ALLELES FROM THE rCRS TO DICTIONARY
    if verbose:
        print('ASSIGNING REFERENCE ALLELES')
        for site in tqdm(range(len(rCRS))):
            ref_sites[str(site + 1)] = rCRS[site]
    else:
        for site in range(len(rCRS)):
            ref_sites[str(site + 1)] = rCRS[site]

    # ASSIGN INFORMATION FROM FASTA FILE TO LISTS
    if verbose:
        print('PULLING INFORMATION FROM FASTA FILE')

    for line in infile:
        if line.startswith('>'):
            line = line.strip() # STRIP WHITE SPACE
            header.append(line.strip('>')) # APPEND SAMPLE LABELS TO HEADERS LIST
            labels.append(line.strip('>')) # APPEND SAMPLE LABELS TO LABELS LIST
        else:
            seqs.append(line.strip()) # APPEND SEQUENCE INFORMATION TO SEQS LIST

    # PULL OUT INFORMATION ABOUT ALTERNATIVE ALLELES
    if verbose:
        print('ASSIGNING ALTERNATIVE ALLELES')
        for site in tqdm(range(len(rCRS))):
            #tmp_set = set()
            tmp_lst = [] # OPEN TEMPORARY LIST TO STORE ALLELES FOR GIVEN SITE
            alleleC = {} # OPEN TEMPORARY DICTIONARY TO STORE ALLELE COUNTS FOR GIVEN SITE
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
    else:
        for site in range(len(rCRS)):
            #tmp_set = set()
            tmp_lst = [] # OPEN TEMPORARY LIST TO STORE ALLELES FOR GIVEN SITE
            alleleC = {} # OPEN TEMPORARY DICTIONARY TO STORE ALLELE COUNTS FOR GIVEN SITE
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
