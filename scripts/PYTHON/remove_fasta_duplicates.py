#!/usr/bin/env python2

import time
import os
import sys
import gzip
import argparse
import fcntl, termios, struct
import itertools
from tqdm import *

def memory_usage_resource():
    # FROM: http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
    import resource
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem

def terminal_size():
    # FROM USER pascal https://stackoverflow.com/users/362947/pascal
    # https://stackoverflow.com/questions/566746/how-to-get-linux-console-window-width-in-python
    h, w, hp, wp = struct.unpack('HHHH',
        fcntl.ioctl(0, termios.TIOCGWINSZ,
        struct.pack('HHHH', 0, 0, 0, 0)))
    return w, h

def getKeysByValue(dictionary, valueToFind):
    matched_keys = []
    items = dictionary.items()
    for item in items:
        if item[1] == valueToFind:
            #matched_items[item[0]] = valueToFind
            matched_keys.append(item[0])
    return matched_keys
            
#
def main():
    
    start_time = time.time()
    start_display = 'Process started at ' + time.strftime("%d/%m/%Y") + " " + time.strftime("%H:%M:%S")
    
    #DEFINE INPUT ARGUMENTS
    parser=argparse.ArgumentParser(description='Remove duplicated sequences from FASTA file',formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    
    parser.add_argument('-f', '--fasta_file', dest='fasta_file', type=str, required=True, help='specify the input FASTA file')
    parser.add_argument('-v', '--verbose', dest='verbose', action="store_true", required=False, help='turn on verbose mode')
    
    args=parser.parse_args()
    
    fasta_file = args.fasta_file
    verbose    = args.verbose
    
    if verbose:
        print
        print(' START REMOVAL OF DUPLICATE SEQUENCES '.center(int(terminal_size()[0]), '='))
        print  
    #
    if fasta_file.endswith(".gz"):
        out_fasta  = os.path.dirname(fasta_file) + "/" + os.path.splitext(os.path.splitext(os.path.basename(fasta_file))[0])[0] + "_duplicatesRemoved.fasta.gz"
        out_csv    = os.path.dirname(fasta_file) + "/" + os.path.splitext(os.path.splitext(os.path.basename(fasta_file))[0])[0] + "_duplicate_list.csv"
    else:
        out_fasta  = os.path.dirname(fasta_file) + "/" + os.path.splitext(os.path.basename(fasta_file))[0] + "_duplicatesRemoved.fasta"
        out_csv    = os.path.dirname(fasta_file) + "/" + os.path.splitext(os.path.basename(fasta_file))[0] + "_duplicate_list.csv"
    
    seq_names = []
    seq_data  = []
    
    unique_seqs    = {}
    duplicate_seqs = {}
    
    unique_seqs_list = []
    
    if verbose:
        print "*\tOPENING FASTA FILE"
    
    tkr = 0
    
    if fasta_file.endswith(".gz"):
        with gzip.open(fasta_file, "r") as ff:
            for line in ff:
                tkr += 1
                if verbose:
                    print "*\tPROCESSING LINE: {0}\r".format(tkr),
                if line.startswith(">"):
                    line = line.strip(">")
                    line = line.strip("\n")
                    seq_names.append(line)
                else:
                    line = line.strip("\n")
                    seq_data.append(line)
    else:
        with open(fasta_file, "r") as ff:
            for line in ff:
                tkr += 1
                if verbose:
                    print "*\tPROCESSING LINE: {0}\r".format(tkr),                
                if line.startswith(">"):
                    line = line.strip(">")
                    line = line.strip("\n")
                    seq_names.append(line)
                else:
                    line = line.strip("\n")
                    seq_data.append(line)    
    
    if verbose:
        print "*\tLOOKING FOR DUPLICATED SEQUENCES"
        for seq in tqdm(range(len(seq_names))):
            if seq_data[seq] in unique_seqs.values():
                tmp_dup = getKeysByValue(unique_seqs, seq_data[seq])[0]
                
                if tmp_dup in duplicate_seqs.keys():
                    duplicate_seqs[tmp_dup].append(seq_names[seq])
                else:
                    duplicate_seqs[tmp_dup] = [seq_names[seq]]
            else:
                unique_seqs[seq_names[seq]] = seq_data[seq]
                unique_seqs_list.append(seq_names[seq])        
    else:
        for seq in range(len(seq_names)):
            if seq_data[seq] in unique_seqs.values():
                tmp_dup = getKeysByValue(unique_seqs, seq_data[seq])[0]
                
                if tmp_dup in duplicate_seqs.keys():
                    duplicate_seqs[tmp_dup].append(seq_names[seq])
                else:
                    duplicate_seqs[tmp_dup] = [seq_names[seq]]
            else:
                unique_seqs[seq_names[seq]] = seq_data[seq]
                unique_seqs_list.append(seq_names[seq])        
    
   
            
    if verbose:
        print "*\tWRITING UNIQUE SEQUENCES TO FILE"
        if fasta_file.endswith(".gz"):
            with gzip.open(out_fasta, "wr") as of:
                for seq in tqdm(range(len(unique_seqs_list))):
                    of.write(">" + unique_seqs_list[seq] + "\n")
                    of.write(unique_seqs[unique_seqs_list[seq]] + "\n")
        else:
            with open(out_fasta, "wr") as of:
                for seq in tqdm(range(len(unique_seqs_list))):
                    of.write(">" + unique_seqs_list[seq] + "\n")
                    of.write(unique_seqs[unique_seqs_list[seq]] + "\n")        
    else:
        if fasta_file.endswith(".gz"):
            with gzip.open(out_fasta, "wr") as of:
                for seq in range(len(unique_seqs_list)):
                    of.write(">" + unique_seqs_list[seq] + "\n")
                    of.write(unique_seqs[unique_seqs_list[seq]] + "\n")
        else:
            with open(out_fasta, "wr") as of:
                for seq in range(len(unique_seqs_list)):
                    of.write(">" + unique_seqs_list[seq] + "\n")
                    of.write(unique_seqs[unique_seqs_list[seq]] + "\n")        
    
    
                
    if verbose:
        print "*\tWRITING DUPLICATE LIST TO FILE"
    
    with open(out_csv, "wr") as oc: 
        oc.write("seq_name,duplicate_seqs\n")
        for key, value in duplicate_seqs.iteritems():
            oc.write(key + "," + ";".join(value) + "\n")    
    #
    stop_time = time.time()
    stop_display = 'Process completed at ' + time.strftime("%d/%m/%Y") + " " + time.strftime("%H:%M:%S")
    mem_used = memory_usage_resource()   
        
    if verbose:
        print(' SUMMARY '.center(int(terminal_size()[0]), ' '))
        print
        print '*\tINPUT FASTA FILE:\t\t%s' % fasta_file
        print '*\tOUTPUT FASTA FILE:\t\t%s' % out_fasta
        print '*\tOUTPUT DUPLICATES LIST:\t\t%s' % out_csv
        print
        print "*\tORIGINAL NUMBER OF SEQUENCES:\t%s" % str(len(seq_names))
        print "*\tOUTPUT NUMBER OF SEQUENCES:\t%s" % str(len(unique_seqs_list))
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
        print
        print(' FINISH HAPLOBLOCK 2 FASTA CONVERSION '.center(int(terminal_size()[0]), '='))
        print    

if __name__=="__main__":
    main()    

'''
fasta_file = "/Users/TimMcInerney/GitCode/MitoImputePrep/DerivedData/MasterAlignment/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing.fasta.gz"

if fasta_file.endswith(".gz"):
    out_fasta  = os.path.dirname(fasta_file) + "/" + os.path.splitext(os.path.splitext(os.path.basename(fasta_file))[0])[0] + "_duplicatesRemoved.fasta.gz"
    out_csv    = os.path.dirname(fasta_file) + "/" + os.path.splitext(os.path.splitext(os.path.basename(fasta_file))[0])[0] + "_duplicate_list.csv"
else:
    out_fasta  = os.path.dirname(fasta_file) + "/" + os.path.splitext(os.path.basename(fasta_file))[0] + "_duplicatesRemoved.fasta"
    out_csv    = os.path.dirname(fasta_file) + "/" + os.path.splitext(os.path.basename(fasta_file))[0] + "_duplicate_list.csv"

seq_names = []
seq_data  = []

unique_seqs    = {}
duplicate_seqs = {}

unique_seqs_list = []

if verbose:
    print
    print "*\tOPENING FASTA FILE"
    print

if fasta_file.endswith(".gz"):
    with gzip.open(fasta_file, "r") as ff:
        for line in ff:
            if line.startswith(">"):
                line = line.strip(">")
                line = line.strip("\n")
                seq_names.append(line)
            else:
                line = line.strip("\n")
                seq_data.append(line)
else:
    with open(fasta_file, "r") as ff:
        for line in ff:
            if line.startswith(">"):
                line = line.strip(">")
                line = line.strip("\n")
                seq_names.append(line)
            else:
                line = line.strip("\n")
                seq_data.append(line)    

if verbose:
    print
    print "*\tLOOKING FOR DUPLICATED SEQUENCES"
    print

for seq in range(len(seq_names)):
    if seq_data[seq] in unique_seqs.values():
        tmp_dup = getKeysByValue(unique_seqs, seq_data[seq])[0]
        
        if tmp_dup in duplicate_seqs.keys():
            duplicate_seqs[tmp_dup].append(seq_names[seq])
        else:
            duplicate_seqs[tmp_dup] = [seq_names[seq]]
    else:
        unique_seqs[seq_names[seq]] = seq_data[seq]
        unique_seqs_list.append(seq_names[seq])
        
if verbose:
    print
    print "*\tWRITING UNIQUE SEQUENCES TO FILE"
    print

if fasta_file.endswith(".gz"):
    with gzip.open(out_fasta, "wr") as of:
        for seq in range(len(unique_seqs_list)):
            of.write(">" + unique_seqs_list[seq] + "\n")
            of.write(unique_seqs[unique_seqs_list[seq]] + "\n")
else:
    with open(out_fasta, "wr") as of:
        for seq in range(len(unique_seqs_list)):
            of.write(">" + unique_seqs_list[seq] + "\n")
            of.write(unique_seqs[unique_seqs_list[seq]] + "\n")
            
if verbose:
    print
    print "*\tWRITING DUPLICATE LIST TO FILE"
    print

with open(out_csv, "wr") as oc: 
    oc.write("seq_name,duplicate_seqs\n")
    for key, value in duplicate_seqs.iteritems():
        oc.write(key + "," + ";".join(value) + "\n")
'''