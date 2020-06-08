#!/usr/bin/env python

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



def main():
    
    start_time = time.time()
    start_display = 'Process started at ' + time.strftime("%d/%m/%Y") + " " + time.strftime("%H:%M:%S")
    
    #DEFINE INPUT ARGUMENTS
    parser=argparse.ArgumentParser(description='Convert HAPLOBLOCKS to FASTA sequences',formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    
    parser.add_argument('-f', '--fasta_file', dest='fasta_file', type=str, required=True, help='specify the input FASTA file')
    parser.add_argument('-s', '--samples_file', dest='samples_file', type=str, required=True, help='specify the input samples file')
    parser.add_argument('-d', '--out_dir', dest='out_dir', type=str, required=False, help='specify the output directory')
    parser.add_argument('-o', '--out_file', dest='out_file', type=str, required=False, help='specify the output file name')
    parser.add_argument('-v', '--verbose', dest='verbose', action="store_true", required=False, help='turn on verbose mode')
    
    args=parser.parse_args()
    
    fasta_file   = args.fasta_file
    samples_file = args.samples_file
    out_dir      = args.out_dir
    out_file     = args.out_file
    verbose      = args.verbose
    
    if out_dir is None and out_file is None:
        out_file = os.path.dirname(fasta_file) + "/" + os.path.splitext(os.path.basename(fasta_file))[0]
    elif out_dir is not None and out_file is not None:
        if out_dir.endswith("/"):
            out_file = out_dir + os.path.basename(out_file)
        else:
            out_file = out_dir + "/" + os.path.basename(out_file)
    else:
        out_file = os.path.dirname(fasta_file) + "/" + os.path.basename(out_file)
        
    if out_file.endswith(".fasta"):
        pass
    else:
        out_file = os.path.dirname(out_file) + "/" + os.path.splitext(os.path.basename(out_file))[0]
        
    if verbose:
        print
        print(' START FASTA FILE SUBSET '.center(int(terminal_size()[0]), '='))
        print       
    
    #
    
    fasta_samples  = []
    fasta_data     = []
    subset_samples = []
    fasta_dict     = {}
    
    if verbose:
        print "*\tOPENING FASTA FILE"
    with open(fasta_file, "r") as ff:
        for line in ff:
            if line.startswith(">"):
                line = line.strip(">")
                line = line.strip("\n")
                fasta_samples.append(line)
            else:
                line = line.strip("\n")
                fasta_data.append(line)
                
    with open(samples_file, "r") as sf:
        for line in sf:
            line = line.strip("\n")
            subset_samples.append(line)
            
    for seq in range(len(fasta_samples)):
        fasta_dict[fasta_samples[seq]] = fasta_data[seq]
    
    if verbose:
        print "*\tWRITING SUBSETTED FASTA FILE"    
    with open(out_file, "wr") as of:
        for seq in range(len(subset_samples)):
            of.write(">" + subset_samples[seq] + "\n")
            of.write(fasta_dict[fasta_samples[seq]] + "\n")    
    #
    stop_time = time.time()
    stop_display = 'Process completed at ' + time.strftime("%d/%m/%Y") + " " + time.strftime("%H:%M:%S")
    mem_used = memory_usage_resource()   
        
    if verbose:
        print(" SUMMARY ".center(int(terminal_size()[0]), " "))
        print
        print "*\tINPUT FASTA FILE:\t%s" % fasta_file
        print "*\tINPUT SAMPLES FILE:\t%s" % samples_file
        print
        print "*\tOUTPUT FASTA FILE:\t%s" % out_file
        print
        print "*\tORIGINAL FASTA SEQUENCES:\t%s" % len(fasta_samples)
        print "*\tSUBSETTED FASTA SEQUENCES:\t%s" % len(subset_samples)
        print
        print "*\t" + str(start_display)
        print "*\t" + str(stop_display)
        time_adj = time.time() - start_time
        if time_adj < 60:
            print("*\tProcess completed in %s seconds" % (round(time_adj, 2)))
        if time_adj >= 60 and time_adj < 3600:
            print("*\tProcess completed in %s minutes" % (round((time_adj / 60), 2)))
        if time_adj >= 3600 and time_adj < 86400:
            print("*\tProcess completed in %s hours" % (round(((time_adj / 60) / 60), 2)))
        if time_adj >= 86400:
            print("*\tProcess completed in %s days" % (round((((time_adj / 60) / 60) / 24), 2)))
        if mem_used < 1024:
            print "*\tProcess used %s MB of memory" % ("%.2f" % (mem_used))
        if mem_used >= 1024:
            print "*\tProcess used %s GB of memory" % ("%.2f" % (mem_used / 1024))
        print
        print(" FINISH FASTA SUBSET ".center(int(terminal_size()[0]), "="))
        print    

if __name__=="__main__":
    main()        


'''

fasta_file   = "/Users/u5015730/Desktop/SANDBOX/MitoImpute/FASTA/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing.fasta"
samples_file = "/Users/u5015730/Desktop/SANDBOX/MitoImpute/FASTA/ReferencePanel_v1_highQual_MAF0.01_filtered_samples.txt"
out_file     = "/Users/u5015730/Desktop/SANDBOX/MitoImpute/FASTA/ReferencePanel_v1_highQual_MAF0.01_filtered.fasta"

fasta_samples  = []
fasta_data     = []
subset_samples = []
fasta_dict     = {}

with open(fasta_file, "r") as ff:
    for line in ff:
        if line.startswith(">"):
            line = line.strip(">")
            line = line.strip("\n")
            fasta_samples.append(line)
        else:
            line = line.strip("\n")
            fasta_data.append(line)
            
with open(samples_file, "r") as sf:
    for line in sf:
        line = line.strip("\n")
        subset_samples.append(line)
        
for seq in range(len(fasta_samples)):
    fasta_dict[fasta_samples[seq]] = fasta_data[seq]
    
with open(out_file, "wr") as of:
    for seq in range(len(subset_samples)):
        of.write(">" + subset_samples[seq] + "\n")
        of.write(fasta_dict[fasta_samples[seq]] + "\n")
'''