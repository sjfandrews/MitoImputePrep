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

fasta_file = "/Users/TimMcInerney/GitCode/MitoImputePrep/DerivedData/MasterAlignment/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing.fasta.gz"
out_file   = ""

seq_names = []
seq_data  = []

unique_seqs = {}

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

for seq in range(len(seq_names)):
    if seq_data[seq] in unique_seqs.values():
        pass
    else:
        unique_seqs[seq_names[seq]] = seq_data[seq]