#!/usr/bin/env python

samples_csv = "/Users/u5015730/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH.csv"

csv_lines = []

with open(samples_csv, 'r') as sc:
    for line in sc:
        if line.startswith("SAMPLE_NUMBER"):
            line = line.strip("\n")
            line = line.split(",")
            csv_header = line