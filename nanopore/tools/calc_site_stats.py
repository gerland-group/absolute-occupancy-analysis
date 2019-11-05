# adapted from calculate_methylation_frequency.py (a Nanopolish script: https://github.com/jts/nanopolish)

import math
import sys
import csv
import argparse
from collections import namedtuple

def make_key(c, s, e):
    return c + ":" + str(s).zfill(8) + ":" + str(e).zfill(8)

class SiteStats:
    def __init__(self, g_size, g_seq):
        self.group_size = g_size
        self.group_sequence = g_seq
	self.hits_total = 0
	self.hits_methylated = 0
	self.hits_ambiguous = 0

def update_call_stats(key, group_size, is_methylated, is_ambiguous, group_sequence):
    if key not in sites:
        sites[key] = SiteStats(group_size, group_sequence)
    sites[key].hits_total += 1
    if is_methylated > 0:
        sites[key].hits_methylated += 1
    if is_ambiguous > 0:
	sites[key].hits_ambiguous += 1 

parser = argparse.ArgumentParser( description='Calculate methylation frequency at genomic CpG sites')
parser.add_argument('-c', '--call-threshold', type=float, required=False, default=2.5)
parser.add_argument('-i', '--input', type=str, required=False)
args = parser.parse_args()
assert(args.call_threshold is not None)

sites = dict()

if args.input:
    in_fh = open(args.input)
else:
    in_fh = sys.stdin

csv_reader = csv.DictReader(in_fh, delimiter='\t')

for record in csv_reader:
    
    group_size = int(record['num_motifs']) 
    llr = float(record['log_lik_ratio'])
    group_sequence = record['sequence']

    is_methylated = llr >= args.call_threshold
    is_ambiguous = llr < args.call_threshold and llr > -args.call_threshold
    
    key = make_key(record['chromosome'], record['start'], record['end'])
    update_call_stats(key, group_size, is_methylated, is_ambiguous, group_sequence)

# header
print("\t".join(["chromosome", "start", "end", "group_size", "hits_total", "hits_methylated", "hits_ambiguous", "group_sequence"]))

for key in sorted(sites.iterkeys()):
    (c, s, e) = key.split(":")
    print("\t".join([str(x) for x in [c, s, e, sites[key].group_size, sites[key].hits_total, sites[key].hits_methylated, sites[key].hits_ambiguous, sites[key].group_sequence]]))

