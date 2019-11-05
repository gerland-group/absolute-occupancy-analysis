# adapted from calculate_methylation_frequency.py (a Nanopolish script)
# assumes unique read_names, thus gives wrong results when fragments are split and aligned to multiple positions, i.e.
# 'supplementary/chimeric alignments', which happens a lot with minimap2!


import math
import sys
import csv
import argparse
from collections import namedtuple

parser = argparse.ArgumentParser( description='Calculate methylation frequency at genomic CpG sites')
parser.add_argument('-c', '--call-threshold', type=float, required=False, default=2.5)
parser.add_argument('-i', '--input', type=str, required=False)
parser.add_argument('-s', '--split-groups', action='store_true')
args = parser.parse_args()
assert(args.call_threshold is not None)

if args.input:
    in_fh = open(args.input)
else:
    in_fh = sys.stdin
csv_reader = csv.DictReader(in_fh, delimiter='\t')

print("\t".join(["read_name", "chromosome", "first_site", "last_site", "sites_total", "sites_methylated", "sites_ambiguous", "groups_total", "groups_methylated", "groups_ambiguous"]))
current_read_name = ""
for record in csv_reader:

	if current_read_name != record['read_name']:
		if current_read_name != "":  # print stats of previous read
			print("\t".join([str(x) for x in [current_read_name, chromosome, first_site, last_site, sites_total, sites_methylated, sites_ambiguous, groups_total, groups_methylated, groups_ambiguous]]))

		current_read_name = record['read_name']
		sites_total = int(record['num_motifs'])
		groups_total = 1
		chromosome = record['chromosome'] 
		first_site = int(record['start'])
		last_site = int(record['end'])  # needed if only one site in the whole read
		sites_methylated = 0
		sites_ambiguous = 0
		groups_methylated = 0 
		groups_ambiguous = 0
	else:
		sites_total = sites_total + int(record['num_motifs'])  # each site in a group counts individually
		groups_total += 1
		last_site = int(record['end'])

	llr = float(record['log_lik_ratio'])
	if llr >= args.call_threshold:
		sites_methylated = sites_methylated + int(record['num_motifs'])
		groups_methylated += 1
	else:
		if llr >= -args.call_threshold:
			sites_ambiguous = sites_ambiguous + int(record['num_motifs'])
			groups_ambiguous += 1
# print last read
print("\t".join([str(x) for x in [current_read_name, chromosome, first_site, last_site, sites_total, sites_methylated, sites_ambiguous, groups_total, groups_methylated, groups_ambiguous]]))

