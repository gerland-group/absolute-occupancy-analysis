import argparse
import re

parser = argparse.ArgumentParser(description='Returns read infos from fastq file')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-c', '--barcode_input', type=str, required=False)
args = parser.parse_args()

if args.input:
    print 'read_name' +'\t'+ 'read_num_in_channel' +'\t'+ 'channel' +'\t'+ 'start_time' +'\t'+ 'barcode' 
    with open(args.input,"r") as fastq_file:
        for ln in fastq_file:
	    if re.search('runid=',ln):
		if re.search('barcode=',ln):
	            print re.search('@(.*) runid=',ln).group(1) +'\t'+ re.search('read=(.*) ch=',ln).group(1) +'\t'+ re.search('ch=(.*) start_time=',ln).group(1) +'\t'+ re.search('start_time=(.*) barcode=',ln).group(1) +'\t'+ re.search('barcode=(.*)',ln).group(1)
		else:
	            print re.search('@(.*) runid=',ln).group(1) +'\t'+ re.search('read=(.*) ch=',ln).group(1) +'\t'+ re.search('ch=(.*) start_time=',ln).group(1) +'\t'+ re.search('start_time=(.*)',ln).group(1) +'\t'+ args.barcode_input
    fastq_file.close()

