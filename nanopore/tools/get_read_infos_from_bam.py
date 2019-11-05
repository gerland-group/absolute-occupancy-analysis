import pysam
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Returns read infos from input bam file')
parser.add_argument('-i', '--input', type=str, required=True)
args = parser.parse_args()

if args.input:
  samfile = pysam.AlignmentFile(args.input, "rb")
  print 'read_name' +'\t'+ 'ref_seq' +'\t'+ 'strand' +'\t'+ 'reference_start' +'\t'+ 'aligned_reference_length' +'\t'+ 'aligned_query_length' +'\t'+ 'length_incl_softclipping' +'\t'+ 'length_incl_hardclipping' +'\t'+ 'mapping_quality' +'\t'+ 'quality_average' +'\t'+ 'quality_std' +'\t'+ 'quality_min' +'\t'+ 'quality_perc_25' +'\t'+ 'quality_median' +'\t'+ 'quality_perc_75' +'\t'+ 'quality_max'
  
  for read in samfile.fetch():
    try:
      s1 = read.query_name +'\t'+ read.reference_name +'\t'+ ('-' if read.is_reverse else '+') +'\t'+ str(read.reference_start) +'\t'+ str(read.reference_length) +'\t'+ str(read.query_alignment_length) +'\t'+ str(read.query_length) +'\t'+ str(read.infer_read_length()) +'\t'+ str(read.mapping_quality)
      
    except TypeError:
      s1 = read.query_name +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN'
      
    try:
      alignment_qualities = np.array(read.query_alignment_qualities)  # contains basecalling qualities for the aligned sequence
      s2 = str(alignment_qualities.mean()) +'\t'+ str(alignment_qualities.std()) +'\t'+ str(alignment_qualities.min()) +'\t'+ str(np.percentile(alignment_qualities,25)) +'\t'+ str(np.median(alignment_qualities)) +'\t'+ str(np.percentile(alignment_qualities,75)) +'\t'+ str(alignment_qualities.max())
      
    except TypeError:
      s2 = 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN' +'\t' + 'NaN'
      
    print s1 +'\t'+ s2
    
  samfile.close()

