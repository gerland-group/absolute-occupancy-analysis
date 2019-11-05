#!/bin/bash

set -x # echo on

NANOPOLISH=$1
RAWDATA=$2
BASECALLS=$3
REFGENOME=$4
MAPPED_READS=$5
NANOPOLISH_OUTPUT=$6
TOOLS=$7
READ_STATS=$8
SITE_STATS=$9
SEQUENCING_SUMMERY=${10}
NUM_THREADS=${11}
PYTHON=${12}
MOTIF=${13}

$NANOPOLISH index -s $SEQUENCING_SUMMERY -d $RAWDATA $BASECALLS  # indexes input files for methylation calling after alignment, later calls to nanopolish need to be started from this folder
echo Indexing run time: $SECONDS
$NANOPOLISH call-methylation -t $NUM_THREADS -r $BASECALLS -g $REFGENOME -b $MAPPED_READS -q $MOTIF > $NANOPOLISH_OUTPUT  # needs the raw data files (provided by nanopolish index)
echo Indexing + methylation calling run time: $SECONDS

$PYTHON $TOOLS/calc_read_stats.py -c 2.5 -i $NANOPOLISH_OUTPUT > $READ_STATS
$PYTHON $TOOLS/calc_site_stats.py -c 2.5 -i $NANOPOLISH_OUTPUT > $SITE_STATS
