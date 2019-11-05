#!/bin/bash

# To process basecalled and barcoded data (downloaded from the Galaxy server and unpacked to local disk)

# Start with "bash process_barcoded_reads.sh |& tee process_barcoded_reads.log" to obtain a log file 
# Start with "bash process_barcoded_reads.sh &> process_barcoded_reads.log &" to obtain a log file and run extern (e.g. via ssh)

set -x # echo on


# START EDIT

ALBACORE_FOLDER=
FASTQ_NAME=
BARCODES="barcode01 barcode02"

REFGENOME=../../reference_genomes/NP_S288C_REL606_pFMP233_whole/S288C_REL606_pFMP233_whole.fa

USE_PORECHOP=false
NANOPOLISH=

# END EDIT

SEQUENCING_SUMMERY=$ALBACORE_FOLDER"/sequencing_summary.txt"
FASTQ_PORECHOP=$FASTQ_NAME"_porechop.fastq"  # only needed when $USE_PORECHOP is true
FASTQ_NAME=$FASTQ_NAME".fastq"

TOOLS=../tools
WD=`pwd`


# Remove existing infos files
rm all_fastq_read_infos.tsv
rm all_bam_read_infos.tsv


# Go through all barcodes and map the reads
for BARCODE in $BARCODES; do

	echo	
	echo $BARCODE

	RAWDATA=$ALBACORE_FOLDER/workspace/$BARCODE
	BASECALLS=$RAWDATA/$FASTQ_NAME

	cd $WD
	mkdir -p $BARCODE
	cd $BARCODE

	mkdir -p alignment/plots_q30 alignment/log_plots_q30
	mkdir -p methylation

	if [ $USE_PORECHOP = true ] ; then
		BASECALLS_PORECHOP=$RAWDATA/$FASTQ_PORECHOP
    	porechop-runner.py -i $BASECALLS -o $BASECALLS_PORECHOP --discard_middle --threads 16
    	BASECALLS=$BASECALLS_PORECHOP
    	fi
    
	# graphmap align -v 1 -t 16 -r $REFGENOME -d $BASECALLS -o alignment/reads.sam
	minimap2 -ax map-ont -t 16 $REFGENOME $BASECALLS > alignment/reads.sam  # 15 times faster, but splits reads by default
	
	# samtools view -F 4 -q 30 alignment/reads.sam -b -o alignment/reads_quality_above_30_temp.bam -U alignment/reads_quality_below_30_temp.bam  # only use reads without the unmapped flag and with mapping quality above 30
	samtools view -F 2052 -q 30 alignment/reads.sam -b -o alignment/reads_quality_above_30_temp.bam -U alignment/reads_quality_below_30_temp.bam  # also exclude supplementary reads
	
	samtools sort -@ 15 -o alignment/reads_quality_above_30.bam -T alignment/reads1.tmp alignment/reads_quality_above_30_temp.bam
	samtools sort -@ 15 -o alignment/reads_quality_below_30.bam -T alignment/reads2.tmp alignment/reads_quality_below_30_temp.bam

	samtools index -@ 16 alignment/reads_quality_above_30.bam
	samtools index -@ 16 alignment/reads_quality_below_30.bam
	
	rm alignment/reads_quality_above_30_temp.bam alignment/reads_quality_below_30_temp.bam alignment/reads.sam

	python $TOOLS/get_read_infos_from_fastq.py -i $BASECALLS > fastq_read_infos.tsv &
	python $TOOLS/get_read_infos_from_bam.py -i alignment/reads_quality_above_30.bam > alignment/bam_read_infos.tsv &
	python $TOOLS/get_read_infos_from_bam.py -i alignment/reads_quality_below_30.bam > alignment/bam_read_infos_unmapped.tsv &
	wait
	
	if [ -e ../all_fastq_read_infos.tsv ]
	then
		tail -n +2 fastq_read_infos.tsv >> ../all_fastq_read_infos.tsv	# this file has a barcode column, so concatenating them doesn't cause a loss of the barcode information
	else
		cp fastq_read_infos.tsv ../all_fastq_read_infos.tsv
	fi
	if [ -e ../all_bam_read_infos.tsv ]
	then
		tail -n +2 alignment/bam_read_infos.tsv >> ../all_bam_read_infos.tsv  # tail -n +2 prints the file without the header line
	else
		cp alignment/bam_read_infos.tsv ../all_bam_read_infos.tsv
	fi

done


cd $WD

echo $SEQUENCING_SUMMERY

# Go through all barcodes and call methylations
for BARCODE in $BARCODES; do

	echo	
	echo $BARCODE

	RAWDATA=$ALBACORE_FOLDER/workspace/$BARCODE
	BASECALLS=$RAWDATA/$FASTQ_NAME
	if [ $USE_PORECHOP = true ] ; then BASECALLS=$RAWDATA/$FASTQ_PORECHOP; fi
	MAPPED_READS=$BARCODE/alignment/reads_quality_above_30.bam
	NANOPOLISH_OUTPUT=$BARCODE/methylation/methylation_q30.tsv
	READ_STATS=$BARCODE/methylation/read_stats_q30.tsv
	SITE_STATS=$BARCODE/methylation/site_stats_q30.tsv
	NUM_THREADS=16  # Running with less threads but all barcodes in parallel (with "&" at the end of the next command) is a bit faster, if there are not too many barcodes
	
	bash $TOOLS/run_nanopolish.sh $NANOPOLISH $RAWDATA $BASECALLS $REFGENOME $MAPPED_READS $NANOPOLISH_OUTPUT $TOOLS $READ_STATS $SITE_STATS $SEQUENCING_SUMMERY $NUM_THREADS |& tee $BARCODE/run_nanopolish.log

done
wait

echo "All done. Run time: $(($SECONDS/3600))h $((($SECONDS%3600)/60))min $(($SECONDS%60))s."

