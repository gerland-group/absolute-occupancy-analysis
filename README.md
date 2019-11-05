# Absolute occupancy analysis code

This repository contains code used in our paper "Absolute nucleosome occupancy map for the Saccharomyces cerevisiae genome" (Oberbeckmann et al., Genome Research, 2019)

Please download the Supplemental Code archive from the journal website to obtain the used reference genomes, N+1/TSS/TTS positions and dyad position files.

## General

Analysis functions for N+1 alignment, comparisons of different samples, dyad occupancy calculation 
and export to bedgraph files can be found in shared_functions.R

## Restriction Enzymes / Bisulfite / EM-seq

Analyses for restriction enzyme (RE) and bisulfite (BS) start with already mapped bam files (see methods section in paper) and use 
the corresponding analysis R scripts running in numbered project folders (with additional project folders for analysing spike-in regions in the case of BS).

## Nanopore

Nanopore analyses start with the fast5 raw data files and the basecalled fastq files needed to run minimap2 and then Nanopolish 
in the script nanopore/NP1/process_barcoded_reads.sh. 

