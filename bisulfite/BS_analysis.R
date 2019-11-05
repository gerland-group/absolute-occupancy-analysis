#### 0 Initialization ####

source("../BS_Rprofile.R")
source("../BS_functions.R")  # do not source like this when using RStudio debugging points or they will be ignored

if(grepl("spike-in_linker", getwd())) {
  read_mapping = "spike-in_linker"
} else if(grepl("spike-in_backbone", getwd())) {
  read_mapping = "spike-in_backbone"
} else if(grepl("REL606", getwd())) {
  read_mapping = "REL606"
} else {
  read_mapping = "yeast"
}

plot_folder <- "plots/"
dir.create(plot_folder)

RData_folder <- "data/RData_files/"
bam_fragments_RData_folder <- paste0(sub(paste0("_", read_mapping), "", getwd()), "/", RData_folder)  # non yeast runs share the bam_fragment .RData files with the yeast run
dir.create(RData_folder, recursive = TRUE)

bam_folder <- "data/bam_files/"
bam_names <- gsub(".bam", "", grep(".bam$", dir(bam_folder), value=TRUE))


#### 1.1 Prepare bisulfite analysis ####

if(read_mapping == "yeast") {
  
  genome <- read_genome_fasta("../../reference_genomes/used_in_BS_and_RE_R_scripts/S288C_reference_genome_R64-1-1_20110203/")
  
  bad_regions <- read.table("../../external_data/bad_regions.tab", header=TRUE, stringsAsFactors=FALSE)
  
} else if(read_mapping == "spike-in_linker") {
  
  # mapping can be done with the sequence of the generic 601 with generic linkers, but for the analysis, 
  # sites where the nucleotide varies due to reads coming from different 601/linker versions need to be excluded!
  genome <- readDNAStringSet("../../reference_genomes/BS_scer_linker_backbone_ecoli/single_versions/pFMP233_linker_with_cut_601s_varied_sites_are_N.fa")
  names(genome)[[1]] <- "pFMP233_linker_with_cut_601s"
  
} else if(read_mapping == "spike-in_backbone") {
  
  genome <- readDNAStringSet("../../reference_genomes/BS_scer_linker_backbone_ecoli/single_versions/pFMP233_backbone.fa")
  names(genome)[[1]] <- "pFMP233_backbone"
  
} else if(read_mapping == "REL606") {
  
  genome <- readDNAStringSet("../../reference_genomes/BS_scer_linker_backbone_ecoli/single_versions/REL606.fa")
  names(genome)[[1]] <- "REL606"
  
} else {
  stop("read_mapping not valid")
}

context_patterns <- c("HCH", "HCG", "GCH", "GCG")

# have to be revcomplement since the search is on the genome forward strand: corresponding sequence on the plus strand
pattern_regexps <- list("plus" = c("HCH"="^[ATC]C[ATC]$", "HCG"="^[ATC]CG$", "GCH"="^GC[ATC]$", "GCG"="^GCG$"),
                        "minus" = c("HCH"="^[ATG]G[ATG]$", "HCG"="^CG[ATG]$", "GCH"="^[ATG]GC$", "GCG"="^CGC$"))

strands <- c("plus"="+", "minus"="-")

genome_nf3 <- lapply(genome, fasta2num, oligo_length=3)
gen_meth_pos <- calc_genomic_methylation_positions(genome_nf3, pattern_regexps)  # to check genomic positions of methylation sites later on

gen_meth_pos_unlisted <- lapply(gen_meth_pos, function(gen_meth_pos_chr) sort(unlist(gen_meth_pos_chr)))
print(lapply(gen_meth_pos_unlisted, head))
print(lapply(gen_meth_pos_unlisted, tail))
rm(gen_meth_pos_unlisted)


#### 1.2 Load bam objects  ####

# This step is independent of the read_mapping and should be done once for yeast only, the other cases then use the data from the yeast run in bam_fragments_RData_folder

if(read_mapping == "yeast") {
  add_samples_with_removed_duplicates <- FALSE
  for(name in bam_names) {
    cat("## ", name, "\n")
    
    bam_fragments_file <- paste0(bam_fragments_RData_folder, "bam_fragments_", name, "__WD.RData")
    if(!file.exists(bam_fragments_file)) {
      bam_fragments <- read_bam_as_GAlignmentPair(bam_folder, name)
      save(bam_fragments, file=bam_fragments_file)
      rm(bam_fragments)
    }
    
    # and also save a version without the read duplicates
    bam_fragments_file <- paste0(bam_fragments_RData_folder, "bam_fragments_", name, "__ND.RData")
    if(add_samples_with_removed_duplicates && !file.exists(bam_fragments_file)) {
      bam_fragments <- read_bam_as_GAlignmentPair(bam_folder, name, remove_duplicates = TRUE)
      save(bam_fragments, file=bam_fragments_file)
      rm(bam_fragments)
    }
  }
}


#### 1.3 Extract read data ####

data_names <- gsub("(bam_fragments_|.RData)", "", grep("bam_fragments_", dir(bam_fragments_RData_folder), value=TRUE))  # continue with files with and without duplicates if present

for(name in data_names ) {
  cat("## ", name, "\n")
  read_data_file <- paste0(RData_folder, "read_data_", name, ".RData")
  
  if(!file.exists(read_data_file)) {
    read_data <- extract_read_data(name, bam_fragments_RData_folder, genome, genome_nf3)
  
    save(read_data, file=read_data_file)
    rm(read_data)
    gc()
  }
}


#### 2.1 Plot conversion ratios along reads ####
# can be ommitted if three_prime_cutoff is already known (will be overwritten with version including pruned reads later on)

cat("\n")
for(name in data_names) {
  cat("## ", name, "\n")

  plot_conversion_ratios_along_reads(name, RData_folder, with_pruned=FALSE, plot_folder)
}


#### 2.2 Filter reads and prune read ends ####

cat("\n")
for(name in data_names) {
  cat("## ", name, "\n")
  if(!file.exists(paste0(RData_folder, "read_data_filtered_", name, ".RData"))) {
    filter_and_prune_reads(name, five_prime_cutoff = 5, three_prime_cutoff = 5, bad_regions, read_mapping)
  }
}


#### 2.3 Plot conversion ratio along reads with filtered and pruned reads ####

cat("\n")
for(name in data_names) {
  cat("## ", name, "\n")
  
  plot_conversion_ratios_along_reads(name, RData_folder, with_pruned=TRUE, plot_folder)
}


#### 3 Plots and calcuations using filtered read_data ####
# loading read_data and conversion_data from file takes quite some time, so we don't want to do that repeatedly

conversion_ratio_tables_file <- paste0(plot_folder, "conversion_ratio_tables.txt")
cat("", file=conversion_ratio_tables_file)  # overwrite from previous runs

cat("\n# Plots and calcuations using filtered read_data\n")
for(name in data_names) {
  cat("\n## ", name, "\n")

  if(file.exists(paste0(RData_folder, "occupancies_raw_", name, ".RData"))) next
  
  
#### 3.1 Calc conversion ratios and quality groups and save coverage ####

  load(paste0(RData_folder, "read_data_filtered_", name, ".RData"))

  temp_list <- calc_conversion_ratios_and_quality_groups(read_data[["converted"]])  # HCH should always be converted!
  list2env(temp_list, .GlobalEnv)
  rm(temp_list)

  save(conversion_summary, conversion_ratio, quality_groups, file=paste0(RData_folder, "conversion_data_", name, ".RData"))

  cov <- coverage(read_data[["bam"]][quality_groups[["good"]]])
  save(cov, file = paste0(RData_folder, "coverage_", name, ".RData"))

#### 3.2 Plot grouped conversion ratio histograms ####

  plot_conversion_ratios_and_print_summary(name, conversion_ratio, quality_groups, plot_folder, output_file=conversion_ratio_tables_file)
  plot_grouped_conversion_ratio_histograms(name, read_data, conversion_ratio, quality_groups, num_pos_per_read_limits=c(1, Inf), plot_folder=plot_folder)


#### 3.3 Plot grouped coverage histograms ####

  plot_coverage_histograms(name, read_data, quality_groups, genome, max_coverage=ifelse(read_mapping=="yeast", 200, Inf), plot_folder)


#### 3.4 Calc genomic occupancies ####

  temp_list <- calc_genomic_occupancies(name, read_data, quality_groups, genome, gen_meth_pos)
  list2env(temp_list, .GlobalEnv)
  rm(temp_list)

  save(occupancies, counts, counts_not_na, converted, file=paste0(RData_folder, "occupancies_raw_", name, ".RData"))

  rm(read_data, conversion_summary, conversion_ratio, quality_groups, occupancies, counts, counts_not_na, converted)
  
}  # end of loop using read_data and conversion_data


#### 4.1 Plot occupancy histograms ####

for(name in data_names) {
  cat("## ", name, "\n")
  
  plot_occupancy_histogram(name, RData_folder, plot_folder)
}


#### 4.2 Plot occupancies in interesting regions ####

if(read_mapping == "yeast") {
  gene_table <- read.table("../../external_data/single_genes.tab", header=TRUE, stringsAsFactors=FALSE)
  gene_table <- rbind(gene_table, "high_nano_coverage"=list("chr11", 16042-1000, 16042+1000))
  
  RE_info <- read.table("../../restriction_enzyme/RE_info.txt", stringsAsFactors=FALSE)[1:2, ]  # only use BamHI and HindIII
  annotations <- sapply(rownames(RE_info), 
                        function(enzyme) lapply(genome, function(chr_seq) RE_info[enzyme, "middle_shift"] + start(matchPattern(RE_info[enzyme, "recognition_site"], chr_seq))), 
                        simplify=FALSE, USE.NAMES=TRUE)
  
} else if(read_mapping == "spike-in_linker") {
  gene_table <- data.frame(row.names=read_mapping, chr=names(genome), xlim_start=1, xlim_end=width(genome))
  start_601 <- list(144)
  names(start_601) <- names(genome)
  end_601 <- list(93)
  names(end_601) <- names(genome)
  annotations <- list("601 start"=start_601, "601 end"=end_601)
  
} else {
  warning("no regions specified for this read_mapping")
}

if(read_mapping %in% c("yeast", "spike-in_linker")) {
  for(name in data_names) {
    cat("## ", name, "\n")
    
    plot_occupancies_in_region(name, RData_folder, gene_table, genome, gen_meth_pos, context_patterns, annotations, plot_folder)
  }
}


#### 5 Save data ####

save(list = setdiff(ls(), lsf.str()), file = "BS_analysis.RData")  # save all data excluding functions
