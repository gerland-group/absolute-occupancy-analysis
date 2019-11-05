fasta2num <- function(fastas, oligo_length, method="lookupTable", simplify=TRUE) {  # checked 2018-03-05
  
  # attention: for oligo_length > 1 the position of the numeric code is the position of the first nucleotide in the subsequence with the defined oligo length
  
  if(class(fastas) == "DNAString") {
    fastas <- DNAStringSet(fastas)
  }
  
  # could use chartr instead of this, but I have to strsplit them anyway, so this probably has similar speed
  base2num <- c("A"=1, "C"=2, "G"=3, "T"=4, "N"=NA)
  
  if(method == "slidingView") {
    seqnum_zero <- sapply(fastas, function(x) colSums(t(letterFrequencyInSlidingView(x, 1, c("A","C","G","T")))*(1:4)))
  } else if (method == "lookupTable") {  # previous line implementation is faster, see benchmark  (on human sized fasta's this seems faster and has less memory footprint)
    # as.integer removes names (ACGT) and saves as integer, with RAM usage decreasing to half for each step
    seqnum_zero <- sapply(fastas, function(fasta) as.integer(base2num[strsplit(as.character(fasta),"")[[1]]]))
  } else if (method == "memoryLimited") {
    seqnum_zero <- sapply(fastas, function(fasta) {
      positions <- c(seq(1,length(fasta), by=10000000), length(fasta)+1)
      tmp <- integer(length(fasta))
      for(i in 1:(length(positions)-1)) {
        tmp[ positions[i]:(positions[i+1]-1) ] <- base2num[ strsplit(as.character( substr(fasta, positions[i], positions[i+1]-1) ) ,"")[[1]] ]
      }
      return(tmp)
    })
  }
  rownames(seqnum_zero) <- NULL  # because they don't make sense
  seqnum_zero[seqnum_zero == 0] <- NA
  
  seqnum_higher <- seqnum_zero
  if(oligo_length > 1) {
    for(i in 2:oligo_length) {
      seqnum_higher <- as.matrix( (seqnum_higher[-nrow(seqnum_higher), ]-1)*4 + seqnum_zero[i:nrow(seqnum_zero), ] )
    }
  }
  
  if(simplify & ncol(seqnum_higher)==1) {
    seqnum_higher <- as.integer(as.vector(seqnum_higher))
  }
  
  return(seqnum_higher)
}


even_fasta2num <- function(fastas) {  # checked 2018-03-05
  
  # could use chartr instead of this, but I have to strsplit them anyway, so this probably has similar speed
  base2num <- c("A"=1, "C"=2, "G"=3, "T"=4, "N"=NA)
  # as.integer removes names (ACGT) and saves as integer, with RAM usage decreasing to half for each step
  seqnum_zero <- lapply(strsplit(as.character(fastas),""), function(x) as.integer(base2num[x]))  

  return(seqnum_zero)
}


read_bam_as_GAlignment <- function(folder, name, remove_duplicates = FALSE) {  # checked 2018-03-05
  # reads bam files to give non-paired aligned reads (the illumina run can be paired-end though)
  
  bam_file <- BamFile(file.path(folder, paste0(name, '.bam')), index=file.path(folder, name))
  
  flag_good <- scanBamFlag(isDuplicate=NA, isNotPassingQualityControls=FALSE)  # mapping tool needs to set these flags, many don't
  bam_fragments <- readGAlignments(bam_file, param=ScanBamParam(what="seq", flag=flag_good))
  
  if(remove_duplicates) {
    SE_duplicates <- duplicated(mcols(bam_fragments)$seq)
    cat(paste0(name, ": removed ", round(100*sum(SE_duplicates)/length(SE_duplicates), digits=1), "% duplicate single-end reads\n"))
    bam_fragments <- bam_fragments[!SE_duplicates]
  }
  
  # we expect only completely matched ("M") cigar values
  fragments_cigar_good <- cigar(bam_fragments) == paste0(max(qwidth(bam_fragments)), "M")
  if(!all(fragments_cigar_good)) {
    bam_fragments <- bam_fragments[fragments_cigar_good]
    warning(paste0(name," : only using the ", round(100*sum(fragments_cigar_good)/length(fragments_cigar_good), digits=1), "% of fragments that have a ", paste0(max(qwidth(bam_fragments)), "M"), " cigar string."))
  }
  
  return(bam_fragments)
}


read_bam_as_GAlignmentPair <- function(folder, name, remove_duplicates = FALSE) {  # new 2018-07-17

  bam_file <- BamFile(file.path(folder, paste0(name, '.bam')), index=file.path(folder, name))
  
  flag_good <- scanBamFlag(isDuplicate=NA, isNotPassingQualityControls=FALSE)  # mapping tool needs to set these flags, many don't
  bam_fragments <- readGAlignmentPairs(bam_file, param=ScanBamParam(what="seq", flag=flag_good))
  
  if(remove_duplicates) {
    PE_duplicates <- duplicated(paste0(as.character(mcols(bam_fragments@first)$seq), as.character(mcols(bam_fragments@last)$seq)))
    bam_fragments_first <- bam_fragments@first[!PE_duplicates]
    bam_fragments_last <- bam_fragments@last[!PE_duplicates]
    cat(paste0(name, ": removed ", round(100*sum(PE_duplicates)/length(PE_duplicates), digits=1), "% duplicate paired-end reads\n"))
  } else {
    bam_fragments_first <- bam_fragments@first
    bam_fragments_last <- bam_fragments@last
  }
  
  mcols(bam_fragments_first)$is_first <- TRUE
  mcols(bam_fragments_last)$is_first <- FALSE
  
  bam_fragments <- c(bam_fragments_first, bam_fragments_last)
  
  fragments_cigar_good <- cigar(bam_fragments) == paste0(max(qwidth(bam_fragments)), "M")
  if(!all(fragments_cigar_good)) {
    bam_fragments <- bam_fragments[fragments_cigar_good]
    warning(paste0(name, ": only using the ", round(100*sum(fragments_cigar_good)/length(fragments_cigar_good), digits=1), "% of fragments that have a ", paste0(max(qwidth(bam_fragments)), "M"), " cigar string."))
  }
  
  return(bam_fragments)
}


extract_read_data <- function(name, RData_folder, genome, genome_nf3, max_fragment_num = 45000000) {  # checked 2018-03-05
  
  strands <- c("plus"="+", "minus"="-")
  
  read_data <- list("bam" = GAlignments(), "converted"=list(), "positions"=list())
  
  converted_list <- list()
  positions_list <- list()
  for(pat in context_patterns) {
    converted_list[[pat]] <- list()
    positions_list[[pat]] <- list()
  }
  
  load(file=paste0(RData_folder, "bam_fragments_", name, ".RData"))
  
  if( any( !(levels(seqnames(bam_fragments)) %in% names(genome)) ) ) {
    bam_fragments <- bam_fragments[as.character(seqnames(bam_fragments)) %in% names(genome)]
  }

  if(length(bam_fragments) > max_fragment_num) {
    warning(paste0("To keep memory footprint below 60GB: reducing the number of used fragments from ", length(bam_fragments), " to ", max_fragment_num, " at random."))
    bam_fragments <- bam_fragments[sort(sample(length(bam_fragments), max_fragment_num, replace = FALSE))]
  }
  
  for(chr in names(genome_nf3)) {
    genome_nf3[[chr]] <- c(NA, genome_nf3[[chr]], NA)  # add NAs such that triple entry position is the middle of the triple and the length is the correct chromosome length
  }
  
  for(st in names(strands)) {
    print(st)
    for(chr in names(genome_nf3)) {
      print(chr)
      
      ## step one: find possible methylation sites on the reads using the reference sequence ('positions')
      # this part is going to do the same calculations multiple times if the coverage is high and reads overlap
      # improvement possible: look for the sites on the reference sequence first and then transfer them to the reads)
      
      bam_fragments_strand <- bam_fragments[ (strand(bam_fragments)==strands[st]) & (as.character(seqnames(bam_fragments)) == chr) ]
      
      reads_genomic_nf3 <- mapply(function(start, end) genome_nf3[[chr]][start:end], start(bam_fragments_strand), end(bam_fragments_strand), SIMPLIFY=FALSE)
      
      # for speed up in trade for RAM usage, replace the inner lapplys with  mclapplys in the following steps
      
      # for each of the four pat (HCH, HCG, GCH, GCG), list of all reads with logic vecs indicating pat match with genomic sequence
      reads_genomic_pat_matches <- lapply(pattern_regexps[[st]], function(regexp) {
        nf3_matched <- grepl(regexp, oligo_names(3))
        lapply(reads_genomic_nf3, function (read_genomic_nf3) {
          read_genomic_pat_matches <- nf3_matched[read_genomic_nf3]
          read_genomic_pat_matches <- read_genomic_pat_matches & !is.na(read_genomic_pat_matches)  # set NAs to FALSE
          })
        })
      rm(reads_genomic_nf3)

      # for each of the four pat (HCH, HCG, GCH, GCG), list of all reads with pattern positions on reference sequence
      positions <- lapply(reads_genomic_pat_matches, function(pattern) lapply(pattern, which))

      ## step two: find converted methylation sites (these are the not methylated sites)
      
      # we expect only completely matched ("M") cigar values
      read_seqs_on_ref_nf1 <- even_fasta2num(mcols(bam_fragments_strand)$seq)
      
      converted <- lapply(reads_genomic_pat_matches, function(reads_genomic_spec_pat_matches) mapply('[', read_seqs_on_ref_nf1, reads_genomic_spec_pat_matches))
      rm(reads_genomic_pat_matches)
      
      if(st == "plus") {  # A C G T (converted Cs show up as Ts)
        conversion <- as.integer(c(NA, 0, NA, 1))
      } else if (st == "minus") {  # on the minus strand converted Cs show up as Ts and are complemented to As on the plus strand sequence
        conversion <- as.integer(c(1, NA, 0, NA))
      }
      converted <- lapply(converted, lapply, function(read_seqs_nf1_at_pat_pos) conversion[read_seqs_nf1_at_pat_pos])
      
      bam_fragments_strand@elementMetadata@listData$seq = NULL  # delete sequence data since we already used it and it takes up too much space in RAM after concatenation (probably a bug)
      
      read_data[["bam"]] <- c(read_data[["bam"]], bam_fragments_strand)
      
      for(pat in context_patterns) {
        converted_list[[pat]][[paste0(st, "_", chr)]] <- converted[[pat]]
        positions_list[[pat]][[paste0(st, "_", chr)]] <- positions[[pat]]
      }
      
      rm(bam_fragments_strand)
      rm(converted)
      rm(positions)
    }
  }
  rm(bam_fragments)
  
  for(pat in context_patterns) {
    read_data[["converted"]][[pat]] <- unlist(converted_list[[pat]], recursive = FALSE, use.names = FALSE)
  }
  rm(converted_list)
  
  for(pat in context_patterns) {
    read_data[["positions"]][[pat]] <- unlist(positions_list[[pat]], recursive = FALSE, use.names = FALSE)
  }
  rm(positions_list)
  
  return(read_data)
}


filter_and_prune_reads <- function(name, five_prime_cutoff = 5, three_prime_cutoff = 10, bad_regions, read_mapping) {  # new 2018-08-17
  # This function combines the code of different old functions for pruning and filtering to avoid copying read_data,
  # as well as loading and saving it often, as this takes much RAM and time.

  # 0) Load read_data
  
  load(paste0(RData_folder, "read_data_", name, ".RData"))
  
  # 1) Prune possibly unmethylated regions due to end-repair
  
  # 1.1) Trim outside three prime ends for paired-end runs
  #
  # Use this for paired-end runs only, since for single-end sequencing the 3' end is always inside the fragment.
  # Use if end repair was done with non-methylated Cs, or if you don't want to use the trim_before_first_methylated_C function and end repair was done with methylated Cs.
  #
  # During sonication double strands are not cut with blunt ends. 
  # If on one side 5'-end is longer than 3'-end, then polymerase fills the 3'-end to the 5'-end and/or enzymes can eat up a longer 5'-end or a longer 3'-end.
  # If non-methylated Cs are used to fill 3'-ends, they are always converted. 
  # So we cut each fragment on the 3'-end by a fixed cutoff.
  #
  # We only need to cut the 3'-ends of the second read in a paired-end sequence read pair, because this was the real 3'-end of the fragment.
  #
  # If we have the is_first information: the read with the three prime end at the fragment end is the second fragment.
  # This is independent of the strand, however the strand determines where on the read it is.
  #
  # sequenced fragment: 5'-----------------------------3'
  # first read:         5'---------3'
  # second read:                          5'-----------3' (paired-end sequencing only)

  if(TRUE) {  # we always paired-end runs
    cat("Pruning outside three prime ends...\n")
    
    first_reads <- mcols(read_data[["bam"]])$is_first
    if(is.null(first_reads)) {
      stop("No first_reads information in read_data[['bam']].")
    }
    
    plus_reads <- as.vector(strand(read_data[["bam"]])=='+')
    read_lengths <- as.vector(width(read_data[["bam"]]))
    
    for(pat in names(read_data[["converted"]])) {
      
      # cut + strand read ends (the 3' end)
      positions_good <- mapply('<=', read_data[["positions"]][[pat]], read_lengths - ifelse(first_reads, 0, three_prime_cutoff))
      
      # cut - strand read starts (the 3' end)
      if(sum(!plus_reads) > 0) {
        positions_good[!plus_reads] <- mapply( '>', read_data[["positions"]][[pat]][!plus_reads], ifelse(first_reads[!plus_reads], 0, three_prime_cutoff))
      }
      
      read_data[["positions"]][[pat]] <- mapply('[', read_data[["positions"]][[pat]], positions_good)
      read_data[["converted"]][[pat]] <- mapply('[', read_data[["converted"]][[pat]], positions_good)
      
      rm(positions_good)
    }
  }
  
  # 1.2) Trim outside five prime ends
  #
  # Use this for paired-end and single-end runs, since the outer 5'-end of the fragment is always the 5'-end of the first read.
  #
  # It is unclear why the 5'-end shows a conversion drop/increase along the average read. Maybe due to lower sequencing quality at the beginning?
  #
  # sequenced fragment: 5'-----------------------------3'
  # first read:         5'---------3'
  # second read:                          5'-----------3' (paired-end sequencing only)
  #
  # In a paired-end run, we only want to cut the 5'-end of the first read.
  #
  # If we don't have the is_first information, this function cuts all 5'-ends.
  
  cat("Pruning outside five prime ends...\n")
  
  first_reads <- mcols(read_data[["bam"]])$is_first
  if(is.null(first_reads)) first_reads <- rep(TRUE, length(read_data[["bam"]]))
  
  plus_reads <- as.vector(strand(read_data[["bam"]])=='+')
  read_lengths <- as.vector(width(read_data[["bam"]]))
  
  for(pat in names(read_data[["converted"]])) {
    
    # cut + strand read ends (the 5' end)
    positions_good <- mapply('>', read_data[["positions"]][[pat]], ifelse(first_reads, five_prime_cutoff, 0))
    
    # cut - strand read starts (the 5' end)
    if(sum(!plus_reads) > 0) {
      positions_good[!plus_reads] <- mapply( '<=', read_data[["positions"]][[pat]][!plus_reads], read_lengths[!plus_reads] - ifelse(first_reads[!plus_reads], five_prime_cutoff, 0))
    }
    
    read_data[["positions"]][[pat]] <- mapply('[', read_data[["positions"]][[pat]], positions_good)
    read_data[["converted"]][[pat]] <- mapply('[', read_data[["converted"]][[pat]], positions_good)
    
    rm(positions_good)
  }
  
  
  # 2) Filter out reads in bad regions
  
  if(read_mapping == "yeast") {
    cat("Filtering out reads in bad regions...\n")
    
    to_rm <- rep(FALSE, length(read_data[["bam"]]))
    for(row in 1:nrow(bad_regions)) {
      to_rm <- to_rm | ( (seqnames(read_data[["bam"]]) == bad_regions$chr[row]) & (start(read_data[["bam"]]) < bad_regions$end[row]) & (end(read_data[["bam"]]) > bad_regions$start[row]) )
    }
    to_rm <- as.vector(to_rm)
    
    read_data[["bam"]] <- read_data[["bam"]][!to_rm]
    for(pat in context_patterns) {
      read_data[["positions"]][[pat]] <- read_data[["positions"]][[pat]][!to_rm]
      read_data[["converted"]][[pat]] <- read_data[["converted"]][[pat]][!to_rm]
    }
  }
  
  
  # 3) Save pruned and filtered read data
  
  cat("Saving read_data...\n")
  save(read_data, file=paste0(RData_folder, "read_data_filtered_", name, ".RData"))
}


plot_conversion_ratios_along_reads <- function(name, RData_folder, with_pruned=TRUE, plot_folder, Illumina_seq = TRUE) {  # checked 2018-03-05
  
  png(paste0(plot_folder, "conversion_ratios_along_reads_", name, ".png"), height = ifelse(with_pruned, 2048, 1024), width = 1440, res = 144)
  par_orig <- par(mfrow = c(ifelse(with_pruned, 4, 2), 2))
  on.exit(par(par_orig))
  
  if(with_pruned) {
    names <- paste0(c("", "filtered_"), name)
  } else {
    names <- name
  }
  for(name in names) {
    load(paste0(RData_folder, "read_data_", name, ".RData"))
    
    con_ratios <- list()
    pos_max <- 50
    for(pat in names(read_data[["converted"]])) {
      read_strand <- strand(read_data[["bam"]])
      con_ratios[[pat]] <- list()
      for(s in c("+", "-")) {
        con_ratios[[pat]][[s]] <- list()
        for(frag_order in c("first", "last")) {  # for !Illumina-seq, "first" and "last" are start and end aligned reads
          is_first <- frag_order=="first"
          read_is_first <- mcols(read_data[["bam"]])$is_first  # equal to NULL for single-end sequencing
          if(is.null(read_is_first)) read_is_first <- rep(TRUE, length(read_data[["bam"]]))
          if(any(as.vector(read_strand == s) & (!Illumina_seq | read_is_first == is_first))) {
            read_data_positions <- read_data[["positions"]][[pat]][as.vector(read_strand == s) & (!Illumina_seq | read_is_first == is_first)]
            if(!Illumina_seq & frag_order == "last") {  # to get the pos wrt the other read end
              read_data_positions <- mapply(function(pos_vec, read_width) read_width - pos_vec, read_data_positions, width(read_data$bam[as.vector(read_strand == s)]))
            }
            pos <- unlist(read_data_positions)
            if(length(pos)>0) {
              pos <- as.ff(pos)
              pos_max <- min(max(pos_max, max(pos)), 2000)
            }
            con <- unlist(read_data[["converted"]][[pat]][as.vector(read_strand == s) & (!Illumina_seq | read_is_first == is_first)])
            if(length(con)>0) con <- as.ff(con)
            if(length(pos) != length(con)) {
              warning("length(pos) != length(con)")
              con_ratios[[pat]][[s]][[frag_order]] <- c()
            } else {
              if(!Illumina_seq & frag_order == "last") {
                
              }
              con_table <- binned_sum(con, pos, nbins=pos_max, na.rm=TRUE)
              con_ratios[[pat]][[s]][[frag_order]] <- con_table[ ,"sum"] / con_table[ ,"count"]
            }
          } else {
            con_ratios[[pat]][[s]][[frag_order]] <- c()
          }
        }
      }
    }
    for(s in c("+", "-")) {
      for(frag_order in c("first", "last")) {
        if(Illumina_seq) {
          main_string <- paste0(frag_order, " read, ", s, " strand: ", ifelse(s=="+", "5' to 3'", "3' to 5'"), ifelse(grepl("pruned", name), ", pruned", ""))
        } else {
          temp <- c("first" = "read start", "last" = "read end")
          main_string <- paste0(temp[frag_order], ", ", s, " strand: ", ifelse(s=="+", "5' to 3'", "3' to 5'"), ifelse(grepl("pruned", name), ", pruned", ""))
        }
        if(!Illumina_seq & frag_order == "last") {
          my_xlim = c(pos_max, 1)
        } else {
          my_xlim = c(1, pos_max)
        }
        plot(NULL, xlim=my_xlim, ylim=c(0, 1), main=main_string, xlab="position", ylab="conversion ratio")
        for(pat in seq_along(con_ratios)) {
          if(length(con_ratios[[pat]][[s]][[frag_order]]) > 0) {
            if(pos_max > 200) {
              con_ratios_temp <- rollmean(con_ratios[[pat]][[s]][[frag_order]], k = 5, fill = NaN)
              names(con_ratios_temp) <- names(con_ratios[[pat]][[s]][[frag_order]])
            } else {
              con_ratios_temp <- con_ratios[[pat]][[s]][[frag_order]]
            }
            lines(names(con_ratios_temp), con_ratios_temp, col=pat+1)
            grid()
          }
        }
        legend("bottom", fill=seq_along(con_ratios)+1, legend=names(con_ratios), ncol=length(con_ratios))
      }
    }
    rm(read_data)
  }
  dev.off()
}


calc_conversion_ratios_and_quality_groups <- function(read_data_converted) {  # checked 2018-02-07
  
  conversion_summary <- list()  # a bit faster than using lapply another time
  for(pat in names(read_data_converted)) {
    print(pat)
    conversion_summary[[pat]] <- sapply(read_data_converted[[pat]], function(conversion_vec) c("unconverted" = sum(conversion_vec == 0, na.rm=TRUE), "converted" = sum(conversion_vec == 1, na.rm=TRUE)))
    
    # first attempt to parallelize the new version was much slower than single threaded version
    # conversion_summary[[pat]] <- do.call(cbind, mclapply(read_data_converted[[pat]], function(conversion_vec) c("unconverted" = sum(conversion_vec == 0, na.rm=TRUE), "converted" = sum(conversion_vec == 1, na.rm=TRUE))))
  }
  
  conversion_ratio <- lapply(conversion_summary, function (read_summary) read_summary['converted', ] / (read_summary['converted', ] + read_summary['unconverted', ]))  # much quicker than using apply to do the calculation for each column

  HCH_conv_threshold <- 1

  quality_groups <- list()
  quality_groups[["good"]]     <- !is.na(conversion_ratio[["HCH"]]) & conversion_ratio[["HCH"]] >= HCH_conv_threshold
  quality_groups[["between"]]  <- !is.na(conversion_ratio[["HCH"]]) & conversion_ratio[["HCH"]] > 0 & conversion_ratio[["HCH"]] < HCH_conv_threshold
  quality_groups[["bad"]]      <- !is.na(conversion_ratio[["HCH"]]) & conversion_ratio[["HCH"]] == 0
  quality_groups[["unknown"]]  <-  is.na(conversion_ratio[["HCH"]])
  
  return(list("conversion_summary" = conversion_summary, "conversion_ratio" = conversion_ratio, "quality_groups" = quality_groups)
  )
}


plot_conversion_ratios_and_print_summary <- function(name, conversion_ratio, quality_groups, plot_folder = "", output_file = "") {  # checked 2018-02-08

  cat("\n## ", name, "\n", file = output_file, append = TRUE)
  
  png(filename = paste0(plot_folder, "conversion_ratio_hists_groups_combined_", name, ".png"), width=1024, height=1024, res=144)
  par(mfrow = c(2, 2), oma=c(0, 0, 2, 0))
  
  for(pat in names(conversion_ratio)) {
    if( sum(!is.na(conversion_ratio[[pat]])) > 0 ) {
      hist(conversion_ratio[[pat]], main = pat)
    } else {
      plot.new()
    }
  }
  mtext(name, side=3, line=0, outer=TRUE, cex=1.5, font=1)
  dev.off()
  
  cat("\n### average conversion (all)\n", file = output_file, append = TRUE)
  cat(kable(as.matrix(t(sapply(conversion_ratio, mean, na.rm=TRUE)))), sep="\n", file = output_file, append = TRUE)

  cat("\n### read group averages\n", file = output_file, append = TRUE)
  cat(kable( cbind("reads" = sapply(quality_groups, sum), t(sapply(quality_groups, function(group) (sapply(conversion_ratio, function (x) mean(x[group], na.rm=TRUE)))))) ), sep="\n", file = output_file, append = TRUE)
}


plot_grouped_conversion_ratio_histograms <- function(name, read_data, conversion_ratio, quality_groups, num_pos_per_read_limits=c(6,6), plot_folder) {  # checked 2018-03-06
  
  png(filename=paste0(plot_folder, "conversion_ratio_hists_for_groups_", name, ".png"), width=1024*2, height=1024*2, res=144)
  par(mfrow=c(length(conversion_ratio), length(quality_groups)), oma=c(0, 0, 2, 0))
  
  for(pat in names(conversion_ratio)) {
    lengths <- unlist(lapply(read_data[["converted"]][[pat]], length))
    for(group in names(quality_groups)) {
      
      reads_to_use <- quality_groups[[group]] & (lengths >= num_pos_per_read_limits[1]) & (lengths <= num_pos_per_read_limits[2])
      
      if( sum(!is.na(conversion_ratio[[pat]][reads_to_use])) > 0 ) {
        hist(conversion_ratio[[pat]][reads_to_use], breaks=seq(0, 1, 0.05), xlim=c(0, 1), main=paste(name, pat, group, "\n # of positions in", paste(num_pos_per_read_limits, collapse=" - ")), xlab="avg. fragment conversion ratio")
      } else {
        plot.new()
      }
    }
  }
  mtext(name, side=3, line=0, outer=TRUE, cex=1.5, font=1)
  dev.off()
}


plot_coverage_histograms <- function(name, read_data, quality_groups, genome, max_coverage=200, plot_folder) {  # checked 2018-02-13
  
  coverage_ff <- list()
  
  for(group in names(quality_groups)) {
    coverage_list <- coverage(read_data[["bam"]][quality_groups[[group]]])[names(genome)]
    coverage_ff[[group]] <- do.call("c", lapply(lapply(coverage_list, as.vector), as.ff))
  }
  
  which_max_ff <- function(x) {
    return( as.ram(ff(1:length(x))[(x == max(x))]) )
  }
  
  png(paste0(plot_folder, "coverage_genomewide_", name, ".png"), width=1024, height=1024, res=144)
  par(mfrow=c(2,2), oma=c(0, 0, 2, 0))
  
  for(group in names(quality_groups)) {
    coverage_temp <- coverage_ff[[group]][coverage_ff[[group]] < max_coverage]
    if(length(coverage_temp) > 0) {
      hist(coverage_temp, main=group, xlab="read coverage (for each genomic position)", breaks = 20,
         sub=paste0("over ", max_coverage, ": ", sum(coverage_ff[[group]] > max_coverage), " max: ", max(coverage_ff[[group]], na.rm=TRUE), " max positions: ", range(which_max_ff(coverage_ff[[group]])) ) )
    } else {
      warning("No data to plot")
      plot.new()
    }
  }
  mtext(name, side=3, line=0, outer=TRUE, cex=1.5, font=1)
  dev.off()
  
  png(paste0(plot_folder, "coverage_group_ratios_genomewide_", name, "_good_vs_", group, ".png"), width=1024, height=1024/2, res=144)
  par(mfrow=c(1,2), oma=c(0, 0, 2, 0))
  
  for(group in c("between", "bad")) {
    ratio <- (coverage_ff[["good"]] / coverage_ff[[group]])[(coverage_ff[["good"]] > 0 & coverage_ff[[group]] > 0)]
    hist(log2(ratio), main=group, xlab=paste0("log2 (coverage_good/coverage_", group, ")"), breaks = 20,
         sub=paste0("good zero: ", sum(coverage_ff[["good"]] == 0 & coverage_ff[[group]] > 0), " ", group, " zero: ", sum(coverage_ff[["good"]] > 0 & coverage_ff[[group]] == 0)) )
  }
  mtext(name, side=3, line=0, outer=TRUE, cex=1.5, font=1)
  dev.off()
}


calc_genomic_methylation_positions <- function(genome_nf3, pattern_regexps) {  # new 2018-03-08
  gen_meth_pos <- list()
  for(chr in names(genome_nf3)) {
    gen_meth_pos[[chr]] <- list()
    for(pat in names(pattern_regexps[[1]])) {
      gen_meth_pos[[chr]][[pat]] <- sapply(names(pattern_regexps), function(st) {
        regexp <- pattern_regexps[[st]][[pat]]
        nf3_matches <- grep(regexp, oligo_names(3))
        return(which(genome_nf3[[chr]] %in% nf3_matches) + 1)
      }, simplify = FALSE, USE.NAMES = TRUE)
    }
  }
  return(gen_meth_pos)
}


calc_genomic_occupancies <- function(name, read_data, quality_groups, genome, gen_meth_pos) {  # improved 2018-03-08
  
  read_shifts <- start(read_data[["bam"]]) - 1
  
  occupancies <- list()
  counts <- list()
  converted <- list()
  counts_not_na <- list()
  
  cat("\ncalc_genomic_occupancies for", name, "\n")
  for(chr in names(genome)) {
    # cat(chr, "\n")
    
    occupancies[[chr]] <- list()
    counts[[chr]] <- list()
    converted[[chr]] <- list()
    counts_not_na[[chr]] <- list()
    
    reads_chr_good <- as.vector(quality_groups[["good"]] & (seqnames(read_data[["bam"]]) == chr))
    
    for(pat in names(read_data[["converted"]])) {

      counts[[chr]][[pat]] <- rep(NA, length(genome[[chr]]))
      counts_not_na[[chr]][[pat]] <- rep(NA, length(genome[[chr]]))
      converted[[chr]][[pat]] <- rep(NA, length(genome[[chr]]))
      occupancies[[chr]][[pat]] <- rep(NA, length(genome[[chr]]))
      
      genomic_pat_pos_on_reads <- mapply('+', read_data[["positions"]][[pat]][reads_chr_good], read_shifts[reads_chr_good])
      
      pos_counts_table <- table(unlist(genomic_pat_pos_on_reads))  # also counts reads that have NA in converted (i.e. not a G or C)
      
      pos_counts_not_na_table <- tapply(unlist(read_data[["converted"]][[pat]][reads_chr_good]), unlist(genomic_pat_pos_on_reads), function(vec) length(vec[!is.na(vec)]))
      
      pos_conversions_table <- tapply(unlist(read_data[["converted"]][[pat]][reads_chr_good]), unlist(genomic_pat_pos_on_reads), sum, na.rm=TRUE)
      
      pat_pos_on_chr <- as.integer(names(pos_conversions_table))
      
      if(!all( pat_pos_on_chr %in% unlist(gen_meth_pos[[chr]][[pat]]) )) {
        stop("Found unexpected methylation sites. This shouldn't happen!")
      }
      
      counts[[chr]][[pat]][pat_pos_on_chr] <- pos_counts_table
      counts_not_na[[chr]][[pat]][pat_pos_on_chr] <- pos_counts_not_na_table
      converted[[chr]][[pat]][pat_pos_on_chr] <- pos_conversions_table
      
      occupancies[[chr]][[pat]] <- converted[[chr]][[pat]] / counts_not_na[[chr]][[pat]]  # can be NA if counts_not_na[[chr]][[pat]] is zero
    }
  }
  return(list("occupancies"=occupancies, "counts"=counts, "counts_not_na"=counts_not_na, "converted"=converted))
}


plot_occupancy_histogram <- function(name, RData_folder, plot_folder="") {  # checked 2018-03-07
  
  load(file=paste0(RData_folder, "occupancies_raw_", name, ".RData"))
  
  png(paste0(plot_folder, "occupancy_hist_", name, ".png"), height = 1024, width = 1024, res = 144)
  par(mfrow=c(2, 2), oma=c(0, 0, 2, 0))
  for(pat in names(occupancies[[1]])) {
    tmp <- unlist( lapply(occupancies, '[[', pat) )
    hist(tmp, breaks = seq(0, 1, 0.05), xlim=c(0, 1),
         xlab="occupancy", main=pat,
         sub=paste0("mean: ", format(mean(tmp, na.rm=TRUE), digits=3)))
  }
  mtext(name, side=3, line=0, outer=TRUE, cex=1.5, font=1)
  dev.off()
}


plot_occupancies_in_region <- function(name, RData_folder, gene_table, genome, gen_meth_pos, context_patterns, annotations=list(), plot_folder="") {  # improved 2018-03-07
  
  load(file=paste0(RData_folder, "occupancies_raw_",name,".RData"))
  
  pats_meth <- c("HCH" = FALSE, "HCG" = grepl("CpG", name), "GCH" = grepl("GpC", name), "GCG" = TRUE)
  
  original_palette <- palette( c( ifelse(!pats_meth, 'gray90', colorschemes$Categorical.12[c(1,7,10,12)]),  # pat colours
                                   'gray90', colorschemes$Categorical.12[c(6,8)], "black") )  # strand colours
  # original_palette <- palette( c( colorschemes$Categorical.12[c(1,7,10,12)],  # pat colours
  #                                 'gray90', colorschemes$Categorical.12[c(6,8)], "black") )  # strand colours
  
  
  on.exit(palette(original_palette))
  
  for(gene_name in rownames(gene_table)) {
    
    chr <- gene_table[gene_name, "chr"]
    chr_seq_num <- even_fasta2num(genome[chr])
    
    png(paste0(plot_folder, "region_occupancies_", gene_name, "_", name, ".png"), height=1024, width=1024*2, pointsize = 18)
    
    par(mfrow=c(1, 2), oma=c(0, 0, 2, 0), mar=c(5, 4, 5.5, 2))
    
    for(colour_scheme in c("pattern", "strand") ) {
      
      for(j in 1:4) {
        
        pos_temp <- gene_table[gene_name, "xlim_start"]:gene_table[gene_name, "xlim_end"]
        pos_not_na <- pos_temp[!is.na( occupancies[[chr]][[j]][pos_temp] )]  # if occupancies is not NA, then counts_not_na is greater than zero
        
        if(!all(pos_not_na %in% unlist(gen_meth_pos[[chr]][[j]]))) {
          stop("Found unexpected methylation positions! This shouldn't happen!")
        }
        
        if(j == 1) {  # HCH
          
          plot(occupancies[[chr]][["HCH"]], col=ifelse(colour_scheme=="strand", 5, 1), pch=20,
            xlim=c(gene_table[gene_name, "xlim_start"], gene_table[gene_name, "xlim_end"]), ylim=c(0, 1), xlab=chr, ylab="occupancy")
          
        } else {  # not HCH

          col_to_use <- j
          
          # -- strand colour -- #
          if(colour_scheme == "strand") {
            col_to_use <- 1 + 4  # + 4 here and below to get to strand palette range
            if(pats_meth[j]) {
              col_to_use <- chr_seq_num[[chr]][pos_not_na] + 4
            }
          }
          
          converted_temp <- converted[[chr]][[j]][pos_not_na]
          counts_not_na_temp <- counts_not_na[[chr]][[j]][pos_not_na]
          occupancies_temp <- occupancies[[chr]][[j]][pos_not_na]
          
          # there has to be a better way to get correct error bars!
          tests <- mapply(binom.test, converted_temp, counts_not_na_temp)
          if(ncol(tests) >= 1) {
            segments(x0=pos_not_na, y0=sapply(1:ncol(tests), function (x) tests[,x][["conf.int"]][1]) ,x1=pos_not_na, y1=sapply( 1:ncol(tests), function (x) tests[,x][["conf.int"]][2]), col=col_to_use)
          }
          
          points(pos_not_na, occupancies_temp, col=col_to_use, pch=20)
        }
      }
      
      if(length(annotations) > 0) {
        for(anno_name in names(annotations)) {
          abline(v=annotations[[anno_name]][[chr]], lty=2)
          text(x=annotations[[anno_name]][[chr]], y=1+0.06, labels=anno_name, xpd=TRUE)
        }
      }
      
      if(colour_scheme == "pattern") {
        legend(x="top", xpd=TRUE, inset=c(0, -0.11), fill=c(1, (1:4)[pats_meth]), horiz=TRUE, legend=c(paste(context_patterns[!pats_meth], collapse=", "), context_patterns[pats_meth]))
        
      } else {
        legend(x="top", xpd=TRUE, inset=c(0, -0.11), fill=5:7, horiz=TRUE, legend=list(paste(context_patterns[!pats_meth], collapse=", "), "plus strand sites", "minus strand sites"))
      }
    }
    
    mtext(paste0(name, "  -  ", gene_name), side=3, line=0, outer=TRUE, cex=1.5, font=1)
    dev.off()
  }
}


get_BS_results_dataframe <- function(BS_run, num_digits = 2, min_site_cov = 0) {  # now using occs_df
  
  data_names <- sub(".RData", "", sub("occupancies_raw_", "", grep("occupancies_raw", dir(paste0("../", BS_run, "/data/RData_files")), value=TRUE)))
  
  results_list <- list()
  for(name in data_names) {
    cat(name, "\n")
    name_addition <- ifelse(grepl("__WD", name), "WD", ifelse(grepl("__ND", name), "ND", ""))
    file_name <- sub(paste0("__", name_addition), "", name)

    load(paste0("../", BS_run, "/data/RData_files/conversion_data_", name, ".RData"))
    conversion_summary <- cbind("reads" = sapply(quality_groups, sum), t(sapply(quality_groups, function(group) (sapply(conversion_ratio, function (x) mean(x[group], na.rm=TRUE))))))
    mio_reads <- sum(conversion_summary[, "reads"]) / 1000000
    mio_reads_HCH_converted <- sum(conversion_summary["good", "reads"]) / 1000000
    good_read_frac <- mio_reads_HCH_converted / mio_reads
    
    conversion_ratios <- sapply(conversion_ratio, mean, na.rm=TRUE)
    
    conversion_ratios_HCH_conv <- sapply(conversion_ratio, function (x) mean(x[quality_groups[["good"]]], na.rm=TRUE))
    
    occs_df <- load_BS_occs_df(run = BS_run, name = file_name)
    mean_counts <- mean(occs_df$counts)
    median_counts <- median(occs_df$counts)
    high_cov_sites <- occs_df$counts > min_site_cov
    high_cov_frac <- sum(high_cov_sites) / nrow(occs_df)
    mean_occ <- mean(occs_df$occ, na.rm=TRUE)
    sd_occ <- sd(occs_df$occ, na.rm=TRUE)
    mean_occ_high_cov <- mean(occs_df$occ[high_cov_sites], na.rm=TRUE)
    sd_occ_high_cov <- sd(occs_df$occ[high_cov_sites], na.rm=TRUE)
    
    occs_GCG_df <- subset(occs_df, pattern == "GCG")
    high_cov_sites <- occs_GCG_df$counts > min_site_cov
    mean_occ_GCG <- mean(occs_GCG_df$occ, na.rm=TRUE)
    sd_occ_GCG <- sd(occs_GCG_df$occ, na.rm=TRUE)
    mean_occ_GCG_high_cov <- mean(occs_GCG_df$occ[high_cov_sites], na.rm=TRUE)
    sd_occ_GCG_high_cov <- sd(occs_GCG_df$occ[high_cov_sites], na.rm=TRUE)

    results_list[[name]] <- data.frame(file_name = file_name, name_addition = name_addition, run = BS_run, mio_reads = mio_reads, mio_reads_HCH_converted = mio_reads_HCH_converted, good_read_frac = good_read_frac, 
                                       HCH_conv = conversion_ratios["HCH"], HCG_conv = conversion_ratios["HCG"], GCH_conv = conversion_ratios["GCH"], GCG_conv = conversion_ratios["GCG"],
                                       HCG_conv_good_reads = conversion_ratios_HCH_conv["HCG"], GCH_conv_good_reads = conversion_ratios_HCH_conv["GCH"], GCG_conv_good_reads = conversion_ratios_HCH_conv["GCG"],
                                       mean_counts = mean_counts, median_counts = median_counts, high_cov_frac = high_cov_frac, 
                                       mean_occ = mean_occ, sd_occ = sd_occ, mean_occ_high_cov = mean_occ_high_cov, sd_occ_high_cov = sd_occ_high_cov, 
                                       mean_occ_GCG = mean_occ_GCG, sd_occ_GCG = sd_occ_GCG, mean_occ_GCG_high_cov = mean_occ_GCG_high_cov, sd_occ_GCG_high_cov = sd_occ_GCG_high_cov, stringsAsFactors = FALSE)
    
    occs_df$chr <- sub("chr0", "chr", occs_df$chr)
    Np1_aligned_info_df <- calc_average_gene_occ_df(occs_df, pos_binning = 10, filter_size = 1, max_promoter_length = 1000, max_ORF_length = 1000)[["info_df"]]
    results_list[[name]] <- cbind(results_list[[name]], Np1_aligned_info_df)
  }
  
  results <- dplyr::bind_rows(results_list)
  
  return(results)
}


get_BS_results_601_dataframe <- function(BS_run, num_digits = 2, sites_601 = c(21:(93-20) , (144+20):(236-20)), sites_linker = (94+8):(143-8)) {
  
  data_names <- sub(".RData", "", sub("occupancies_raw_", "", grep("occupancies_raw", dir(paste0("../", BS_run, "/data/RData_files")), value=TRUE)))
  
  results_list <- list()
  for(name in data_names) {
    cat(name, "\n")
    name_addition <- ifelse(grepl("__WD", name), "WD", ifelse(grepl("__ND", name), "ND", ""))
    file_name <- sub(paste0("__", name_addition), "", name)
    
    load(paste0("../", BS_run, "/data/RData_files/occupancies_raw_", name, ".RData"))
    pats_meth <- c("HCH" = FALSE, "HCG" = grepl("CpG", name), "GCH" = grepl("GpC", name), "GCG" = TRUE)
    occ_values_601 <- unlist(lapply(lapply(occupancies, "[", pats_meth)[[1]], function(occs) occs[sites_601]))
    occ_values_linker <- unlist(lapply(lapply(occupancies, "[", pats_meth)[[1]], function(occs) occs[sites_linker]))
    mean_occ_601 <- mean(occ_values_601, na.rm=TRUE)
    mean_occ_linker <- mean(occ_values_linker, na.rm=TRUE)
    sd_occ_601 <- sd(occ_values_601, na.rm=TRUE)
    sd_occ_linker <- sd(occ_values_linker, na.rm=TRUE)
    
    occ_values_GCG_601 <- unlist(lapply(lapply(occupancies, "[", c(FALSE, FALSE, FALSE, TRUE))[[1]], function(occs) occs[sites_601]))
    occ_values_GCG_linker <- unlist(lapply(lapply(occupancies, "[", c(FALSE, FALSE, FALSE, TRUE))[[1]], function(occs) occs[sites_linker]))
    mean_occ_GCG_601 <- mean(occ_values_GCG_601, na.rm=TRUE)
    mean_occ_GCG_linker <- mean(occ_values_GCG_linker, na.rm=TRUE)
    sd_occ_GCG_601 <- sd(occ_values_GCG_601, na.rm=TRUE)
    sd_occ_GCG_linker <- sd(occ_values_GCG_linker, na.rm=TRUE)
    
    load(paste0("../", BS_run, "/data/RData_files/conversion_data_", name, ".RData"))
    load(paste0("../", BS_run, "/data/RData_files/read_data_filtered_", name, ".RData"))
    median_cov_601_linker <- median(unlist(coverage(read_data[["bam"]][quality_groups[["good"]]])["pFMP233_linker_with_cut_601s"]))
    
    results_list[[name]] <- data.frame(file_name = file_name, name_addition = name_addition, run = sub("_spike-in_linker", "", BS_run), 
                                       median_cov_601_linker = median_cov_601_linker, mean_occ_601 = mean_occ_601, mean_occ_linker = mean_occ_linker, sd_occ_601 = sd_occ_601, sd_occ_linker = sd_occ_linker,
                                       mean_occ_GCG_601 = mean_occ_GCG_601, mean_occ_GCG_linker = mean_occ_GCG_linker, sd_occ_GCG_601 = sd_occ_GCG_601, sd_occ_GCG_linker = sd_occ_GCG_linker, stringsAsFactors = FALSE)
  }
  
  results <- dplyr::bind_rows(results_list)
  
  return(results)
}


get_BS_results_BB_dataframe <- function(BS_run, num_digits = 2) {
  
  data_names <- sub(".RData", "", sub("occupancies_raw_", "", grep("occupancies_raw", dir(paste0("../", BS_run, "/data/RData_files")), value=TRUE)))
  
  results_list <- list()
  for(name in data_names) {
    cat(name, "\n")
    name_addition <- ifelse(grepl("__WD", name), "WD", ifelse(grepl("__ND", name), "ND", ""))
    file_name <- sub(paste0("__", name_addition), "", name)
    
    load(paste0("../", BS_run, "/data/RData_files/occupancies_raw_", name, ".RData"))
    pats_meth <- c("HCH" = FALSE, "HCG" = grepl("CpG", name), "GCH" = grepl("GpC", name), "GCG" = TRUE)
    occ_values <- unlist(lapply(occupancies, "[", pats_meth)[[1]])
    mean_occ <- mean(occ_values, na.rm=TRUE)
    sd_occ <- sd(occ_values, na.rm=TRUE)
    
    occ_values_GCG <- unlist(lapply(occupancies, "[", c(FALSE, FALSE, FALSE, TRUE))[[1]])
    mean_occ_GCG <- mean(occ_values_GCG, na.rm=TRUE)
    sd_occ_GCG <- sd(occ_values_GCG, na.rm=TRUE)
    
    load(paste0("../", BS_run, "/data/RData_files/conversion_data_", name, ".RData"))
    load(paste0("../", BS_run, "/data/RData_files/read_data_filtered_", name, ".RData"))
    median_cov_backbone <- median(unlist(coverage(read_data[["bam"]][quality_groups[["good"]]])["pFMP233_backbone"]))

    results_list[[name]] <- data.frame(file_name = file_name, name_addition = name_addition, run = sub("_spike-in_backbone", "", BS_run), 
                                       median_cov_backbone = median_cov_backbone, mean_occ_backbone = mean_occ, sd_occ_backbone = sd_occ,
                                       mean_occ_GCG_backbone = mean_occ_GCG, sd_occ_GCG_backbone = sd_occ_GCG, stringsAsFactors = FALSE)
  }
  
  results <- dplyr::bind_rows(results_list)
  
  return(results)
}


get_BS_results_REL606_dataframe <- function(BS_run, num_digits = 2) {
  
  data_names <- sub(".RData", "", sub("occupancies_raw_", "", grep("occupancies_raw", dir(paste0("../", BS_run, "/data/RData_files")), value=TRUE)))
  
  results_list <- list()
  for(name in data_names) {
    cat(name, "\n")
    name_addition <- ifelse(grepl("__WD", name), "WD", ifelse(grepl("__ND", name), "ND", ""))
    file_name <- sub(paste0("__", name_addition), "", name)
    
    load(paste0("../", BS_run, "/data/RData_files/occupancies_raw_", name, ".RData"))
    pats_meth <- c("HCH" = FALSE, "HCG" = grepl("CpG", name), "GCH" = grepl("GpC", name), "GCG" = TRUE)
    occ_values <- unlist(lapply(occupancies, "[", pats_meth)[[1]])
    mean_occ <- mean(occ_values, na.rm=TRUE)
    sd_occ <- sd(occ_values, na.rm=TRUE)
    
    occ_values_GCG <- unlist(lapply(occupancies, "[", c(FALSE, FALSE, FALSE, TRUE))[[1]])
    mean_occ_GCG <- mean(occ_values_GCG, na.rm=TRUE)
    sd_occ_GCG <- sd(occ_values_GCG, na.rm=TRUE)
    
    load(paste0("../", BS_run, "/data/RData_files/conversion_data_", name, ".RData"))
    load(paste0("../", BS_run, "/data/RData_files/read_data_filtered_", name, ".RData"))
    median_cov_REL606 <- median(unlist(coverage(read_data[["bam"]][quality_groups[["good"]]])["REL606"]))
    
    results_list[[name]] <- data.frame(file_name = file_name, name_addition = name_addition, run = sub("_REL606", "", BS_run), 
                                       median_cov_REL606 = median_cov_REL606, mean_occ_REL606 = mean_occ, sd_occ_REL606 = sd_occ,
                                       mean_occ_GCG_REL606 = mean_occ_GCG, sd_occ_GCG_REL606 = sd_occ_GCG, stringsAsFactors = FALSE)
  }
  
  results <- dplyr::bind_rows(results_list)
  
  return(results)
}


calc_occs_and_sites_chr_lists <- function(occs_df_list, occ_type = "occ_cut_uncut_2") {  # for RE occs_df
  
  occs = list()
  sites = list()
  
  for(name in names(occs_df_list)) {
    occs_df <- occs_df_list[[name]]
    occs[[name]] <- list()
    sites[[name]] <- list()
    
    chrs <- levels(occs_df$chr)[1:16]
    occs_df <- subset(occs_df, chr %in% chrs)
    occs_df_chr_list <- split(occs_df, as.character(occs_df$chr))
    
    for(chr in chrs) {
      occs[[name]][[chr]] <- occs_df_chr_list[[chr]][, occ_type]
      sites[[name]][[chr]] <- occs_df_chr_list[[chr]]$pos
    }
    
  }
  return(list(occs = occs, sites = sites))
}


calc_and_save_occs_df_from_chr_list <- function(run, name, WD_or_ND = "WD",
                                                genome = readDNAStringSet("../../reference_genomes/BS_scer_linker_backbone_ecoli/scer_mt_pFMP233_linker_backbone_REL606.fa")
                                                ) {  # saving the results probably should have been done this way from the beginning
  
  RData_file <- paste0("../", run, "/data/RData_files/occupancies_raw_", name, "__", WD_or_ND, ".RData")
  load(RData_file)
  
  occs_df_list <- list()
  
  pats_meth <- c("HCH" = FALSE, "HCG" = grepl("CpG", name), "GCH" = grepl("GpC", name), "GCG" = TRUE)
  pats_meth_to_use <- names(pats_meth)[pats_meth]
  
  for(chr in names(occupancies)) {
    for(pat in pats_meth_to_use) {
      sites_not_na <- which(!is.na(occupancies[[chr]][[pat]]))
      nucleotide_num <- even_fasta2num(genome[chr])[[1]][sites_not_na]  # nucleotides on plus strand at site positions with "A"=1, "C"=2, "G"=3, "T"=4, "N"=NA
      if(!all(nucleotide_num %in% c(2, 3, NA))) stop("Found occ entry at non-C site!")
      occs_df_list[[paste0(chr, "_", pat)]] <- data.frame(chr = chr, pos = sites_not_na, pattern = pat,
                                                          strand = ifelse(nucleotide_num == 2, "+", "-"),
                                                          counts = counts_not_na[[chr]][[pat]][sites_not_na],
                                                          converted = converted[[chr]][[pat]][sites_not_na],
                                                          occ = occupancies[[chr]][[pat]][sites_not_na], stringsAsFactors = FALSE)
    }
  }
  
  occs_df <- bind_rows(occs_df_list)
  occs_df <- occs_df[order(occs_df$chr, occs_df$pos), ]
  row.names(occs_df) <- NULL
  
  # calc statistical error of occupancies (this still has to be checked properly)
  # occs_df <- calc_occ_conf_intervals(occs_df)

  save(occs_df, file = paste0("../", run, "/data/RData_files/occs_df_", name, "__", WD_or_ND, ".RData"))
  
  return(occs_df)
}


calc_occ_conf_intervals <- function(occs_df) {
  
  tests <- mapply(binom.test, occs_df$converted, occs_df$counts)
  
  occs_df$occ_lower_bound <- sapply(1:ncol(tests), function (x) tests[,x][["conf.int"]][1])
  occs_df$occ_upper_bound <- sapply(1:ncol(tests), function (x) tests[,x][["conf.int"]][2])
  
  return(occs_df)
}


load_BS_occs_df <- function(run, name, WD_or_ND = "WD") {

  RData_file <- paste0("../", run, "/data/RData_files/occs_df_", name, "__", WD_or_ND, ".RData")
  if(file.exists(RData_file)) {
    load(RData_file)
  } else {
    occs_df <- calc_and_save_occs_df_from_chr_list(run, name, WD_or_ND = WD_or_ND)
  }
  return(occs_df)
}

