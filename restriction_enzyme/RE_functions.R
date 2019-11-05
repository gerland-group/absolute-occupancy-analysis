

#' save_PE_bam_as_chr_df_list
#'
#' @param name name of file (without extension)
#' @param bam_folder location of file (default: '')
#' @param extension file extension (default: '.bam')
#' @param bam_duplicates_flag should the duplicate reads be kept instead of the unique. NA keeps both. BWA does not flag duplicates, so this option doesn't do anything!
#'
#' Saves the bam read positions in a .RData file of the same name. The object is called bam.
#' 
#' @export
#'
#' @examples
save_PE_bam_as_chr_df_list <- function (name, bam_folder = "", extension = ".bam", RData_folder = "", bam_duplicates_flag = NA) {  # use this for paired end sequencing
  
  bam_file <- BamFile(file.path(bam_folder, paste0(name, extension)), index=file.path(bam_folder, name))
  open(bam_file)
  
  flag_sel <- scanBamFlag(isDuplicate=bam_duplicates_flag)
  fragments <- readGAlignmentPairs(bam_file, param=ScanBamParam(what="seq", flag=flag_sel))
  close(bam_file)
  
  if(isTRUE(bam_duplicates_flag)) {
    name <- paste0(name, "_bam_duplicates_only")
  } else if(isTRUE(!bam_duplicates_flag)) {
    name <- paste0(name, "_no_bam_duplicates")
  }
  
  # Filter out pairs that aren't on the same chromosome
  fragments <- fragments[!is.na(seqnames(fragments))]
  
  # Create Granges objects
  gen_rang <- granges(fragments)
  
  chr_list_table <- list()
  
  for(chr_name in levels(seqnames(gen_rang))) {
    tmp <- subset(gen_rang, seqnames(.(gen_rang)) == .(chr_name))
    chr_list_table[[chr_name]] <- data.frame("start"=start(tmp), "end"=end(tmp)+1, "strand"=strand(tmp))
  }
  
  bam <- chr_list_table
  save(bam, file = paste0(RData_folder, gsub("%", "", name), ".RData"), compress = TRUE)  # remove % from name because it causes problems later
}


#' count_fragment_starts_at_cut_sites
#'
#' @param chr_fragment_starts vector of all fragment starts
#' @param cut_sites RE recognition site
#' @param chr_length chromosome length
#' @param cut_site_shift shift between RE recognition site and where the DNA is cut
#' @param count_window count the fragment ends in this region
#'
#' fragment starts are defined as the fragment end with the lower genomic position
#'
#' @return vector of the fragment counts per RE recognition site
#' @export
#'
#' @examples
count_fragment_starts_at_cut_sites <- function(chr_fragment_starts, cut_sites, chr_length, cut_site_shift, count_window) {  # checked 2018-02-25
  
  starts_vec <- calc_genomic_count_vec(chr_fragment_starts, chr_length)
  
  return(sapply(cut_sites - cut_site_shift, function (i) sum(starts_vec[i + count_window])))
}


#' count_fragment_ends_at_cut_sites
#'
#' @seealso count_fragment_starts_at_cut_sites
#'
count_fragment_ends_at_cut_sites <- function(chr_fragment_ends, cut_sites, chr_length, cut_site_shift, count_window) {  # checked 2018-02-25
  
  return(count_fragment_starts_at_cut_sites(chr_fragment_ends, cut_sites, chr_length, -cut_site_shift, -count_window))
}


count_uncut_fragments_at_cut_sites <- function(bam_chr, cut_sites, chr_length, cut_site_shift) {
  
  # count_covering_fragments_old <- function(cut_site) {
  #   return(sum(bam_chr$start < cut_site - cut_site_shift & bam_chr$end > cut_site + cut_site_shift))
  # }
  # 
  # return(sapply(cut_sites, count_covering_fragments_old))
  
  bam_chr <- bam_chr[order(bam_chr$start, bam_chr$end), ]
  
  max_indeces_to_test <- findInterval(cut_sites - cut_site_shift - 1, bam_chr$start)  # test all starts for all sites
  max_length <- max(bam_chr$end - bam_chr$start) + 1
  
  min_indeces_to_test <- findInterval(cut_sites - cut_site_shift - 1 - max_length, bam_chr$start)  # gives all starts that could have valid end locations
  
  count_covering_fragments <- function(i) {
    indeces <- min_indeces_to_test[i]:max_indeces_to_test[i]
    return(sum(bam_chr$end[indeces] > cut_sites[i] + cut_site_shift))  # now also test all ends for each site
  }
  
  return(sapply(1:length(cut_sites), count_covering_fragments))
}


calc_genomic_count_vec <- function(location_vec, chr_length) {
  
  location_table <- table(location_vec)
  genomic_count_vec <- rep(0, chr_length)
  genomic_count_vec[as.numeric(names(location_table))] <- location_table
  return(genomic_count_vec)
}


count_fragments_in_background <- function(chr_fragment_starts, bg_bool_vec) {
  
  starts_vec <- calc_genomic_count_vec(chr_fragment_starts, length(bg_bool_vec))
  
  return(sum(starts_vec[bg_bool_vec], na.rm=TRUE))
}


calc_coverage_in_background <- function(bam_chr, bg_bool_vec) {  # counting fragments ending at, starting at and spanning over the boundary between base x-1 and x
  
  starts_vec <- calc_genomic_count_vec(bam_chr$start, length(bg_bool_vec)+1)
  ends_vec <- calc_genomic_count_vec(bam_chr$end, length(bg_bool_vec)+1)  # "+1" needed due to end shift of +1
  coverage_vec <- cumsum(starts_vec) - cumsum(ends_vec)
  
  return(sum(coverage_vec[bg_bool_vec]))
}


truncate_to_valid_index <- function(x, vec_length) {  # checked 2018-02-25
  return(remove_negative_numbers(x[x <= vec_length]))
}


remove_negative_numbers <- function (x) {  # checked 2018-02-25
  return(x[x >= 0])
}


calc_cut_sites <- function(data_names, genome, sample_to_RE_table, RE_info, bad_regions, large_window_limit=200, dist_ign_site=0, dist_ign_half=0, filter_type = 'either') {  # new 2018-02-24
  
  # filter_type can be 'either', 'start' or 'end':
  # 'either' filters out sites that have a close neighbor (d < max(dist_ign_site, dist_ign_half)) in any direction
  # 'start'  filters out sites that have a close neighbor (d < dist_ign_site) in any direction 
  #          AND filters out sites that have a close downstream neighbor (d < dist_ign_half)
  # 'end'    filters out sites that have a close neighbor (d < dist_ign_site) in any direction 
  #          AND filters out sites that have a close upstream neighbor (d < dist_ign_half)
  
  if(!(filter_type %in% c("either", "start", "end"))) {
    stop("invalid filter_type")
  }
  
  sites <- list()
  
  for(name in data_names) {
    
    sites[[name]] <- list()
    
    all_enzymes <- rownames(RE_info)[sapply(rownames(RE_info), grepl, name)]
    
    enzymes <- strsplit(sample_to_RE_table[name, 1], "+", fixed = TRUE)[[1]]
    if(!all(enzymes %in% rownames(RE_info))) {
      stop("enzyme not listed")
    }
    
    sites_temp_list <- list()
    for(chr in names(genome)) {
      
      sites_temp <- lapply(all_enzymes, function(enzyme) RE_info[enzyme, "middle_shift"] + start(matchPattern(RE_info[enzyme, "recognition_site"], genome[[chr]])))
      names(sites_temp) <- paste0(all_enzymes, ".")
      
      sites_temp <- sort(unlist(sites_temp))
      
      if(length(sites_temp) > 0) {
        
        # filter out sites that have a close neighbor (d < dist_ign_site) in any direction
        close_sites_1 <- calc_close_RE_sites(list(chr=sites_temp), close_distance=dist_ign_site, list(chr=width(genome[chr])))
        
        # filter out sites that have a close downstream/upstream neighbor (d < dist_ign_half)
        close_sites_2 <- calc_close_RE_sites(list(chr=sites_temp), close_distance=dist_ign_half, list(chr=width(genome[chr])))
        
        sites_temp <- sites_temp[!(close_sites_1[["either"]][['chr']] | close_sites_2[[filter_type]][['chr']])]
        
        # now exclude enzymes that occur in name, but are not in sample_to_RE_table[name, 1]
        sites_temp <- sites_temp[sapply(strsplit(names(sites_temp), ".", fixed = TRUE), "[", 1) %in% enzymes]
        
        # filter those too close to the chromsome start/end
        sites_temp <- sites_temp[ sites_temp > large_window_limit & sites_temp < (width(genome[chr]) - large_window_limit) ]
        
        # filter those in bad regions
        for(bad_index in which(bad_regions$chr == chr)) {
          sites_temp <- sites_temp[ (sites_temp < bad_regions$start[bad_index] - large_window_limit) | (sites_temp > bad_regions$end[bad_index] + large_window_limit) ]
        }
        sites_temp_list[[chr]] <- sites_temp
        
      }
    }
    
    for(enzyme in enzymes) {
      sites[[name]][[enzyme]] <- lapply(sites_temp_list, function(sites_temp) sites_temp[grepl(enzyme, names(sites_temp))])
    }
  }
  
  return(sites)
}


count_fragments_at_cut_sites <- function(data_names, genome, chr_name_list, sample_to_RE_table, RE_info, count_window_df, max_length, bad_regions, window_limit_factor, 
                                         use_both_strands = TRUE, bg_exclusion_window = -300:300) {
  
  cut_counts_df_list <- list()  # will contain dataframes with chr, site position and count columns per site
  bg_counts_df_list <- list()  # will contain dataframes with background counts and sizes per genome
  fragment_counts_list <- list()
  count_window_list <- list()
  
  for(name in data_names) {
    
    load(paste0(RData_folder, name, ".RData"))
    
    # use the count region of the X sample and the sample genome for both samples and genomes
    name_X <- sub("1$", "X", name)
    
    # ALL sites
    sites <- calc_cut_sites(list(name), genome, sample_to_RE_table, RE_info, bad_regions, 
                            large_window_limit=max(subset(count_window_df, genome_type == "main" & name == name_X)[, "count_window_limit"]))[[1]]
    
    cut_counts_df_enzyme_list <- list()
    
    fragment_counts_list[[name]] <- list()
    count_window_list[[name]] <- list()
    
    for(enzyme in names(sites)) {
      
      count_window <- 0:round(window_limit_factor * subset(count_window_df[count_window_df$enzyme == enzyme, ], genome_type == "main" & name == name_X)[, "count_window_limit"])
      count_window_list[[name]][[enzyme]] <- count_window
      cat(name, ", enzyme = ", enzyme, ", count window: ", min(count_window), " to ", max(count_window), "\n", sep = "")
      
      chrs_vec <- do.call(c, lapply(names(sites[[enzyme]]), function(chr) rep(chr, times=length(sites[[enzyme]][[chr]]))))
      
      cut_counts_df_enzyme_list[[enzyme]] <- data.frame(chr = chrs_vec, enzyme = enzyme, pos = unlist(sites[[enzyme]]), row.names = names(unlist(sites[[enzyme]])))
      
      cut_starts_plus_list <- list()
      cut_starts_minus_list <- list()
      cut_ends_plus_list <- list()
      cut_ends_minus_list <- list()
      
      uncut_plus_list <- list()
      uncut_minus_list <- list()
      
      for(chr in names(sites[[enzyme]])) {
        
        sites_chr <- sites[[enzyme]][[chr]]
        
        chr_bam <- subset(bam[[chr]], end - start + 1 <= max_length)
        
        fragment_counts_list[[name]][[chr]] <- nrow(chr_bam)  # independent of the enzyme
        
        if(!use_both_strands) {
          cut_starts_plus_list[[chr]] <- count_fragment_starts_at_cut_sites(chr_bam$start[chr_bam$strand=="+"], sites_chr, width(genome[chr]), RE_info[enzyme, "site_shift"], count_window)
          cut_ends_minus_list[[chr]] <- count_fragment_ends_at_cut_sites(chr_bam$end[chr_bam$strand=="-"], sites_chr, width(genome[chr]), RE_info[enzyme, "site_shift"], count_window)
        } else {  # default case
          cut_starts_plus_list[[chr]] <- count_fragment_starts_at_cut_sites(chr_bam$start[chr_bam$strand=="+"], sites_chr, width(genome[chr]), RE_info[enzyme, "site_shift"], count_window)
          cut_starts_minus_list[[chr]] <- count_fragment_starts_at_cut_sites(chr_bam$start[chr_bam$strand=="-"], sites_chr, width(genome[chr]), RE_info[enzyme, "site_shift"], count_window)
          cut_ends_plus_list[[chr]] <- count_fragment_ends_at_cut_sites(chr_bam$end[chr_bam$strand=="+"], sites_chr, width(genome[chr]), RE_info[enzyme, "site_shift"], count_window)
          cut_ends_minus_list[[chr]] <- count_fragment_ends_at_cut_sites(chr_bam$end[chr_bam$strand=="-"], sites_chr, width(genome[chr]), RE_info[enzyme, "site_shift"], count_window)
        }
        
        uncut_plus_list[[chr]] <- count_uncut_fragments_at_cut_sites(chr_bam[chr_bam$strand=="+", ], sites_chr, width(genome[chr]), RE_info[enzyme, "site_shift"])
        uncut_minus_list[[chr]] <- count_uncut_fragments_at_cut_sites(chr_bam[chr_bam$strand=="-", ], sites_chr, width(genome[chr]), RE_info[enzyme, "site_shift"])
      }
      
      if(!use_both_strands) {
        cut_counts_df_enzyme_list[[enzyme]] <- cbind(cut_counts_df_enzyme_list[[enzyme]], data.frame(starts_plus = unlist(cut_starts_plus_list), ends_minus = unlist(cut_ends_minus_list), 
                                                                                                     uncut_plus = unlist(uncut_plus_list), uncut_minus = unlist(uncut_minus_list)))
      } else {
        cut_counts_df_enzyme_list[[enzyme]] <- cbind(cut_counts_df_enzyme_list[[enzyme]], data.frame(starts_plus = unlist(cut_starts_plus_list), starts_minus = unlist(cut_starts_minus_list), 
                                                                                                     ends_plus = unlist(cut_ends_plus_list), ends_minus = unlist(cut_ends_minus_list),
                                                                                                     uncut_plus = unlist(uncut_plus_list), uncut_minus = unlist(uncut_minus_list)))
      }
    }
    
    cut_counts_df_list[[name]] <- do.call(rbind, cut_counts_df_enzyme_list)
    row.names(cut_counts_df_list[[name]]) <- do.call(c, lapply(cut_counts_df_enzyme_list, rownames))
    cut_counts_df_list[[name]] <- cut_counts_df_list[[name]][order(cut_counts_df_list[[name]]$chr, cut_counts_df_list[[name]]$pos), ]
    
    
    # Preparation for background treatment 
    # Take the read starts and ends from regions outside the windows around cut sites
    
    bg_starts_plus_list <- list()
    bg_starts_minus_list <- list()
    bg_ends_plus_list <- list()
    bg_ends_minus_list <- list()
    bg_coverage_plus_list <- list()
    bg_coverage_minus_list <- list()
    
    bg_starts_size <- list()
    bg_ends_size <- list()
    bg_coverage_size <- list()
    
    for(genome_type in names(chr_name_list)) {
      
      bg_starts_plus_list[[genome_type]] <- 0
      bg_starts_minus_list[[genome_type]] <- 0
      bg_ends_plus_list[[genome_type]] <- 0
      bg_ends_minus_list[[genome_type]] <- 0
      bg_coverage_plus_list[[genome_type]] <- 0
      bg_coverage_minus_list[[genome_type]] <- 0
      
      bg_starts_size[[genome_type]] <- 0
      bg_ends_size[[genome_type]] <- 0
      bg_coverage_size[[genome_type]] <- 0
      
      for(chr in chr_name_list[[genome_type]]) {
        
        chr_bam <- bam[[chr]]
        
        bg_bool_vec <- rep(TRUE, width(genome[chr]))
        
        # bad regions
        for(bad_index in which(bad_regions$chr == chr)) {
          bg_bool_vec[ bad_regions$start[bad_index]:bad_regions$end[bad_index] ] <- FALSE
        }
        bg_bool_vec_starts <- bg_bool_vec  # for starts
        bg_bool_vec_ends <- bg_bool_vec  # for ends
        
        # unfiltered cutsites
        # exclude the neighbourhood of all sites of all enzymes that occur in name (even though not all sites of all enzymes might have been counted before)
        all_enzymes <- rownames(RE_info)[sapply(rownames(RE_info), grepl, name)]
        sites_chr_all <- lapply(all_enzymes, function(enzyme) RE_info[enzyme, "middle_shift"] + start(matchPattern(RE_info[enzyme, "recognition_site"], genome[[chr]])))
        sites_chr_all <- sort(unlist(sites_chr_all))
        sites_chr_all_starts <- lapply(all_enzymes, function(enzyme) RE_info[enzyme, "middle_shift"] - RE_info[enzyme,"site_shift"] + start(matchPattern(RE_info[enzyme, "recognition_site"], genome[[chr]])))
        sites_chr_all_starts <- sort(unlist(sites_chr_all_starts))
        sites_chr_all_ends <- lapply(all_enzymes, function(enzyme) RE_info[enzyme, "middle_shift"] + RE_info[enzyme,"site_shift"] + start(matchPattern(RE_info[enzyme, "recognition_site"], genome[[chr]])))
        sites_chr_all_ends <- sort(unlist(sites_chr_all_ends))
        
        # do not use positions within bg_exclusion_window of any site
        for(x in bg_exclusion_window) {
          bg_bool_vec_starts[ truncate_to_valid_index(sites_chr_all_starts + x, length(bg_bool_vec)) ] <- FALSE
          bg_bool_vec_ends[ truncate_to_valid_index(sites_chr_all_ends - x, length(bg_bool_vec)) ] <- FALSE
          bg_bool_vec[ truncate_to_valid_index(sites_chr_all + x, length(bg_bool_vec)) ] <- FALSE
        }
        
        if(!use_both_strands) {
          bg_starts_plus_list[[genome_type]] <- bg_starts_plus_list[[genome_type]] + count_fragments_in_background(chr_bam$start[chr_bam$strand=="+"], bg_bool_vec_starts)
          bg_ends_minus_list[[genome_type]] <- bg_ends_minus_list[[genome_type]] + count_fragments_in_background(chr_bam$end[chr_bam$strand=="-"], bg_bool_vec_ends)
        } else {
          bg_starts_plus_list[[genome_type]] <- bg_starts_plus_list[[genome_type]] + count_fragments_in_background(chr_bam$start[chr_bam$strand=="+"], bg_bool_vec_starts)
          bg_starts_minus_list[[genome_type]] <- bg_starts_minus_list[[genome_type]] + count_fragments_in_background(chr_bam$start[chr_bam$strand=="-"], bg_bool_vec_starts)
          bg_ends_plus_list[[genome_type]] <- bg_ends_plus_list[[genome_type]] + count_fragments_in_background(chr_bam$end[chr_bam$strand=="+"], bg_bool_vec_ends)
          bg_ends_minus_list[[genome_type]] <- bg_ends_minus_list[[genome_type]] + count_fragments_in_background(chr_bam$end[chr_bam$strand=="-"], bg_bool_vec_ends)
        }
        
        bg_coverage_plus_list[[genome_type]] <- bg_coverage_plus_list[[genome_type]] + calc_coverage_in_background(chr_bam[chr_bam$strand=="+", ], bg_bool_vec)
        bg_coverage_minus_list[[genome_type]] <- bg_coverage_minus_list[[genome_type]] + calc_coverage_in_background(chr_bam[chr_bam$strand=="-", ], bg_bool_vec)
        
        bg_starts_size[[genome_type]] <- bg_starts_size[[genome_type]] + sum(bg_bool_vec_starts, na.rm=TRUE)
        bg_ends_size[[genome_type]] <- bg_ends_size[[genome_type]] + sum(bg_bool_vec_ends, na.rm=TRUE)
        bg_coverage_size[[genome_type]] <- bg_coverage_size[[genome_type]] + sum(bg_bool_vec, na.rm=TRUE)
      }
    }
    
    bg_df <- data.frame(genome = names(chr_name_list), starts_size = unlist(bg_starts_size), ends_size = unlist(bg_ends_size), coverage_size = unlist(bg_coverage_size),
                        starts_plus = unlist(bg_starts_plus_list), starts_minus = unlist(bg_starts_minus_list), 
                        ends_plus = unlist(bg_ends_plus_list), ends_minus = unlist(bg_ends_minus_list), 
                        coverage_sum_plus = unlist(bg_coverage_plus_list), coverage_sum_minus = unlist(bg_coverage_minus_list), row.names = NULL)
    
    bg_df$starts_plus_per_bp <- with(bg_df, starts_plus / starts_size)
    bg_df$starts_minus_per_bp <- with(bg_df, starts_minus / starts_size)
    bg_df$ends_plus_per_bp <- with(bg_df, ends_plus / ends_size)
    bg_df$ends_minus_per_bp <- with(bg_df, ends_minus / ends_size)
    bg_df$coverage_plus <- with(bg_df, coverage_sum_plus / coverage_size)  # here a mean "per window" doesn't make sence
    bg_df$coverage_minus <- with(bg_df, coverage_sum_minus / coverage_size)
    genome_sizes <- c(sum(width(genome[chr_name_list[["main"]]])), sum(width(genome[chr_name_list[["norm"]]])))
    bg_df$size_fraction <- with(bg_df, (starts_size + ends_size) / 2 / genome_sizes)
    
    if(bg_df$size_fraction[1] < 0.1) {
      warning(paste0(name, ": background calculated from less than 25% of the sample genome: ", bg_df$size_fraction[1]))
    }
    
    bg_counts_df_list[[name]] <- bg_df
    
  }
  
  return( list( "cut_counts_nf_df_list" = cut_counts_df_list, "bg_counts_df_list" = bg_counts_df_list, "fragment_counts_list" = fragment_counts_list, "count_window_list" = count_window_list) )
}


set_close_site_cut_counts_NA <- function(cut_counts_df_list, genome, sample_to_RE_table, RE_info, bad_regions, close_distances) {
  
  dist_ign_site <- close_distances[["dist_ign_site"]]
  dist_ign_half <- close_distances[["dist_ign_half"]]
  
  for(name in names(cut_counts_df_list)) {
    df <- cut_counts_df_list[[name]]
    
    sites_starts <- calc_cut_sites(list(name), genome, sample_to_RE_table, RE_info, bad_regions, dist_ign_site = dist_ign_site, dist_ign_half = dist_ign_half, filter_type = 'start')
    sites_starts <- lapply(sites_starts[[1]], unlist)
    names(sites_starts) <- NULL
    sites_starts <- unlist(sites_starts)
    
    sites_ends <- calc_cut_sites(list(name), genome, sample_to_RE_table, RE_info, bad_regions, dist_ign_site = dist_ign_site, dist_ign_half = dist_ign_half, filter_type = 'end')
    sites_ends <- lapply(sites_ends[[1]], unlist)
    names(sites_ends) <- NULL
    sites_ends <- unlist(sites_ends)
    
    sites_uncut <- calc_cut_sites(list(name), genome, sample_to_RE_table, RE_info, bad_regions, dist_ign_site = dist_ign_site, dist_ign_half = 0)
    sites_uncut <- lapply(sites_uncut[[1]], unlist)
    names(sites_uncut) <- NULL
    sites_uncut <- unlist(sites_uncut)
    
    df[setdiff(row.names(df), names(sites_starts)), c("starts_plus", "starts_minus")] <- NA
    df[setdiff(row.names(df), names(sites_ends)), c("ends_plus", "ends_minus")] <- NA
    df[setdiff(row.names(df), names(sites_uncut)), c("uncut_plus", "uncut_minus")] <- NA
    
    cut_counts_df_list[[name]] <- df
  }
  
  return(cut_counts_df_list)
}


helper_print_table <- function(table, file="", append=TRUE) {
  
  require(knitr)
  
  if(!append) {
    cat( "", file=file, append=FALSE)
  }
  
  cat(knitr::kable(table), sep="\n", file=file, append=TRUE)
  cat("\n", file=file, append=TRUE)
  write.table(table, file=file, sep="\t", quote=FALSE, append=TRUE)
}


remove_duplicate_fragments <- function(name, RData_folder) {
  load(paste0(RData_folder, name, ".RData"))
  
  for(chr in names(bam)) {
    bam[[chr]] <- plyr::count(bam[[chr]])[, 1:3 ]
  }
  
  save(bam, file=paste0(RData_folder, name, "_no_duplicates", ".RData"))
}


calc_close_RE_sites <- function(sites_list, close_distance, chr_sizes) {  # checked 22-05-2018
  # gives logical vectors for each chromosome, saying whether or not a site has a close neighbor for starting/ending reads or either
  
  close_gaps <- sapply(names(sites_list), function(s) diff(c(0, sites_list[[s]], chr_sizes[[s]])) <= close_distance, simplify=FALSE)
  
  close_start <- lapply(close_gaps, function(x) x[-1])
  close_end <- lapply(close_gaps, function(x) x[-length(x)])
  
  close_either <- mapply('|', close_start, close_end, SIMPLIFY=FALSE)  # close neighbour site in either direction
  
  return(list('either' = close_either, 'start' = close_start, 'end' = close_end))
}


set_high_outliers_to_quantile <- function(vec, outlier_quantile) {
  quantile_value <- quantile(vec, outlier_quantile, na.rm = TRUE)
  vec[vec > quantile_value] <- quantile_value
  
  return(vec)
}


plot_count_comparisons <- function(cut_counts_df_list, chr_name_list, genome_names, plot_folder, outlier_quantile = 0.995, name_add = "") {
  
  helper_plot <- function(...) {
    plot(..., pch=20, col=adjustcolor(palette()[1], alpha.f=min(0.2, 0.2*5000/length(list(...)[[1]]))))
  }
  
  plot_folder = paste0(plot_folder, "count_comparisons", name_add, "/")
  dir.create(plot_folder)
  
  for(name_X in grep("_X$", names(cut_counts_df_list), value = TRUE)) {
    
    name <- sub("_X$", "", name_X)
    name_1 <- paste0(name, "_1")
    
    for(genome_type in names(chr_name_list)) {
      
      png(paste0(plot_folder, "count_comparisons_", name, "_", genome_type, ".png"), res=144, width=3500, height=2000)
      par(mfrow=c(4, 7), oma=c(2, 0, 3, 0))
      
      chrs <- chr_name_list[[genome_type]]
      
      df_X <- subset(cut_counts_df_list[[name_X]], chr %in% chrs)
      df_1 <- subset(cut_counts_df_list[[name_1]], chr %in% chrs)
      
      for(A in 1:8) {
        for(B in A:8) {
          if(A != B) {
            
            name_A <- c(rep(name_1, 4), rep(name_X, 4))[A]
            name_B <- c(rep(name_1, 4), rep(name_X, 4))[B]
            
            type_A <- rep(c("starts_plus", "starts_minus", "ends_plus", "ends_minus"), 2)[A]
            type_B <- rep(c("starts_plus", "starts_minus", "ends_plus", "ends_minus"), 2)[B]
            
            name_add_A <- ifelse(name_A == name_X, "X", "1")
            name_add_B <- ifelse(name_B == name_X, "X", "1")
            
            if(name_add_A == "X") {
              counts_A <- df_X[, type_A]
            } else {
              counts_A <- df_1[, type_A]
            }
            
            if(name_add_B == "X") {
              counts_B <- df_X[, type_B]
            } else {
              counts_B <- df_1[, type_B]
            }
            
            counts_A <- set_high_outliers_to_quantile(counts_A, outlier_quantile)
            counts_B <- set_high_outliers_to_quantile(counts_B, outlier_quantile)
            
            counts_not_na <- !is.na(counts_A) & !is.na(counts_B)
            counts_A <- counts_A[counts_not_na]
            counts_B <- counts_B[counts_not_na]
            
            helper_plot(counts_A, counts_B, xlab=paste0(type_A, "_", name_add_A), ylab=paste0(type_B, "_", name_add_B))
          }
        }
      }
      
      mtext(paste0(name, ", ", genome_names[genome_type]), side=3, line=0, outer=TRUE, cex=1.5, font=1)
      dev.off()
    }
  }
}


plot_cut_uncut_comparisons <- function(cut_counts_df_list, chr_name_list, genome_names, plot_folder, outlier_quantile = 0.995, name_add = "") {
  
  helper_plot <- function(...) {
    plot(..., pch=20, col=adjustcolor(palette()[1], alpha.f=min(0.2, 0.2*5000/length(list(...)[[1]]))))
  }
  
  plot_folder = paste0(plot_folder, "cut_uncut_comparisons", name_add, "/")
  dir.create(plot_folder)
  
  for(name_X in grep("_X$", names(cut_counts_df_list), value = TRUE)) {
    
    name <- sub("_X$", "", name_X)
    name_1 <- paste0(name, "_1")
    
    for(genome_type in names(chr_name_list)) {
      
      png(paste0(plot_folder, "cut_uncut_comparisons_", name, "_", genome_type, ".png"), res=144, width=3500, height=2000)
      par(mfrow=c(4, 8), oma=c(2, 0, 3, 0))
      
      chrs <- chr_name_list[[genome_type]]
      
      df_X <- subset(cut_counts_df_list[[name_X]], chr %in% chrs)
      df_1 <- subset(cut_counts_df_list[[name_1]], chr %in% chrs)
      
      for(A in 1:8) {
        for(B in 1:4) {
          
          name_A <- c(rep(name_1, 4), rep(name_X, 4))[A]
          name_B <- c(rep(name_1, 2), rep(name_X, 2))[B]
          
          type_A <- rep(c("starts_plus", "starts_minus", "ends_plus", "ends_minus"), 2)[A]
          type_B <- rep(c("uncut_plus", "uncut_minus"), 2)[B]
          
          name_add_A <- ifelse(name_A == name_X, "X", "1")
          name_add_B <- ifelse(name_B == name_X, "X", "1")
          
          if(name_add_A == "X") {
            counts_A <- df_X[, type_A]
          } else {
            counts_A <- df_1[, type_A]
          }
          
          if(name_add_B == "X") {
            counts_B <- df_X[, type_B]
          } else {
            counts_B <- df_1[, type_B]
          }
          
          counts_A <- set_high_outliers_to_quantile(counts_A, outlier_quantile)
          counts_B <- set_high_outliers_to_quantile(counts_B, outlier_quantile)
          
          counts_not_na <- !is.na(counts_A) & !is.na(counts_B)
          counts_A <- counts_A[counts_not_na]
          counts_B <- counts_B[counts_not_na]
          
          helper_plot(counts_A, counts_B, xlab=paste0(type_A, "_", name_add_A), ylab=paste0(type_B, "_", name_add_B))
        }
      }
      
      mtext(paste0(name, ", ", genome_names[genome_type]), side=3, line=0, outer=TRUE, cex=1.5, font=1)
      dev.off()
    }
  }
}


plot_counts_vs_nn_distance <- function(cut_counts_df_list, chr_name_list, chr_sizes, genome_names, plot_folder = window_folder, outlier_quantile=0.99, max_dist=600) {
  
  plot_folder = paste0(plot_folder, "counts_vs_nn_distance/")
  dir.create(plot_folder)
  
  for(name_X in grep("_X$", names(cut_counts_df_list), value = TRUE)) {
    
    name <- sub("_X$", "", name_X)
    name_1 <- paste0(name, "_1")
    
    for(genome_type in names(chr_name_list)) {
      
      png(paste0(plot_folder, "cut_counts_vs_nn_distance_", name, "_", genome_type, ".png"), res=100, width=2000, height=1400)
      par(mfrow=c(4, 6), oma=c(2, 0, 3, 0))
      
      chrs <- chr_name_list[[genome_type]]
      
      df_X <- subset(cut_counts_df_list[[name_X]], chr %in% chrs)
      df_1 <- subset(cut_counts_df_list[[name_1]], chr %in% chrs)
      
      nn_dist_upstream <- unlist(sapply(chrs, function(chr) diff(c(0, df_X$pos[df_X$chr == chr]))))  # these are the same for all following type as well as X and 1 sample
      nn_dist_downstream <- unlist(sapply(chrs, function(chr) diff(c(df_X$pos[df_X$chr == chr], chr_sizes[[chr]]))))
      nn_dist_pmin <- pmin(nn_dist_upstream, nn_dist_downstream)
      
      nn_dist_upstream[nn_dist_upstream > max_dist] <- max_dist
      nn_dist_downstream[nn_dist_downstream > max_dist] <- max_dist
      nn_dist_pmin[nn_dist_pmin > max_dist] <- max_dist
      
      rm(df_X, df_1)
      
      helper_plot <- function(...) {
        plot(..., pch=20, col=adjustcolor(palette()[1], alpha.f=min(0.1, 0.1*5000/length(nn_dist_upstream))), xlim=c(0, max_dist))
      }
      helper_plot_binned_mean <- function(x, y) {
        x_binned <- round(x/20)*20
        temp_df <- aggregate(y, list(x_binned), mean, na.rm = TRUE)
        colnames(temp_df) <- c("binned_distance", "mean_count")
        lines(temp_df$binned_distance, temp_df$mean_count, col="green", lwd=2.5)
      }
      
      for(name_X_1 in c(name_X, name_1)) {
        name_add <- ifelse(name_X_1 == name_X, "X", "1")
        for(type in c("starts_plus", "ends_plus", "uncut_plus", "starts_minus", "ends_minus", "uncut_minus")) {
          df <- subset(cut_counts_df_list[[name_X_1]], chr %in% chrs)
          counts <- df[, type]
          counts <- set_high_outliers_to_quantile(counts, outlier_quantile)
          
          helper_plot(nn_dist_upstream, counts, ylab=paste0(type, "_", name_add))
          helper_plot_binned_mean(nn_dist_upstream, counts)
          
          helper_plot(nn_dist_downstream, counts, ylab=paste0(type, "_", name_add))
          helper_plot_binned_mean(nn_dist_downstream, counts)
        }
      }
      
      mtext(paste0(name, ", ", genome_names[genome_type]), side=3, line=0, outer=TRUE, cex=1.5, font=1)
      dev.off()
    }
  }
}


plot_count_histograms <- function(cut_counts_df_list, chr_name_list, genome_names, plot_folder = window_folder, outlier_quantile=0.99, name_add="") {
  
  plot_folder = paste0(plot_folder, "count_histograms", name_add, "/")
  dir.create(plot_folder)
  
  for(name_X in grep("_X$", names(cut_counts_df_list), value = TRUE)) {
    
    name <- sub("_X$", "", name_X)
    name_1 <- paste0(name, "_1")
    
    for(genome_type in names(chr_name_list)) {
      
      png(paste0(plot_folder, "counts_histogram_", name, "_", genome_type, ".png"), res=144, width=2048, height=1048)
      par(mfrow=c(2, 6), oma=c(2, 0, 3, 0))
      
      chrs <- chr_name_list[[genome_type]]
      
      df_X <- subset(cut_counts_df_list[[name_X]], chr %in% chrs)
      df_1 <- subset(cut_counts_df_list[[name_1]], chr %in% chrs)
      
      count_limit <- quantile(c(unlist(df_X[, c("starts_plus", "starts_minus", "ends_plus", "ends_minus")]), unlist(df_1[, c("starts_plus", "starts_minus", "ends_plus", "ends_minus")])), outlier_quantile, na.rm = TRUE)
      
      rm(df_X, df_1)
      
      helper_hist <- function(...) {
        num_bars_goal <- 40
        x <- list(...)[[1]]
        x_max <- ceiling(max(x))
        seqs <- sapply(1:10, function(s) seq(-0.5, x_max+0.5+s, s))
        lengths <- sapply(seqs, length)
        deviations <- abs(lengths-num_bars_goal)
        i <- which(deviations == min(deviations))[1]
        hist(..., breaks = seqs[[i]] )
      }
      
      for(name_X_1 in c(name_X, name_1)) {
        name_add <- ifelse(name_X_1 == name_X, "X", "1")
        for(type in c("starts_plus", "starts_minus", "ends_plus", "ends_minus", "uncut_plus", "uncut_minus")) {
          
          df <- subset(cut_counts_df_list[[name_X_1]], chr %in% chrs)
          counts <- df[, type]
          
          counts[counts > count_limit] <- count_limit
          counts <- counts[!is.na(counts)]
          
          helper_hist(counts, xlab=paste0(type, "_", name_add), main=paste0(type, "_", name_add))
        }
      }
      
      mtext(paste0(name, ", ", genome_names[genome_type]), side=3, line=0, outer=TRUE, cex=1.5, font=1)
      dev.off()
    }
  }
}


append_cut_count_vectors <- function(cut_count_df) {
  return(c(cut_count_df[, "starts_plus"], cut_count_df[, "starts_minus"], cut_count_df[, "ends_plus"], cut_count_df[, "ends_minus"]))
}


print_cut_count_statistics <- function(cut_counts_df_list, bg_counts_df_list, fragment_counts_list, count_window_list, chr_name_list, output_file="") {  # checked 2018-02-25
  
  
  calc_valid_site_number <- function(cut_counts_df) {  # these are sites with valid count values (starts or ends)
    return(with(cut_counts_df, sum( !( is.na(starts_plus) & is.na(ends_plus) ) )))
  }
  
  calc_mean_cut_fragments <- function(cut_counts_df) {  
    # sums up counts on plus and minus strand, averages over sites with valid start/end count and sums up starts and ends
    # i.e. gives the mean cuts for a site without close neighbours (using the counts from sites with 0 or 1 close neighbour)
    return(with(cut_counts_df, mean((starts_plus + starts_minus)[!is.na(starts_plus)], na.rm = TRUE) + 
                  mean((ends_plus + ends_minus)[!is.na(ends_plus)], na.rm = TRUE)))
  }
  
  calc_mean_uncut_fragments <- function(cut_counts_df) {  
    # sums up counts on plus and minus strand and averages over sites with valid start/end counts and takes the average of mean values for starts and ends
    # i.e. gives the mean uncut fragment count for a site without close neighbours (using the counts from sites with 0 or 1 close neighbour)
    return(with(cut_counts_df, ( mean((uncut_minus + uncut_plus)[!is.na(starts_plus)], na.rm = TRUE) + 
                                   mean((uncut_minus + uncut_plus)[!is.na(ends_plus)], na.rm = TRUE) ) / 2))
  }
  
  data_names <- names(cut_counts_df_list)
  
  main_df_list <- lapply(cut_counts_df_list, function(df) subset(df, chr %in% chr_name_list[[1]]))
  norm_df_list <- lapply(cut_counts_df_list, function(df) subset(df, chr %in% chr_name_list[[2]]))  # both have the same name order as in data_names
  
  cut_statistics <- data.frame(row.names = data_names)
  
  cut_statistics$reads_main <- sapply(data_names, function(name) sum(unlist(fragment_counts_list[[name]][chr_name_list[["main"]]])))
  cut_statistics$cuts_main <- sapply(main_df_list, function(df) sum(append_cut_count_vectors(df), na.rm = TRUE))
  cut_statistics$cuts_reads_ratio_main <- with(cut_statistics, cuts_main / reads_main)
  
  cut_statistics$sites_main <- sapply(main_df_list, calc_valid_site_number)
  cut_statistics$cut_mean_main <-  sapply(main_df_list, calc_mean_cut_fragments)
  cut_statistics$uncut_mean_main <- sapply(main_df_list, calc_mean_uncut_fragments)
  
  cut_statistics$bg_main_per_bp <- sapply(data_names, function(name) sum(bg_counts_df_list[[name]][1, c("starts_plus_per_bp", "starts_minus_per_bp", "ends_plus_per_bp", "ends_minus_per_bp")]))
  cut_statistics$bg_cov_main <- sapply(data_names, function(name) sum(bg_counts_df_list[[name]][1, c("coverage_plus", "coverage_minus")]))
  cut_statistics$bg_cut_cov_ratio_main <- sapply(data_names, function(name) with(bg_counts_df_list[[name]][1, ], (starts_plus_per_bp/coverage_plus + starts_minus_per_bp/coverage_minus + ends_plus_per_bp/coverage_plus + ends_minus_per_bp/coverage_minus) / 4))
  
  
  cut_statistics$reads_norm <- sapply(data_names, function(name) sum(unlist(fragment_counts_list[[name]][chr_name_list[["norm"]]])))
  cut_statistics$cuts_norm <- sapply(norm_df_list, function(df) sum(append_cut_count_vectors(df), na.rm = TRUE))
  cut_statistics$cuts_reads_ratio_norm <- with(cut_statistics, cuts_norm / reads_norm)
  
  cut_statistics$sites_norm <- sapply(norm_df_list, calc_valid_site_number)
  cut_statistics$cut_mean_norm <-  sapply(norm_df_list, calc_mean_cut_fragments)
  cut_statistics$uncut_mean_norm <- sapply(norm_df_list, calc_mean_uncut_fragments)
  
  cut_statistics$bg_norm_per_bp <- sapply(data_names, function(name) sum(bg_counts_df_list[[name]][2, c("starts_plus_per_bp", "starts_minus_per_bp", "ends_plus_per_bp", "ends_minus_per_bp")]))
  cut_statistics$bg_cov_norm <- sapply(data_names, function(name) sum(bg_counts_df_list[[name]][2, c("coverage_plus", "coverage_minus")]))
  cut_statistics$bg_cut_cov_ratio_norm <- sapply(data_names, function(name) with(bg_counts_df_list[[name]][2, ], (starts_plus_per_bp/coverage_plus + starts_minus_per_bp/coverage_minus + ends_plus_per_bp/coverage_plus + ends_minus_per_bp/coverage_minus) / 4))
  
  cut_statistics$count_limit <- sapply(data_names, function(name) mean(unlist(lapply(count_window_list[[name]], max))))  # mean if several REs in the same sample
  
  cut_statistics$bg_main_per_bp_per_read_mio <- with(cut_statistics, bg_main_per_bp/reads_main*1000000)
  cut_statistics$bg_norm_per_bp_per_read_mio <- with(cut_statistics, bg_norm_per_bp/reads_norm*1000000)
  
  if(output_file != "none") {
    cat(kable(cut_statistics, digits = 4), sep = "\n", file = output_file)
    print(cut_statistics, digits = 4)
  }
  
  return(cut_statistics)
}


print_occs_from_mean_counts <- function(cut_statistics, normalization_method = "norm_cuts", output_file = "") {  # new 2018-02-27
  
  # calculate some preliminary results: accessibilities using the counts averaged over all sites
  # (this is not the same as the average accessibility over all sites, which is the real average accessibility)
  
  names_X <- grep("_X$", rownames(cut_statistics), value = TRUE)
  names_1 <- grep("_1$", rownames(cut_statistics), value = TRUE)
  cuts_X <- cut_statistics[names_X, "cut_mean_main"]
  uncuts_X <- cut_statistics[names_X, "uncut_mean_main"]
  if(normalization_method == "norm_reads") {
    cuts_1 <- with(cut_statistics[names_1, ], cut_mean_main / reads_norm) * cut_statistics[names_X, "reads_norm"]
  } else if(normalization_method == "norm_cuts") {
    cuts_1 <- with(cut_statistics[names_1, ], cut_mean_main / cuts_norm) * cut_statistics[names_X, "cuts_norm"]
  } else if(normalization_method == "main_reads") {
    cuts_1 <- with(cut_statistics[names_1, ], cut_mean_main / reads_main) * cut_statistics[names_X, "reads_main"]
  }
  bg <- cut_statistics[names_X, "bg_main_per_bp"] * (cut_statistics[names_X, "count_limit"] + 1)  # this is not correct if there are two different REs in the same sample
  
  accs_occs_main <- data.frame(row.names = names_X)
  
  accs_occs_main$acc_bg_untreated <- cuts_X / cuts_1
  accs_occs_main$acc_bg_michael <- (cuts_X - bg) / (cuts_1 - bg)
  accs_occs_main$acc_cut_uncut_untreated <- cuts_X / (cuts_X + uncuts_X)
  accs_occs_main$acc_cut_uncut_bg_michael <- (cuts_X - bg) / (cuts_X - bg + uncuts_X)
  
  accs_occs_main$occ_bg_untreated <- 1 - accs_occs_main$acc_bg_untreated
  accs_occs_main$occ_bg_michael <- 1 - accs_occs_main$acc_bg_michael
  accs_occs_main$occ_cut_uncut_untreated <-  1- accs_occs_main$acc_cut_uncut_untreated
  accs_occs_main$occ_cut_uncut_bg_michael <- with(accs_occs_main, 1-acc_cut_uncut_bg_michael)
  
  accs_occs_main$shear_prob_per_window_in_pc <- 100 * bg / cuts_1
  accs_occs_main$bg_cut_cov_ratio_in_pc <- 100 * cut_statistics[names_X, "bg_cut_cov_ratio_main"]
  
  accs_occs_main$shear_prob_per_bp_in_pc <- accs_occs_main$shear_prob_per_window_in_pc / (cut_statistics[names_X, "count_limit"] + 1)
  accs_occs_main$bg_cut_cov_ratio_per_bp_in_pc <- accs_occs_main$bg_cut_cov_ratio_in_pc / (cut_statistics[names_X, "count_limit"] + 1)
  
  accs_occs_main <- accs_occs_main[sort(rownames(accs_occs_main)), ]
  
  cat(kable(accs_occs_main, digits = 3, format.args = list(scientific = FALSE)), sep = "\n", file = output_file)
  
  return(accs_occs_main)
}


calc_genomic_counts <- function(positions_vec, chr_length) {  # positions meaning starts or ends
  count_table <- table(positions_vec)
  counts_genomic <- rep(0, chr_length)
  counts_genomic[as.numeric(names(count_table))] <- count_table
  
  return(counts_genomic)
}


calc_cuts_of_average_site_and_plot_genomic_cuts <- function(data_names, genome, chr_name_list, sample_to_RE_table, RE_info, 
                                                            bad_regions, large_window_limit, small_window_limits, plot_folder, close_distances, RData_folder="data/RData_files/", max_length=500) {  # checked 2018-02-23
  
  plot_folder <- paste0(plot_folder, "genomic_cut_distribution/")
  
  large_window_region <- -large_window_limit:large_window_limit
  
  reads_in_window_list <- list()
  
  if(Sys.info()['sysname'] == "Linux") {
    warning("Genomic plots are too large for Linux png and aren't useful when plotted smaller, they are omitted.")
    plot_genomic <- FALSE
  } else {
    plot_genomic <- TRUE
    dir.create(plot_folder)
  }
  
  count_window_df_list <- list()
  
  for(name in data_names) {
    
    reads_in_window_list[[name]] <- list()
    
    load(paste0(RData_folder, name, ".RData"))
    
    for(i in 1:length(close_distances[[1]])) {
      
      dist_ign_site <- close_distances[["dist_ign_site"]][i]
      dist_ign_half <- close_distances[["dist_ign_half"]][i]
      close_dist_str <- paste0(dist_ign_site, "_", dist_ign_half)
      
      reads_in_window_list[[name]][[close_dist_str]] <- list()
      
      sites_either <- calc_cut_sites(name, genome, sample_to_RE_table, RE_info, bad_regions, large_window_limit, dist_ign_site, dist_ign_half, filter_type = "either")
      sites_start <- calc_cut_sites(name, genome, sample_to_RE_table, RE_info, bad_regions, large_window_limit, dist_ign_site, dist_ign_half, filter_type = "start")
      sites_end <- calc_cut_sites(name, genome, sample_to_RE_table, RE_info, bad_regions, large_window_limit, dist_ign_site, dist_ign_half, filter_type = "end")
      
      for(enzyme in names(sites_either[[name]])) {
        
        if(grepl("BamHI", enzyme) & grepl("RE16", getwd())) {  # remove sites that are at the start or end of the plasmid inserts
          load(file="../../reference_genomes/pGP564/bam_sites_near_insert_ends.RData")
          for(i in 1:16) {
            bam_sites_near_insert_ends$chr[bam_sites_near_insert_ends$chr == i] <- chr_name_list[["main"]][i]
          }
          cat(paste0(name, ": removing sites near plasmid inserts.\n"))
          for(this_chr in names(sites_start[[1]][[enzyme]])) {
            sites_start[[1]][[enzyme]][[this_chr]] <- sites_start[[1]][[enzyme]][[this_chr]][!sites_start[[1]][[enzyme]][[this_chr]] %in% subset(bam_sites_near_insert_ends, chr == this_chr)$pos]
            sites_end[[1]][[enzyme]][[this_chr]] <- sites_end[[1]][[enzyme]][[this_chr]][!sites_end[[1]][[enzyme]][[this_chr]] %in% subset(bam_sites_near_insert_ends, chr == this_chr)$pos]
            sites_either[[1]][[enzyme]][[this_chr]] <- sites_either[[1]][[enzyme]][[this_chr]][!sites_either[[1]][[enzyme]][[this_chr]] %in% subset(bam_sites_near_insert_ends, chr == this_chr)$pos]
          }
        }
        
        reads_in_window_list[[name]][[close_dist_str]][[enzyme]] <- list()
        
        cat(name, ", enzyme = ", enzyme, ", close_distances = ", close_dist_str, "\n", sep = "")
        
        # genomic distribution plots
        for(genome_type in names(chr_name_list)) {
          
          chrs <- chr_name_list[[genome_type]]
          if(genome_type == "main") {
            chrs_with_sites <- chrs[mapply(length, sites_either[[name]][[enzyme]][chr_name_list$main]) > 0]
            # pdf(paste0(plot_folder, "genomic_cut_distribution_", gsub("%", "%%", name), "_", genome_type, "_close_distances_", close_dist_str, "_", enzyme, "_sites.png"), 
            #     height = 1.5*length(chrs_with_sites), width = 20, pointsize = 0.5)
            if(plot_genomic) {
              png(paste0(plot_folder, "genomic_cut_distribution_", gsub("%", "%%", name), "_", genome_type, "_close_distances_", close_dist_str, "_", enzyme, "_sites.png"), 
                  height = 400*length(chrs_with_sites), width = 20000, res = 2*72)
            }
          } else if(genome_type == "norm") {
            chrs_with_sites <- chrs[mapply(length, sites_either[[name]][[enzyme]][chr_name_list$norm]) > 0]
            # pdf(paste0(plot_folder, "genomic_cut_distribution_", gsub("%", "%%", name), "_", genome_type, "_close_distances_", close_dist_str, "_", enzyme, "_sites.png"), 
            #     height = 1.5*length(chrs_with_sites), width = 3.6*10, pointsize = 0.5)
            if(plot_genomic) {
              png(paste0(plot_folder, "genomic_cut_distribution_", gsub("%", "%%", name), "_", genome_type, "_close_distances_", close_dist_str, "_", enzyme, "_sites.png"), 
                  height = 400*length(chrs_with_sites), width = 3.6*20000, res = 2*72)
            }
          }
          
          par(mfrow = c(length(chrs_with_sites), 1))
          
          read_starts_p_in_window_list <- list()
          read_starts_m_in_window_list <- list()
          read_ends_m_in_window_list <- list()
          read_ends_p_in_window_list <- list()
          
          for(chr in chrs_with_sites) {
            
            sites_start_chr <- sites_start[[name]][[enzyme]][[chr]]
            sites_end_chr <- sites_end[[name]][[enzyme]][[chr]]
            
            sites_start_chr <- sites_start_chr[sites_start_chr > large_window_limit + RE_info[enzyme, "site_shift"]]  # to avoid negative indeces in read_starts[i+large_window_region] 
            sites_end_chr <- sites_end_chr[sites_end_chr > large_window_limit + RE_info[enzyme, "site_shift"]]
            
            
            chr_bam <- bam[[chr]]
            chr_bam <- chr_bam[chr_bam$end - chr_bam$start <= max_length, ]
            
            chr_length <- width(genome[chr]) + 1  # a lot of reads extend the chromosome by one base (this is due to the +1 shift of read ends when reading in the bam file)
            
            reads_p_strand <- chr_bam$strand == "+"
            reads_m_strand <- chr_bam$strand == "-"
            
            read_starts_p <- calc_genomic_counts(chr_bam$start[reads_p_strand], chr_length)
            read_starts_m <- calc_genomic_counts(chr_bam$start[reads_m_strand], chr_length)
            read_ends_m <- calc_genomic_counts(chr_bam$end[reads_m_strand], chr_length)
            read_ends_p <- calc_genomic_counts(chr_bam$end[reads_p_strand], chr_length)
            
            # Explanation of site_shift:
            #
            # middle_shift =  middle of the recognition sequence
            # site_shift = how far left the cut site is from the middle for the upper strand
            # 
            # example HindIII
            # 5' - A|A G C T T - 3'  (plus strand)
            # 3' - T T C G A|A - 5'  (minus strand)
            # 
            # If the 3' end is shorter than the 5' end after cutting, the 3' end is elongated to match the 5' end (by polymerase).
            # If the 3' end is longer  that the 5' end after cutting, the 3' end is digested  to match the 5' end.
            # 
            # Thus site_shift is also how many bases to add from the middle of the sequence when looking for reads that were cut.
            # 
            # sites_chr = recognition_site start on the plus strand + middle_shift
            # 
            # cut_fragment_position = sites_chr - site_shift (for reads starting at the cut site using the genomic position)
            # cut_fragment_position = sites_chr + site_shift (for reads ending at the cut site using the genomic position)
            
            read_starts_p_in_window_list[[chr]] <- sapply(large_window_region, function (i) read_starts_p[i + sites_start_chr - RE_info[enzyme, "site_shift"]])
            read_starts_m_in_window_list[[chr]] <- sapply(large_window_region, function (i) read_starts_m[i + sites_start_chr - RE_info[enzyme, "site_shift"]])
            read_ends_m_in_window_list[[chr]] <- sapply(large_window_region, function (i) read_ends_m[i + sites_end_chr + RE_info[enzyme, "site_shift"]])
            read_ends_p_in_window_list[[chr]] <- sapply(large_window_region, function (i) read_ends_p[i + sites_end_chr + RE_info[enzyme, "site_shift"]])
            
            if(plot_genomic) {
              read_starts <- read_starts_p + read_starts_m
              read_ends <- read_ends_m + read_ends_p
              
              for(bad_index in which(bad_regions$chr == chr)) {
                read_starts[bad_regions$start[bad_index]:bad_regions$end[bad_index]] <- NA
                read_ends[bad_regions$start[bad_index]:bad_regions$end[bad_index]] <- NA
              }
              
              pos_binning <- 10
              binned_positions <- seq(1, length(read_starts)-pos_binning, pos_binning)
              read_starts_binned <- do.call(rbind, lapply(1:pos_binning, function(i) read_starts[binned_positions+i] ))
              read_starts_binned_max <- apply(read_starts_binned, 2, max)
              read_starts_binned_mean <- colSums(read_starts_binned) / pos_binning  # can also be used
              
              x_max <- max(width(genome[chr_name_list[[genome_type]]]))
              
              plot(x = binned_positions, y = read_starts_binned_max, type="l", lwd = 0.1, xlim = c(0, x_max), xaxt = "n", 
                   main = paste(name, chr, "read starts with pos_binning =", pos_binning, "and close_distances =", close_dist_str), panel.first = abline(v = sites_start_chr, col = "lightblue", lwd = 0.1),
                   xaxs="i")
              axis(1, at = seq(0, x_max, 10000))
              grid(col = "grey", lwd = 0.2)
              for(bad_index in which(bad_regions$chr == chr)) {
                rect(bad_regions$start[bad_index], 0, bad_regions$end[bad_index], max(read_starts), col = adjustcolor("red", alpha.f = 0.5), border = NA)
              }
            }
          }
          if(plot_genomic) {
            dev.off()
          }
          
          # cutsite distribution plots
          
          read_starts_p_in_window <- do.call(rbind, read_starts_p_in_window_list)
          read_starts_m_in_window <- do.call(rbind, read_starts_m_in_window_list)
          read_ends_m_in_window <- do.call(rbind, read_ends_m_in_window_list)
          read_ends_p_in_window <- do.call(rbind, read_ends_p_in_window_list)
          
          read_starts_in_window <- read_starts_p_in_window + read_starts_m_in_window
          read_ends_in_window <- read_ends_p_in_window + read_ends_m_in_window
          
          read_all_in_window <- rbind(read_starts_in_window, read_ends_in_window[, seq(ncol(read_ends_in_window), 1, -1)])
          
          reads_in_window <- data.frame("position" = large_window_region)
          
          reads_in_window$starts_p_mean <- colMeans(read_starts_p_in_window)
          reads_in_window$starts_m_mean <- colMeans(read_starts_m_in_window)
          reads_in_window$ends_m_mean <- colMeans(read_ends_m_in_window)
          reads_in_window$ends_p_mean <- colMeans(read_ends_p_in_window)
          
          reads_in_window$starts_mean <- colMeans(read_starts_in_window)
          reads_in_window$ends_mean <- colMeans(read_ends_in_window)
          reads_in_window$all_mean <- 2*colMeans(read_all_in_window)
          
          reads_in_window$starts_p_sd <- apply(read_starts_p_in_window, 2, sd)
          reads_in_window$starts_m_sd <- apply(read_starts_m_in_window, 2, sd)
          reads_in_window$ends_m_sd <- apply(read_ends_m_in_window, 2, sd)
          reads_in_window$ends_p_sd <- apply(read_ends_p_in_window, 2, sd)
          
          reads_in_window$starts_sd <- apply(read_starts_in_window, 2, sd)
          reads_in_window$ends_sd <- apply(read_ends_in_window, 2, sd)
          # reads_in_window$all_sd <- apply(read_all_in_window, 2, sd)  # we should add up starts and ends at each site instead (problem: start sites and ends sites are not always the same)
          
          reads_in_window$num_sites_start <- nrow(read_starts_p_in_window)
          reads_in_window$num_sites_end <- nrow(read_ends_p_in_window)
          
          reads_in_window_list[[name]][[close_dist_str]][[enzyme]][[genome_type]] <- reads_in_window
        }
      }
    }
  }
  
  return(reads_in_window_list)
}


calc_count_limits_and_plot_average_site <- function(cuts_near_average_site_list, chr_name_list, genome_names, large_window_limit, small_window_limits, plot_folder, RData_folder="data/RData_files/") {
  
  plot_folder <- paste0(plot_folder, "cut_distribution_near_sites/")
  
  dir.create(plot_folder)
  
  large_window_region <- -large_window_limit:large_window_limit
  
  count_window_df_list <- list()
  
  for(name in names(cuts_near_average_site_list)) {
    
    for(close_dist_str in names(cuts_near_average_site_list[[name]])) {
      
      for(enzyme in names(cuts_near_average_site_list[[name]][[close_dist_str]])) {
        
        cat(name, ", close_distances = ", close_dist_str, ", enzyme = ", enzyme, "\n", sep = "")
        
        for(genome_type in names(chr_name_list)) {
          
          reads_in_window <- cuts_near_average_site_list[[name]][[close_dist_str]][[enzyme]][[genome_type]]
          
          png(paste0(plot_folder,  "cut_distribution_near_sites_", gsub("%", "%%" , name), "_", genome_type, "_close_distances_", close_dist_str, "_", enzyme, "_sites.png"), res = 144, width = 1440, height = 1280)
          
          par(mfrow = c(4, 3), oma = c(0, 0, 5, 0))
          
          helper_plot_windows <- function(mean_vec, sd_vec, large_window_region, small_window_limits, main_title, count_window_limit, y_max) {
            if(missing(y_max)) {
              y_max <- 1.2*max(mean_vec + sd_vec)
            }
            if(count_window_limit >= 0) {
              rect_pos <- c(-0.5, 0, count_window_limit+0.5, y_max)
            } else {
              rect_pos <- c(-0.5+count_window_limit, 0, 0.5, y_max)
            }
            plot(large_window_region, mean_vec + sd_vec, col = "grey", type = 'l', 
                 xlab = "distance to expected RE cut site", ylab = "count", main = paste0(main_title, ", large window"), ylim = c(0, y_max))
            rect(rect_pos[1], rect_pos[2], rect_pos[3], rect_pos[4], density = NULL, angle = 45, col = "greenyellow", border = NA)
            lines(large_window_region, mean_vec + sd_vec, col = "grey", type = 'l')
            lines(large_window_region, mean_vec, type = 'l')
            legend("topright", legend = c("mean", "mean + std"), lty = c(1, 1), col = c("black", "grey"), cex = 0.75)
            
            plot(large_window_region, mean_vec + sd_vec, col = "grey", type = 'b',
                 xlim = small_window_limits, xlab = "distance to expected RE cut site", 
                 ylab = "count", main = paste0(main_title, ", small window"), ylim = c(0, y_max))
            rect(rect_pos[1], rect_pos[2], rect_pos[3], rect_pos[4], density = NULL, angle = 45, col = "greenyellow", border = NA)
            lines(large_window_region, mean_vec + sd_vec, col = "grey", type = 'b')
            lines(large_window_region, mean_vec, type = 'b')
            
            plot(large_window_region, mean_vec + sd_vec, col = "grey", type = 'b',
                 xlim = small_window_limits, xlab = "distance to expected RE cut site", 
                 ylab = "count", main = paste0(main_title, ", small window, zoomed"), ylim = c(0, y_max/5))
            rect(rect_pos[1], rect_pos[2], rect_pos[3], y_max/5, density = NULL, angle = 45, col = "greenyellow", border = NA)
            lines(large_window_region, mean_vec + sd_vec, col = "grey", type = 'b')
            lines(large_window_region, mean_vec, type = 'b')
            
            return(y_max)
          }
          
          coord_starts <- seq(large_window_limit+1, 2*large_window_limit+1, 1)
          coord_ends <- seq(large_window_limit+1, 1, -1)
          mean_counts <- (reads_in_window$starts_p_mean[coord_starts] + reads_in_window$starts_m_mean[coord_starts] + reads_in_window$ends_p_mean[coord_ends] + reads_in_window$ends_m_mean[coord_ends]) / 4
          
          ## Finding a good count_window_limit
          
          if(grepl("_0_X$", name) & genome_type == "main") {
            cat("Found ..._0_X sample, setting count_window_limit to 5 for main genome\n")
            count_window_limit <- 5
            mean_resection_length <- NA
          } else {
            
            # simple background treatment (just for finding the count_window_limit and mean resection length): correct first half by the mean of the second half
            mean_counts_corrected <- mean_counts[1:100] - mean(mean_counts[101:201], na.rm = TRUE)
            counts_so_far <- cumsum(mean_counts_corrected)
            
            N <- 5
            f <- 0.01
            count_window_limits <- rep(0, N)
            for(n in 1:N) {
              # count_window_limits(n): largest distance from cutsite where the count sum of the next n sites is greater or equal a fraction f of the counts so far
              #                        i. e. also counting the next n sites after count_window_limit would add less than a fraction of f to the counts up to count_window_limit
              mean_of_next_n_sites <- sapply(1:(100-n), function(x) sum(mean_counts_corrected[x+(1:n)]))
              logic_vec <- !(f * pmax(0, counts_so_far[1:(100-n)]) <= mean_of_next_n_sites)  # the first TRUE is the position we are looking for
              if(any(logic_vec)) {
                count_window_limits[n] <- which(logic_vec)[1] - 1  # "- 1" because the first entry is the position 0
              } else {
                count_window_limits[n] <- 0
              }
            }
            print(count_window_limits)
            count_window_limit <- max(count_window_limits)
            
            resection_counts <- mean_counts_corrected[1+(0:count_window_limit)]
            mean_resection_length <- sum((0:count_window_limit) * resection_counts) / sum(resection_counts)
          }
          
          mean_counts_within_limit <- sum(mean_counts[1:(count_window_limit+1)], na.rm = TRUE)
          
          count_window_df_list[[length(count_window_df_list)+1]] <- data.frame(name = name, genome_type = genome_type, close_dist_str = close_dist_str, enzyme = enzyme, count_window_limit = count_window_limit, mean_counts_within_limit = mean_counts_within_limit, mean_resection_length = mean_resection_length)
          
          y_max <- helper_plot_windows(reads_in_window$starts_p_mean, reads_in_window$starts_p_sd, large_window_region, small_window_limits, "starts on + strand", count_window_limit)
          helper_plot_windows(reads_in_window$starts_m_mean, reads_in_window$starts_m_sd, large_window_region, small_window_limits, "starts on - strand", count_window_limit, y_max)
          helper_plot_windows(reads_in_window$ends_m_mean, reads_in_window$ends_m_sd, large_window_region, -small_window_limits[c(2, 1)], "ends on - strand", -count_window_limit, y_max)
          helper_plot_windows(reads_in_window$ends_p_mean, reads_in_window$ends_p_sd, large_window_region, -small_window_limits[c(2, 1)], "ends on + strand", -count_window_limit, y_max)
          
          mtext(paste0(name, ", ", genome_names[genome_type], ", ", enzyme), side = 3, line = 3, outer = TRUE, cex = 1.2, font = 1)
          
          mtext(paste0("count_window_limit = ", format(count_window_limit, digits = 2), ", mean_counts_within_limit = ", format(mean_counts_within_limit, digits = 2), ", mean_resection_length = ", format(mean_resection_length, digits = 2), ", averaged over ", reads_in_window$num_sites_start[1], ", ", reads_in_window$num_sites_end[1], " sites (starts, ends, respectively)"), side = 3, line = 0, outer = TRUE, cex = 0.8, font = 1)
          
          dev.off()
        }
      }
    }
  }
  
  return(do.call(rbind, count_window_df_list))
}


calc_accs_from_mus <- function(mus, shear_prob_site_vec) {  # valid for the method with X and 1 sample
  accs <- (mus - shear_prob_site_vec) / (1 - shear_prob_site_vec)
  return(accs)
}


pmean <- function(x, y) {  # alternative to rowMeans
  return(ifelse(is.na(x), y, ifelse(is.na(y), x, (x + y) / 2)))
}


plogmean <- function(x, y) {
  x <- log(x)
  y <- log(y)
  return(exp(ifelse(is.na(x), y, ifelse(is.na(y), x, (x + y) / 2))))
}


calc_accessibilities <- function(cut_counts_df_list, bg_counts_df_list, count_window_list, background_method, chr_name_list, acc_max = 1.5, 
                                 uncut_loss_corrections, uncut_loss_corrections_per_enzyme) {
  
  if(missing(uncut_loss_corrections) & missing(uncut_loss_corrections_per_enzyme)) stop("No uncut correction provided!")
  
  accs_df_list <- list()
  mus_df_list <- list()
  shear_prob_per_bp_list <- list()
  bg_cut_per_bp_cov_ratio_list <- list()
  
  for(name_X in grep("_X$", names(cut_counts_df_list), value = TRUE)) {
    
    if(!missing(uncut_loss_corrections_per_enzyme)) {
      enzyme <- "combined"
      for(this_enzyme in names(uncut_loss_corrections_per_enzyme)[-1]) {
        if(grepl(this_enzyme, name_X)) {
          enzyme <- this_enzyme
        }
      }
      uncut_loss_corrections <- uncut_loss_corrections_per_enzyme[[enzyme]]
    }
    
    name <- sub("_X$", "", name_X)
    name_1 <- paste0(name, "_1")
    
    chrs <- chr_name_list[["main"]]
    df_X <- subset(cut_counts_df_list[[name_X]], chr %in% chrs)
    df_1 <- subset(cut_counts_df_list[[name_1]], chr %in% chrs)
    
    chrs <- chr_name_list[["norm"]]
    df_X_norm <- subset(cut_counts_df_list[[name_X]], chr %in% chrs)
    df_1_norm <- subset(cut_counts_df_list[[name_1]], chr %in% chrs)
    
    scale_1_to_X <- sum(append_cut_count_vectors(df_X_norm), na.rm=TRUE) / sum(append_cut_count_vectors(df_1_norm), na.rm=TRUE)
    
    window_length_enzyme_vec <- unlist(lapply(count_window_list[[name_X]], function(count_window) max(count_window) - min(count_window) + 1))
    
    bg_cut_main_per_bp <- with(bg_counts_df_list[[name_X]][1, ], starts_plus_per_bp + starts_minus_per_bp + ends_plus_per_bp + ends_minus_per_bp) / 4
    bg_cut_main_enzyme_vec <- window_length_enzyme_vec * bg_cut_main_per_bp
    
    if(background_method == "michael") {
      shear_prob_enzyme_vec <- bg_cut_main_enzyme_vec / mean(append_cut_count_vectors(df_1), na.rm=TRUE) / scale_1_to_X
    } else if(background_method == "untreated") {
      shear_prob_enzyme_vec <- window_length_enzyme_vec * 0
    } else {
      stop("Unknown background method")
    }
    
    if(any(shear_prob_enzyme_vec > 1) || any(shear_prob_enzyme_vec < 0)) {
      stop("Invalid shear_prob!")
    }
    
    shear_prob_per_bp_list[[name]] <- bg_cut_main_per_bp / mean(append_cut_count_vectors(df_1), na.rm=TRUE) / scale_1_to_X
    bg_cut_per_bp_cov_ratio_list[[name]] <- with(bg_counts_df_list[[name_X]][1, ], (starts_plus_per_bp/coverage_plus + starts_minus_per_bp/coverage_minus + ends_plus_per_bp/coverage_plus + ends_minus_per_bp/coverage_minus) / 4)
    
    mu_max_enzyme_vec <- acc_max * (1-shear_prob_enzyme_vec) + shear_prob_enzyme_vec
    
    calc_mus <- function(counts_X, counts_1, enzyme_vec, scale_1_to_X) {
      mus <- counts_X / counts_1 / scale_1_to_X
      mus <- ifelse(counts_1 > 0, mus, ifelse(counts_X > 0, NA, NA))  # We set all sites with counts_1==0 to NA
      over_max <- !is.na(mus) & mus > mu_max_enzyme_vec[enzyme_vec]
      mus[over_max] <- mu_max_enzyme_vec[enzyme_vec[over_max]]
      return(mus)
    }
    
    mus_df <- list()
    for(type in c("starts_plus", "starts_minus", "ends_plus", "ends_minus")) {
      mus_df[[type]] <- calc_mus(df_X[, type], df_1[, type], df_X[, "enzyme"], scale_1_to_X)
    }
    
    mus_df$starts_mean <- with(mus_df, pmean(starts_plus, starts_minus))
    mus_df$starts_logmean <- with(mus_df, plogmean(starts_plus, starts_minus))
    mus_df$starts_countmean <- calc_mus( with(df_X, pmean(starts_plus, starts_minus)), with(df_1, pmean(starts_plus, starts_minus)), df_X[, "enzyme"], scale_1_to_X )
    
    mus_df$ends_mean <- with(mus_df, pmean(ends_plus, ends_minus))
    mus_df$ends_logmean <- with(mus_df, plogmean(ends_plus, ends_minus))
    mus_df$ends_countmean <- calc_mus( with(df_X, pmean(ends_plus, ends_minus)), with(df_1, pmean(ends_plus, ends_minus)), df_X[, "enzyme"], scale_1_to_X )
    
    mus_df$all_mean <- with(mus_df, pmean(starts_mean, ends_mean))
    mus_df$all_logmean <- with(mus_df, plogmean(starts_logmean, ends_logmean))
    mus_df$all_countmean_mean <- with(mus_df, pmean(starts_countmean, ends_countmean))
    mus_df$all_countmean_logmean <- with(mus_df, plogmean(starts_countmean, ends_countmean))
    
    mus_df_list[[name]] <-  cbind(df_X[, c(1, 2, 3)], mus_df)
    
    accs_df <- lapply(mus_df[!grepl("cut_uncut", names(mus_df))], calc_accs_from_mus, shear_prob_site_vec = shear_prob_enzyme_vec[df_X$enzyme])
    
    
    accs_cut_uncut_df <- list()
    
    calc_acc_from_cut_uncut_ratio <- function(cut_uncut_ratio, cut_coverage_ratio_bg, window_length, uncut_loss_correction) {
      if(!is.finite(cut_coverage_ratio_bg)) {
        error("cut_coverage_ratio_bg not finite!")
      } else {
        cut_uncut_ratio_bg <- cut_coverage_ratio_bg / (1 - cut_coverage_ratio_bg) / uncut_loss_correction
      }
      cut_uncut_ratio <- cut_uncut_ratio/uncut_loss_correction
      accs <- 1 - (1 + cut_uncut_ratio_bg) / (cut_uncut_ratio + 1 - cut_uncut_ratio_bg * (window_length-1))  # gives 1 if cut_uncut_ratio is Inf (when everything is cut)
      return(accs)
    }
    
    cut_coverage_ratio_bg <- with(bg_counts_df_list[[name_X]][1, ], (starts_plus_per_bp + ends_plus_per_bp + starts_minus_per_bp + ends_minus_per_bp) / 2 / (coverage_plus + coverage_minus))
    
    # properly modelled background correction (uncut loss not corrected)
    
    accs_cut_uncut_df$cut_uncut_starts_plus_1 <- with(df_X, calc_acc_from_cut_uncut_ratio(starts_plus / uncut_plus, cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_correction = 1))
    
    accs_cut_uncut_df$cut_uncut_ends_plus_1 <- with(df_X, calc_acc_from_cut_uncut_ratio(ends_plus / uncut_plus, cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_correction = 1))
    
    accs_cut_uncut_df$cut_uncut_starts_minus_1 <- with(df_X, calc_acc_from_cut_uncut_ratio(starts_minus / uncut_minus, cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_correction = 1))
    
    accs_cut_uncut_df$cut_uncut_ends_minus_1 <- with(df_X, calc_acc_from_cut_uncut_ratio(ends_minus / uncut_minus, cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_correction = 1))
    
    accs_cut_uncut_df$cut_uncut_all_1 <- with(accs_cut_uncut_df, pmean(pmean(cut_uncut_starts_plus_1, cut_uncut_ends_plus_1), pmean(cut_uncut_starts_minus_1, cut_uncut_ends_minus_1)))
    
    
    # properly modelled background correction (uncut loss is corrected)
    
    accs_cut_uncut_df$cut_uncut_starts_plus_2 <- with(df_X, calc_acc_from_cut_uncut_ratio(starts_plus / uncut_plus, cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_corrections[1]))
    
    accs_cut_uncut_df$cut_uncut_ends_plus_2 <- with(df_X, calc_acc_from_cut_uncut_ratio(ends_plus / uncut_plus, cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_corrections[1]))
    
    accs_cut_uncut_df$cut_uncut_starts_minus_2 <- with(df_X, calc_acc_from_cut_uncut_ratio(starts_minus / uncut_minus, cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_corrections[1]))
    
    accs_cut_uncut_df$cut_uncut_ends_minus_2 <- with(df_X, calc_acc_from_cut_uncut_ratio(ends_minus / uncut_minus, cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_corrections[1]))
    
    accs_cut_uncut_df$cut_uncut_starts_2 <- with(accs_cut_uncut_df, pmean(cut_uncut_starts_plus_2, cut_uncut_starts_minus_2))
    
    accs_cut_uncut_df$cut_uncut_ends_2 <- with(accs_cut_uncut_df, pmean(cut_uncut_ends_plus_2, cut_uncut_ends_minus_2))
    
    accs_cut_uncut_df$cut_uncut_all_2 <- with(accs_cut_uncut_df, pmean(cut_uncut_starts_2, cut_uncut_ends_2))
    
    
    plus_with_na_rm <- function(x, y) ifelse(!is.na(x), ifelse(!is.na(y), x + y, x), y)
    
    calc_eff_cuts <- function(cuts, uncuts, cut_coverage_ratio_bg, window_length, uncut_loss_correction) {
      if(!is.finite(cut_coverage_ratio_bg)) {
        error("cut_coverage_ratio_bg not finite!")
      } else {
        cut_uncut_ratio_bg <- cut_coverage_ratio_bg / (1 - cut_coverage_ratio_bg) / uncut_loss_correction
      }
      return(cuts - window_length * cut_uncut_ratio_bg * uncuts * uncut_loss_correction)
    }
    
    calc_eff_uncuts <- function(uncuts, cut_coverage_ratio_bg, uncut_loss_correction) {
      if(!is.finite(cut_coverage_ratio_bg)) {
        error("cut_coverage_ratio_bg not finite!")
      } else {
        cut_uncut_ratio_bg <- cut_coverage_ratio_bg / (1 - cut_coverage_ratio_bg) / uncut_loss_correction
      }
      return(uncuts * (1 + cut_uncut_ratio_bg) * uncut_loss_correction)
    }
    
    # averaging counts first (uncut loss not corrected)
    
    accs_cut_uncut_df$cut_uncut_starts_3 <- with(df_X, calc_acc_from_cut_uncut_ratio((starts_plus + starts_minus) / (uncut_plus + uncut_minus), cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_correction = 1))
    
    accs_cut_uncut_df$cut_uncut_ends_3 <- with(df_X, calc_acc_from_cut_uncut_ratio((ends_plus + ends_minus) / (uncut_plus + uncut_minus), cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_correction = 1))
    
    accs_cut_uncut_df$cut_uncut_all_3 <- with(df_X, calc_acc_from_cut_uncut_ratio( pmean(starts_plus + starts_minus, ends_plus + ends_minus) / (uncut_plus + uncut_minus), cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_correction = 1))
    
    
    # averaging counts first (uncut loss is corrected)
    
    accs_cut_uncut_df$eff_start_cuts <- with(df_X, calc_eff_cuts(starts_plus + starts_minus, uncut_plus + uncut_minus, cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_corrections[2]))
    accs_cut_uncut_df$eff_end_cuts <- with(df_X, calc_eff_cuts(ends_plus + ends_minus, uncut_plus + uncut_minus, cut_coverage_ratio_bg, window_length_enzyme_vec[enzyme], uncut_loss_corrections[2]))
    accs_cut_uncut_df$eff_uncuts <- with(df_X, calc_eff_uncuts(uncut_plus + uncut_minus, cut_coverage_ratio_bg, uncut_loss_corrections[2]))
    
    accs_cut_uncut_df$eff_cuts <- with(accs_cut_uncut_df, plus_with_na_rm(eff_start_cuts , eff_end_cuts))
    accs_cut_uncut_df$eff_uncuts <- with(accs_cut_uncut_df, ifelse(is.na(eff_cuts), NA, eff_uncuts))
    
    accs_cut_uncut_df$eff_coverage <- with(accs_cut_uncut_df, eff_cuts + ((!is.na(eff_start_cuts)) + (!is.na(eff_end_cuts))) * eff_uncuts)
    
    accs_cut_uncut_df$cut_uncut_starts_4 <- with(accs_cut_uncut_df, eff_start_cuts / (eff_start_cuts + eff_uncuts))
    accs_cut_uncut_df$cut_uncut_ends_4 <- with(accs_cut_uncut_df, eff_end_cuts / (eff_end_cuts + eff_uncuts))
    accs_cut_uncut_df$cut_uncut_all_4 <- with(accs_cut_uncut_df, eff_cuts / eff_coverage)
    
    
    accs_df_list[[name]] <- cbind(df_X[, c(1, 2, 3)], accs_df, accs_cut_uncut_df)
  }
  
  return(list(accs_df_list = accs_df_list, mus_df_list = mus_df_list, shear_prob_per_bp_list = shear_prob_per_bp_list, bg_cut_per_bp_cov_ratio_list = bg_cut_per_bp_cov_ratio_list))
}


calc_occs_df_list <- function(accs_df_list, min_coverage) {
  occs_df_list <- list()
  for(name in names(accs_df_list)) {
    occs_df <- cbind(accs_df_list[[name]][, c("chr", "enzyme", "pos", "eff_coverage", "eff_cuts")], 1-accs_df_list[[name]][, c("all_mean", "cut_uncut_all_2", "cut_uncut_all_4")])
    occs_df$chr <- as.character(occs_df$chr)
    occs_df$enzyme <- as.character(occs_df$enzyme)
    colnames(occs_df) <- c("chr", "enzyme", "pos", "eff_coverage", "eff_cuts", "occ_X_1", "occ_cut_uncut_2", "occ_cut_uncut_4")
    occs_df <- occs_df[with(occs_df, !is.na(eff_coverage) & eff_coverage >= min_coverage), ]
    occs_df_list[[name]] <- occs_df
  }
  return(occs_df_list)
}



calc_occs_df_list_with_not_corrected_version <- function(accs_df_list, min_coverage = 40) {
  occs_df_list <- list()
  for(name in names(accs_df_list)) {
    occs_df <- cbind(accs_df_list[[name]][, c("chr", "enzyme", "pos", "eff_coverage", "eff_cuts")], 1-accs_df_list[[name]][, c("cut_uncut_all_3", "all_mean", "cut_uncut_all_4")])
    occs_df$chr <- as.character(occs_df$chr)
    occs_df$enzyme <- as.character(occs_df$enzyme)
    colnames(occs_df) <- c("chr", "enzyme", "pos", "eff_coverage", "eff_cuts", "occ_cut_uncut_not_corrected", "occ_cut_all_cut", "occ_cut_uncut_corrected")
    occs_df <- occs_df[with(occs_df, !is.na(eff_coverage) & eff_coverage >= min_coverage), ]
    occs_df_list[[name]] <- occs_df
  }
  return(occs_df_list)
}


plot_accessibility_histograms <- function(accs_df_list, plot_folder, acc_limits = c(-0.5, 2), name_add = "", min_coverage = 0) {
  
  plot_folder = paste0(plot_folder, "accessibility_histograms", name_add, "/")
  dir.create(plot_folder)
  
  for(name in names(accs_df_list)) {
    
    png(paste0(plot_folder, "accessibility_histogram_", name, "_min_cov_", min_coverage, ".png"), res=144, width=2000, height=1800)
    par(mfrow=c(5, 5), oma=c(2, 0, 3, 0))
    
    for(type in c("starts_plus", "starts_minus", "ends_plus", "ends_minus", "all_mean", 
                  "cut_uncut_starts_plus_1", "cut_uncut_ends_plus_1", "cut_uncut_starts_minus_1", "cut_uncut_ends_minus_1", "cut_uncut_all_1",
                  "cut_uncut_starts_plus_2", "cut_uncut_ends_plus_2", "cut_uncut_starts_minus_2", "cut_uncut_ends_minus_2", "cut_uncut_all_2",
                  "cut_uncut_starts_3", "cut_uncut_starts_3", "cut_uncut_ends_3", "cut_uncut_ends_3", "cut_uncut_all_3",
                  "cut_uncut_starts_4", "cut_uncut_starts_4", "cut_uncut_ends_4", "cut_uncut_ends_4", "cut_uncut_all_4")) {
      
      if(type == "empty") {
        plot.new()
      } else {
        accs_df <- accs_df_list[[name]]
        accs <- accs_df[accs_df$eff_coverage >= min_coverage, type]
        accs <- accs[!is.na(accs)]
        accs <- set_outliers_to_limits(accs, acc_limits)
        
        hist(accs, main = type, xlim = acc_limits, breaks = seq(acc_limits[1], acc_limits[2], 0.05))
      }
    }
    
    mtext(name, side=3, line=0, outer=TRUE, cex=1.5, font=1)
    dev.off()
  }
}


plot_accessibility_comparisons <- function(accs_df_list, plot_folder, acc_limits = c(-0.5, 2), name_add = "", min_coverage = 0, cov_plot_lim = c(0, 400)) {
  
  plot_folder = paste0(plot_folder, "accessibility_comparisons", name_add, "/")
  dir.create(plot_folder)
  
  helper_plot <- function(...) {
    plot(..., pch=20, col=adjustcolor(palette()[1], alpha.f=min(0.2, 0.2*5000/length(list(...)[[1]]))))
  }
  
  helper_plot_2 <- function(types, xlim = acc_limits, ylim = acc_limits) {
    for(A in 1:length(types)) {
      for(B in A:length(types)) {
        if(A != B) {
          type_A <- types[A]
          type_B <- types[B]
          
          accs_df <- accs_df_list[[name]]
          
          accs_A <- accs_df[, type_A]
          accs_A[accs_df$eff_coverage < min_coverage] <- NA
          accs_A <- set_outliers_to_limits(accs_A, xlim)
          
          accs_B <- accs_df[, type_B]
          accs_B[accs_df$eff_coverage < min_coverage] <- NA
          accs_B <- set_outliers_to_limits(accs_B, ylim)
          
          accs_not_na <- !is.na(accs_A) & !is.na(accs_B)
          accs_A <- accs_A[accs_not_na]
          accs_B <- accs_B[accs_not_na]
          
          helper_plot(accs_A, accs_B, xlim = xlim, ylim = ylim, xlab = type_A, ylab = type_B)
        }
      }
    }
  }
  
  for(name in names(accs_df_list)) {
    
    png(paste0(plot_folder, "accessibility_comparison_", name, "_min_cov_", min_coverage, ".png"), res=100, width=2300, height=1500)
    par(mfrow=c(4, 7), oma=c(2, 0, 3, 0))
    
    helper_plot_2(c("starts_plus", "starts_minus", "ends_plus", "ends_minus"))
    helper_plot_2(c("eff_coverage", "all_mean"), xlim = cov_plot_lim)
    
    helper_plot_2(c("cut_uncut_starts_plus_2", "cut_uncut_starts_minus_2", "cut_uncut_ends_plus_2", "cut_uncut_ends_minus_2"))
    helper_plot_2(c("eff_coverage", "cut_uncut_all_2"), xlim = cov_plot_lim)
    
    helper_plot_2(c("cut_uncut_starts_2", "cut_uncut_ends_2", "cut_uncut_starts_4", "cut_uncut_ends_4"))
    helper_plot_2(c("eff_coverage", "cut_uncut_all_4"), xlim = cov_plot_lim)
    
    helper_plot_2(c("all_mean", "cut_uncut_all_1", "cut_uncut_all_2"))
    helper_plot_2(c("all_mean", "cut_uncut_all_3", "cut_uncut_all_4"))
    helper_plot_2(c("cut_uncut_all_2", "cut_uncut_all_4"))
    
    mtext(name, side=3, line=0, outer=TRUE, cex=1.5, font=1)
    dev.off()
  }
}


readTable2SparseGenome <- function(filename, chromosomes, col_chr=1, col_pos=2, col_data=3, ...) {
  genome_table <- read.table(filename, ...)
  
  genome_list <- list()
  for(chr in levels(as.factor(genome_table[,col_chr]))) {
    genome_list[[chr]] <- as.matrix(genome_table[ (genome_table[,col_chr]==chr), c(col_pos, col_data)])
  }
  invisible(genome_list)
}


plot_single_genes <- function(single_genes, occs_df_list, plot_folder, occ_limits = c(-0.25, 1.25)) {
  
  plot_folder <- paste0(plot_folder, "single_genes_and_promoters/")
  dir.create(plot_folder)
  
  types <- c("occ_X_1", "occ_cut_uncut_2", "occ_cut_uncut_4")
  
  for(gene in rownames(single_genes)) {
    
    for(name in names(occs_df_list)) {
      
      row_number <- ceiling(length(types))
      
      sites <- subset(occs_df_list[[name]], chr == single_genes[gene, "chr"])[, "pos"]
      gene_start <- single_genes[gene, 'xlim_start']
      gene_end <- single_genes[gene, 'xlim_end']
      
      if(sum(sites >= gene_start & sites <= gene_end) > 0) {
        png(file=paste0(plot_folder, "single_gene_occupancy_", gene, "_", name,".png"), res=144, width=2000, height=600*row_number)
        
        par(mfrow=c(row_number, 1))
        
        for(type in types) {  
          occs <- subset(occs_df_list[[name]], chr == single_genes[gene, "chr"])[, type]
          occs <- set_outliers_to_limits(occs, occ_limits)
          
          plot(sites, occs, main=paste(gene, name, type), xlim=c(gene_start-100, gene_end+100), ylim=occ_limits, xlab=single_genes[gene, 'chr'], ylab="occupancy")
          text(occs~sites, labels = round(occs, digits = 3), pos = 3)
          grid()
        }
        dev.off()
      } else {
        cat("Plot omitted:", name, "has no sites on", gene, "\n")
      }
    }
  }
}


load_mnase_occs_genomic <- function(chrs_main) {  # new 2018-02-25
  
  mnase_tab <- readTable2SparseGenome("../../external_data/Kaplan_inVivo.train.tab", chromosomes=chrs_main)
  names(mnase_tab) <- chrs_main
  
  mnase_occs_genomic <- lapply( lapply(convertSparse2Complete_ff(mnase_tab, chr_sizes[chrs_main]), function(x) {x[, 2]}), smear_ff, from=-73, to=73)
  mnase_occs_genomic <- lapply( mnase_occs_genomic, cut, breaks=c(seq(0, 0.8*max(sapply(mnase_occs_genomic, max)), length.out=254), max(sapply(mnase_occs_genomic, max))+1) )
  
  return(mnase_occs_genomic)
}


load_aligned_mnase_occs <- function(genome, genes_to_align, chrs_main, alignment_window) {  # new 2018-02-25
  
  mnase_occs_genomic <- load_mnase_occs_genomic(chrs_main)
  
  mnase_occ_aligned_chr_list <- list()
  for(chr in chrs_main) {
    
    forward_selection <- (genes_to_align[,1]==chr) & (genes_to_align[,4]=='+') & (genes_to_align[,2] > alignment_window) & (genes_to_align[,2] < length(genome[[chr]])-alignment_window)
    reverse_selection <- (genes_to_align[,1]==chr) & (genes_to_align[,4]=='-') & (genes_to_align[,3] > alignment_window) & (genes_to_align[,3] < length(genome[[chr]])-alignment_window)
    
    forward_mnase <- sapply(genes_to_align[forward_selection,2], function (x) as.numeric(mnase_occs_genomic[[chr]][x+(-alignment_window:alignment_window)]) )
    reverse_mnase <- sapply(genes_to_align[reverse_selection,3], function (x) rev(as.numeric(mnase_occs_genomic[[chr]][x+(-alignment_window:alignment_window)])) )
    
    mnase_occ_aligned_chr_list[[chr]] <- cbind(forward_mnase, reverse_mnase)
  }
  mnase_occs_aligned <- do.call(cbind, mnase_occ_aligned_chr_list)
  
  return(mnase_occs_aligned)
}


align_occupancies <- function(occs_df_list, genome, genes_to_align, chr_sizes, alignment_window) {  # checked 2018-02-25
  
  occs_aligned <- list()
  
  for(name in names(occs_df_list)) {
    
    occs_aligned[[name]] <- list()
    
    for(occ_type in grep("occ", colnames(occs_df_list[[name]]), value = TRUE)) {
      
      occs_df_chr_list <- split(occs_df_list[[name]], as.character(occs_df_list[[name]]$chr))
      
      occs_genomic <- list()
      for(chr in names(occs_df_chr_list)) {
        occs_genomic[[chr]] <- rep(NA, length(genome[[chr]]))
        occs_genomic[[chr]][occs_df_chr_list[[chr]]$pos] <- occs_df_chr_list[[chr]][, occ_type]
      }
      
      occs_aligned_chr_list <- list()
      
      for(chr in names(occs_df_chr_list)) {
        
        forward_selection <- (genes_to_align[,1]==chr) & (genes_to_align[,4]=='+') & (genes_to_align[,2] > alignment_window) & (genes_to_align[,2] < length(occs_genomic[[chr]])-alignment_window)
        reverse_selection <- (genes_to_align[,1]==chr) & (genes_to_align[,4]=='-') & (genes_to_align[,3] > alignment_window) & (genes_to_align[,3] < length(occs_genomic[[chr]])-alignment_window)
        
        forward <- sapply(genes_to_align[forward_selection,2], function (x) occs_genomic[[chr]][x+(-alignment_window:alignment_window)])
        reverse <- sapply(genes_to_align[reverse_selection,3], function (x) rev(occs_genomic[[chr]][x+(-alignment_window:alignment_window)]) )
        
        occs_aligned_chr_list[[chr]] <- cbind(forward, reverse)
      }
      occs_aligned[[name]][[occ_type]] <- do.call(cbind, occs_aligned_chr_list)
    }
  }
  
  return(occs_aligned)
}


running_mean_ignore_na <- function(vec, smooth, ...) {  # checked 2018-02-14
  smooth_minus_one <- smooth-1
  
  na_pos <- is.na(vec)
  vec[na_pos] <- 0
  
  if(smooth_minus_one > 0) {
    if(smooth_minus_one < length(vec)) {
      return( smear(vec, from=-ceiling(smooth_minus_one/2), to=floor(smooth_minus_one/2)) / smear(as.numeric(!na_pos), from=-ceiling(smooth_minus_one/2), to=floor(smooth_minus_one/2)) )
      #return(smear(vec, from=-ceiling(smooth_minus_one/2), to=floor(smooth_minus_one/2)) / (1+c(ceiling(smooth_minus_one/2):(smooth_minus_one), rep((smooth_minus_one), length(vec)-smooth_minus_one-2) , (smooth_minus_one):floor(smooth_minus_one/2))) )
    } else {
      return(rep(mean(vec, ...), length(vec)))
    }
  } else {
    return(vec)
  }
}


plot_promoter_regions <- function(occs_aligned, mnase_occs_aligned, alignment_window, plot_folder, occ_limits = c(-0.25, 1.25)) {  # checked 2018-02-14
  
  plot_folder <- paste0(plot_folder, "promoter_regions/")
  dir.create(plot_folder)
  
  require(dichromat)
  
  smoothing_size <- 41
  
  x_window <- -alignment_window:alignment_window
  
  avg_gene_occs <- list()
  avg_gene_occs_smoothed <- list()
  
  for(name in names(occs_aligned)) {
    
    avg_gene_occs[[name]] <- list()
    avg_gene_occs_smoothed[[name]] <- list()
    
    png(file=paste0(plot_folder, "promoter_region_posav+smoothed", smoothing_size, "_", name, ".png"), res=216, width=2000, height=3000)
    
    occ_types <- names(occs_aligned[[name]])
    
    par(mfrow=c(length(occ_types), 1))
    for(occ_type in occ_types) {
      average_occ <- rowMeans(occs_aligned[[name]][[occ_type]], na.rm = TRUE)
      average_occ_smoothed <- running_mean_ignore_na(average_occ, smoothing_size)
      plot(x_window,average_occ,main=paste0(name, ":  ", occ_type), typ='p', ylim=occ_limits, ylab="occupancy", xlab="distance to TSS", pch=20, col=adjustcolor('black', alpha.f=0.2))
      # lines(smooth.spline(x_window[!is.na(average_occ)], average_occ[!is.na(average_occ)]), col="blue", lwd=4)
      lines(x_window, average_occ_smoothed, col="blue", lwd=4)
      avg_gene_occs[[name]][[occ_type]] <- average_occ
      avg_gene_occs_smoothed[[name]][[occ_type]] <- average_occ_smoothed
    }
    dev.off()
    
    png(paste0(plot_folder, "promoter_region_allgenes+polyfit_", name, ".png"), res=216, width=2000, height=3000)
    
    par(mfrow=c(length(occ_types), 1))
    for(occ_type in occ_types) {
      gene_table <- occs_aligned[[name]][[occ_type]]
      plot(rep(x_window,ncol(gene_table)), gene_table, main=paste0(name, ":  ", occ_type), typ='p', ylim=occ_limits, ylab="occupancy", xlab="distance to TSS", pch=20, col=adjustcolor('black', alpha.f=0.2))
      # lines(smooth.spline(rep(x_window,ncol(gene_table))[!is.na(gene_table)], gene_table[!is.na(gene_table)]), col="blue", lwd=4)
      if(sum(!is.na(gene_table)) >= 100) {
        span <- 0.1
        lines(x_window, predict(loess(gene_table[!is.na(gene_table)]~rep(x_window,ncol(gene_table))[!is.na(gene_table)], span=span), x_window), col="blue", lwd=4)  # strange fitting
      }
    }
    dev.off()
    
    # colour_range <- colorRampPalette(dichromat::colorschemes$BrowntoBlue.12, space = "Lab")(256)
    #
    # png(paste0(plot_folder, "promoter_region_allgenes_mnase_coloured_", sites_name, ".png"), res=216, width=2048, height=1024*length(occ_names))
    # par(mfrow=c(length(occ_names), 1))
    # for(name in occ_names) {
    #   gene_table <- occs_aligned[[name]]
    #   plot(rep(x_window,ncol(gene_table))[!is.na(gene_table)], gene_table[!is.na(gene_table)], main=name, typ='p', ylim=occ_limits, ylab="occupancy", xlab="distance to TSS", col=adjustcolor(colour_range[mnase_occs_aligned[!is.na(gene_table)]], alpha.f=0.5), pch=20)
    # }
    # dev.off()
  }
  
  # comparison of averages
  
  group_names_all <- gsub("_0_U$|_100_U$|_300_U$|_400_U$|_0u$|_100u$|_300u$|_400u$|_0$|_100$|_200$|_400$|_800$|_1000$|_1600$|_4000$|.0.1$|.0.X$|.20.1$|.20.X$|.5.1$|.5.X", "zzzz", names(occs_aligned))  # replaces inconsistent nomenclature with abcd which is cut off in the next step
  group_names_all <- gsub("zzzz", "", group_names_all)  # groups together BY_AluI_0, BY_AluI_1, BY_AluI_3
  group_names <- unique(group_names_all)
  
  for(group_name in group_names) {
    
    this_group <- names(occs_aligned)[group_name == group_names_all]
    
    # col_scheme <- dichromat::colorschemes$Categorical.12[1:length(this_group)]
    col_scheme <- c("yellow", "orange", "red", "violet")[1:length(this_group)]
    names(col_scheme) <- this_group
    
    png(file=paste0(plot_folder, "promoter_region_comparison_posav_smoothed", smoothing_size, "_", group_name, ".png"), res=216, width=2000, height=3000)
    
    par(mfrow=c(length(occ_types), 1))
    for(occ_type in occ_types) {
      
      plot(NA, xlim=range(x_window), main=paste0(group_name, ":  ", occ_type), typ='p', ylim=occ_limits, ylab="occupancy", xlab="distance to TSS", pch=20, col=adjustcolor('black', alpha.f=0.2))
      
      polygon(c(min(x_window), x_window,max(x_window)), y=c(0, rowMeans(mnase_occs_aligned, na.rm=TRUE) / max(mnase_occs_aligned, na.rm=TRUE)*2, 0), col="lightgray", border=NA)
      
      for(name in this_group) {
        average_occ <- rowMeans(occs_aligned[[name]][[occ_type]], na.rm = TRUE)
        lines(x_window, running_mean_ignore_na(average_occ, smoothing_size), col=col_scheme[name], lwd=2)
        # tmp <- loess(y~x, data.frame(x=x_window, y=average_occ), span=0.25)
        # lines(x_window, predict(tmp, x_window), col=col_scheme[name], lwd=2)
      }
      legend("bottomleft", legend=gsub("_filtered_close", "", this_group), fill=col_scheme, ncol=ceiling(length(this_group)/3))
      
    }
    dev.off()
  }
  
  return(list(avg_gene_occs = avg_gene_occs, avg_gene_occs_smoothed = avg_gene_occs_smoothed))
}


plot_PE_fragment_length_histograms <- function(data_names, chr_name_list, max_fragment_size = 500, 
                                               low_coverage_cutoff = 10000, plot_folder = "", RData_folder = "data/RData_files/") {  # checked 2018-02-23
  
  plot_folder <- paste0(plot_folder, "read_length_histograms/")
  dir.create(plot_folder)
  
  data_names_low_coverage <- c()
  fragment_num_list <- list()
  mean_fragment_length_list <- list()
  
  for(name in unique(sub("(_X$|_1$)", "", data_names))) {
    
    png(paste0(plot_folder, "fragment_length_histogram_", name, ".png"), res=144, width=1024*2, height=2*1024)
    par(mfrow=c(4,2))
    
    for(name_X_1 in paste0(name, c("_X", "_1"))) {
      
      if(!(name_X_1 %in% data_names)) stop(paste0("Didn't find ", name_X_1, " in data_names!"))
      
      load(paste0(RData_folder, name_X_1, ".RData"))
      
      for(genome_type in names(chr_name_list)) {
        
        frag_lengths <- c()
        
        for(chr in chr_name_list[[genome_type]]) {
          frag_lengths <- c(frag_lengths, bam[[chr]]$end - bam[[chr]]$start)  # ends were shifted by + 1 in save_PE_bam_as_chr_df_list(...)
        }
        
        if(length(frag_lengths) <= low_coverage_cutoff) {  # do not continue to analyze samples with too few reads in either genome
          cat("Low read count of", length(frag_lengths), "for", name_X_1, genome_type, "\n")
          data_names_low_coverage <- c(data_names_low_coverage, name_X_1)
        }
        
        fragment_num_list[[paste0(name_X_1, "|", genome_type)]] <- length(frag_lengths[frag_lengths <= max_fragment_size])
        mean_fragment_length_list[[paste0(name_X_1, "|", genome_type)]] <- mean(frag_lengths[frag_lengths <= max_fragment_size])
        
        hist(frag_lengths[frag_lengths <= max_fragment_size], xlim = c(0, max_fragment_size), breaks = seq(0, max_fragment_size, 5), main = paste(name_X_1, genome_names[genome_type]), xlab = "Fragment lengths")
        if(any(frag_lengths > max_fragment_size)) {
          hist(log10(frag_lengths[frag_lengths > max_fragment_size]), xlim = c(2, 8), breaks = seq(2, 8, 0.05), main = paste(name_X_1, genome_names[genome_type]), xlab = paste0("log10 of fragment lengths that are > ", max_fragment_size))
        } else {
          plot.new()
        }
        
      }
    }
    dev.off()
  }
  
  l <- length(mean_fragment_length_list) / 4
  
  names <- unique(sub("(_X$|_1$)", "", sapply(strsplit(names(mean_fragment_length_list), "|", fixed = TRUE), "[", 1)))
  
  mean_fragment_length_df <- data.frame()
  
  for(name in names) {
    name_X <- paste0(name, "_X")
    name_1 <- paste0(name, "_1")
    
    mean_fragment_length_df <- rbind(mean_fragment_length_df, data.frame(name = name, 
                                                                         num_main_X = fragment_num_list[[paste0(name_X, "|", names(chr_name_list)[1])]],
                                                                         num_main_1 = fragment_num_list[[paste0(name_1, "|", names(chr_name_list)[1])]],
                                                                         num_norm_X = fragment_num_list[[paste0(name_X, "|", names(chr_name_list)[2])]],
                                                                         num_norm_1 = fragment_num_list[[paste0(name_1, "|", names(chr_name_list)[2])]],
                                                                         length_main_X = mean_fragment_length_list[[paste0(name_X, "|", names(chr_name_list)[1])]], 
                                                                         length_main_1 = mean_fragment_length_list[[paste0(name_1, "|", names(chr_name_list)[1])]],
                                                                         length_norm_X = mean_fragment_length_list[[paste0(name_X, "|", names(chr_name_list)[2])]],
                                                                         length_norm_1 = mean_fragment_length_list[[paste0(name_1, "|", names(chr_name_list)[2])]]))
  }
  
  write.table(mean_fragment_length_df, file = paste0(plot_folder, "../mean_fragment_length.tsv"), sep = "\t", row.names = FALSE)
  
  return(list(unique(data_names_low_coverage), mean_fragment_length_df))
}


set_outliers_to_limits <- function(vec, limits) {
  vec[vec < limits[1]] <- limits[1]
  vec[vec > limits[2]] <- limits[2]
  return(vec)
}


plot_genomic_occs <- function(occs_df_list, plot_folder = "", occ_limits = c(-0.25, 1.25)) {
  
  plot_folder <- paste0(plot_folder, "genomic_occupancies/")
  dir.create(plot_folder)
  
  types <- c("occ_X_1", "occ_cut_uncut_2", "occ_cut_uncut_4")
  
  for(name in names(occs_df_list)) {
    
    row_number <- ceiling(length(types))
    png(file=paste0(plot_folder, "genomic_occupancy_", name, ".png"), res=144, width=2000, height=600*row_number)
    
    par(mfrow=c(row_number, 1))
    for(type in types) {  
      occs_df <- occs_df_list[[name]]
      occs <- occs_df[, type]
      occs <- set_outliers_to_limits(occs, occ_limits)
      chr_infos <- as.character(occs_df[, "chr"])
      chr_site_pops <- table(chr_infos)
      chr_number <- length(chr_site_pops)
      occ_means <- sapply(split(occs, chr_infos), mean, na.rm=TRUE)
      if(length(occs) > 0) {
        plot(occs, pch=20, xaxt='n', xlab='', ylab="Occupancy", main=paste0(name, " using ", type), ylim=occ_limits,
             col=rep(rep(c("black", "darkgray"), chr_number)[1:chr_number], chr_site_pops),
             panel.last = {
               axis(side=1, tick = FALSE, at=cumsum(chr_site_pops)-chr_site_pops/2, labels=names(chr_site_pops))
               axis(side=1, at=cumsum(c(0, chr_site_pops)), labels=FALSE)
               segments( y0=occ_means,
                         x0=cumsum(c(0, chr_site_pops[-length(chr_site_pops)])),
                         x1=cumsum(chr_site_pops),
                         col=adjustcolor(c("blue", "darkblue"), alpha.f=0.5), lwd=4 )
             }
        )
      }
    }
    dev.off()
  }
  
}


get_results_dataframe <- function(RE_run, method_name, min_coverage, num_digits = 2) {
  
  load(paste0("../", RE_run, "/", method_name, ".RData"))  # occ_df_list here has low coverage sites excluded
  if(method_name == "window_limit_times_1_max_length_500_close_distances_200_300_background_michael") {
    load(paste0("../", RE_run, "/analysis_results/window_limit_times_1_max_length_500/close_distances_200_300/background_michael/occs_df_list.RData"))  # reload occs_df_list with low coverage sites
  } else {
    stop("implement loading of occs_df_list for this method!")
  }
  samples <- names(occs_df_list)
  names_X <- paste0(samples, "_X")
  names_1 <- paste0(samples, "_1")
  
  # update cut_statistics to changes from 2019-03-10
  cut_statistics_fc <- print_cut_count_statistics(cut_counts_df_list, bg_counts_df_list, fragment_counts_list, count_window_list, chr_name_list, output_file="none") 
  
  results <- data.frame(file_name = samples, run = RE_run, method_name = method_name, stringsAsFactors = FALSE)
  results$sample_problems <- ""
  results$critical_problems_X_1 <- ""
  results$critical_problems_cut_uncut <- ""
  
  results$mio_main_reads_X <- signif(cut_statistics_fc[names_X, "reads_main"]/1000000, digits = num_digits)
  results$mio_main_reads_1 <- signif(cut_statistics_fc[names_1, "reads_main"]/1000000, digits = num_digits)
  
  results$mio_norm_reads_X <- signif(cut_statistics_fc[names_X, "reads_norm"]/1000000, digits = num_digits)
  results$mio_norm_reads_1 <- signif(cut_statistics_fc[names_1, "reads_norm"]/1000000, digits = num_digits)
  
  results$min_spike_in_fraction <- signif(with(results, pmin(mio_norm_reads_1/mio_main_reads_1, mio_norm_reads_X/mio_main_reads_X)), digits = num_digits)
  
  mean_fragment_length_df <- mean_fragment_length_df[, c(1, 6, 7, 8, 9)]
  mean_fragment_length_df[, 2:5] <- signif(mean_fragment_length_df[, 2:5], digits = num_digits+1)
  results <- merge(results, mean_fragment_length_df, by.x = "file_name", by.y = "name", sort = FALSE)  # merge sorts the results with respect to the 'by' column by default!
  
  results$mean_main_cuts_X <- signif(cut_statistics_fc[names_X, "cut_mean_main"], digits = num_digits+1)
  results$mean_main_cuts_1 <- signif(cut_statistics_fc[names_1, "cut_mean_main"], digits = num_digits+1)
  
  results$mean_main_uncut_fragments_X <- signif(cut_statistics_fc[names_X, "uncut_mean_main"], digits = num_digits+1)
  results$mean_main_uncut_fragments_1 <- signif(cut_statistics_fc[names_1, "uncut_mean_main"], digits = num_digits+1)
  
  results$mean_norm_cuts_X <- signif(cut_statistics_fc[names_X, "cut_mean_norm"], digits = num_digits+1)
  results$mean_norm_cuts_1 <- signif(cut_statistics_fc[names_1, "cut_mean_norm"], digits = num_digits+1)
  
  results$mean_norm_uncut_fragments_X <- signif(cut_statistics_fc[names_X, "uncut_mean_norm"], digits = num_digits+1)
  results$mean_norm_uncut_fragments_1 <- signif(cut_statistics_fc[names_1, "uncut_mean_norm"], digits = num_digits+1)
  
  results$mean_main_coverage_X <- signif(with(results, mean_main_cuts_X + 2*mean_main_uncut_fragments_X), digits = num_digits+1)
  results$mean_main_coverage_1 <- signif(with(results, mean_main_cuts_1 + 2*mean_main_uncut_fragments_1), digits = num_digits+1)
  results$mean_norm_coverage_X <- signif(with(results, mean_norm_cuts_X + 2*mean_norm_uncut_fragments_X), digits = num_digits+1)
  results$mean_norm_coverage_1 <- signif(with(results, mean_norm_cuts_1 + 2*mean_norm_uncut_fragments_1), digits = num_digits+1)
  
  results$main_cut_fraction_X <- signif(with(results, mean_main_cuts_X/(mean_main_cuts_X + 2*mean_main_uncut_fragments_X)), digits = num_digits+1)
  results$main_cut_fraction_1 <- signif(with(results, mean_main_cuts_1/(mean_main_cuts_1 + 2*mean_main_uncut_fragments_1)), digits = num_digits+1)
  results$norm_cut_fraction_X <- signif(with(results, mean_norm_cuts_X/(mean_norm_cuts_X + 2*mean_norm_uncut_fragments_X)), digits = num_digits+1)
  results$norm_cut_fraction_1 <- signif(with(results, mean_norm_cuts_1/(mean_norm_cuts_1 + 2*mean_norm_uncut_fragments_1)), digits = num_digits+1)
  
  temp <- mclapply(names_X, function(name) {
    this_file <- paste0("../", RE_run, "/data/RData_files/cut_uncut_fragment_lengths_", name, ".RData")
    if(file.exists(this_file)) {
      load(this_file)
      cut_reads_list <- lapply(cut_reads_list, "[", 1:16)
      uncut_reads_list <- lapply(uncut_reads_list, "[", 1:16)
      mean_cut_reads <- round(mean(unlist(cut_reads_list)))
      mean_uncut_reads <- round(mean(unlist(uncut_reads_list)))
    } else {
      warning(paste0("Didn't find: ", this_file))
      mean_cut_reads <- NA
      mean_uncut_reads <- NA
    }
    return(data.frame(mean_cut_fragment_length_X = mean_cut_reads, mean_uncut_fragment_length_X = mean_uncut_reads))
  }, mc.cores = 8)
  results <- cbind(results, bind_rows(temp))
  
  temp <- mclapply(names_1, function(name) {
    this_file <- paste0("../", RE_run, "/data/RData_files/cut_uncut_fragment_lengths_", name, ".RData")
    if(file.exists(this_file)) {
      load(this_file)
      cut_reads_list <- lapply(cut_reads_list, "[", 1:16)
      uncut_reads_list <- lapply(uncut_reads_list, "[", 1:16)
      mean_cut_reads <- round(mean(unlist(cut_reads_list)))
      mean_uncut_reads <- round(mean(unlist(uncut_reads_list)))
    } else {
      warning(paste0("Didn't find: ", this_file))
      mean_cut_reads <- NA
      mean_uncut_reads <- NA
    }
    return(data.frame(mean_cut_fragment_length_1 = mean_cut_reads, mean_uncut_fragment_length_1 = mean_uncut_reads))
  }, mc.cores = 8)
  results <- cbind(results, bind_rows(temp))
  
  set_null_to_na <- function(x) ifelse(is.null(x), NA, x)
  set_vec_to_na <- function(x) ifelse(length(x)>1, NA, x)
  results$window_limit <- sapply(samples, function(name) set_null_to_na(set_vec_to_na(count_window_df[count_window_df$name == paste0(name, "_X") & count_window_df$close_dist_str== "200_300" & count_window_df$genome_type == "main", "count_window_limit"])))
  
  results$shear_prob_per_bp <- signif(sapply(samples, function(name) set_null_to_na(shear_prob_per_bp_list[[name]])), digits = num_digits+1)
  results$bg_cut_per_bp_cov_ratio <- signif(sapply(samples, function(name) set_null_to_na(bg_cut_per_bp_cov_ratio_list[[name]])), digits = num_digits+1)
  
  results$valid_sites_X_1 <- sapply(samples, function(name) sum(!is.na(occs_df_list[[name]][, "occ_X_1"])))
  results$valid_sites_cut_uncut <- sapply(samples, function(name) sum(!is.na(occs_df_list[[name]][, "occ_cut_uncut_4"]) & occs_df_list[[name]][, "eff_coverage"] >= min_coverage))
  
  results$median_eff_coverage_cut_uncut <- signif(sapply(samples, function(name) median(occs_df_list[[name]][, "eff_coverage"], na.rm = TRUE)), digits = num_digits+1)
  results$mean_eff_coverage_cut_uncut <- signif(sapply(samples, function(name) mean(occs_df_list[[name]][, "eff_coverage"], na.rm = TRUE)), digits = num_digits+1)
  results$mean_eff_cuts <- signif(sapply(samples, function(name) mean(occs_df_list[[name]][, "eff_cuts"], na.rm = TRUE)), digits = num_digits+1)
  
  results$high_cov_frac <- sapply(samples, function(name) sum(occs_df_list[[name]][, "eff_coverage"] >= min_coverage, na.rm = TRUE))
  results$high_cov_frac <- signif(results$high_cov_frac / results$valid_sites_cut_uncut, digits = num_digits+1)
  
  occ_columns <- grep("occ", colnames(occs_df_list[[1]]))
  
  occ_site_means <- data.frame(sapply(occ_columns, function(col) sapply(occs_df_list, function(occs_df) mean(occs_df[, col], na.rm = TRUE))))
  occ_site_sds <- data.frame(sapply(occ_columns, function(col) sapply(occs_df_list, function(occs_df) sd(occs_df[, col], na.rm = TRUE))))
  colnames(occ_site_means) <- colnames(occs_df_list[[1]][occ_columns])
  colnames(occ_site_sds) <- colnames(occs_df_list[[1]][occ_columns])
  
  occ_site_means_high_cov <- data.frame(sapply(occ_columns, function(col) sapply(occs_df_list, function(occs_df) mean(occs_df[occs_df$eff_coverage >= min_coverage, col], na.rm = TRUE))))
  occ_site_sds_high_cov <- data.frame(sapply(occ_columns, function(col) sapply(occs_df_list, function(occs_df) sd(occs_df[occs_df$eff_coverage >= min_coverage, col], na.rm = TRUE))))
  colnames(occ_site_means_high_cov) <- colnames(occs_df_list[[1]][occ_columns])
  colnames(occ_site_sds_high_cov) <- colnames(occs_df_list[[1]][occ_columns])
  
  results$mean_occ_X_1 <- signif(occ_site_means[samples, "occ_X_1"], digits = num_digits+1)
  results$mean_occ_cut_uncut <- signif(occ_site_means_high_cov[samples, "occ_cut_uncut_4"], digits = num_digits+1)
  results$mean_occ_cut_uncut_all_sites <- signif(occ_site_means[samples, "occ_cut_uncut_4"], digits = num_digits+1)
  results$mean_occ_cut_uncut_unfitted <- signif(1-acc_site_means[samples, "cut_uncut_all_3"], digits = num_digits+1)
  results$mean_occ_cut_uncut_alternative <- signif(occ_site_means_high_cov[samples, "occ_cut_uncut_2"], digits = num_digits+1)
  results$mean_occ_cut_uncut_alternative_all_sites <- signif(occ_site_means[samples, "occ_cut_uncut_2"], digits = num_digits+1)
  results$mean_occ_cut_uncut_alternative_unfitted <- signif(1-acc_site_means[samples, "cut_uncut_all_1"], digits = num_digits+1)
  
  results$sd_occ_X_1 <- signif(occ_site_sds[samples, "occ_X_1"], digits = num_digits+1)
  results$sd_occ_cut_uncut <- signif(occ_site_sds_high_cov[samples, "occ_cut_uncut_4"], digits = num_digits)
  results$sd_occ_cut_uncut_all_sites <- signif(occ_site_sds[samples, "occ_cut_uncut_4"], digits = num_digits)
  results$sd_occ_cut_uncut_unfitted <- signif(acc_site_sds[samples, "cut_uncut_all_3"], digits = num_digits)
  results$sd_occ_cut_uncut_alternative <- signif(occ_site_sds_high_cov[samples, "occ_cut_uncut_2"], digits = num_digits)
  results$sd_occ_cut_uncut_alternative_all_sites <- signif(occ_site_sds[samples, "occ_cut_uncut_2"], digits = num_digits)
  results$sd_occ_cut_uncut_alternative_unfitted <- signif(acc_site_sds[samples, "cut_uncut_all_1"], digits = num_digits)
  
  calc_outlier_fraction <- function(occ_vec, acc_min, acc_max) {
    acc_vec <- 1 - occ_vec
    acc_vec <- acc_vec[!is.na(acc_vec)]
    return( sum(acc_vec < acc_min | acc_vec > acc_max) / length(acc_vec))
  }
  
  results$outlier_fraction_X_1 <- round(sapply(samples, function(name) {
    occs_df <- occs_df_list[[name]]
    return(calc_outlier_fraction(occs_df[occs_df$eff_coverage >= min_coverage, "occ_X_1"], -0.05, 1.05))
  }), digits = num_digits+1)
  
  results$outlier_fraction_cut_uncut <- round(sapply(samples, function(name) {
    occs_df <- occs_df_list[[name]]
    return(calc_outlier_fraction(occs_df[occs_df$eff_coverage >= min_coverage, "occ_cut_uncut_4"], -0.05, 1.05))
  }), digits = num_digits+1)
  
  results$resection_X <- signif(sapply(samples, function(name) set_null_to_na(set_vec_to_na(count_window_df[count_window_df$name == paste0(name, "_X") & count_window_df$close_dist_str== "200_300" & count_window_df$genome_type == "main", "mean_resection_length"]))), digits = num_digits+1)
  results$resection_1 <- signif(sapply(samples, function(name) set_null_to_na(set_vec_to_na(count_window_df[count_window_df$name == paste0(name, "_1") & count_window_df$close_dist_str== "200_300" & count_window_df$genome_type == "main", "mean_resection_length"]))), digits = num_digits+1)
  
  calc_CI <- function(cuts, coverage) {
    if(length(cuts) == 0) return(NA_real_)
    
    test <- binom.test(round(coverage - max(0, cuts)), round(coverage), conf.level = 0.95)
    
    return(c(test[["conf.int"]][1], test[["conf.int"]][2]))
  }
  
  calc_single_site_results <- function(site_name, locus_occs_df_list) {
    
    df <- data.frame("occ_X_1" =       sapply(locus_occs_df_list, function(occ_df) ifelse(nrow(occ_df)>0, occ_df$occ_X_1, NA_real_)),
                     "occ_cut_uncut" = sapply(locus_occs_df_list, function(occ_df) ifelse(nrow(occ_df)>0, occ_df$occ_cut_uncut_4, NA_real_)),
                     "CI_low" =  sapply(locus_occs_df_list, function(occ_df) calc_CI(occ_df$eff_cuts, occ_df$eff_coverage)[1]),
                     "CI_high" = sapply(locus_occs_df_list, function(occ_df) calc_CI(occ_df$eff_cuts, occ_df$eff_coverage)[2]))
    df <- signif(df, digits = num_digits+1)
    names(df) <- paste0(site_name, "_", names(df))
    return(df)
  }
  
  locus_occs_df_list <- lapply(occs_df_list, function(occs_df) subset(occs_df, chr == "chr02" & pos == 431491))
  results <- cbind(results, calc_single_site_results("PHO5_prom_Bam", locus_occs_df_list))
  
  locus_occs_df_list <- lapply(occs_df_list, function(occs_df) subset(occs_df, chr == "chr02" & pos == 431591))
  results <- cbind(results, calc_single_site_results("PHO5_prom_Alu", locus_occs_df_list))
  
  locus_occs_df_list <- lapply(occs_df_list, function(occs_df) subset(occs_df, chr == "chr02" & pos == 431707))
  results <- cbind(results, calc_single_site_results("PHO5_prom_Alu_Nm5", locus_occs_df_list))
  
  locus_occs_df_list <- lapply(occs_df_list, function(occs_df) subset(occs_df, chr == "chr04" & pos == 1420409))
  results <- cbind(results, calc_single_site_results("PHO8_prom_Hind", locus_occs_df_list))
  
  Np1_alignment_info_df_list <- mclapply(occs_df_list[!grepl("Asi", occs_df_list)], function(occs_df) {
    occs_df$chr <- sub("chr0", "chr", occs_df$chr)
    occs_df$occ <- occs_df$occ_cut_uncut_4
    occs_df$occ[occs_df$occ < 0] <- 0
    occs_df$occ[occs_df$occ > 1] <- 1
    info_df_c_u <- calc_average_gene_occ_df(occs_df, pos_binning = 10, filter_size = 1, max_promoter_length = 1000, max_ORF_length = 1000)[[3]]
    
    occs_df$occ <- occs_df$occ_X_1
    occs_df$occ[occs_df$occ < 0] <- 0
    occs_df$occ[occs_df$occ > 1] <- 1
    info_df_X_1 <- calc_average_gene_occ_df(occs_df, pos_binning = 10, filter_size = 1, max_promoter_length = 1000, max_ORF_length = 1000)[[3]]
    colnames(info_df_X_1) <- paste0(colnames(info_df_X_1), "_X_1")
    
    return(cbind(info_df_c_u, info_df_X_1[, -1]))
  }, mc.cores = 8)
  
  Np1_alignment_info_df <- signif(bind_rows(Np1_alignment_info_df_list), digits = num_digits+1)
  Np1_alignment_info_df$file_name <- names(occs_df_list[!grepl("Asi", occs_df_list)])
  
  results <- merge(results, Np1_alignment_info_df, all = TRUE)
  
  return(results)
}


calc_deviation_from_test_runs <- function(accs_df_list, background_folder, name_add = "", min_coverage, enzyme = "combined") {
  
  calibration_samples <- c(1, 4, 7, 10:27)
  accs_df_list <- accs_df_list[calibration_samples]
  
  if(enzyme == "combined") {
    calibration_samples_accs <- c(0, 0, 0, rep(c(10, 100, 30, 50, 70, 90) / 100, 3))
  } else if(enzyme == "AluI") {
    accs_df_list <- accs_df_list[grepl("AluI", names(accs_df_list))]
    calibration_samples_accs <- c(0, c(10, 100, 30, 50, 70, 90) / 100)
  } else if(enzyme == "BamHI") {
    accs_df_list <- accs_df_list[grepl("BamHI", names(accs_df_list))]
    calibration_samples_accs <- c(0, c(10, 100, 30, 50, 70, 90) / 100)
  } else if(enzyme == "HindIII") {
    accs_df_list <- accs_df_list[grepl("HindIII", names(accs_df_list))]
    calibration_samples_accs <- c(0, c(10, 100, 30, 50, 70, 90) / 100)
  } else {
    stop("enzyme not valid!")
  }
  
  occ_columns <- grep("(all_mean|cut_uncut_all)", colnames(accs_df_list[[1]]), value = TRUE)
  
  high_cov_fracs <- sapply(accs_df_list, function(accs_df) sum(accs_df$eff_coverage >= min_coverage, na.rm = TRUE) / sum(!is.na(accs_df$eff_coverage)))
  if(any(high_cov_fracs < 0.5)) {
    warning("Some high_cov_fracs < 0.5")
    if(any(high_cov_fracs < 0.25)) stop("Some high_cov_fracs < 0.25")
  }
  
  accs_df_list <- lapply(accs_df_list, function(accs_df) accs_df[accs_df$eff_coverage >= min_coverage, occ_columns])
  
  # mean over samples of mean over sites of diffs
  deviation_matrix_sites <- sapply(occ_columns, function(col) rowMeans(sapply(1:length(accs_df_list), function(sample) {
    diff_sites <- accs_df_list[[sample]][, col] - calibration_samples_accs[sample]
    return(c( mean(diff_sites, na.rm = TRUE), mean(abs(diff_sites), na.rm = TRUE), sqrt(mean((diff_sites)^2, na.rm = TRUE)), mean((diff_sites)^2, na.rm = TRUE) ))
  })))
  rownames(deviation_matrix_sites) <- c("sample_mean_site_mean_diff", "sample_mean_site_mean_abs_diff", "sample_mean_sqrt_site_mean_squ_diff", "sample_mean_site_mean_squ_diff")
  deviation_matrix_sites <- rbind(deviation_matrix_sites, sqrt(deviation_matrix_sites["sample_mean_site_mean_squ_diff", ]))
  rownames(deviation_matrix_sites)[5] <- "sqrt_sample_mean_site_mean_squ_diff"
  deviation_matrix_sites <- deviation_matrix_sites[-which("sample_mean_site_mean_squ_diff" == rownames(deviation_matrix_sites)), ]
  
  # fixed samples, diffs of global mean
  acc_site_means <- data.frame(sapply(occ_columns, function(col) sapply(accs_df_list, function(accs_df) mean(accs_df[, col], na.rm = TRUE))))
  acc_site_sds <- data.frame(sapply(occ_columns, function(col) sapply(accs_df_list, function(accs_df) sd(accs_df[, col], na.rm = TRUE))))
  colnames(acc_site_sds) <- paste0("sd_", colnames(acc_site_sds))
  
  deviation_matrix_samples <- sapply(occ_columns, function(col) sapply(1:nrow(acc_site_means), function(sample) acc_site_means[sample, col] - calibration_samples_accs[sample]))
  rownames(deviation_matrix_samples) <- names(accs_df_list)   
  
  # mean over samples of diffs of global mean
  deviation_matrix_global <- sapply(occ_columns, function(col) {
    diff_samples <- deviation_matrix_samples[, col]
    return(c( mean(abs(diff_samples)), sqrt(mean(diff_samples^2)) ))
  })
  rownames(deviation_matrix_global) <- c("sample_mean_abs_diff_site_mean", "sqrt_sample_mean_squ_diff_site_mean")
  
  deviation_matrix <- rbind(deviation_matrix_sites, deviation_matrix_global)
  deviation_df <- data.frame("error_method" = rownames(deviation_matrix), deviation_matrix, "enzyme" = enzyme, stringsAsFactors = FALSE)
  
  sample_df <- data.frame("file_name" = names(accs_df_list), "enzyme" = enzyme, "calibration_samples_accs" = calibration_samples_accs, acc_site_means, acc_site_sds, stringsAsFactors = FALSE)
  
  #cat(kable(deviation_df), file = paste0(background_folder, "deviation_df_", enzyme, name_add, ".txt"), sep = "\n")
  
  return(list("deviation_df" = deviation_df, "sample_df" = sample_df))
}


run_function_parallel <- function(data_names, f, core_factor = 1) {
  num_cores <- round(detectCores(all.tests = FALSE, logical = TRUE) * core_factor)  # too many cores use too much RAM and then some results files just go missing or are not readable
  
  if(Sys.info()['sysname'] == "Linux") {
    result_list <- mclapply(data_names, f, mc.cores = num_cores)
  } else {  # needs to be tested on Windows
    cl <- makeCluster(num_cores)
    result_list <- parLapply(cl, data_names, f)
    stopCluster(cl)
  }
  
  if(length(result_list) < length(data_names)) {
    warning("not all cores finished properly!")
  }
  
  return(unlist(result_list, recursive = FALSE))
}


load_RE_occs_df <- function(run, names, method_folder = "window_limit_times_1_max_length_500/close_distances_200_300/background_michael/", as_list = FALSE) {
  load(paste0("../", run, "/analysis_results/", method_folder, "occs_df_list.RData"))
  if(as_list == FALSE) {
    return(occs_df_list[[names]])
  } else {
    if(!missing(names)) return(occs_df_list[names])
    else return(occs_df_list)
  }
}


add_ <- function(char_vec, sep = "_") {
  return(ifelse(char_vec == "", "", paste0(sep, char_vec)))
}


calc_fragment_strand_pos_length_lists <- function(bam) {
  fragment_start_length_lists <- list("+" = list(), "-" = list())
  fragment_end_length_lists <- list("+" = list(), "-" = list())
  for(chr in names(bam)) {
    this_bam <- subset(bam[[chr]], strand == "+")
    fragment_start_length_lists[["+"]][[chr]] <- split(as.integer(this_bam$end - this_bam$start), this_bam$start)  # ends were shifted by + 1 in save_PE_bam_as_chr_df_list(...)
    fragment_end_length_lists[["+"]][[chr]] <- split(as.integer(this_bam$end - this_bam$start), this_bam$end)
    
    this_bam <- subset(bam[[chr]], strand == "-")
    fragment_start_length_lists[["-"]][[chr]] <- split(as.integer(this_bam$end - this_bam$start), this_bam$start)
    fragment_end_length_lists[["-"]][[chr]] <- split(as.integer(this_bam$end - this_bam$start), this_bam$end)
  }
  return(list(fragment_start_length_lists = fragment_start_length_lists, fragment_end_length_lists = fragment_end_length_lists))
}


calc_fragment_pos_length_lists <- function(bam) {
  fragment_start_length_lists <- list()
  fragment_end_length_lists <- list()
  for(chr in names(bam)) {
    this_bam <- subset(bam[[chr]])
    fragment_start_length_lists[[chr]] <- split(as.integer(this_bam$end - this_bam$start), this_bam$start)
    fragment_end_length_lists[[chr]] <- split(as.integer(this_bam$end - this_bam$start), this_bam$end)
  }
  return(list(fragment_start_length_lists = fragment_start_length_lists, fragment_end_length_lists = fragment_end_length_lists))
}


plot_cut_uncut_length_distributions <- function(name, cut_counts_df_list, count_window_list, RE_info, max_length = 500, max_plot_length = 300, plot_folder, name_add = "") {
  
  library("viridis")
  
  plot_folder <- paste0(plot_folder, "cut_uncut_lengths_distribution", name_add, "/")
  dir.create(plot_folder)
  
  load(paste0("data/RData_files/", name, ".RData"))
  bam <- lapply(bam, function(bam_chr) bam_chr[bam_chr$end - bam_chr$start <= max_length, ])
  temp <- calc_fragment_pos_length_lists(bam)
  fragment_start_length_lists <- temp[["fragment_start_length_lists"]]
  fragment_end_length_lists <- temp[["fragment_end_length_lists"]]
  rm(temp)
  
  count_window <- count_window_list[[name]]
  cut_reads_list <- list()
  uncut_reads_list <- list()
  for(this_enzyme in unique(as.character(cut_counts_df_list[[name]]$enzyme))) {
    site_shift <- RE_info[this_enzyme, "site_shift"]
    cut_reads_list[[this_enzyme]] <- list()
    uncut_reads_list[[this_enzyme]] <- list()
    cut_counts_df <- subset(cut_counts_df_list[[name]], enzyme = this_enzyme)
    for(this_chr in unique(cut_counts_df$chr)) {
      
      # cut fragments:
      
      cut_reads_list[[this_enzyme]][[this_chr]] <- list()
      
      cut_counts_chr_df <- subset(cut_counts_df, chr == this_chr & enzyme == this_enzyme & !is.na(starts_plus))
      starts <- as.integer(names(fragment_start_length_lists[[this_chr]]))
      
      cut_reads_list[[this_enzyme]][[this_chr]][["starts"]] <- lapply(count_window[[this_enzyme]], function(w) {
        list_indeces <- starts %in% (cut_counts_chr_df$pos - site_shift + w)
        temp <- fragment_start_length_lists[[this_chr]][list_indeces]
        if(length(temp) > 0) names(temp) <- paste0(as.character(as.numeric(names(temp)) + site_shift - w), "_")
        return(unlist(temp))
      })
      names(cut_reads_list[[this_enzyme]][[this_chr]][["starts"]]) <- count_window[[this_enzyme]]
      
      cut_counts_chr_df <- subset(cut_counts_df, chr == this_chr & enzyme == this_enzyme & !is.na(ends_plus))
      ends <- as.integer(names(fragment_end_length_lists[[this_chr]]))
      cut_reads_list[[this_enzyme]][[this_chr]][["ends"]] <- lapply(count_window[[this_enzyme]], function(w) {
        list_indeces <- ends %in% (cut_counts_chr_df$pos + site_shift - w)
        temp <- fragment_end_length_lists[[this_chr]][list_indeces]
        if(length(temp) > 0) names(temp) <- paste0(as.character(as.numeric(names(temp)) - site_shift + w), "_")
        return(unlist(temp))
      })
      names(cut_reads_list[[this_enzyme]][[this_chr]][["ends"]]) <- count_window[[this_enzyme]]
      
      # uncut fragments:
      
      uncut_reads_list[[this_enzyme]][[this_chr]] <- list()
      
      cut_counts_chr_df <- subset(cut_counts_df, chr == this_chr & enzyme == this_enzyme & !(is.na(ends_plus) & is.na(starts_plus)))
      
      bam_chr <- bam[[this_chr]]
      bam_chr <- bam_chr[order(bam_chr$start, bam_chr$end), ]
      
      max_indeces_to_test <- findInterval(cut_counts_chr_df$pos - site_shift - 1, bam_chr$start)  # test all starts for all sites
      
      min_indeces_to_test <- findInterval(cut_counts_chr_df$pos - site_shift - 1 - max_length + 1, bam_chr$start)  # gives all starts that could have valid end locations
      
      calc_covering_fragments_lengths <- function(i) {
        indeces <- min_indeces_to_test[i]:max_indeces_to_test[i]
        good_indeces <- indeces[bam_chr$end[indeces] > cut_counts_chr_df$pos[i] + site_shift]
        return(as.integer(bam_chr$end[good_indeces] - bam_chr$start[good_indeces]))  # now also test all ends for each site
      }
      
      uncut_reads_list[[this_enzyme]][[this_chr]] <- lapply(1:length(cut_counts_chr_df$pos), calc_covering_fragments_lengths)
      names(uncut_reads_list[[this_enzyme]][[this_chr]]) <- paste0(as.character(cut_counts_chr_df$pos), "_")
      uncut_reads_list[[this_enzyme]][[this_chr]] <- unlist(uncut_reads_list[[this_enzyme]][[this_chr]])
    }
    
    cut_reads <- unlist(cut_reads_list[[this_enzyme]][1:16])
    uncut_reads <- unlist(uncut_reads_list[[this_enzyme]][1:16])
    
    png(file = paste0(plot_folder, "cut_uncut_fragment_lengths_", name, "_", this_enzyme, ".png"), width = 1400, height = 1400, res = 150)
    
    par(mfrow = c(3, 1), oma=c(1, 0, 2, 0))
    
    hist(uncut_reads[uncut_reads <= max_plot_length], breaks = seq(0, max_plot_length, 1), main = paste0("uncut reads: mean length = ", round(mean(uncut_reads))))
    
    hist(cut_reads[cut_reads <= max_plot_length], breaks = seq(0, max_plot_length, 1), main = paste0("cut reads: mean length = ", round(mean(cut_reads))))
    
    cut_reads_w_list <- lapply(count_window[[this_enzyme]], function(w) cut_reads[grep(paste0(".", w, "."), names(cut_reads), fixed = TRUE)])
    
    good_w <- sapply(cut_reads_w_list, function(x) length(x) >= 10)  # only plot if there are at least 10 fragments at w (density need at least 2)
    density_w_list <- lapply(cut_reads_w_list[good_w], density)
    
    colors <- viridis(length(density_w_list))
    
    plot(density_w_list[[1]], col = adjustcolor(col = colors[1], alpha.f = 0.5), main = paste0("uncut reads: densities for w = ", min(count_window[[this_enzyme]]), ":", max(count_window[[this_enzyme]]), " (increasing with brightness)"), xlim = c(0, max_plot_length), ylim = c(0, max(sapply(density_w_list, function(dens) max(dens$y)))))
    if(length(density_w_list) > 1) {
      for(i in seq(2, length(density_w_list), 1)) {
        lines(density_w_list[[i]], col = adjustcolor(col = colors[i], alpha.f = 0.5))
      }
    }
    
    mtext(paste0(name), side=3, line=0, outer=TRUE, cex=1, font=1)
    dev.off()
  }
  
  save(cut_reads_list, uncut_reads_list, file = paste0("data/RData_files/cut_uncut_fragment_lengths", name_add, "_", name, ".RData"))
}

