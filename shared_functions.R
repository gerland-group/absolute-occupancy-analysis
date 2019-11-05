combine_cols_for_plot <- function(df, new_col_names, col_types) {  # restructure df to plot multiple cols of different types marked by for example different color
  
  # example: df with cols x_1, x_2 to be restructured to col x and col type with entries "type_1" for x_1 values and "type_2" for x_2 values. 
  # Then use new_col_names = c("x", "type") and col_types = c("x_1" = "type_1", "x_2" = "type_2").
  
  if(length(new_col_names) != 2) stop("Need two new column names!")
  
  df_list <- list()
  for(col in names(col_types)) {
    df_new_cols <- data.frame(df[, col], col_types[col], row.names = NULL, stringsAsFactors = FALSE)
    colnames(df_new_cols) <- new_col_names
    df_list[[col]] <- cbind(df, df_new_cols)
  }
  
  df_new <- bind_rows(df_list)
  df_new[, names(col_types)] <- NULL
  
  return(df_new)
}


helper_plot <- function(... , col=1, alpha.f=0.2, pch=19) {
  if(class(col) == "character") {
    plot(... ,  pch=pch, col=adjustcolor(col, alpha.f=alpha.f))
  } else {
    plot(... ,  pch=pch, col=adjustcolor(palette()[col], alpha.f=alpha.f))
  }
}


compare_two_occs_dfs <- function(method, sample_1, df_1, sample_2, df_2, min_coverage = NA, plot_folder = NA, 
                                 no_title = FALSE, my_alpha = NA, use_pdf = FALSE, full_info = !no_title, flip_axes = FALSE, name_add = "") {
  
  # Sample properties:
  
  # RE sites are at least 200bp apart from each other. For this comparison we define a RE site area as the 10bp where the cut happens between 5bp and 6bp.
  
  # NP sites (called by nanopolish) are at least 11bp apart, as methylation sites that are closer are grouped together to one site. 
  # For nanopolish the sites on the + and - strand in the "CG" motif are always grouped together. 
  # NP dfs should have a start and end column such that end - start is at least 10bp
  
  # BS sites are arbitrary close to each other ("CG" already has two neighbouring sites with distance 1, the C on the + strand and the C on the - strand)
  
  # methods:
  
  # RE_vs_BS_5: BS sites within a fixed RE site are averaged (i.e. 10bp bins at RE cut positions)
  # RE_vs_RE_5: RE_1 and RE_2 sites are compared if their area overlap is at least 5 (i.e. 10bp bins at RE_1 cut positions)
  # BS_vs_BS_0: Compare identical sites 
  # BS_vs_BS_5: Compare average in bins of 10bp (genomic bins)
  
  # RE_vs_NP_5: RE and NP sites if the RE site touches the NP site
  # NP_vs_BS_5: Average all BS sites within the NP site (i.e. bins at NP sites with boundarys start and end)
  # NP_vs_NP_0: Compare identical sites (since so far only CpG NP data is available)
  
  # All sites need to have minimal coverage of min_coverage if count data is provided
  
  colnames(df_1)[grepl("mean_occ", colnames(df_1))] <- "occ"  # RE maps have a mean_occ column instead of an occ column
  colnames(df_2)[grepl("mean_occ", colnames(df_2))] <- "occ"
  
  colnames(df_1)[grepl("mean_eff_coverage", colnames(df_1))] <- "counts"  # RE maps have a mean_eff_coverage column instead of an counts column
  colnames(df_2)[grepl("mean_eff_coverage", colnames(df_2))] <- "counts"
  
  if(is.finite(min_coverage)) {
    if("counts" %in% colnames(df_1)) {
      cat("Removed", round(100*sum(df_1$counts < min_coverage)/nrow(df_1)), "% low coverage sites in sample 1.\n")
      df_1 <- subset(df_1, counts >= min_coverage)
    } else {
      cat("No counts column in df_1, min_coverage not checked.\n")
    }
    if("counts" %in% colnames(df_2)) {
      cat("Removed", round(100*sum(df_2$counts < min_coverage)/nrow(df_2)), "% low coverage sites in sample 2.\n")
      df_2 <- subset(df_2, counts >= min_coverage)
    } else {
      cat("No counts column in df_2, min_coverage not checked.\n")
    }
  }
  
  radius <- as.numeric(sub("[A-Za-z_]*", "", method))
  method_str <- sub("_[0-9]*$", "", method)
  
  if(grepl("(RE_vs_NP|NP_vs_BS)", method) & radius != 5) {  # NP methods have fixed radius given by the NP calling resolution
    stop("NP method with radius != 5 not implemented")
  }
  
  if(grepl("^RE", method) | method == "BS_vs_BS_0" | method == "NP_vs_BS_5" | method == "NP_vs_NP_0") {
    if(method_str == "RE_vs_BS") {
      bin <- -radius:(radius-1)
    } else if(method_str == "RE_vs_RE" | method == "BS_vs_BS_0" | method == "NP_vs_NP_0") {
      bin <- -radius:radius
    } else if (method == "NP_vs_BS_5") {
    } else {
      print(method)
      stop("Invalid method!")
    }
    
    pos_valid_list <- list()
    occs_1_list <- list()
    occs_2_bin_mean_list <- list()
    for(this_chr in unique(df_1$chr)) {
      
      df_1_chr <- subset(df_1, chr == this_chr & !is.na(occ))
      df_2_chr <- subset(df_2, chr == this_chr & !is.na(occ))
      
      pos_1_good <- df_1_chr$pos
      pos_2_good <- df_2_chr$pos
      
      occs_1 <- df_1_chr$occ
      
      if(method != "NP_vs_BS_5") {
        occs_2_vec <- rep(NA, max(c(pos_1_good + radius, pos_2_good)))
        occs_2_vec[pos_2_good] <- df_2_chr$occ
        
        occs_2_bin_mean <- sapply(pos_1_good, function(x) mean(occs_2_vec[x+bin], na.rm = TRUE))
      } else {
        occs_2_vec <- rep(NA, max(c(df_1_chr$end, pos_2_good)))
        occs_2_vec[pos_2_good] <- df_2_chr$occ
        
        occs_2_bin_mean <- mapply(function(s, e) mean(occs_2_vec[s:e], na.rm = TRUE), df_1_chr$start, df_1_chr$end)
      }
      
      pos_valid <- pos_1_good[is.finite(occs_2_bin_mean)]
      occs_2_bin_mean <- occs_2_bin_mean[is.finite(occs_2_bin_mean)]
      
      pos_valid_list[[this_chr]] <- pos_valid
      occs_1_list[[this_chr]] <- subset(df_1_chr, pos %in% pos_valid)$occ
      occs_2_bin_mean_list[[this_chr]] <- occs_2_bin_mean
    }
    occs_1 <- unlist(occs_1_list)
    occs_2 <- unlist(occs_2_bin_mean_list)
    
  } else if(grepl("^BS_vs_BS", method) & radius > 0) {  # bin BS sites of all samples (e.g radius=5: pos_binned = 10 <==> pos = 5,6,...,13,14)
    
    pos_binning <- 2*radius
    
    calc_binned_occs_df <- function(df) {
      df$pos_binned <- round(df$pos/pos_binning) * pos_binning
      df$chr_pos <- with(df, paste0(chr, "_", pos_binned))
      
      df_binned <- aggregate(df$occ, list(df$chr_pos), mean, na.rm = TRUE)
      colnames(df_binned) <- c("chr_pos", "occ")
      df_binned <- df_binned[!is.na(df_binned$occ), ]
      df_binned <- df_binned[order(df_binned$chr_pos), ]
      return(df_binned)
    }
    
    df_1 <- calc_binned_occs_df(df_1)
    df_2 <- calc_binned_occs_df(df_2)
    
    valid_chr_pos <- intersect(df_1$chr_pos, df_2$chr_pos)
    
    occs_1 <- with(df_1, occ[chr_pos %in% valid_chr_pos])
    occs_2 <- with(df_2, occ[chr_pos %in% valid_chr_pos])
    
  } else {
    print(method)
    stop("Method not valid or to be implemented!")
  }
  
  delta = mean(abs(occs_1 - occs_2))
  delta_2 = mean((occs_1 - occs_2)^2)
  
  r_value = cor(occs_1, occs_2)
  d_value = delta + abs(mean(occs_1) - mean(occs_2))
  
  num_comparisons <- length(occs_1)
  
  # plotting
  
  if(!is.na(plot_folder)) {
    
    if(is.na(my_alpha)) my_alpha <- min(0.3, 2500/num_comparisons)
    
    my_alpha <- max(0.01, my_alpha)
    
    this_title <- paste0(sample_1, " vs ", sample_2)
    
    occs_df <- data.frame(occs_1 = occs_1, occs_2 = occs_2)
    
    geom_text_size <- theme_paper[[1]]$axis.text$size
    
    occ_limits <- c(0, 1)
    p <- ggplot(occs_df, aes(x=occs_1, y=occs_2)) + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",  color = "blue", size = 0.4) +
      geom_point(alpha = my_alpha, size = 0.4, shape = 16, stroke = 0) + 
      scale_x_continuous(labels = percent_no_sign, limits = (occ_limits)) + scale_y_continuous(labels = percent_no_sign, limits = (occ_limits)) +
      xlab(paste0("abs. occupancy ", sample_1, " / %")) + ylab(paste0("abs. occupancy ", sample_2, " / %")) + theme_paper
    
    if(!flip_axes) {
      plot_infos_1 <- paste0("<|x-y|> = ", sprintf("%.3f", round(delta, 3)), 
                             ",  <x> = ", sprintf("%.3f", round(mean(occs_df$occs_1), 3)), 
                             ",  <y> = ", sprintf("%.3f", round(mean(occs_df$occs_2), 3)), 
                             ",  d = ", sprintf("%.3f", round(d_value, 3)),
                             ",\nradius = ", radius, ",  ", num_comparisons, " sites,  r = ", sprintf("%.3f", round(r_value, 3)))
      
      plot_infos_2 <- paste0(num_comparisons, " sites",
                             "\n<x> = ", sprintf("%.3f", round(mean(occs_df$occs_1), 3)), 
                             ",  <y> = ", sprintf("%.3f", round(mean(occs_df$occs_2), 3)),
                             ",  <|x-y|> = ", sprintf("%.3f", round(delta, 3)),
                             ",  cor = ", sprintf("%.3f", round(r_value, 3)))
    } else {
      p <- p + coord_flip()
      plot_infos_1 <- paste0("<|x-y|> = ", sprintf("%.3f", round(delta, 3)), 
                             ",  <x> = ", sprintf("%.3f", round(mean(occs_df$occs_2), 3)), 
                             ",  <y> = ", sprintf("%.3f", round(mean(occs_df$occs_1), 3)), 
                             ",  d = ", sprintf("%.3f", round(d_value, 3)),
                             ",\nradius = ", radius, ",  ", num_comparisons, " sites,  r = ", sprintf("%.3f", round(r_value, 3)))
      
      plot_infos_2 <- paste0(num_comparisons, " sites",
                             "\n<x> = ", sprintf("%.3f", round(mean(occs_df$occs_2), 3)), 
                             ",  <y> = ", sprintf("%.3f", round(mean(occs_df$occs_1), 3)),
                             ",  <|x-y|> = ", sprintf("%.3f", round(delta, 3)),
                             ",  cor = ", sprintf("%.3f", round(r_value, 3)))
        
    }

    if(full_info) {
      plot_infos <- plot_infos_1
    } else {
      plot_infos <- plot_infos_2
    }
    
    if(!no_title) {
      p <- p + labs(title = this_title, subtitle = plot_infos)
    } else {
      p <- p + labs(subtitle = plot_infos)
    }
    
    if(use_pdf) {
      ggsave(plot = p, filename = paste0(plot_folder, "scatterplot_", name_add, method, "_" , sample_1, "_vs_", sample_2, "_min_cov_", min_coverage, ".pdf"),
             width = 6, height = 6, units = "cm")
    } else {
      ggsave(plot = p, filename = paste0(plot_folder, "scatterplot_", name_add, method, "_" , sample_1, "_vs_", sample_2, "_min_cov_", min_coverage, ".png"),
             width = 6, height = 6, units = "cm", type = "cairo")
    }
  }
  
  comparison_table <- data.frame(method=method, sample_1=sample_1, sample_2=sample_2, mean_occ_1=mean(occs_1), mean_occ_2=mean(occs_2), radius = radius, min_coverage = min_coverage,
                                 delta=delta, delta_2=delta_2, r_value=r_value, d_value = d_value, num_comparisons=num_comparisons,
                                 stringsAsFactors = FALSE)
  return(comparison_table)
}


add_ <- function(char_vec, sep = "_") {
  return(ifelse(char_vec == "", "", paste0(sep, char_vec)))
}


calc_tss_Np1_distance <- function(vec, pos, s) {
  if(s == "+") {
    d <- vec - pos
  } else if(s == "-") {
    d <- pos - vec
  } else {
    stop("strand not valid!")
  }
  return(min(d[d>=0]))
}


load_scer_gene_data <- function(data_set = "BY0_CL_pm9") {
  start_path <- "../../"
  
  
  if(data_set == "transcripts_steinmetz") {  # TSS alignment
    gene_data <- read.delim(paste0(start_path, "external_data/transcripts_R64-1-1_steinmetz.txt"))
    gene_data$chr <- paste0("chr", as.numeric(as.roman(sub("chr", "", gene_data$chr))))
    gene_data$aligned_pos <- with(gene_data, ifelse(strand == "+", start, end))
  
 
  } else {

    gene_data <- read.delim(paste0(start_path, "external_data/+1coordiantesETC_tirosh_32U.tab"))
    gene_data$BY0_all <- gene_data$plus1NK
    gene_data$plus1NK <- NULL
    
    gene_data$chr <- paste0("chr", as.numeric(as.roman(sub("chr", "", gene_data$chr))))
    
    Np1_files <- grep(".rda", dir(paste0(start_path, "external_data/called_dyads_and_Np1")), value = TRUE)
    for(f in Np1_files) {
      load(paste0(start_path, "external_data/called_dyads_and_Np1/", f))
      col <- sub(".rda", "", f)
      df <- centers[, c("ID", "plus1")]
      colnames(df)[2] <- col
      gene_data <- merge(gene_data, df[, c("ID", col)], all.x = TRUE, by = "ID")
    }
    
   
    gene_data$aligned_pos <- gene_data[, data_set]
    gene_data$TSS <- gene_data$tss
    gene_data$tss <- NULL
    gene_data$TTS <- gene_data$tts
    gene_data$tts <- NULL
  }
  
  if(any(duplicated(gene_data$ID))) stop("Found duplicated gene IDs!")
  
  gene_data$name <- as.character(gene_data$name)
  
  if("TSS" %in% colnames(gene_data) & "TTS" %in% colnames(gene_data)) {
    gene_data <- gene_data[!is.na(gene_data$aligned_pos), c("ID", "chr", "start", "end", "strand", "name", "aligned_pos", "TSS", "TTS")]
  } else {
    gene_data <- gene_data[!is.na(gene_data$aligned_pos), c("ID", "chr", "start", "end", "strand", "name", "aligned_pos")]
  }
  gene_data <- gene_data[order(gene_data$chr, gene_data$start), ]
  rownames(gene_data) <- NULL
  
  return(gene_data)
}


load_dyad_positions <- function(data_set) {
  if(data_set == "BY0_CL_pm9") {
    dyad_df <- read.delim("../../external_data/called_dyads_and_Np1/BY0_CL_pm9.csv", sep = ";", header = TRUE)
    dyad_df$chr <- sub("chr0", "chr", dyad_df$chr)
    dyad_df <- subset(dyad_df, grepl("chr[0-9]*$", chr))
    dyad_df <- dyad_df[, c("chr", "smt_pos")]
    colnames(dyad_df) <- c("chr", "aligned_pos")
    dyad_df$name <- paste0("nucl", 1:nrow(dyad_df))
  } else if(data_set == "chereji_email") {
    dyad_df <- read.delim("../../external_data/called_dyads_and_Np1/Chereji.bed", sep = "\t", header = TRUE)
    dyad_df$start <- dyad_df$start_bed + 1
    dyad_df$aligned_pos <- with(dyad_df, floor((start + end) / 2))
    dyad_df$chr <- paste0("chr", as.numeric(as.roman(sub("chr", "",  dyad_df$chr))))
    dyad_df <- dyad_df[, c("chr", "aligned_pos", "name")]
    dyad_df$name <- as.character(dyad_df$name)
  } else if(data_set == "zhang") {
    dyad_df <- read.delim("../../external_data/called_dyads_and_Np1/Zhang.csv", sep = "\t", header = TRUE)
    dyad_df$chr <- paste0("chr", as.numeric(as.roman(sub("chr", "",  dyad_df$chr))))
    dyad_df <- subset(dyad_df, grepl("chr[0-9]*$", chr))
    dyad_df <- dyad_df[, c("chr", "smt_pos")]
    colnames(dyad_df) <- c("chr", "aligned_pos")
    dyad_df$name <- paste0("nucl", 1:nrow(dyad_df))
  }
  dyad_df <- dyad_df[order(dyad_df$chr, dyad_df$aligned_pos), ]
  dyad_df$ID <- dyad_df$name
  dyad_df$strand <- "+"
  
  unique(dyad_df$chr)
  sapply(unique(dyad_df$chr), function(this_chr) max(subset(dyad_df, chr == this_chr)[, "aligned_pos"]))
  any(duplicated(dyad_df$name))
  
  return(dyad_df)
}


calc_aligned_positions <- function(occs_df, max_promoter_length = 1000, max_ORF_length = 1000, cut_at_gene_ends = FALSE, data_set = "BY0_CL_pm9", 
                                   dyads_instead_of_genes = FALSE) {
  
  occs_df <- occs_df[!is.na(occs_df$occ), grepl("^(chr|pos|occ|sd_occ|counts)$", colnames(occs_df))]  
  
  chrs <- sort(unique(as.character(occs_df$chr)))
  
  expected_labels <- sort(c(paste0("chr", 1:9), paste0("chr", 10:16)))  # NO LEADING ZERO
  if(length(chrs)==0 || !all(diag(sapply(expected_labels, grepl, chrs)))) {
    stop("Check chr labels in input occs_df! Not all 16 chromosomes were found!\n")
  }
  
  occs_df_chr_list <- split(occs_df, list(occs_df$chr))
  if(nrow(occs_df) > 1500000) {
    method <- "access_extended_df"
    occs_df_chr_list <- lapply(occs_df_chr_list, function(occs_df) merge(occs_df, data.frame(pos = 1:max(occs_df$pos)), by = "pos", all = TRUE))
  } else {
    method <- "subset_normal_df"
  }
  occs_df_chr_list <- lapply(occs_df_chr_list, function(occs_df) occs_df[order(occs_df$pos), ])
  rm(occs_df)
  
  if(!dyads_instead_of_genes) {
    gene_data <- load_scer_gene_data(data_set)
  } else {
    gene_data <- load_dyad_positions(data_set)
  }

  aligned_genes_list <- list()
  for (i in 1:nrow(gene_data)) {
    chr_temp <- gene_data$chr[i]
    aligned_pos <- gene_data[i, "aligned_pos"]
    if (gene_data[i, "strand"] == '+') {
      genome_pos_min <- aligned_pos - max_promoter_length
      if(cut_at_gene_ends) {
        genome_pos_max <- min(max_ORF_length + aligned_pos, gene_data[i, "end"])
      } else {
        genome_pos_max <- max_ORF_length + aligned_pos
      }
      if(method == "subset_normal_df") {
        new_gene <- subset(occs_df_chr_list[[chr_temp]], pos >= genome_pos_min & pos <= genome_pos_max)
      } else if(method == "access_extended_df") {
        new_gene <- occs_df_chr_list[[chr_temp]][genome_pos_min:genome_pos_max, ]
        new_gene <- new_gene[!is.na(new_gene$chr), ]
      } else {
        stop("invalid method!")
      }
      new_gene$pos = new_gene$pos - aligned_pos
    } else {
      if(cut_at_gene_ends) {
        genome_pos_min <- max(gene_data[i, "start"], aligned_pos - max_ORF_length)
      } else {
        genome_pos_min <- aligned_pos - max_ORF_length
      }
      genome_pos_max <- aligned_pos + max_promoter_length
      if(method == "subset_normal_df") {
        new_gene <- subset(occs_df_chr_list[[chr_temp]], pos >= genome_pos_min & pos <= genome_pos_max)
      } else if(method == "access_extended_df") {
        new_gene <- occs_df_chr_list[[chr_temp]][genome_pos_min:genome_pos_max, ]
        new_gene <- new_gene[!is.na(new_gene$chr), ]
      } else {
        stop("invalid method!")
      }
      new_gene$pos = aligned_pos - new_gene$pos
      new_gene = new_gene[nrow(new_gene):1,]  # revert pos order (not needed)
    }
    if(nrow(new_gene)>0) {
      new_gene$gene_number = i
      new_gene$gene_name = gene_data[i, "name"]
      new_gene$gene_id = gene_data[i, "ID"]
      new_gene$chr <- NULL
      aligned_genes_list[[i]] <- new_gene
    }
  }
  aligned_genes <- bind_rows(aligned_genes_list)  # makes one data.frame from the list
  
  aligned_genes
}


calc_average_gene_occ_df_from_genes_df <- function(genes_df) {
  calc_average_gene_occ_df(genes_df = genes_df)
}


calc_average_gene_occ_df <- function(occs_df, genes_df, pos_binning = 10, filter_size = 1, max_promoter_length = 1000, max_ORF_length = 1000, 
                                     cut_at_gene_ends = FALSE, normalize_cov_for_each_gene = FALSE, data_set = "BY0_CL_pm9", dyads_instead_of_genes = FALSE) {
  if(missing(genes_df) & !missing(occs_df)) {
    genes_df_given <- FALSE
    if(!("occ" %in% colnames(occs_df)) & ("cov" %in% colnames(occs_df))) {
      coverage_given <- TRUE
      occs_df$occ <- occs_df$cov  # treat coverage as occupancy
      occs_df$cov <- NULL
    } else {
      coverage_given <- FALSE
    }
    genes_df <- calc_aligned_positions(occs_df, max_promoter_length, max_ORF_length, cut_at_gene_ends, data_set, dyads_instead_of_genes = dyads_instead_of_genes)
   
  } else {
    genes_df_given <- TRUE
    if(!("occ" %in% colnames(genes_df)) & ("cov" %in% colnames(genes_df))) {
      coverage_given <- TRUE
      genes_df$occ <- genes_df$cov  # treat coverage as occupancy
      genes_df$cov <- NULL
    } else {
      coverage_given <- FALSE
    }
  }
  
  if(coverage_given) {
    if(normalize_cov_for_each_gene) {
      genes_df_list <- split(genes_df, list(genes_df$gene_id))
      genes_df_list <- lapply(genes_df_list, function(df) {
        df$occ <- df$occ / mean(df$occ)
        return(df)
      })
      genes_df <- bind_rows(genes_df_list)
    } else {  # normalize with the global mean
      genes_df$occ <- genes_df$occ / mean(genes_df$occ, na.rm = TRUE)
    }
  }
  
  genes_df$pos_binned = round(genes_df$pos/pos_binning, 0) * pos_binning
  
  av_gene_df = aggregate(genes_df[, c('occ')], list(genes_df$pos_binned), mean, na.rm = TRUE)
  colnames(av_gene_df) <- c("pos", "occ")

  av_gene_df_temp = aggregate(genes_df[, c('occ')], list(genes_df$pos_binned), length)
  colnames(av_gene_df_temp) <- c("pos", "num_sites_in_bin")
  av_gene_df <- merge(av_gene_df, av_gene_df_temp)
  
  av_gene_df_temp = aggregate(genes_df[, c('occ')], list(genes_df$pos_binned), sd, na.rm = TRUE)
  colnames(av_gene_df_temp) <- c("pos", "sd_of_binned_occs")
  av_gene_df <- merge(av_gene_df, av_gene_df_temp)
  
  if("counts" %in% colnames(genes_df)) {
    av_gene_df_temp = aggregate(genes_df[, c('counts')], list(genes_df$pos_binned), mean, na.rm = TRUE)
    colnames(av_gene_df_temp) <- c("pos", "mean_counts_per_site")
    av_gene_df <- merge(av_gene_df, av_gene_df_temp)
  }
  
  if("sd_occ" %in% colnames(genes_df)) {
    av_gene_df_temp = aggregate(genes_df[, c('sd_occ')], list(genes_df$pos_binned), mean, na.rm = TRUE)
    colnames(av_gene_df_temp) <- c("pos", "mean_sd_occ")
    av_gene_df <- merge(av_gene_df, av_gene_df_temp)
    
    av_gene_df_temp = aggregate(genes_df[, c('sd_occ')], list(genes_df$pos_binned), function(x) sqrt(sum(x^2, na.rm = TRUE))/sum(!is.na(x)))
    colnames(av_gene_df_temp) <- c("pos", "sd_of_av_occ")
    av_gene_df <- merge(av_gene_df, av_gene_df_temp)
  }
  
  av_gene_df$occ_smoothed <- NA
  av_gene_df$num_sites_smoothed <- NA
  n <- nrow(av_gene_df)
  av_gene_df$occ_smoothed[ceiling(filter_size/2):(n-floor(filter_size/2))] <- rollmean(av_gene_df$occ, filter_size)
  av_gene_df$num_sites_smoothed[ceiling(filter_size/2):(n-floor(filter_size/2))] <- rollmean(av_gene_df$num_sites_in_bin, filter_size)
  
  if(!coverage_given & !dyads_instead_of_genes & data_set != "transcripts_steinmetz") {
  
    lag_min_max <- round(c(100/pos_binning, 250/pos_binning))
    
    autocor <- acf(av_gene_df$occ[av_gene_df$pos >= 0 & av_gene_df$pos <= min(1000, max_ORF_length)], na.action=na.pass, lag.max = lag_min_max[2], plot = FALSE)[["acf"]][lag_min_max[1]:lag_min_max[2]]
    autocor_max <- max(autocor)
    autocor_lag <- lag_min_max[1]-1 + mean(which(autocor == max(autocor)))
    
    # max min differences
    #period <- lag_min_max[1]-1 + which(autocor == max(autocor))[1]
    period <- round(180 / pos_binning)
    if(is.finite(period)) {
      maxima <- rep(NA, 3)
      minima <- rep(NA, 3)
      for(i in 1:3) {
        max_area <- which(av_gene_df$pos == 0) + period * (i-1) + seq(-round(period/2), round(period/2), 1)
        min_area <- which(av_gene_df$pos == 0) + period * (i-1) + seq(0, period, 1)
        maxima[i] <- max(av_gene_df$occ[max_area])
        minima[i] <- min(av_gene_df$occ[min_area])
      }
      mean_max_min_diff <- mean(maxima - minima)
    } else {
      mean_max_min_diff <- NA
    }
  } else {
    autocor_lag <- NA
    autocor_max <- NA
    mean_max_min_diff <- NA
  }
  if(coverage_given) {  # bin coverage to reduce memory needed
    temp_df <- aggregate(genes_df[, c('occ')], list(paste0(genes_df$gene_id, "__", genes_df$gene_name, "__", genes_df$pos_binned)), mean, na.rm = TRUE)
    names(temp_df) <- c("info", "cov")
    temp_info <- str_split(temp_df$info, pattern = "__")
    temp_df$gene_id <- sapply(temp_info, "[", 1)
    temp_df$gene_name <- sapply(temp_info, "[", 2)
    temp_df$pos <- as.integer(sapply(temp_info, "[", 3))
    temp_df$info <- NULL
    genes_df <- temp_df
    
    av_gene_df$cov <- av_gene_df$occ
    av_gene_df$occ <- NULL
  } else {
    genes_df$pos_binned <- NULL
  }
  
  if(!genes_df_given) {
    genes_df$gene_id <- as.character(genes_df$gene_id)
    return(list("av_gene_df" = av_gene_df, "genes_df" = genes_df, "info_df" = data.frame(pos_binning, autocor_lag, autocor_max, mean_max_min_diff)))
  } else {
    return(list("av_gene_df" = av_gene_df, "info_df" = data.frame(pos_binning, autocor_lag, autocor_max, mean_max_min_diff)))
  }
}


add_sample_names_to_aligned_dfs_list <- function(aligned_dfs_list, name_add = "") {
  
  saved_names <- names(aligned_dfs_list)
  aligned_dfs_list <- lapply(names(aligned_dfs_list), function(name) {
    aligned_dfs_list[[name]][[1]]$sample <- paste0(name_add, name)
    aligned_dfs_list[[name]][[2]]$sample <- paste0(name_add, name)
    aligned_dfs_list[[name]][[3]] <- cbind(data.frame(sample = paste0(name_add, name)), aligned_dfs_list[[name]][[3]])
    return(aligned_dfs_list[[name]])
  })
  names(aligned_dfs_list) <- saved_names
  return(aligned_dfs_list)
}


plot_sorted_genes <- function(genes_df, sorted_genes_df, list_name, map_name, plot_folder, pos_binnings = c(5, 10, 20), file_type = "png_heatmap", 
                              window = NA, plot_quintiles_heatmap = TRUE) {
  
  occ_given <- "occ" %in% colnames(genes_df)
  if(!occ_given) {
    genes_df$occ <- genes_df$cov
  }
  
  if(missing(sorted_genes_df)) {
    if(all(!is.na(window))) {
      genes_df_to_sort <- genes_df[genes_df$pos %in% window, ]
    } else {
      genes_df_to_sort <- genes_df
    }
    sorted_genes_df <- aggregate(genes_df_to_sort$occ, list(genes_df_to_sort$gene_id), mean, na.rm = TRUE)
    colnames(sorted_genes_df) <- c("gene_id", "mean_signal")
    sorted_genes_df <- sorted_genes_df[order(sorted_genes_df$mean_signal), ]
    sorted_genes_df$sorted_number <- 1:nrow(sorted_genes_df)
    sorted_genes_df$sample <- map_name
    if(missing(list_name)) list_name <- paste0("self_sorted_from_", min(window), "_to_", max(window))
    
    min_cov <- 1.5*quantile(sorted_genes_df$mean_signal, 0.01, na.rm = TRUE)
    max_cov <- 1.5*quantile(sorted_genes_df$mean_signal, 0.99, na.rm = TRUE)
    
    if(!is.na(plot_folder)) {
      p <- ggplot(sorted_genes_df, aes(x=mean_signal)) + geom_histogram(bins = 40) + theme_paper +
        xlab(paste0(ifelse(occ_given, "mean abs. gene occupancy / %", "mean norm. gene coverage"), 
                    ifelse(all(!is.na(window)), paste0(" (", min(window), " to ", max(window), ")"), ""))) +
        ggtitle(map_name)
      if(occ_given) {
        p <- p + scale_x_continuous(labels = percent_no_sign, limit = c(0, 1), breaks = seq(0, 1, 0.25))
      } else {
        p <- p + scale_x_continuous(limit = c(ifelse(min_cov < 0, min_cov, 0), max_cov))
      }
      if(file_type == "png") {
        ggsave(plot = p, paste0(plot_folder, "signal_histogram_", list_name, "_using_", map_name, ".png"), width = 7, height = 5, units = "cm", type = "cairo")
      } else {
        ggsave(plot = p, paste0(plot_folder, "signal_histogram_", list_name, "_using_", map_name, ".pdf"), width = 7, height = 5, units = "cm")
      }
    }
  }
  
  if(plot_quintiles_heatmap & !is.na(plot_folder)) {
  
    genes_df <- merge(genes_df, sorted_genes_df)
  
    for(pos_binning in pos_binnings) {
      genes_df$pos_binned <- round(genes_df$pos/pos_binning) * pos_binning
      binned_df <- unique(genes_df[, c("gene_id", "sorted_number", "pos_binned")])
      binned_df$ID <- with(binned_df, paste0(gene_id, "_", pos_binned))
      temp_df <- aggregate(genes_df$occ, list(with(genes_df, paste0(gene_id, "_", pos_binned))), mean)
      colnames(temp_df) <- c("ID", "occ")
      binned_df <- merge(binned_df, temp_df, by = "ID")
      rm(temp_df)
      max_cov <- quantile(binned_df$occ, probs = 0.95, na.rm = TRUE)
  
      ggplot(binned_df, aes(x=pos_binned, y=sorted_number, fill=occ)) + geom_raster() +
        scale_fill_viridis(direction = -1, limits = c(0, ifelse(occ_given, 1, max_cov))) +
        theme_paper + labs(fill=ifelse(occ_given, "abs. occupancy", "norm. coverage"))
      if(file_type %in% c("png", "png_heatmap")) {
        ggsave(paste0(plot_folder, "heatmap_", list_name, "_using_", map_name, "_pos_binning_", pos_binning, ".png"), width = 15, height = 25, units = "cm", type = "cairo")
      } else {
        ggsave(paste0(plot_folder, "heatmap_", list_name, "_using_", map_name, "_pos_binning_", pos_binning, ".pdf"), width = 15, height = 25, units = "cm")
      }
    }
  
    genes_df$quintile <- with(genes_df, ceiling(5 * sorted_number / max(sorted_number)))
    genes_df_list <- split(genes_df, list(genes_df$quintile))
    temp_list <- lapply(genes_df_list, calc_average_gene_occ_df_from_genes_df)
    quintiles_df <- bind_rows(lapply(names(temp_list), function(name) {
      df <- temp_list[[name]][[1]]
      df$quintile <- name
      return(df)
    }))
    rm(temp_list)
    max_cov <- max(quintiles_df$occ, na.rm = TRUE)
  
    p <- ggplot(quintiles_df, aes(x=pos, y=occ, color=quintile)) +
      geom_line(size = 0.35, alpha = 0.8) + scale_color_viridis(direction = -1, discrete = TRUE, end = 0.9) + theme_paper +
      ylab(ifelse(occ_given, "gene averaged abs. occupancy / %", "gene averaged norm. coverage")) + xlab("distance to in vivo +1 nucleosome / bp")
    if(occ_given) {
      p <- p + scale_y_continuous(labels = percent_no_sign, limit = c(0, 1), breaks = seq(0, 1, 0.25))
    } else {
      p <- p + scale_y_continuous(limit = c(0, max_cov), breaks = seq(0, max_cov, 0.5))
    }
    if(file_type == "png") {
      ggsave(plot = p, paste0(plot_folder, "averaged_quintiles_", list_name, "_using_", map_name, ".png"), width = 4, height = 4, units = "cm", type = "cairo")
    } else {
      ggsave(plot = p, paste0(plot_folder, "averaged_quintiles_", list_name, "_using_", map_name, ".pdf"), width = 4, height = 4, units = "cm")
    }
  }
  
  return(sorted_genes_df)
}


save_RE_occs_df_as_table_and_bedgraph <- function(df, file_name, output_path = getwd()) {
  
  df[, c("eff_coverage", "eff_cuts", "occ_cut_uncut_not_corrected", "occ_cut_all_cut", "occ_cut_uncut_corrected")] <- round(df[, c("eff_coverage", "eff_cuts", "occ_cut_uncut_not_corrected", "occ_cut_all_cut", "occ_cut_uncut_corrected")], digits = 4)
  
  write.table(x = df, file = paste0(output_path, file_name, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  df_2 <- data.frame(seqnames = paste0("chr", as.roman(gsub("chr", "", df$chr))), start=df$pos, end=df$pos, score = round(df$occ_cut_uncut_corrected*100, digits = 2))
  
  df_2$score[df_2$score > 100] <- 100  # for RE samples, BS samples are always between 0 and 100%
  df_2$score[df_2$score < 0] <- 0
  
  gr <- as(df_2, "GRanges")
  export.bedGraph(gr, paste0(output_path, file_name, ".bedgraph"))
  
  return(NULL)
}


save_BS_occs_df_as_table_and_bedgraph <- function(df, min_cov = 20, file_name, output_path = getwd()) {
  if("counts" %in% colnames(df)) {  # for BS samples
    df <- df[df$counts >= min_cov & is.finite(df$occ), ]
  } else {
    df <- df[is.finite(df$occ), ]
  }
  
  write.table(x = df, file = paste0(output_path, file_name, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  df_2 <- data.frame(seqnames = paste0("chr", as.roman(gsub("chr","" ,df$chr))), start=df$pos, end=df$pos, strand = df$strand, score = round(df$occ*100, digits = 2))
  
  gr <- as(df_2, "GRanges")
  export.bedGraph(gr, paste0(output_path, file_name, ".bedgraph"))
  
  return(NULL)
}

