load_site_stats <- function(site_stats_file) {
  site_stats <- read.table(file = site_stats_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE);

  names(site_stats)[1] <- "ref_seq_detailed"
  site_stats$ref_seq <- sub(".*\\|", "", site_stats$ref_seq_detailed)
  site_stats$ref_seq_detailed <- as.factor(site_stats$ref_seq_detailed)
  site_stats$ref_seq <- as.factor(site_stats$ref_seq)
  
  site_stats$hits_called <- site_stats$hits_total - site_stats$hits_ambiguous
  site_stats$hit_fraction_methylated <- site_stats$hits_methylated / site_stats$hits_total
  site_stats$hit_call_fraction_methylated <- site_stats$hits_methylated / site_stats$hits_called
  site_stats$hit_fraction_ambiguous <- site_stats$hits_ambiguous / site_stats$hits_total
  
  # Shifting start and end positions by +1 to give correct methylation cytosine position
  site_stats$start <- site_stats$start + 1
  site_stats$end <- site_stats$end + 1
  
  site_stats$position <- (site_stats$start + site_stats$end) / 2
  
  site_stats
}


load_barcoded_site_stats <- function(site_stats_file, barcode_file, project_folder = "../") {
  barcode_infos <- read.table(paste0(project_folder, barcode_file), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  barcode_folders <- barcode_infos$barcode

  site_stats_list <- vector(mode = "list", length = length(barcode_folders))
  for(i in 1:length(barcode_folders)) {
    barcode_folder <- barcode_folders[i]
    site_stats <- load_site_stats(site_stats_file = paste0(project_folder, barcode_folder, "/", site_stats_file))
    site_stats$ref_seq_detailed <- as.character(site_stats$ref_seq_detailed)
    site_stats$ref_seq <- as.character(site_stats$ref_seq)
    site_stats$barcode <- sub("barcode", "", barcode_folder)
    site_stats$sample <- barcode_infos$sample[i]
    if("methylation_time" %in% colnames(barcode_infos)) {
      site_stats$methylation_time <- barcode_infos$methylation_time[i]
    } else {
      site_stats$methylation_time <- NA
    }
    site_stats_list[[i]] <- site_stats
  }
  
  site_stats_combined <- bind_rows(site_stats_list) 
  site_stats_combined$barcode <- factor(site_stats_combined$barcode, sub("barcode", "", barcode_folders))
  site_stats_combined$sample <- factor(site_stats_combined$sample, barcode_infos$sample)
  if("methylation_time" %in% colnames(barcode_infos)) {
    site_stats_combined$methylation_time  <-  factor(site_stats_combined$methylation_time, unique(barcode_infos$methylation_time))
  }
  site_stats_combined$ref_seq_detailed <- as.factor(site_stats_combined$ref_seq_detailed)
  site_stats_combined$ref_seq <- as.factor(site_stats_combined$ref_seq)
  
  site_stats_combined
}
