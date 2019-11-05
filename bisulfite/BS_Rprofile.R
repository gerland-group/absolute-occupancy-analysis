#### Install packages

# install.packages(c("ff", "ffbase", "knitr", "dichromat", "parallel", "zoo", "viridis", "cowplot", "tictoc", "stringi", "profvis"), dependencies = TRUE)

# installing Mark Heron's packages
# install.packages("devtools")
# devtools::install_bitbucket(repo="markheron/maRs", keep_source=TRUE)
# devtools::install_bitbucket(repo="markheron/pRon", keep_source=TRUE)
# devtools::install_bitbucket(repo="markheron/nucular", keep_source=TRUE)

# installing bioconductor packages
# (try http:// if https:// URLs are not supported)
source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("Biostrings")
# biocLite("rtracklayer")
# biocLite("GenomicRanges")
# biocLite("GenomicAlignments")
# biocLite("GenomeInfoDbData")
# biocLite("DelayedArray")


#### Some useful options ####

Sys.setenv("R_HISTSIZE"=10000)
options(max.print=100)


#### Load packages ####

library(ff)
library(ffbase)
library(knitr)
library(dichromat)
library(parallel)
library(zoo)
library(viridis)
library(cowplot)
library(tictoc)
library(profvis)
library(dplyr)

# Marks Herons packages
library(maRs)
library(pRon)

# bioconductor packages
library(Biostrings)
library(rtracklayer)
library(GenomicRanges)
library(GenomicAlignments)


#### Plot setup ####

theme_set(theme_cowplot(font_size = 10, line_size = 1) + background_grid(major = "xy", minor = "none"))
update_geom_defaults("bar", list(fill="grey", colour=I("black"), size=0.25))
cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot_colours <-  list(scale_colour_manual(values = cbPalette), scale_fill_manual(values = cbPalette))

