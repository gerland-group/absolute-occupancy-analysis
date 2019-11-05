#### Some useful options ####

Sys.setenv("R_HISTSIZE"=10000)
options(max.print=100)


#### Load packages ####

library(ff)
library(ffbase)
library(knitr)
library(parallel)
library(plyr)
library(profvis)
library(dplyr)

# installing Mark Heron's packages
#install.packages("devtools")
#devtools::install_bitbucket(repo="markheron/maRs", keep_source=TRUE)
#devtools::install_bitbucket(repo="markheron/pRon", keep_source=TRUE)

# Marks Herons packages
library(maRs)
library(pRon)

# installing bioconductor packages
# (try http:// if https:// URLs are not supported)
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")

# bioconductor packages
library(Biostrings)
library(rtracklayer)
library(GenomicRanges)
library(GenomicAlignments)

