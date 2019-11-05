#### Some useful options ####

Sys.setenv("R_HISTSIZE"=10000)
options(max.print=100)


#### Install packages

# install.packages(c("ff", "ffbase", "knitr", "dichromat", "parallel", "zoo", "viridis", "cowplot", "tictoc", "stringi", "profvis"), dependencies = TRUE)  # these are probably not all needed common packages

# installing Mark Heron's packages
# install.packages("devtools")
# devtools::install_bitbucket(repo="markheron/maRs", keep_source=TRUE)
# devtools::install_bitbucket(repo="markheron/pRon", keep_source=TRUE)
# devtools::install_bitbucket(repo="markheron/nucular", keep_source=TRUE)

# installing bioconductor packages
# (try http:// if https:// URLs are not supported)
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("Biostrings")
# biocLite("rtracklayer")
# biocLite("GenomicRanges")
# biocLite("GenomicAlignments")
# biocLite("GenomeInfoDbData")
# biocLite("DelayedArray")


#### Load packages ####

library(ff)
library(ffbase)
library(knitr)
library(parallel)
library(plyr)
library(profvis)
library(pryr)
library(dplyr)
library(stringr)
library(reshape2)
library(zoo)
library(tictoc)
library(dichromat)

# Mark Herons packages
library(maRs)
library(pRon)

# bioconductor packages
library(Biostrings)
library(rtracklayer)
library(GenomicRanges)
library(GenomicAlignments)

# plotting
library(cowplot)
library(viridis)
library(scales)

percent_no_sign <- function(x) {  # to use with scale_y_continuous(labels = percent_no_sign) instead of scale_y_continuous(labels = percent)
  if(length(x) == 0) return(character())
  x * 100
}


#### Set plot defaults and colors ####

theme_set(theme_cowplot(font_size = 10, line_size = 1) + background_grid(major = "xy", minor = "none"))
update_geom_defaults("bar", list(fill="grey", colour=I("black"), size=0.25))

#color_palette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")  # my old colors
color_palette <- c("#e31a1c", "#253494", "#fecc5c", "#56B4E9", "#888888", "#fd8d3c", "#000000")

ggplot_colours <-  list(scale_colour_manual(values = color_palette), scale_fill_manual(values = color_palette))

diag_line <- geom_abline(slope = 1, intercept = 0, linetype = "dashed",  color = "black", size = 0.3, alpha = 0.3)
my_gap <- theme(legend.box.margin=margin(0, 8, 0, 0, "pt"))

if(Sys.info()['sysname'] == "Linux") {  # for some reason this need different values to work out the same way
  axis_title_x_top_margin <- 2.5
  strip_text_x <- 2.4 
} else {
  axis_title_x_top_margin <- 0
  strip_text_x <- 1
}

theme_paper <- list(theme_cowplot(font_size = 6, line_size = 0.5) + background_grid(major = "xy", minor = "none", size.major = 0.2, size.minor = 0.4),
                    guides(shape=guide_legend(keywidth=0.3, keyheight=0.3, default.unit="cm"),
                           color=guide_legend(keywidth=0.3, keyheight=0.3, default.unit="cm"),
                           fill=guide_legend(keywidth=0.3, keyheight=0.3, default.unit="cm")
                           ),
                    theme(legend.margin = margin(0, 0, 0, 0, "pt"), 
                          legend.box.margin = margin(0, 0, 0, 0, "pt"),
                          legend.box.spacing = margin(0, 3, 0, 3, "pt"),
                          strip.text.x = element_text(margin = margin(strip_text_x, 0, strip_text_x, 0, "pt")),  # facet_wrap headers
                          strip.text.y = element_text(margin = margin(1, 0, 1, 0, "pt")),  # facet_wrap headers
                          axis.title.y = element_text(margin = margin(0, 2.5, 0, 0.5, "pt")),
                          axis.title.y.right = element_text(margin = margin(0, 0, 0, 1.5, "pt")),
                          axis.title.x = element_text(margin = margin(axis_title_x_top_margin, 0, 0, 0, "pt")),
                          axis.title.x.top = element_text(margin = margin(0, 0, 1.5, 0, "pt")),
                          plot.margin = margin(2, 2, 1, 0, "pt"),
                          strip.switch.pad.grid = unit(1.5, "pt"),
                          strip.switch.pad.wrap = unit(1.5, "pt"),
                          plot.subtitle = element_text(size = 5.5)
                          ))
rm(strip_text_x)

source("../../shared_functions.R")


