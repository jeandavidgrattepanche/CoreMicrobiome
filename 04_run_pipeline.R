#==========
# LIBRARY
#=========
library(vegan)
library(dplyr)


source("00_parameters.R")
source("01_core_functions.R")
source("02_sensitivity_tests.R")

pipeline_summary()

# 01 Import
otu <- read.delim(OTU_FILE, row.names = 1, check.names = FALSE)
meta <- read.delim(META_FILE)

dat <- subset_data(
    otu,
    meta,
    filters = FILTERS,
    min_reads = MIN_READS
)

x <- run_core(dat, group = GROUP)

stratified <- x$stratified
rare <- x$rare
occ <- x$occ
classified <- x$core

sumcore <- summary_core(classified)
print(sumcore$classification)
print(sumcore$pattern)

cat("\nAnalysis completed successfully.\n")
