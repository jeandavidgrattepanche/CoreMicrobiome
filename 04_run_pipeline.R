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
otu <- read.delim("OTUtable.txt", row.names = 1, check.names = FALSE) 
# your OTU table
#OTU sample1 sample2 ... 
#OTU1 number number ...
#OTU2 number number ...
#OTU3 number number ...
# ...
meta <- read.delim("metadata.txt")
# metadata
#Samplename	city	site	molecule	time	Sampleid	dataset
# sample1	cityA	site1	16Sc	t1	A_1_t1	exp2
# sample2	cityB	site1	16Sc	t4	B_1_t4	exp2
# sample3	cityB	site2	16SD	t10	B_2_t10	exp4 
dat <- subset_data(
    otu,
    meta,
    city %in% c("LMDC","FHNY"),
    dataset == "Environmental",
    molecule == "16Sc",
    min_reads = MIN_READS
)


# stratified <- stratified_sampling(
#     dat,
#     group = "city",
#     n = NULL,
# ##    stratify = c("site","time"), ### for LATER
#     stratify = "site",
#     seed = RANDOM_SEED
# )
# 
# 
# rare <- repeat_rarefy(
#     otu = stratified$otu,
#     depth = RARE_DEPTH,
#     nrep = N_RAREFACTIONS
# )
# 
# 
# occ <- estimate_occupancy(rare,
#     meta = stratified$meta,
#     group = "city", 
#     min_count = MIN_COUNT, 
#     min_rel = MIN_REL_ABUND
# )
# 
# 
# core <- classify_core(
#     occ,
#     max_abs = MAX_ABSENCES
# )
x <- run_core(dat)

stratified <- x$stratified
rare <- x$rare
occ <- x$occ
core <- x$core

sumcore <- summary_core(core)
print(sumcore$classification)
print(sumcore$pattern)

cat("\nAnalysis completed successfully.\n")
