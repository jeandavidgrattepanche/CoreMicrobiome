############################################################
## CoreMicrobiome Validation
##
## Run this script after completing the pipeline to:
##   1. verify that all intermediate objects are consistent,
##   2. assess the robustness of the core-community definition
##      to key analytical parameters.
############################################################


############################################################
## 1. Pipeline diagnostics
############################################################

check_pipeline(
    dat,
    stratified,
    rare,
    occ,
    classified
)


############################################################
## 2. Sensitivity analyses
############################################################

## Number of rarefactions
rarefaction_test <- test_rarefactions(dat)

## Rarefaction depth
depth_test <- test_depth(dat)

## occupancy thresholds
occ_threshold_test <- test_core_threshold(occ)

## Allowed absences
absence_test <- test_absences(occ)

## Sample size
samplesize_test <- test_sample_size(dat, occ$group)

## Minimum read count
count_test <- test_min_count(stratified, occ$group)

## Minimum relative abundance
relab_test <- test_relative_abundance(stratified, occ$group)


cat("\n")
cat("========================================\n")
cat("Validation completed successfully.\n")
cat("Sensitivity results stored in:\n")
cat("  rarefaction_test\n")
cat("  depth_test\n")
cat("  occ_threshold_test\n")
cat("  absence_test\n")
cat("  samplesize_test\n")
cat("  count_test\n")
cat("  relab_test\n")
cat("========================================\n")