run_core <- function(
    dat,
    group = GROUP,
    n = NULL,
    depth = RARE_DEPTH,
    nrep = N_RAREFACTIONS,
    max_abs = MAX_ABSENCES,
    threshold = NULL,
    seed = RANDOM_SEED,
    min_count = MIN_COUNT,
    min_rel = MIN_REL_ABUND
){
    stratified <- stratified_sampling(dat,group = group,n = n,seed = seed)
    rare <- repeat_rarefy(stratified$otu,depth = depth,nrep = nrep,seed = seed)
    occ <- estimate_occupancy(rare,meta = stratified$meta,group = group, min_count = min_count, min_rel = min_rel)
    core <- classify_core(occ, max_abs = max_abs, threshold = threshold)
    return(list(
        stratified = stratified,
        rare = rare,
        occ = occ,
        core = core$core,
        groups = core$groups
))
}

count_core <- function(core, groups){
    tab <- table(core$Classification)
    expected <- c("Global Core",paste0(groups, " Specific"))
    out <- setNames(rep(0L, length(expected)), sub(" Core$", "", sub(" Specific$", "", expected)))

    for(i in seq_along(expected)){
        if(expected[i] %in% names(tab)){
            out[i] <- tab[[expected[i]]]
        }
    }
    as.data.frame(as.list(out), check.names = FALSE)
}
# ======================
# Check pipeline
# ======================
#
# Validate that the pipeline completed successfully and that
# intermediate objects are consistent before interpreting the
# results or running sensitivity analyses.
#
# The following objects must already exist:
#   dat         : output of subset_data()
#   stratified  : output of stratified_sampling()
#   rare        : output of repeat_rarefy()
#   occ         : output of estimate_occupancy()
#   core        : output of classify_core()
#
# Usage:
#
# check_pipeline(
#     dat,
#     stratified,
#     rare,
#     occ,
#     core
# )
#
# The function reports:
#   - sample distribution before and after stratification
#   - site representation
#   - sequencing depth summary
#   - OTU/metadata synchronization
#   - rarefaction dimensions
#   - occupancy groups
#   - core-community summary
#
# It stops with an error if OTU and metadata are not synchronized.

check_pipeline <- function(
    dat,
    stratified,
    rare,
    occ,
    classified
){
    group <- occ$group
    cat("\n")
    cat("========================================\n")
    cat(" Pipeline diagnostics\n")
    cat("========================================\n")

    ## ----------------------------
    ## Original dataset
    ## ----------------------------
    cat("\nOriginal dataset\n")
    print(table(dat$meta[[group]]))

    cat("\nOriginal sites\n")
    print(table(dat$meta[[group]], dat$meta$site))

    cat("\nSequencing depth\n")
    print(summary(colSums(dat$otu)))

    ## ----------------------------
    ## Stratified subset
    ## ----------------------------
    cat("\nStratified dataset\n")
    print(table(stratified$meta[[group]]))

    cat("\nSelected sites\n")
    print(table(stratified$meta[[group]],
                stratified$meta$site))

    ## ----------------------------
    ## OTU / metadata synchronization
    ## ----------------------------
    stopifnot(
        ncol(stratified$otu) == nrow(stratified$meta),
        identical(
            colnames(stratified$otu),
            stratified$meta$Samplename
        ),
        !anyDuplicated(stratified$meta$Samplename)
    )
    cat("\nDetection summary\n")
    print(table(classified$nDetected))

    cat("\n✓ OTU table and metadata synchronized\n")

    ## ----------------------------
    ## Rarefaction
    ## ----------------------------
    cat("\nRarefaction\n")
    cat("Replicates :", length(rare), "\n")
    cat("Dimensions :", dim(rare[[1]]), "\n")

    ## ----------------------------
    ## Occupancy
    ## ----------------------------
    cat("\nOccupancy\n")
    cat("Groups :", paste(colnames(occ$mean), collapse = ", "), "\n")
    cat("Group sizes :", paste(occ$group_sizes, collapse = ", "), "\n")

    ## ----------------------------
    ## Core classification
    ## ----------------------------
    cat("\nCore classification\n")
    print(table(classified$Classification))

    cat("\nPattern summary\n")
    print(table(classified$Pattern))

    cat("\n========================================\n")
    cat(" Pipeline OK\n")
    cat("========================================\n")
}


# ==========================
# for testing Core threshold 
# ==========================
test_core_threshold <- function (
    occ,
    thresholds = c(0.75, 0.80, 0.90, 0.95, 0.99))
    {
    results <- data.frame(Threshold = integer(),Global = integer(),FHNY = integer(),LMDC = integer())
    for(th in thresholds){
        x <- classify_core(occ,threshold = th)
        results <- rbind(results,cbind(Threshold = th ,count_core(x$core, x$groups)))
    }
    cat("Running test_core_threshold...\n")
    return(results)
}

# ==========================
# for testing number of rarefaction 
# ==========================

test_rarefactions <- function(
    dat,
    group,
    n_rare = c(1,2,5,10) #,20,50,100)
){
    cat("Running test_rarefactions...\n")
    ## Empty results table
    results <- data.frame(Rarefactions = integer(),Global = integer(),FHNY = integer(),LMDC = integer())
    ## Loop over the number of rarefactions
    for(nr in n_rare){
        x <- run_core(dat, nrep = nr)        
        results <- rbind(results,cbind(Rarefactions = nr,count_core(x$core, x$groups)))
    cat("---------------------------------\n")
    cat("Rarefactions:", nr, "\n")
    }
    return(results)
}


# ==========================
# Test depth of rarefaction
# ==========================

test_depth <- function(
    dat,
    depths = c(100,500,1000,2500,5000,6419), #,10000),
    nrep = nrep,
    max_abs = max_abs
){
    cat("Running test_depth...\n")
    results <- data.frame()
    for(d in depths){
        x <- run_core(dat, depth = d)
        results <- rbind(results,cbind(Depth = d,count_core(x$core, x$groups)))
    }
    return(results)
}

# ======================
# Test missing in sample
# ======================

test_absences <- function(
    occ,
    max_abs = 0:5
){
    cat("Running test_absences...\n")
    results <- data.frame()
    for(a in max_abs){
        x <- classify_core(occ,max_abs = a)
        threshold <- threshold_from_absence(occ$group_sizes, a)
        results <- rbind(results,data.frame(MaxAbs = a,Threshold = round(min(threshold),3),count_core(x$core, x$groups)))
    }

    return(results)
}
# ===========================
# Test minimum samples needed
# ===========================

test_sample_size <- function(
    dat,
    group,
    sample_sizes = NULL
){
    cat("Running test_sample_sizes...\n")
    if (is.null(sample_sizes))
        sample_sizes <- 2:min(table(dat$meta[[group]]))
    results <- data.frame()
    for(n in sample_sizes){
        x <- run_core(dat, n = n)
#         tmp <- count_core(x$core, x$groups)
#         print(names(tmp))
        results <- rbind(results,cbind(Samples = n,count_core(x$core, x$groups)))
    }
    return(results)
}

# ===================================
# Test filtering by minimum abundance
# ===================================

test_min_count <- function(
    stratified, group,
    min_counts = 1:5
){
    cat("Running test_min_count...\n")
    results <- data.frame()
    rare <- repeat_rarefy(stratified$otu)
    for(mc in min_counts){
        occ <- estimate_occupancy(rare,meta = stratified$meta,group = group ,min_count = mc,min_rel = 0)
        x <- classify_core(occ, max_abs = MAX_ABSENCES)
        results <- rbind(results,cbind(MinCount = mc,count_core(x$core, x$groups)))
    }
    return(results)
}

# ============================================
# Test filtering by minimum relative abundance
# =============================================

test_relative_abundance <- function(
    stratified, group,
    min_rel = c(0, 0.0001, 0.0005, 0.001, 0.002)
){
    cat("Running test_relative_abundance...\n")
    results <- data.frame()
    rare <- repeat_rarefy(stratified$otu)
    for(r in min_rel){
        occ <- estimate_occupancy(rare,meta = stratified$meta,group = group ,min_count = 1,min_rel = r)
        x <- classify_core(occ, max_abs = MAX_ABSENCES)
        results <- rbind(results,cbind(MinRel = r,count_core(x$core, x$groups)))
    }
    return(results)
}
