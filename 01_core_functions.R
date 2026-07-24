
# "We estimated occupancy using repeated rarefaction, quantified its uncertainty, and classified core communities while allowing sensitivity analyses for sequencing depth, sample size, and occupancy threshold."
# 
# SAB_core_pipeline_v3/01_core_functions.R
# │
# ├── 01 Data handling
# |      subset_data()
# |      stratified_sampling()
# |          allocate_sites()
# |          sample_sites()
# │
# ├── 02 Rarefaction
# |     repeat_rarefy()
# ├── 03 Occupancy
# |     occupancy()
# |     estimate_occupancy()
# │
# ├── 04 Classification
# |     threshold_from_absence()
# |     classify_core()
# |     summary_core()
# |
# ├── 05 Utilities
# |     pipeline_summary()

# 01 Data handling
subset_data <- function(otu, meta, filters = NULL, min_reads = 0) {
#   browser()
  ## Filter metadata using dplyr syntax
    meta_sub <- meta
    if (!is.null(filters)) {
        for (v in names(filters)) {
            meta_sub <- meta_sub[meta_sub[[v]] %in% filters[[v]],,drop = FALSE]
        }}  
  sample_depth <- colSums(otu)[meta_sub$Samplename]
  good_samples <- names(sample_depth)[sample_depth >= min_reads]

  ## Keep only samples present in the OTU table
  samples <- intersect(meta_sub$Samplename, good_samples)
  if(length(samples) == 0)
     stop("No matching samples found.")

  ## Subset OTU table
  otu_sub <- otu[, samples, drop = FALSE]
  
  ## Remove OTUs absent from the selected samples
  keep_otus <- rowSums(otu_sub) > 0
  otu_sub <- otu_sub[keep_otus, , drop = FALSE]
  
  ## Reorder metadata to match OTU table
  meta_sub <- meta_sub[match(samples, meta_sub$Samplename), ]
  missing <- setdiff(meta_sub$Samplename, colnames(otu))

  if(length(missing) > 0){

      warning(length(missing),
              " samples found in metadata but not OTU table.")
  
  }
  cat("---------------------------------\n")
  cat("Subset summary\n")
  cat("---------------------------------\n")
  cat("Samples:", ncol(otu_sub), "\n")
  cat("OTUs before filtering   :", nrow(otu), "\n")
  cat("OTUs after filtering   :", nrow(otu_sub), "\n")
  cat("OTUs removed (absent):", sum(!keep_otus), "\n")
  cat("Minimum reads    :", min_reads, "\n")
  cat("Samples removed by read depth:", sum(sample_depth < min_reads), "\n")
  cat("---------------------------------\n")

  ## Return both objects
  return(list(
    otu = otu_sub,
    meta = meta_sub
 #   parameters = list(filter = deparse(substitute(...)))
 #    history = "subset_data" 
 ))
}

# 1b: stratified_sampling
## To account for unequal sampling effort among cities while preserving spatial representation, samples were selected using stratified random sampling in which each site was represented before additional samples were allocated randomly among sites with remaining observations.

stratified_sampling <- function(dat, group, n = NULL, stratify = "site", seed = RANDOM_SEED) {

  ## target sample size
  if(is.null(n))
       n <- min(table(dat$meta[[group]]))
  if(n <= 0)
       stop("'n' must be > 0.")  
  set.seed(seed)
  selected_samples <- character()

  ## select in group
  groups <- unique(dat$meta[[group]])
    for(g in groups){
        meta_g <- dat$meta[dat$meta[[group]] == g,]
        allocation <- allocate_sites(meta_g, n, stratify)
        selected <- sample_sites(meta_g, allocation, stratify = stratify) ##, id = "Samplename")
        selected_samples <- c(selected_samples,selected)
    }	
  otu_new <- dat$otu[, selected_samples]
  meta_new <- dat$meta[match(selected_samples,dat$meta$Samplename),]
#   sampling_summary <- character()
# 
#   cat("---------------------------------\n")
#   cat("Subset summary\n")
#   cat("---------------------------------\n")
#   cat("Samples:", ncol(otu_new), "\n")
#   cat("OTUs   :", nrow(otu_new), "\n")
#   cat("---------------------------------\n")
# 
  ## Return both objects
  return(list(
    otu = otu_new,
    meta = meta_new
#     sampling = sampling_summary 
 ))
}

allocate_sites <- function(
    meta_group,
    n,
    stratify = "site"
){
  if (!is.numeric(n) || length(n) != 1 || is.na(n))
       stop("'n' must be a single numeric value.")
  if (n <= 0)
       stop("'n' must be greater than 0.")
  if (n != as.integer(n))
       stop("'n' must be an integer.")
  if (n > nrow(meta_group))
    stop("Requested sample size (n=", n,
         ") exceeds available samples (",
         nrow(meta_group), ").")
  site_sizes <- table(meta_group[[stratify]])
  if(length(site_sizes) == 0)
    stop("No strata found.")
  allocation <- setNames(
    rep(0, length(site_sizes)),
    names(site_sizes)
)
  for(site in names(site_sizes)){
    if(sum(allocation) < n)
	allocation[site] <- allocation[site] + 1
}

  while(sum(allocation) < n){
    available <- names(
      site_sizes[allocation < site_sizes])
    chosen <- sample(available,1)
    allocation[chosen] <-
        allocation[chosen] + 1
}
  return(allocation)  
}

sample_sites <- function(meta_group, allocation, stratify = "site", id = "Samplename", seed = RANDOM_SEED){
  set.seed(seed)
  selected <- character()
  for(site in names(allocation)){
    ids <- meta_group[[id]][meta_group[[stratify]] == site]
    if(allocation[site] >0){
      pick <- sample(ids,allocation[site])
      selected <- c(selected, pick)}
   
}
  return(selected)
}

# 02 Rarefaction
repeat_rarefy <- function(otu,
                          depth = RARE_DEPTH,
                          nrep = N_RAREFACTIONS,
                          seed = RANDOM_SEED){
  set.seed(seed)
  out <- vector("list", nrep)
  for(i in seq_len(nrep)){
    out[[i]] <- t(
      rrarefy(
        t(otu),
        sample = depth
      )
    )
  }
  names(out) <- paste0("Rep", seq_len(nrep))
  return(out)
}

# 03 Occupancy
occupancy <- function(otu, meta, group,  min_count = MIN_COUNT, min_rel = MIN_REL_ABUND){
#     cat("run occupancy\n")
    TotalReads <- rowSums(otu)
    ## make sure sample order is identical
    meta <- meta[match(colnames(otu), meta$Samplename), ]
    ## relative abundance
    rel <- sweep(otu, 2, colSums(otu), "/")
    ## define presence
    present <- (otu >= min_count) &
               (rel >= min_rel)
    groups <- sort(unique(meta[[group]]))
    occ <- lapply(groups, function(g){
        samples <- meta$Samplename[meta[[group]] == g]
#         rowMeans(otu[, samples, drop = FALSE] > 0)
        rowMeans(present[, samples, drop = FALSE])
    })
    occ <- do.call(cbind, occ)
    colnames(occ) <- groups
    list(occ = occ,
    TotalReads = TotalReads
    )
}


estimate_occupancy <- function(rare.list, meta, group, min_count = MIN_COUNT, min_rel = MIN_REL_ABUND){
    occ.list <- lapply(rare.list,occupancy,meta = meta,group = group, min_count = min_count,min_rel = min_rel)
    occ_mats <- lapply(occ.list, `[[`, "occ")
    mean_occ <- Reduce("+", occ_mats) / length(occ_mats)
    sd_occ <- sqrt(Reduce("+",lapply(occ_mats,function(x) (x-mean_occ)^2)) / (length(occ_mats)-1))
    read_list <- lapply(occ.list, `[[`, "TotalReads")
    mean_reads <- Reduce("+", read_list) / length(read_list)
    cat("estimate occupancy\n")
    list(
      mean = mean_occ,
      sd = sd_occ,
      all = occ_mats,
      TotalReads = mean_reads,
      groups = colnames(mean_occ),
      group  = group,
      group_sizes = table(meta[[group]]))		

}

# 04 Classification
threshold_from_absence <- function(n, max_absences = 2){
    (n - max_absences) / n
}

max_abs_from_threshold <- function(n, threshold){
    ceiling(n * (1 - threshold))
}

classify_core <- function(occ, max_abs = NULL, threshold = NULL){
    cat("run classifier\n")
    occ_mean <- occ$mean
    group_sizes <- occ$group_sizes
    detected <- occ_mean > 0
    ## One (and only one) definition must be provided
    if (is.null(max_abs) && is.null(threshold))
        stop("Provide either 'max_abs' or 'threshold'.")

    if (!is.null(max_abs) && !is.null(threshold))
        stop("Provide only one of 'max_abs' or 'threshold'.")

    ## Convert max_abs to thresholds
    if (!is.null(max_abs))
        threshold <- threshold_from_absence(group_sizes, max_abs)
    if(length(threshold) == 1)
        threshold <- rep(threshold, ncol(occ_mean))
    names(threshold) <- colnames(occ_mean)
#     cat("Group        :",  names(group_sizes),"\n")
#     cat("# samples    :", group_sizes,"\n")
#     cat("Threshold    :", threshold,"\n")
#     core <- occ >= threshold
#     print(colnames(occ))
#     print(names(threshold))
    core <- matrix(FALSE, nrow = nrow(occ_mean), ncol = ncol(occ_mean), dimnames = dimnames(occ_mean))

    for(g in colnames(occ_mean)){
        core[, g] <- occ_mean[, g] >= threshold[g]
    }
#     cat("\nCore OTUs by group\n")
#     print(colSums(core))
#     cat("\nCore combinations\n")
#     print(table(core[,1], core[,2]))
    detected_pattern <- apply(detected * 1L, 1, paste0, collapse = "")
    core_pattern     <- apply(core * 1L, 1, paste0, collapse = "")
    nDetected <- rowSums(detected)
    nCore <- rowSums(core)
    CoreType <- character(nrow(occ_mean))
    for(i in seq_len(nrow(occ_mean))){
        groups_core <- names(group_sizes)[core[i,]]
        if(nDetected[i] == 0){
           CoreType[i]<- "Absent"
        } else if (nCore[i] == 0){
            CoreType[i] <- "Not Core"
        } else if(nCore[i] == ncol(occ_mean)){
            CoreType[i] <- "Global Core"
        } else if(nCore[i] == 1){
            CoreType[i] <- paste(groups_core, "Specific")
        } else{
            CoreType[i] <- paste(
                paste(groups_core, collapse="-"),
                "Shared"
            )
        }
    }

    out <- data.frame(
        OTU = rownames(occ_mean),
        Pattern =  core_pattern,
        detectedPattern = detected_pattern,
        nDetected = nDetected,
        nCore = nCore,
        CoreType = CoreType,
        TotalReads = occ$TotalReads,
       # occ,
        row.names = NULL
    )
    return(list(core = out, groups = occ$groups))
}

list_core <- function(core, coretype = NULL){
    x <- core
    if(!is.null(coretype))
        x <- subset(x, CoreType %in% coretype)
    x
}

split_core <- function(core){
    split(core, core$CoreType)
}

core_table <- function(classified,
                       occ = NULL,
                       otu = NULL,
                       file = "CoreOTUs.tsv",
                       remove = c("Absent", "Not Core")) { #, "Not Core"
    x <- classified
    x <- subset(x, !(CoreType %in% remove))
    if (!is.null(occ))
        x <- cbind(x, occ$mean[match(x$OTU, rownames(occ$mean)), ])
    if (!is.null(otu))
        x$TotalReads <- rowSums(otu[x$OTU, , drop = FALSE])
 ## Order for convenience
    x <- x[order(x$CoreType, -x$TotalReads), ]
    write.table(x,file = file,sep = "\t",quote = FALSE,row.names = FALSE)
    invisible(x)
}

summary_core <- function(core){
    list(
        coretype =
            dplyr::count(
                core,
                CoreType,
                sort = TRUE
            ),
        pattern =
            dplyr::count(
                core,
                Pattern,
                sort = TRUE
            )
    )
}

# 05 Utilities
pipeline_summary <- function(){

    cat("\n")
    cat("========================================\n")
    cat(" Core Community Pipeline Parameters\n")
    cat("========================================\n")
    cat("Minimum reads               :", MIN_READS, "\n")
    cat("Minimum count               :", MIN_COUNT, "\n")
    cat("Minimum relative abundance  :", MIN_REL_ABUND, "\n")
    cat("Rarefaction depth           :", RARE_DEPTH, "\n")
    cat("Rarefactions                :", N_RAREFACTIONS, "\n")
    cat("missing in max samples      :", MAX_ABSENCES, "\n")
    cat("Random seed                 :", RANDOM_SEED, "\n")
    cat("group investigated          :", GROUP, "\n")
    cat("========================================\n\n")

}

check_alignment <- function(otu, meta){

    stopifnot(
        identical(
            colnames(otu),
            meta$Samplename
        )
    )
}