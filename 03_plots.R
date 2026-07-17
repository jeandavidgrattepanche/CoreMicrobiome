plot_sensitivity <- function(
    results,
    x ,
    xlab,
    title = NULL,
    cols = NULL
){
    groups <- names(results)[!(names(results) %in% x)]
if (is.null(cols)) {
    cols <- c(
        "black",
        "steelblue",
        "firebrick",
        "darkgreen",
        "orange",
        "purple",
        "brown"
    )[seq_along(groups)]
}
    matplot(
        x = results[[x]],
        y= as.matrix(results[, groups]),
        type = "b",
        pch = 19,
        lty = 1,
        cex = 1,
        lwd = 2,
        cex.lab = 1.2,
        cex.axis = 1,
        col = cols,
        xlab = xlab,
        ylab = "Number of core OTUs",
        main = title
    )
    legend(
        "topright",
        legend = groups,
        col = cols,
        lty = 1,
        pch = 19,
        bty = "n"
    )
    invisible(NULL)
}


plot_depth <- function(results, title = NULL){
    plot_sensitivity(results, x= "Depth", xlab="Rarefaction depth", title = title)
}

plot_rarefactions <- function(results, title = NULL){
    plot_sensitivity(results, x= "Rarefactions", xlab="Number of Rarefactions", title = title)
}

plot_absences <- function(results, title = NULL){
    plot_sensitivity(results, x= "MaxAbs", xlab="Maximum allowed missing samples", title = title)
}

plot_thocc <- function(results, title = NULL){
    plot_sensitivity(results, x= "Threshold", xlab="Occupancy threshold", title = title)
}

plot_sample_size <- function(results, title = NULL){
    plot_sensitivity(results, x= "Samples", xlab="Samples per group", title = title)
}

plot_min_count <- function(results, title = NULL){
    plot_sensitivity(results, x= "MinCount", xlab="Minimum Number of read per OTUs", title = title)
}

plot_relative_abundance <- function(results, title = NULL){
    plot_sensitivity(results, x= "MinRel", xlab="Minimum relative abundance", title = title)
}


plot_validation <- function(
    depth_test,
    rarefaction_test,
    absence_test,
    occ_threshold_test,
    samplesize_test,
    count_test,
    relab_test
){

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    par(mfrow = c(3, 3), mar = c(4,4,2,1))

    plot_depth(depth_test)
    plot_rarefactions(rarefaction_test)
    plot_absences(absence_test)
    plot_thocc(occ_threshold_test)
    plot_sample_size(samplesize_test)
    plot_min_count(count_test)
    plot_relative_abundance(relab_test)
}

plot_validation(
    depth_test,
    rarefaction_test,
    absence_test,
    occ_threshold_test,
    samplesize_test,
    count_test,
    relab_test
)


