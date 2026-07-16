CoreMicrobiome is an R framework for identifying microbial core communities from amplicon sequencing data. It implements stratified sampling, repeated rarefaction, occupancy estimation, flexible core classification (occupancy threshold or allowed absences), and a suite of sensitivity analyses to evaluate the robustness of core-community definitions across sequencing depth, sample size, rarefaction, and detection thresholds.

#----------------------------------------------------------
# Function: estimate_occupancy
#
# Description:
#   Estimate OTU occupancy across groups using repeated
#   rarefaction.
#
# Arguments:
#   rare : list of rarefied OTU tables
#   meta      : metadata table
#   group     : grouping variable
#
# Returns:
#   List containing:
#     - mean
#     - sd
#     - all
#     - group_sizes
#----------------------------------------------------------
