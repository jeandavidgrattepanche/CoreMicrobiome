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