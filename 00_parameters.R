# ===========================
# Input files
# ===========================

OTU_FILE  <- "/OTUtable.txt"
META_FILE <- "/metadata.txt"

# ===========================
# Dataset selection
# ===========================
FILTERS <- list(
    city     = c("NY","DC","PH")
    dataset  = c("Environmental","Experimental"),
    molecule = c("16Sc","16SD","18Sc"),
)

GROUP      <- "dataset"

# ===========================
# Pipeline parameters
# ===========================

MIN_READS        <- 5000
RARE_DEPTH       <- 5000
N_RAREFACTIONS   <- 20
MAX_ABSENCES     <- 3
RANDOM_SEED      <- 1234
MIN_COUNT        <- 1
MIN_REL_ABUND    <- 0.0005