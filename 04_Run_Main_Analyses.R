# Main analyses

input_file <- "data/dat_classed.rda"
if (!file.exists(input_file)) stop("Missing input: ", input_file)

run_level <- function(level) {
  SPATIAL_LEVEL <- level
  PREY_FAMILY <- "_1"
  RUN_TAG <- "_T"

  N_MIN <- 5
  N_STRICT <- 25
  MIN_SETS <- 3
  R_PERM <- 999
  ALPHA <- 0.05

  BS_TEST <- "welch"
  COMP_RELATIVE <- TRUE
  OCC_BINARY <- FALSE
  ADD_NICHE_PART <- TRUE
  ADD_IND_NICHE <- FALSE

  source("03_Statistical_Pipeline.R", local = TRUE)
  run_all_scenarios(data_path = input_file, tag = RUN_TAG)
}

# Spatial levels
run_level("all_gulf")
run_level("ecoregion")
run_level("stratum")
