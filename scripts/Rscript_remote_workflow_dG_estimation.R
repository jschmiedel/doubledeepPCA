# remotely execute workflow_dG_estimation.R

require(optparse)

# read in varlist file deposited in the copied folder structure
option_list <- list(
  make_option(
    opt_str = c("-f", "--first_stage"),
    type = "integer",
    default = 1,
    help = "first stage to be performed"
  ),
  make_option(
    opt_str = c("-l", "--last_stage"),
    type = "integer",
    default = 0,
    help = "final stage to be performed"
  ),
  make_option(
    opt_str = c("-d", "--dataset_name"),
    type = "character",
    help = "dataset name"
  ),
  make_option(
    opt_str = c("-c", "--n_cores"),
    type = "integer",
    default = 15,
    help = "# cores to be used"
  ),
  make_option(
    opt_str = c("-b", "--n_bootstraps"),
    type = "integer",
    default = 15,
    help = "# bootstraps per model fit"
  ),
  make_option(
    opt_str = c("-s", "--dataset_suffix"),
    type = "character",
    default = ""
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## set working directory
setwd("/nfs/users/blehner/jschmiedel/doubledeepPCA/")

## source functions
filelist <- list.files("functions/")
sapply(paste0("functions/", filelist), source, .GlobalEnv)

source("scripts/workflow_dG_estimation.R")

workflow_dG_estimation(
  first_stage = opt$first_stage,
  last_stage = opt$last_stage,
  dataset_name = opt$dataset_name,
  dataset_suffix = opt$dataset_suffix,
  n_cores = opt$n_cores,
  n_bootstraps = opt$n_bootstraps
)