# remotely execute function_dG_estimation.R

require(optparse)

# read in varlist file deposited in the copied folder structure
option_list <- list(
  make_option(
    opt_str = c("-d", "--dataset_name"),
    type = "character",
    help = "dataset name"
  ),
  make_option(
    opt_str = c("-m", "--method"),
    type = "integer",
    help = "which optimization method to use"
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

function_dG_estimation(
  dataset_name = opt$dataset_name,
  dataset_suffix = opt$dataset_suffix,
  method = opt$method
)