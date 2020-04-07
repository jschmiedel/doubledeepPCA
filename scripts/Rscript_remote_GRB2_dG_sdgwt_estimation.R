# remotely execute function_dG_estimation.R

require(optparse)

# read in varlist file deposited in the copied folder structure
option_list <- list(
  make_option(
    opt_str = c("-d", "--dataset_name"),
    default = "GRB2_dG_dataset",
    type = "character",
    help = "dataset name"
  ),
  make_option(
    opt_str = c("-m", "--method"),
    type = "integer",
    help = "which optimization method to use"
  ),
  make_option(
    opt_str = c("-b", "--predict_binding0"),
    default = 0,
    type = "integer"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## set working directory
setwd("/nfs/users/blehner/jschmiedel/doubledeepPCA/")

## source functions
filelist <- list.files("functions/")
sapply(paste0("functions/", filelist), source, .GlobalEnv)

function_GRB2_dG_sdgwt_estimation(
  dataset_name = opt$dataset_name,
  method = opt$method,
  predict_binding0 = as.logical(opt$predict_binding0)
)