# remotely execute function_dG_estimation.R

require(optparse)

# read in varlist file deposited in the copied folder structure
option_list <- list(
  make_option(
    opt_str = c("-d", "--dataset_name"),
    default = "GB1_dG_dataset",
    type = "character",
    help = "dataset name"
  ),
  make_option(
    opt_str = c("-s", "--dataset_size"),
    type = "double",
    default = 1,
    help = "what fraction of doubles to use"
  ),
  make_option(
    opt_str = c("-b", "--predict_binding0"),
    default = 0,
    type = "integer"
  ),
  make_option(
    opt_str = c("-a", "--approach"),
    default = 1,
    type = "integer"
  ),
  make_option(
    opt_str = c("-i", "--iteration"),
    default = 0,
    type = "integer"
  ),
  make_option(
    opt_str = c("-t", "--task"),
    default = 1,
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

function_GB1_dG_estimation(
  dataset_name = opt$dataset_name,
  dataset_size = opt$dataset_size,
  iteration = opt$iteration,
  task = opt$task,
  approach = opt$approach,
  predict_binding0 = as.logical(opt$predict_binding0)
)