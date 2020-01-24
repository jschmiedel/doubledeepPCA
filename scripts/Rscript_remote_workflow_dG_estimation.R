#remotely execute workflow_dG_estimation.R

require(optparse)

#read in varlist file deposited in the copied folder structure
option_list = list(
  make_option(opt_str = c("-f", "--first_stage"), type="integer", default=1,
              help="first stage to be performed", metavar = "character"),
  make_option(opt_str = c("-l", "--last_stage"), type="integer", default=0,
              help="final stage to be performed", metavar = "character"),
  make_option(opt_str = c("-d", "--dataset_name"), type="character", default="GRB2_epPCR",
              help="dataset name", metavar = "character"),
  make_option(opt_str = c("-c", "--Ncores"), type="integer", default=15,
              help="# cores to be used", metavar = "character"),
  make_option(opt_str = c("-b", "--Nbootstraps"), type="integer", default=15,
              help="# bootstraps per model fit", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## set working directory
setwd("/nfs/users/blehner/jschmiedel/doubledeepPCA/")

## source functions
filelist = list.files('functions/')
sapply(paste0('functions/',filelist),source,.GlobalEnv)

source("scripts/workflow_dG_estimation.R")

workflow_dG_estimation(first_stage = opt$first_stage,
  last_stage = opt$last_stage,
  dataset_name = opt$dataset_name,
  Ncores = opt$Ncores,
  Nbootstraps = opt$Nbootstraps)