## source functions
filelist <- list.files("functions/")
sapply(paste0("functions/", filelist), source, .GlobalEnv)

# create directory structure
dir.create("processed_data/", showWarnings = FALSE)
dir.create("processed_data/dG/", showWarnings = FALSE)
dir.create("processed_data/tmp/", showWarnings = FALSE)
dir.create("results/", showWarnings = FALSE)
dir.create("results/dG/", showWarnings = FALSE)

##################
###### GRB2 ######
##################

### prepare dataset
function_GRB2_dG_prepare_dataset(
  dataset_name = "GRB2",
  DMS_file_list = c(
    "processed_data/GRB2_wildtype_alldata.txt",
    "processed_data/GRB2_singles_alldata.txt",
    "processed_data/GRB2_doubles_alldata.txt"
  ),
  PDB_interaction_file = "dataset/GRB2/PDB_contactmap_2vwf_AB.txt"
)

#####################
### dG estimation
### method 1: using data from both assays only for single variants
### method 2: using data from both assays for all variants
### method 3: using only data from binding assay for all variants
### method 4: using only doubles from both assays
# this is run on the cluster via bash_remote_dG_estimation.sh
function_dG_estimation(
  dataset_name = "GRB2_dG_dataset", # or with sufix '_allvars'
  dataset_suffix = "", # used to run iterations
  method = 1, # any of 1 to 4
)

#################
###### PDZ ######
#################

### prepare dataset
function_PDZ_dG_prepare_dataset(
  dataset_name = "PDZ",
  DMS_file_list = c(
    "processed_data/PDZ_wildtype_alldata.txt",
    "processed_data/PDZ_singles_alldata.txt",
    "dataset/PDZ/McLaughlin2012_DLG4_singlemutants.txt"
  )
)

#####################
### dG estimation
### using data single mutants from both ddPCA assays and McLaughlin2012
# this is run on the cluster via bash_remote_dG_estimation.sh
function_PDZ_dG_estimation(
  dataset_name = "PDZ_dG_dataset", # or with sufix '_allvars'
  dataset_suffix = "", # used to run iterations
)