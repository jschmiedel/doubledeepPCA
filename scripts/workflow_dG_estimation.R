
## source functions
filelist = list.files('functions/')
sapply(paste0('functions/',filelist),source,.GlobalEnv)

dir.create("processed_data/", showWarnings = FALSE)
dir.create("results/", showWarnings = FALSE)
dir.create("results/dG/", showWarnings = FALSE)

##################
###### GRB2 ######
##################

### prepare dataset
#epPCR datasets
function_dGestimation_prepare_dataset(
  name = "GRB2_epPCR",
  DMS_file_list = c("dataset/DiMSum_GRB2/GRB2_epPCR_GPD_fitness_singles.txt",
                    "dataset/DiMSum_GRB2/GRB2_epPCR_GPD_fitness_doubles.txt",
                    "dataset/DiMSum_GRB2/GRB2_epPCR_CYC_fitness_singles.txt",
                    "dataset/DiMSum_GRB2/GRB2_epPCR_CYC_fitness_doubles.txt"),
  PDB_interaction_file = "dataset/PDB_contactmap_2vwf_AB.txt",
  RSA_file = "dataset/2vwf_A.rsa"
)

#####################
### dG estimation ###
#####################

### method 1: using only single data from both assays
function_dGestimation_method1_singles_bothassays(
  name = "GRB2_epPCR",
  dataset_file = "processed_data/GRB2_epPCR_dG_dataset.txt",
  Ncores = 15,
  Nbootstraps = 300
)
#minimal error if scaling factors are 0
#whats wrong here?

#with background growth set
function_dGestimation_method1_singles_bothassays(
  name = "GRB2_epPCR_bgrset",
  dataset_file = "processed_data/GRB2_epPCR_dG_dataset.txt",
  Ncores = 15,
  Nbootstraps = 300,
  bgr_set = rep(0.55,2)
)
