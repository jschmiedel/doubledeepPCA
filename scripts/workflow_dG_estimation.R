
workflow_dG_estimation <- function(
                                   first_stage = 1,
                                   last_stage = 0,
                                   dataset_name = "GRB2",
                                   dataset_suffix = "",
                                   n_cores = 15,
                                   n_bootstraps = 10) {

  ## source functions
  filelist <- list.files("functions/")
  sapply(paste0("functions/", filelist), source, .GlobalEnv)



  # create directory structure
  dir.create("processed_data/", showWarnings = FALSE)
  dir.create("processed_data/tmp/", showWarnings = FALSE)
  dir.create("results/", showWarnings = FALSE)
  dir.create("results/dG/", showWarnings = FALSE)

  ##################
  ###### GRB2 ######
  ##################

  ### prepare dataset
  stagenum <- 1
  function_dG_prepare_dataset(
    dataset_name = dataset_name,
    DMS_file_list = c(
      "processed_data/GRB2_wildtype_alldata.txt",
      "processed_data/GRB2_singles_alldata.txt",
      "processed_data/GRB2_doubles_alldata.txt"
    ),
    PDB_interaction_file = "dataset/GRB2/PDB_contactmap_2vwf_AB.txt",
    RSA_file = "dataset/GRB2/2vwf_A.rsa",
    execute = (first_stage <= stagenum &
      (last_stage == 0 | last_stage >= stagenum))
  )

  #####################
  ### dG estimation ###
  #####################

  #####################
  ### method 1: using data from both assays only for single variants

  # 1) with all parameters free to fit
  stagenum <- 2
  function_dG_method1_singles_bothassays(
    dataset_name = dataset_name,
    dataset_suffix = dataset_suffix,
    n_cores = n_cores,
    n_bootstraps = n_bootstraps,
    execute = (first_stage <= stagenum &
      (last_stage == 0 | last_stage >= stagenum))
  )
  
  #####################
  ### method 2: using data from both assays for all variants
  stagenum <- 3
  function_dG_method2_allvars_bothassays(
    dataset_name = dataset_name,
    dataset_suffix = dataset_suffix,
    n_cores = n_cores,
    n_bootstraps = n_bootstraps,
    execute = (first_stage <= stagenum &
      (last_stage == 0 | last_stage >= stagenum))
  )


  #####################
  ### method 3: using only data from binding assay for all variants
  stagenum <- 4
  function_dG_method3_allvars_bindingassay(
    dataset_name = dataset_name,
    dataset_suffix = dataset_suffix,
    n_cores = n_cores,
    n_bootstraps = n_bootstraps,
    execute = (first_stage <= stagenum &
      (last_stage == 0 | last_stage >= stagenum))
  )
}
