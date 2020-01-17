
workflow_dG_estimation = function(
  first_stage = 1,
  last_stage = 0,
  dataset_name = "GRB2_epPCR",
  Ncores = 15,
  Nbootstraps = 15
  ) 
{

  ## source functions
  filelist = list.files('functions/')
  sapply(paste0('functions/',filelist),source,.GlobalEnv)



  # create directory structure
  dir.create("processed_data/", showWarnings = FALSE)
  dir.create("results/", showWarnings = FALSE)
  dir.create("results/dG/", showWarnings = FALSE)

  ##################
  ###### GRB2 ######
  ##################

  ### prepare dataset
  #epPCR datasets
  stagenum = 1
  function_dG_prepare_dataset(
    name = dataset_name,
    DMS_file_list = c("dataset/DiMSum_GRB2/GRB2_epPCR_GPD_fitness_singles.txt",
                      "dataset/DiMSum_GRB2/GRB2_epPCR_GPD_fitness_doubles.txt",
                      "dataset/DiMSum_GRB2/GRB2_epPCR_CYC_fitness_singles.txt",
                      "dataset/DiMSum_GRB2/GRB2_epPCR_CYC_fitness_doubles.txt"),
    PDB_interaction_file = "dataset/PDB_contactmap_2vwf_AB.txt",
    RSA_file = "dataset/2vwf_A.rsa",
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum))
  )

  #####################
  ### dG estimation ###
  #####################

  #####################
  ### method 1: using data from both assays only for single variants

  # 1) with all parameters free to fit
  stagenum = 2
  function_dG_method1_singles_bothassays(
    name = dataset_name,
    dataset_file = paste0("processed_data/",dataset_name,"_dG_dataset.txt"),
    Ncores = Ncores,
    Nbootstraps = Nbootstraps*10,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum))
  )
  #minimal error if scaling factors are 0
  #whats wrong here?

  # 2) with background parameters set
  stagenum = 3
  function_dG_method1_singles_bothassays(
    name = dataset_name,
    dataset_file = paste0("processed_data/",dataset_name,"_dG_dataset.txt"),
    Ncores = Ncores,
    Nbootstraps = Nbootstraps*10,
    bgr_set = rep(0.55,2), #peak of fitness values for both assays in GRB2_epPCR data
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum))
  )

  #####################
  ### method 2: using data from both assays for all variants
  stagenum = 4
  function_dG_method2_allvars_bothassays(
    name = dataset_name,
    dataset_file = paste0("processed_data/",dataset_name,"_dG_dataset.txt"),
    Ncores = Ncores,
    Nbootstraps = Nbootstraps,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum))
  )


  #####################
  ### method 3: using only data from binding assay for all variants
  stagenum = 5
  function_dG_method3_allvars_bindingassay(
    name = dataset_name,
    dataset_file = paste0("processed_data/",dataset_name,"_dG_dataset.txt"),
    Ncores = Ncores,
    Nbootstraps = Nbootstraps,
    execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum))
  )
}

