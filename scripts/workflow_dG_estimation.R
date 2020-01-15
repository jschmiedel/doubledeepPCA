

##################
###### GRB2 ######
##################

### prepare dataset
#epPCR datasets
function_prepare_dG_dataset(
  name = "GRB2_epPCR",
  DMS_file_list = c("dataset/DiMSum_GRB2/GRB2_epPCR_GPD_fitness_singles.txt",
                    "dataset/DiMSum_GRB2/GRB2_epPCR_GPD_fitness_doubles.txt",
                    "dataset/DiMSum_GRB2/GRB2_epPCR_CYC_fitness_singles.txt",
                    "dataset/DiMSum_GRB2/GRB2_epPCR_CYC_fitness_doubles.txt"),
  PDB_interaction_file = "dataset/PDB_contactmap_2vwf_AB.txt",
  RSA_file = "dataset/2vwf_A.rsa"
)