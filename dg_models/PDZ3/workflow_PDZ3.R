# fit 3 and 4 state thermodynamic models to SH3 doubledeepPCA data
require(tempura)

dataset_folder = "/nfs/users/blehner/jschmiedel/doubledeepPCA/dg_models/PDZ3"
model_name3 = "three_state"
model_name4 = "four_state"
model_name4_l0 = "four_state_l0"

setwd(dataset_folder)


## 1) preprocess DiMSum outputs using script "doubledeepPCA/dg_models/SH3/dataset_preprocessing"

## 2) prepare datasets
dg_prepare_datasets(
    dataset_folder = dataset_folder,
    abundancepca_files = c(
        paste0(dataset_folder, "/data/01e-PDZ3_NM2_abundancePCA_thresholded.RData"),
        paste0(dataset_folder, "/data/01c-PDZ3_NM1_abundancePCA_thresholded.RData"),
        paste0(dataset_folder, "/data/01a-PDZ3_epPCR_abundancePCA_thresholded.RData")),
    bindingpca_files = c(
        paste0(dataset_folder, "/data/01f-PDZ3_NM2_bindingPCA_thresholded.RData"),
        paste0(dataset_folder, "/data/01d-PDZ3_NM1_bindingPCA_thresholded.RData"),
        paste0(dataset_folder, "/data/01b-PDZ3_epPCR_bindingPCA_thresholded.RData")),
    wt_seq = "PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTIIAQYKP*"
)


################# three state model
## 3a) prepare dg model parameters for three state model
dg_prepare_model(
    dataset_folder = dataset_folder,
    model_name = model_name3
)

## 4a) estimate model parameters
# qsub doubledeepPCA/dg_models/PDZ3/bash_3state_model_estimation.sh

## 5a) collect models
dg_collect_models(
    dataset_folder = dataset_folder,
    model_name = model_name3,
    model_averaging = "mean"
)

## 6a) plot basic analysis plots
dg_basic_analyses(
    dataset_folder = dataset_folder,
    model_name = model_name3,
    color_type = "type",
    stage = "model",
    datasets_ab = c(1,1)
)

## 7a) bootstrap model parameters
# qsub doubledeepPCA/dg_models/PDZ3/bash_3state_bootstrap.sh

## 8a) collect bootstrapped models
dg_collect_models(
    dataset_folder = dataset_folder,
    model_name = model_name3,
    stage = "bootstrap"
)
dg_basic_analyses(
    dataset_folder = dataset_folder,
    model_name = model_name3,
    color_type = "type",
    datasets_ab = c(1,1)
)


################# four state model
## 3b) prepare dg model parameters for three state model
dg_prepare_model(
    dataset_folder = dataset_folder,
    model_name = model_name4,
    no_folded_states = 2
)

## 4b) estimate model parameters
# qsub doubledeepPCA/dg_models/PDZ3/bash_4state_model_estimation.sh

## 5b) collect models
dg_collect_models(
    dataset_folder = dataset_folder,
    model_name = model_name4,
    model_averaging = "mean"
)

## 6b) plot basic analysis plots
dg_basic_analyses(
    dataset_folder = dataset_folder,
    model_name = model_name4,
    color_type = "type",
    datasets_ab = c(1,1)
)

## 4 state model has issues, with weird fA<>fB distribution etc
## check whether this might be because of lambda = 1e-1 or because of the older datasets


################# four state model with lambda = 0
## 3c) prepare dg model parameters for three state model
dg_prepare_model(
    dataset_folder = dataset_folder,
    model_name = model_name4_l0,
    no_folded_states = 2,
    lambda = 0
)

## 4c) estimate model parameters
# qsub doubledeepPCA/dg_models/PDZ3/bash_4state_model_estimation.sh

## 5c) collect models
dg_collect_models(
    dataset_folder = dataset_folder,
    model_name = model_name4_l0,
    model_averaging = "mean"
)

## 6c) plot basic analysis plots
dg_basic_analyses(
    dataset_folder = dataset_folder,
    model_name = model_name4_l0,
    color_type = "type",
    datasets_ab = c(1,1)
)


################# three/four state model with only the newest dataset
dataset_folder_new = "/nfs/users/blehner/jschmiedel/doubledeepPCA/dg_models/PDZ3_newonly"
model_name3_new = "three_state"
model_name4_new = "four_state"
model_name4_new_lm2 = "four_state_lm2"
setwd(dataset_folder_new)
## 2d) prepare datasets
dg_prepare_datasets(
    dataset_folder = dataset_folder_new,
    abundancepca_files = c(
        paste0(dataset_folder_new, "/data/01e-PDZ3_NM2_abundancePCA_thresholded.RData")),
    bindingpca_files = c(
        paste0(dataset_folder_new, "/data/01f-PDZ3_NM2_bindingPCA_thresholded.RData")),
    wt_seq = "PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTIIAQYKP*"
)


############# three state
## 3d) prepare dg model parameters for three state model
dg_prepare_model(
    dataset_folder = dataset_folder_new,
    model_name = model_name3_new,
    no_folded_states = 1
)

## 4d) estimate model parameters
# qsub doubledeepPCA/dg_models/PDZ3_newonly/bash_3state_model_estimation.sh

## 5d) collect models
dg_collect_models(
    dataset_folder = dataset_folder_new,
    model_name = model_name3_new,
    model_averaging = "mean"
)

## 6d) plot basic analysis plots
dg_basic_analyses(
    dataset_folder = dataset_folder_new,
    model_name = model_name3_new,
    color_type = "type",
    datasets_ab = c(1,1)
)
dg_basic_analyses(
    dataset_folder = dataset_folder_new,
    model_name = model_name3_new,
    color_type = "HAmin_ligand",
    datasets_ab = c(1,1)
)

############# four state
## 3e) prepare dg model parameters for three state model
dg_prepare_model(
    dataset_folder = dataset_folder_new,
    model_name = model_name4_new,
    no_folded_states = 2
)

## 4e) estimate model parameters
# qsub doubledeepPCA/dg_models/PDZ3_newonly/bash_4state_model_estimation.sh

## 5e) collect models
dg_collect_models(
    dataset_folder = dataset_folder_new,
    model_name = model_name4_new,
    model_averaging = "mean"
)

## 6e) plot basic analysis plots
dg_basic_analyses(
    dataset_folder = dataset_folder_new,
    model_name = model_name4_new,
    color_type = "type",
    datasets_ab = c(1,1)
)


############# four state with lambda = 1e-2
## 3f) prepare dg model parameters for three state model
dg_prepare_model(
    dataset_folder = dataset_folder_new,
    model_name = model_name4_new_lm2,
    no_folded_states = 2,
    lambda = 1e-2
)

## 4f) estimate model parameters
# qsub doubledeepPCA/dg_models/PDZ3_newonly/bash_4state_model_estimation.sh

## 5f) collect models
dg_collect_models(
    dataset_folder = dataset_folder_new,
    model_name = model_name4_new_lm2,
    model_averaging = "mean"
)

## 6f) plot basic analysis plots
dg_basic_analyses(
    dataset_folder = dataset_folder_new,
    model_name = model_name4_new_lm2,
    color_type = "type",
    datasets_ab = c(1,1)
)
