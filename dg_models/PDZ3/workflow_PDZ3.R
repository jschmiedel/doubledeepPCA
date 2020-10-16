# fit 3 and 4 state thermodynamic models to SH3 doubledeepPCA data
require(tempura)

dataset_folder = "/nfs/users/blehner/jschmiedel/doubledeepPCA/dg_models/PDZ3"
model_name3 = "three_state"
model_name4 = "four_state"

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
    model_name = model_name3
)

## 6a) bootstrap model parameters
# qsub doubledeepPCA/dg_models/PDZ3/bash_3state_bootstrap.sh

## 7a) collect bootstrapped models
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
    model_name = model_name4
)
