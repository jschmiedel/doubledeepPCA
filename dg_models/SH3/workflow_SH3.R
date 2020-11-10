# fit 3 and 4 state thermodynamic models to SH3 doubledeepPCA data
require(tempura)

dataset_folder = "dg_models/SH3"
model_name3 = "three_state"
model_name4 = "four_state"

## 1) preprocess DiMSum outputs using script "doubledeepPCA/dg_models/SH3/dataset_preprocessing"

## 2) prepare datasets
dg_prepare_datasets(dataset_folder = dataset_folder,
    abundancepca_files = c(
        paste0(dataset_folder, "/data/01a-GRB2_epPCR_stabilityPCA_aa012_thresholded.RData"),
        paste0(dataset_folder, "/data/01c-GRB2_NM2_stabilityPCA_aa012_thresholded.RData")),
    bindingpca_files = paste0(dataset_folder, "/data/01b-GRB2_epPCR_bindingPCA_aa012_thresholded.RData"),
    wt_seq = "TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYVTPVN*"
)


################# three state model
## 3a) prepare dg model parameters for three state model
dg_prepare_model(
    dataset_folder = dataset_folder,
    model_name = model_name3,
    fix_f_dgwt = TRUE,
    fix_b_dgwt = TRUE
)

## 4a) estimate model parameters
# qsub -t 1:500 doubledeepPCA/dg_models/bash_model_estimation.sh SH3 three_state

## 5a) collect models
dg_collect_models(
    dataset_folder = dataset_folder,
    model_name = model_name3,
    model_averaging = "mean"
)

## 6a) basic plots
dg_basic_analyses(
    dataset_folder = dataset_folder,
    model_name = model_name3,
    color_type = "type",
    stage = "model",
    datasets_ab = c(1,1)
)

## 7a) bootstrap model parameters
# qsub -t 1:100 doubledeepPCA/dg_models/bash_bootstrap.sh SH3 three_state

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
    no_folded_states = 2,
    fix_f_dgwt = TRUE,
    fix_b_dgwt = TRUE
)

## 4b) estimate model parameters
# qsub -t 1:500 doubledeepPCA/dg_models/bash_model_estimation.sh SH3 four_state

## 5b) collect models
dg_collect_models(
    dataset_folder = dataset_folder,
    model_name = model_name4,
    model_averaging = "mean"
)

## 6b) basic plots
dg_basic_analyses(
    dataset_folder = dataset_folder,
    model_name = model_name4,
    color_type = "type",
    stage = "model",
    datasets_ab = c(1,1)
)

## 7b) bootstrap model parameters
# qsub -t 1:100 doubledeepPCA/dg_models/bash_bootstrap.sh SH3 four_state


## 8b) collect bootstrapped models
dg_collect_models(
    dataset_folder = dataset_folder,
    model_name = model_name4,
    stage = "bootstrap"
)
dg_basic_analyses(
    dataset_folder = dataset_folder,
    model_name = model_name4,
    color_type = "type",
    datasets_ab = c(1,1)
)
