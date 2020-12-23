# fit 3 and 4 state thermodynamic models to SH3 doubledeepPCA data
require(tempura)

dataset_folder = "dg_models/PDZ3"
model_name3 = "three_state"
model_name3_fix = "three_state_fixdgwt"
# model_name4 = "four_state"
# model_name4_fix = "four_state_fixdgwt"

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


################# three state model with fixed dG wildtype parameters
## 3a) prepare dg model parameters for three state model
dg_prepare_model(
    dataset_folder = dataset_folder,
    model_name = model_name3_fix,
    fix_f_dgwt = TRUE,
    fix_b_dgwt = TRUE
)

## 4a) estimate model parameters
# qsub -t 1:500 doubledeepPCA/dg_models/bash_model_estimation.sh PDZ3 model_name3_fix

## 5a) collect models
dg_collect_models(
    dataset_folder = dataset_folder,
    model_name = model_name3_fix,
    model_averaging = "mean"
)

## 6a) plot basic analysis plots
dg_basic_analyses(
    dataset_folder = dataset_folder,
    model_name = model_name3_fix,
    color_type = "type",
    stage = "model",
    datasets_ab = c(1,1)
)

## 7a) bootstrap model parameters
# qsub -t 1:100 doubledeepPCA/dg_models/bash_bootstrap.sh PDZ3 model_name3_fix

## 8a) collect bootstrapped models
dg_collect_models(
    dataset_folder = dataset_folder,
    model_name = model_name3_fix,
    stage = "bootstrap"
)
dg_basic_analyses(
    dataset_folder = dataset_folder,
    model_name = model_name3_fix,
    color_type = "type",
    datasets_ab = c(1,1)
)


################# three state model with variable dG wildtype parameters
## 3b) prepare dg model parameters for three state model
dg_prepare_model(
    dataset_folder = dataset_folder,
    model_name = model_name3,
    fix_f_dgwt = FALSE,
    fix_b_dgwt = FALSE
)

## 4b) estimate model parameters
# qsub -t 1:500 doubledeepPCA/dg_models/bash_model_estimation.sh PDZ3 three_state

## 5b) collect models
dg_collect_models(
    dataset_folder = dataset_folder,
    model_name = model_name3,
    model_averaging = "mean"
)

## 6b) plot basic analysis plots
dg_basic_analyses(
    dataset_folder = dataset_folder,
    model_name = model_name3,
    color_type = "type",
    stage = "model",
    datasets_ab = c(1,1)
)

## 7b) bootstrap model parameters
# qsub -t 1:100 doubledeepPCA/dg_models/bash_bootstrap.sh PDZ3 three_state

## 8b) collect bootstrapped models
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


# ################# four state model
# ## 3c) prepare dg model parameters for four state model
# dg_prepare_model(
#     dataset_folder = dataset_folder,
#     model_name = model_name4_fix,
#     no_folded_states = 2,
#     fix_f_dgwt = TRUE,
#     fix_b_dgwt = TRUE
# )

# ## 4c) estimate model parameters
# # qsub -t 1:500 doubledeepPCA/dg_models/bash_model_estimation.sh PDZ3 four_state

# ## 5c) collect models
# dg_collect_models(
#     dataset_folder = dataset_folder,
#     model_name = model_name4,
#     model_averaging = "median"
# )
# ## two parameter subsets:
# dg_collect_models(
#     dataset_folder = dataset_folder,
#     model_name = model_name4,
#     which_test_set = "[1578]",
#     model_averaging = "mean"
# )
# dg_collect_models(
#     dataset_folder = dataset_folder,
#     model_name = model_name4,
#     which_test_set = "[23469]",
#     model_averaging = "mean"
# )

# ## 6c) plot basic analysis plots
# dg_basic_analyses(
#     dataset_folder = dataset_folder,
#     model_name = model_name4,
#     color_type = "type",
#     datasets_ab = c(1,1)
# )

# ## 4 state model has issues, with weird fA<>fB distribution etc
# ## check whether this might be because of lambda = 1e-1 or because of the older datasets
