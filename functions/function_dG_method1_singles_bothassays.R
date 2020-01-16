######## dG estimation method 1: calculate free energies of folding and binding from stabilityPCA and deepPCA data of single mutants ########
#### version: 3.0
#### created: 2019/11/27 by J??rn
#### last modified: 2020/01/15 by J??rn
#v3: parallelized v2

function_dG_method1_singles_bothassays =  function(
  name = "",
  dataset_file = "processed_data/GRB2_epPCR_dG_dataset.txt",
  Ncores = 15,
  Nbootstraps = 100,
  bgr_set = c()
) {
  require(data.table)
  require(foreach)
  require(doMC)
  require(ggplot2)
  require(GGally)
  theme_set(theme_bw(base_size=9))
  
  registerDoMC(cores=Ncores)

  #read dataset
  all_data = fread(dataset_file)
  setkey(all_data,Nmut,id1)
  
  #use only singles
  list_fitness_sigma = list()
  list_fitness_sigma[[1]] = matrix(all_data[Nmut==1,s_fitness])
  list_fitness_sigma[[2]] = matrix(all_data[Nmut==1,s_sigma])
  list_fitness_sigma[[3]] = matrix(all_data[Nmut==1,b_fitness])
  list_fitness_sigma[[4]] = matrix(all_data[Nmut==1,b_sigma])
  
  #which/how many variants?
  id_L = all_data[Nmut==1,.N]
  id_var = all_data[Nmut==1,id]
  
  models = foreach(m=1:Nbootstraps) %dopar% {
    ####### first, fit global relationship by assuming binding dG are all zero
    if (is.null(bgr_set)) {
      start_par = c(rep(0,id_L),
                    1e-5+(1-2e-5)*runif(1),
                    1e-5+(1-2e-5)*runif(1),
                    1e-5+(1-2e-5)*runif(1),
                    1e-5+(1-2e-5)*runif(1))
      global_par = c()
    } else {
      start_par = c(rep(0,id_L),
        1e-5+(1-2e-5)*runif(1),
        1e-5+(1-2e-5)*runif(1))
      global_par = bgr_set
    }
    fit_global_relationship = optim(par = start_par,
                                    fn = function_dG_method1_fitting,
                                    method = "L-BFGS-B",
                                    lower = c(rep(-10,id_L),rep(1e-5,ifelse(is.null(bgr_set),4,2))),
                                    upper = c(rep(10,id_L),rep(1-1e-5,ifelse(is.null(bgr_set),4,2))),
                                    control = list(),
                                    id_L = id_L,global_par = global_par,list_fitness_sigma=list_fitness_sigma)
    
    global_par = c(fit_global_relationship$par[(id_L+1) : (id_L+ifelse(is.null(bgr_set),4,2))],bgr_set)

    ####### second, estimate binding dGs (and again stability dG)
    fit_dGs = optim(par = c(fit_global_relationship$par[1:id_L],
                            rep(0,id_L)),
                    fn = function_dG_method1_fitting,
                    method = "L-BFGS-B",
                    lower = c(rep(-10,2*id_L)),
                    upper = c(rep(10,2*id_L)),
                    control = list(),
                    id_L = id_L,global_par=global_par,list_fitness_sigma=list_fitness_sigma)
    
    #gather model parameters
    fit_dGs$par = c(fit_dGs$par,
                    global_par,
                    function_folding_F2dG(1,global_par[3],global_par[1]),
                    function_binding_F2dG(1,function_folding_F2dG(1,global_par[3],global_par[1]),global_par[4],global_par[2]))
    names(fit_dGs$par) = c(id_var,id_var,"s_scale","b_scale","s_bgr","b_bgr","s_dGwt","b_dGwt")

    return(fit_dGs)
  }
  #save
  save(models,file=paste0("processed_data/dG_method1_",name,"_",Nbootstraps,"models.Rdata"))
  
  #plot models
  function_dG_plot_models(models_name = paste0("dG_method1_",name,"_",Nbootstraps,"models"),
                          dataset_file = dataset_file)
}




