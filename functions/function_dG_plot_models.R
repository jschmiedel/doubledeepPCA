function_dG_plot_models = function(
	models_name,
	dataset_file = "processed_data/GRB2_epPCR_dG_dataset.txt",
	model_dir = "processed_data/",
	print_dir = "results/dG/",
	which_model = "best") {

  require(ggplot2)
  require(GGally)
  theme_set(theme_bw(base_size=9))

	#load Rdata file with models
	load(file=file.path(model_dir,paste0(models_name,".Rdata")))

	## extract parameters and objective from models
	parameters = t(sapply(X=1:length(models),FUN = function(X){models[[X]]$par}))
  # global parameters
  global_par_idx = grep("[sb]_",colnames(parameters))
  global_par = data.table(parameters[,global_par_idx])
  # which/how many variants?
  Nvar = (ncol(parameters) - length(global_par_idx))/2
  id_var = names(models[[1]]$par[1:Nvar])
  # optim objective
  objective = sapply(X=1:length(models),FUN = function(X){models[[X]]$value})

  #add wild-type dG values if missing
  if (sum(names(global_par) %in% "s_dGwt") == 0) {
    global_par[,s_dGwt = function_folding_F2dG(1,s_bgr = s_bgr,s_scale = s_scale)]
  }
  if (sum(names(global_par) %in% "b_dGwt") == 0) {
    global_par[,b_dGwt = function_binding_F2dG(1,s_dGwt,b_bgr=b_bgr,b_scale = b_scale)]
  }

	
	#### first plot relationship between global parameters and objective
	dt = data.table(objective = log10(objective),global_par)
  p=ggpairs(dt)
  ggsave(plot = p,filename = file.path(print_dir,paste0(models_name,"_globalpar.pdf")),width=10,height=10)

  	#### next plot dG relationship for select model
	if (which_model == "best") {
		sel_model = which.min(objective)
	} else { #if which_model is integer
		sel_model = which_model
	}
	sel_parameters = parameters[sel_model,]

	global_par_sel = global_par[sel_model,]
	table_s_dG = sel_parameters[1:Nvar]
	table_b_dG = sel_parameters[(Nvar+1) : (2*Nvar)]

  
  	## plot global relationship between stability and binding versus fitness of variants
	all_data = fread(dataset_file)[id %in% id_var]
	setkey(all_data,id)
	all_data[names(table_s_dG),s_ddG := table_s_dG]
	all_data[names(table_b_dG),b_ddG := table_b_dG]

  if (sum(grep("s_scale",names(global_par))) > 0 ) {
  	all_data[,s_fitness_pred := function_folding_dG2F(s_ddG,
                          														global_par_sel$s_dGwt,
                                                      global_par_sel$s_bgr,
                                                      global_par_sel$s_scale),s_ddG]
    ggpairs_cols = c("s_fitness","s_fitness_pred","s_ddG","b_fitness","b_fitness_pred","b_ddG")
  } else {ggpairs_cols = c("s_fitness","s_ddG","b_fitness","b_fitness_pred","b_ddG")}
	all_data[,b_fitness_pred_bdG0 := function_binding_dG2F(0,
                            														s_ddG,
                            														global_par_sel$b_dGwt,
                            														global_par_sel$s_dGwt,
                                                        global_par_sel$b_bgr,
                                                        global_par_sel$b_scale),s_ddG]
  all_data[,b_fitness_pred := function_binding_dG2F(b_ddG,
                        														s_ddG,
                        														global_par_sel$b_dGwt,
                        														global_par_sel$s_dGwt,
                                                    global_par_sel$b_bgr,
                                                    global_par_sel$b_scale),.(b_ddG,s_ddG)]


    	## relationship between parameters and fitness

	p=ggpairs(all_data,
		columns = ggpairs_cols,
		aes(color=res_type1,alpha=0.3),title = which_model)
	ggsave(plot=p,file = file.path(print_dir,paste0(models_name,"_parameters_fitness.pdf")),width=10,height=10)
  

	## global relationship in bindingPCA and stabilityPCA space
  if (sum(grep("s_scale",names(global_par))) > 0 ) {
  	p=ggplot(all_data) +
    	geom_point(aes(b_fitness,s_fitness)) +
    	geom_point(aes(b_fitness_pred_bdG0,s_fitness_pred),color="red",size=2) +
    	geom_segment(aes(x=b_fitness_pred_bdG0,y=s_fitness_pred,xend = b_fitness,yend=s_fitness),alpha=0.3) +
    	labs(x="fitness bindingPCA",y="fitness stabilityPCA",color="used in fitting",title = which_model)
  	ggsave(plot=p,file = file.path(print_dir,paste0(models_name,"_fitness_globalrelationship.pdf")),width=6,height=6)
  }
}