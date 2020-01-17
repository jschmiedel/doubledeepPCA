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

	#extract parameters and objective from models
	parameters = t(sapply(X=1:length(models),FUN = function(X){models[[X]]$par}))
	Lpar = ncol(parameters)
	id_var = names(models[[1]]$par[1:((Lpar-6)/2)])
	objective = sapply(X=1:length(models),FUN = function(X){models[[X]]$value})

	#### first plot relationship between global parameters and objective
	dt = data.table(log10(objective),parameters[,(Lpar-5) : (Lpar-2)])
  	names(dt) = c("objective","s_scale","b_scale","s_bgr","b_bgr")
  	p=ggpairs(dt)
  	ggsave(plot = p,filename = file.path(print_dir,paste0(models_name,"_globalpar.pdf")),width=10,height=10)

  	#### next plot dG relationship for select model
	if (which_model == "best") {
		sel_model = which.min(objective)
	} else { #if which_model is integer
		sel_model = which_model
	}
	sel_parameters = parameters[sel_model,]

	table_globalpar = data.table(t(sel_parameters[(Lpar-5) : (Lpar)]))
	table_s_dG = sel_parameters[1:((Lpar-6)/2)]
	table_b_dG = sel_parameters[((Lpar-6)/2+1) : (Lpar-6)]

	## plot global relationship between stability and binding versus fitness of variants
	all_data = fread(dataset_file)[id %in% id_var]
	setkey(all_data,id)
	all_data[names(table_s_dG),s_ddG := table_s_dG]
  	all_data[names(table_b_dG),b_ddG := table_b_dG]
  	all_data[,s_fitness_pred := function_folding_dG2F(s_ddG,
  														table_globalpar$s_dGwt,
                                                        table_globalpar$s_bgr,
                                                        table_globalpar$s_scale),s_ddG]
  	all_data[,b_fitness_pred_bdG0 := function_binding_dG2F(0,
  														s_ddG,
  														table_globalpar$b_dGwt,
  														table_globalpar$s_dGwt,
                                                        table_globalpar$b_bgr,
                                                        table_globalpar$b_scale),s_ddG]
	all_data[,b_fitness_pred := function_binding_dG2F(b_ddG,
  														s_ddG,
  														table_globalpar$b_dGwt,
  														table_globalpar$s_dGwt,
                                                        table_globalpar$b_bgr,
                                                        table_globalpar$b_scale),.(b_ddG,s_ddG)]


  	## relationship between parameters and fitness
  	p=ggpairs(all_data,
  		columns = c("s_fitness","s_fitness_pred","s_ddG","b_fitness","b_fitness_pred","b_ddG"),
  		aes(color=res_type1,alpha=0.3),title = which_model)
  	ggsave(plot=p,file = file.path(print_dir,paste0(models_name,"_parameters_fitness.pdf")),width=10,height=10)

	## global relationship in bindingPCA and stabilityPCA space
  	p=ggplot(all_data) +
    	geom_point(aes(b_fitness,s_fitness)) +
    	geom_point(aes(b_fitness_pred_bdG0,s_fitness_pred),color="red",size=2) +
    	geom_segment(aes(x=b_fitness_pred_bdG0,y=s_fitness_pred,xend = b_fitness,yend=s_fitness),alpha=0.3) +
    	labs(x="fitness bindingPCA",y="fitness stabilityPCA",color="used in fitting",title = which_model)
  	ggsave(plot=p,file = file.path(print_dir,paste0(models_name,"_fitness_globalrelationship.pdf")),width=6,height=6)
}