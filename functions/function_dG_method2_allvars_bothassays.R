function_dG_method2_allvars_bothassays = function(
	name = "",
	dataset_file = "processed_data/GRB2_epPCR_dG_dataset.txt",
	Ncores = 15,
	Nbootstraps = 30)
{
	require(data.table)
	require(foreach)
	require(doMC)
	require(ggplot2)
	require(GGally)
	theme_set(theme_bw(base_size=9))

	registerDoMC(cores=Ncores)

  	#read dataset
  	all_data = fread(dataset_file)
  	setkey(all_data,Nmut,id1,id2)

  	#which variants are present?
  	id_var = unique(all_data[!is.na(s_fitness) & !is.na(b_fitness),id1])
  	all_data[,id1_key := which(id_var == id1),id1]
  	all_data[,id2_key := which(id_var == id2),id2]
  	id_L = length(id_var)

  	#10-fold cross validation
  	for (Xval in 1:10) {
  		print(Xval)
  		
		#grab singles and doubles
		singles = all_data[Nmut==1,.(id1_key,s_fitness,b_fitness,s_sigma,b_sigma)]
		doubles = all_data[Nmut==2 & tenfold_Xval != Xval,.(id1_key,id2_key,s_fitness,b_fitness,s_sigma,b_sigma)]
		#key them
		list_keys = list(singles_key = singles[,id1_key],
			doubles_key1 = doubles[,id1_key],
			doubles_key2 = doubles[,id2_key])

		#create fitness & error vectors
		list_fs = list(s_fitness = c(singles[,s_fitness],doubles[,s_fitness]),
			s_sigma = c(singles[,s_sigma],doubles[,s_sigma]),
			b_fitness = c(singles[,b_fitness],doubles[,b_fitness]),
			b_sigma = c(singles[,b_sigma],doubles[,b_sigma]))

		# fit parameters
		registerDoMC(cores=Ncores)
		models = foreach(m=1:Nbootstraps) %dopar% {

			global_model = optim(par = c(rep(0,id_L),
				rep(0,id_L),
				rep(1e-5,4)),
			fn = function_dG_method2_fitting,
			method = "L-BFGS-B",
			lower = c(rep(-10,2*id_L),rep(1e-5,4)),
			upper = c(rep(10,2*id_L),rep(1-1e-5,4)),
			control = list(),
			id_L = id_L,list_fs=list_fs,list_keys=list_keys)

			#gather model parameters
			global_model$par = c(global_model$par,
				function_folding_F2dG(1,global_model$par[2*id_L+3],global_model$par[2*id_L+1]),
				function_binding_F2dG(1,function_folding_F2dG(1,global_model$par[2*id_L+3],global_model$par[2*id_L+1]),global_model$par[2*id_L+4],global_model$par[2*id_L+2]))
			names(global_model$par) = c(id_var,id_var,"s_scale","b_scale","s_bgr","b_bgr","s_dGwt","b_dGwt")

			global_model$Xval = Xval

			return(global_model)
		}
	#save
	save(models,file=paste0("processed_data/dG_method2_",name,"_",Nbootstraps,"models_",Xval,"Xval.Rdata"))

	#plot models
	function_dG_plot_models(models_name = paste0("dG_method2_",name,"_",Nbootstraps,"models_",Xval,"Xval"),
							dataset_file = dataset_file)
	}
}
