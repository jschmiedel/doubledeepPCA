function_dG_plot_globalparameter_relationships = function(models_file,print_file) {
  load(file=models_file)
  parameters = t(sapply(X=1:length(models),FUN = function(X){models[[X]]$par}))
  Lpar = ncol(parameters)
  objective = sapply(X=1:length(models),FUN = function(X){models[[X]]$value})
  dt = data.table(log10(objective),parameters[,(Lpar-3) : (Lpar)])
  names(dt) = c("objective","s_scale","b_scale","s_bgr","b_bgr")
  p=ggpairs(dt)
  ggsave(plot = p,filename = print_file,width=10,height=10)
  
}
