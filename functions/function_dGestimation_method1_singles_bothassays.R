######## dG estimation method 1: calculate free energies of folding and binding from stabilityPCA and deepPCA data of single mutants ########
#### version: 3.0
#### created: 2019/11/27 by J??rn
#### last modified: 2020/01/15 by J??rn
#v3: parallelized v2

function_dGestimation_method1_singles_bothassays =  function(
  name = "GRB2_epPCR",
  dataset_file = "processed_data/GRB2_epPCR_dG_dataset.txt",
  Ncores = 15,
  Nbootstraps = 100
) {
  require(data.table)
  require(foreach)
  require(doMC)
  require(ggplot2)
  require(GGally)
  theme_set(theme_bw(base_size=9))
  
  #parameters
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  
  # DT[,is.global := T]
  # old data
  DT[,is.global := binding > 0.75  | binding < 0.75 & stability < 0.55]
  #new DiMSum data
  # DT[,is.global := binding > -1  | binding < -1 & stability < -1.75 & stability > -3.5]
  
  binding_fitness = DT[is.global==T,binding]
  stability_fitness = DT[is.global==T,stability]
  binding_error = DT[is.global==T,sigma_binding]
  stability_error = DT[is.global==T,sigma_stab]
  
  L = DT[is.global==T,.N]
  
  # fix background
  b_bgr = 0.3
  f_bgr = 0.35
  
  # use regularization
  # lambda = 1e-8
  
  #new DiMSum data
  # models = list()
  # objective = c()
  # for (i in 1:10) {
  
  registerDoMC(cores=Ncores)
  models = foreach(m=1:10) %dopar% {
    # print(i)
    # fit_global_relationship = nlminb(start = c(rep(0,L),
    #                                            0.1+0.3*runif(2),
    #                                            0.9+0.1*runif(1),
    #                                            0.9+0.1*runif(1),
    #                                            rep(0,L)),
    #                                  objective = function_fit_global_relationship,
    #                                  lower = c(rep(-10,L),1e-5,1e-5,1e-5,1e-5,rep(-10,L)),
    #                                  upper = c(rep(10,L),1-1e-5,1-1e-5,1-1e-5,1-1e-5,rep(10,L)),
    #                                  control = list())
    # objective[i] = fit_global_relationship$objective
    # fit_global_relationship = optim(par = c(rep(0,L),
    #                                         0.1+0.3*runif(2),
    #                                         0.9+0.1*runif(1),
    #                                         0.9+0.1*runif(1)),
    #                                 fn = function_fit_global_relationship,
    #                                 method = "L-BFGS-B",
    #                                 lower = c(rep(-10,L),1e-5,1e-5,1e-5,1e-5),
    #                                 upper = c(rep(10,L),1-1e-5,1-1e-5,1-1e-5,1-1e-5),
    #                                 control = list())
    fit_global_relationship = optim(par = c(rep(0,L),
                                            0.8+0.2*runif(1),
                                            0.8+0.2*runif(1)),
                                    fn = function_fit_global_relationship,
                                    method = "L-BFGS-B",
                                    lower = c(rep(-10,L),1e-5,1e-5),
                                    upper = c(rep(10,L),1-1e-5,1-1e-5),
                                    control = list())
    # fit_global_relationship = optim(par = c(rep(0,L),
    #                                         0.1+0.3*runif(2),
    #                                         0.9+0.1*runif(1),
    #                                         0.9+0.1*runif(1),
    #                                         rep(0,L)),
    #                                 fn = function_fit_global_relationship,
    #                                 method = "L-BFGS-B",
    #                                 lower = c(rep(-10,L),1e-5,1e-5,1e-5,1e-5,rep(-10,L)),
    #                                 upper = c(rep(10,L),1-1e-5,1-1e-5,1-1e-5,1-1e-5,rep(10,L)),
    #                                 control = list())
    # objective[i] = fit_global_relationship$value
    # models[[i]] = fit_global_relationship
    return(fit_global_relationship)
    
  }
  
  parameters = t(sapply(X=1:length(models),FUN = function(X){models[[X]]$par}))
  objective = sapply(X=1:length(models),FUN = function(X){models[[X]]$value})
  # dt = data.table(log10(objective),parameters[,(L+1) : (L+4)])
  # names(dt) = c("objective","b_bgr","f_bgr","b_scale","f_scale")
  dt = data.table(log10(objective),parameters[,(L+1) : (L+2)])
  names(dt) = c("objective","b_scale","f_scale")
  ggpairs(dt)
  
  fit_global_relationship = models[[which(rank(objective)==1)]]
  # table_global_parameters = data.table(b_bgr = fit_global_relationship$par[L+1],
  # f_bgr = fit_global_relationship$par[L+2],
  # b_scale = fit_global_relationship$par[L+3],
  # f_scale = fit_global_relationship$par[L+4])
  table_global_parameters = data.table(b_bgr = 0.3,
                                       f_bgr = 0.35,
                                       b_scale = fit_global_relationship$par[L+1],
                                       f_scale = fit_global_relationship$par[L+2])
  
  #calcualte again dG wildtype values from final parameters
  #old
  f_dGwt = function_folding_F2dG(1,table_global_parameters$f_bgr,table_global_parameters$f_scale)
  b_dGwt = function_binding_F2dG(1,f_dGwt,table_global_parameters$b_bgr,table_global_parameters$b_scale)
  #new DiMSum data
  # f_dGwt = function_folding_F2dG(0,table_global_parameters$f_bgr,table_global_parameters$f_s1,table_global_parameters$f_scale)
  # b_dGwt = function_binding_F2dG(0,f_dGwt,table_global_parameters$b_bgr,table_global_parameters$b_s1,table_global_parameters$b_scale)
  
  table_global_parameters[,b_dGwt := b_dGwt]
  table_global_parameters[,f_dGwt := f_dGwt]
  print(table_global_parameters)
  write.table(x = table_global_parameters,file = "dataset/dG/method1_globalpar.txt",row.names=F,col.names=T,quote=F)
  
  
  #plot results
  DT[is.global==T,f_ddG_global := fit_global_relationship$par[1:L]]
  # DT[is.global==T,b_ddG_global := fit_global_relationship$par[(L+5):(2*L+4)]]
  DT[is.global==T,f_fitness_global := function_folding_dG2F(f_ddG_global,f_dGwt,
                                                            table_global_parameters$f_bgr,table_global_parameters$f_scale),f_ddG_global]
  DT[is.global==T,b_fitness_global := function_binding_dG2F(0,f_ddG_global,b_dGwt,f_dGwt,
                                                            table_global_parameters$b_bgr,table_global_parameters$b_scale),f_ddG_global]
  
  #parameter/fitness scatter
  # ggpairs(DT[is.global==T],columns = c("stability","f_fitness_global","f_ddG_global","binding","b_fitness_global","b_ddG_global"),
  # columnLabels = c("Fstab","pred. Fstab","ddGstab","Fbind","pred. Fbind","ddGbind"))
  ggpairs(DT[is.global==T],columns = c("stability","f_fitness_global","f_ddG_global","binding","b_fitness_global"),
          columnLabels = c("Fstab","pred. Fstab","ddGstab","Fbind","pred. Fbind"))
  ggsave("figures/dG/method1_singles_dG_fromglobalparameters.pdf",width=10,height=7)
  
  #global relationship in bindingPCA and stabilityPCA space
  ggplot(DT) +
    geom_point(aes(binding,stability,color=is.global)) +
    geom_point(aes(b_fitness_global,f_fitness_global),color="red",size=2) +
    geom_segment(aes(x=b_fitness_global,y=f_fitness_global,xend = binding,yend=stability),alpha=0.3) +
    labs(x="fitness bindingPCA",y="fitness stabilityPCA",color="used in fitting")
  ggsave("figures/dG/method1_singles_binding_stability_relationship.pdf",width=6,height=6)
  
  ggplot(DT,aes(f_dGwt+f_ddG_global)) + 
    stat_ecdf() + 
    geom_vline(xintercept = f_dGwt,color="red") +
    geom_vline(xintercept = 0,color="black")
  
  
  ## now fit b_ddG and f_ddG
  function_fit_dGfromglobal = function(parameters) {
    bf = function_binding_dG2F(parameters[1:L],parameters[(L+1):(L*2)],
                               table_global_parameters$b_dGwt,
                               table_global_parameters$f_dGwt,
                               table_global_parameters$b_bgr,
                               # table_global_parameters$b_s1,
                               table_global_parameters$b_scale)
    ff = function_folding_dG2F(parameters[(L+1):(L*2)],
                               table_global_parameters$f_dGwt,
                               table_global_parameters$f_bgr,
                               # table_global_parameters$f_s1,
                               table_global_parameters$f_scale)
    binding_fitness2 = binding_fitness
    binding_fitness2[is.na(bf)] = NA
    stability_fitness2 = stability_fitness
    stability_fitness2[is.na(ff)] = NA
    binding_error2 = binding_error
    binding_error2[is.na(bf)] = NA
    stability_error2 = stability_error
    stability_error2[is.na(ff)] = NA
    
    #mean square deviation; divide by weights to correct for variants being NA
    MSD = sum((binding_fitness2 - bf)^2 / binding_error2^2 + (stability_fitness2 - ff)^2 / stability_error2^2,na.rm=T) / sum(binding_error2^-2 + stability_error2^-2,na.rm=T)
    return(MSD)
  }
  
  # DT2 = copy(GRB2_singles)
  binding_fitness = DT[,binding]
  stability_fitness = DT[,stability]
  binding_error = DT[,sigma_binding]
  stability_error = DT[,sigma_stab]
  L = nrow(DT)
  
  fit_dGfromglobal = nlminb(start = c(rep(0,2*L)),
                            # start = c(rnorm(n = L+2),0.3,0.3,0,1,0,1),
                            objective = function_fit_dGfromglobal,
                            control = list(iter.max = 500))
  print(fit_dGfromglobal)
  
  DT[,b_ddG_local := fit_dGfromglobal$par[1:L]]
  DT[,f_ddG_local := fit_dGfromglobal$par[(L+1):(L*2)]]
  DT[,f_fitness_local := function_folding_dG2F(f_ddG_local,f_dGwt,
                                               table_global_parameters$f_bgr,table_global_parameters$f_scale),f_ddG_local]
  DT[,b_fitness_local := function_binding_dG2F(0,f_ddG_local,b_dGwt,f_dGwt,
                                               table_global_parameters$b_bgr,table_global_parameters$b_scale),f_ddG_local]
  
  #save
  write.table(x=DT,file="dataset/dG/method1_bestmodel.txt",row.names=F,quote=F,col.names=T)
  
  helper=melt(DT,measure.vars = grep("ddG",names(DT),value=T))
  helper[,variable := factor(variable,levels=c("f_ddG_global","f_ddG_local","b_ddG_local"))]
  ggplot(helper) +
    geom_point(aes(binding,stability,color=value)) +
    geom_line(aes(b_fitness_local,f_fitness_local),color="red",size=2) +
    geom_segment(aes(x=b_fitness_local,y=f_fitness_local,xend = binding,yend=stability),alpha=0.3) +
    scale_color_gradient2(midpoint=0,mid="grey",na.value = "black") +
    facet_wrap(variable~.,nrow = 2) +
    labs(x="fitness bindingPCA",y="fitness stabilityPCA",color="ddG")
  ggsave("figures/dG/method1_singles_ddG_global_local.pdf",width=6,height=6)
  
  
  ggpairs(DT,columns = c("binding","stability","b_ddG_local","f_ddG_global","f_ddG_local"),
          aes(color=type,alpha=0.5))
  ggsave("figures/dG/method1_singles_fitness_ddG.pdf",width=10,height=8)
  
  #pleiotropy
  DT[,pleiotropy := "none"]
  DT[b_ddG_local > 0.5 & f_ddG_local > 1.5,pleiotropy := "b_ddG > 0.5 & f_ddG > 1.5"]
  DT[b_ddG_local > 1 & f_ddG_local > 3,pleiotropy := "b_ddG > 1 & f_ddG > 3"]
  p1=ggplot(DT) + 
    geom_density2d(aes(b_ddG_local,f_ddG_local,color=type)) + 
    geom_point(aes(b_ddG_local,f_ddG_local,color=type,shape=pleiotropy))
  
  p2=ggplot(DT) +
    geom_point(aes(binding,stability,color=pleiotropy)) +
    geom_line(aes(b_fitness_local,f_fitness_local),color="red",size=2) +
    geom_segment(inherit.aes = F,data=DT[pleiotropy!="none"],aes(x=b_fitness_local,y=f_fitness_local,xend = binding,yend=stability),alpha=0.3)
  p=grid.arrange(p1,p2,nrow=2)
  ggsave(plot=p,file="figures/dG/method1_singles_ddG_binding_stability.pdf",width=6,height=7)
}




