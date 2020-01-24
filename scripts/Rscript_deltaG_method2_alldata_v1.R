######## calculate free energies of folding and binding from stabilityPCA and deepPCA data ########
#### version: 1.0
#### created: 2019/12/04 by J??rn
#### last modified: 2019/12/04 by J??rn

### method 2: use single & double mutant data from both assays


require(data.table)
require(ggplot2)
require(GGally)
require(foreach)
require(doMC)
require(Matrix)
theme_set(theme_bw(base_size=9))


#julia
setwd("deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/") #This is the only path common to any computer
#joern
setwd("/Users/jschmidel/GoogleDrive/deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/")
Ncores=3
## Joern new Macbook
# setwd("/Users/jschmiedel/GoogleDrive/deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/")
# Ncores=2

####################
### prepare data ###
####################

#constants/parameters
R= 1.98*10^(-3) # kcal/mol
Temp= 303.15



#####################
### dG estimation ###
#####################

## set up indicies for mutant:variant matching
setkey(all_data,id1,id2)
id_unique = unique(c(all_data_wide[!is.na(s_fitness) & !is.na(b_fitness),id1],all_data_wide[!is.na(s_fitness) & !is.na(b_fitness) & !is.na(id2),id2]))
id_unique = id_unique[!is.na(id_unique)]
# id_unique = names(sort(table(c(all_data$id1,all_data$id2)),decreasing = T)[1:100])
id_L = length(id_unique)

all_data[,id1_key := which(id_unique == id1),id1]
all_data[,id2_key := which(id_unique == id2),id2]

all_data[,has_keys := !is.na(id1_key) & (grepl("1",type) | !is.na(id2_key))]
ggplot(all_data,aes(fitness,color=type,linetype=has_keys)) +
  geom_density()
all_data[,.N,.(type,has_keys)]

#stability data
s_singles = all_data[has_keys == T & type=="s_F1"]
s_singles_key = s_singles[,id1_key]
s_doubles = all_data[has_keys == T & type=="s_F2"]
s_doubles_key1 = s_doubles[,id1_key]
s_doubles_key2 = s_doubles[,id2_key]

#binding data
b_singles = all_data[has_keys == T & type=="b_F1"]
b_singles_key = b_singles[,id1_key]
b_doubles = all_data[has_keys == T & type=="b_F2"]
b_doubles_key1 = b_doubles[,id1_key]
b_doubles_key2 = b_doubles[,id2_key]

#fitness & error
stability_fitness = c(s_singles[,fitness],s_doubles[,fitness])
stability_error = c(s_singles[,sigma],s_doubles[,sigma])
binding_fitness = c(b_singles[,fitness],b_doubles[,fitness])
binding_error = c(b_singles[,sigma],b_doubles[,sigma])


###
# fit parameters
registerDoMC(cores=Ncores)
global_models = foreach(m=1:25) %dopar% {
  
  global_model = optim(par = c(rep(0,id_L),
                               rep(0,id_L),
                               0.9+0.1*runif(1),
                               0.9+0.1*runif(1),
                               0.1+0.3*runif(2)),
                       fn = function_fit_allvars_bothassays,
                       method = "L-BFGS-B",
                       lower = c(rep(-10,id_L),rep(-10,id_L),1e-2,1e-2,1e-2,1e-2),
                       upper = c(rep(10,id_L),rep(10,id_L),1-1e-2,1-1e-2,1-1e-2,1-1e-2),
                       control = list(trace=3))
  return(global_model)
}

#save
save(global_models,file = "dataset/dG/method2_Xvars_Xmodel.Rdata")

# evaluate parameters versus value of objective function
parameters = t(sapply(X=1:length(global_models),FUN = function(X){global_models[[X]]$par}))
objective = sapply(X=1:length(global_models),FUN = function(X){global_models[[X]]$value})
dt = data.table(log10(objective),parameters[,(2*id_L+1) : (2*id_L+4)])
names(dt) = c("objective","b_scale","s_scale","b_bgr","s_bgr")
ggpairs(dt)
ggsave("figures/dG/method2_parameter_dependencies_top100.pdf")
# best model
best_global_model = global_models[[which(rank(objective)==1)]]
best_global_par = data.table(b_scale = best_global_model$par[2*id_L+1],
                             s_scale = best_global_model$par[2*id_L+2],
                             b_bgr = best_global_model$par[2*id_L+3],
                             s_bgr = best_global_model$par[2*id_L+4])

#calcualte again dG wildtype values from final parameters
s_dGwt = function_folding_F2dG(1,best_global_par$s_bgr,best_global_par$s_scale)
b_dGwt = function_binding_F2dG(1,s_dGwt,best_global_par$b_bgr,best_global_par$b_scale)
best_global_par[,b_dGwt := b_dGwt]
best_global_par[,s_dGwt := s_dGwt]
print(best_global_par)
write.table(x = best_global_par,file = "dataset/dG/method2_top100vars_bestpar.txt",row.names=F,col.names=T,quote=F)

#transfer results to all_data table
ddG_global_pars = data.table(id = id_unique,s_ddG = best_global_model$par[1:id_L],b_ddG = best_global_model$par[(id_L+1):(2*id_L)])
all_data_wide[,s_ddG_global1 := ddG_global_pars[id == id1,s_ddG],id1]
all_data_wide[,s_ddG_global2 := ddG_global_pars[id == id2,s_ddG],id2]
all_data_wide[,b_ddG_global1 := ddG_global_pars[id == id1,b_ddG],id1]
all_data_wide[,b_ddG_global2 := ddG_global_pars[id == id2,b_ddG],id2]

all_data_wide[Nmut==1,s_fitness_global := function_folding_dG2F(s_ddG_global1,best_global_par$s_dGwt,
                                                                best_global_par$s_bgr,best_global_par$s_scale),.(id1,id2)]
all_data_wide[Nmut==2,s_fitness_global := function_folding_dG2F(s_ddG_global1+s_ddG_global2,best_global_par$s_dGwt,
                                                                best_global_par$s_bgr,best_global_par$s_scale),.(id1,id2)]
all_data_wide[Nmut==1,b_fitness_global := function_binding_dG2F(b_ddG_global1,s_ddG_global1,best_global_par$b_dGwt,best_global_par$s_dGwt,
                                                                best_global_par$b_bgr,best_global_par$b_scale),.(id1,id2)]
all_data_wide[Nmut==2,b_fitness_global := function_binding_dG2F(b_ddG_global1+b_ddG_global2,s_ddG_global1+s_ddG_global2,
                                                                best_global_par$b_dGwt,best_global_par$s_dGwt,best_global_par$b_bgr,best_global_par$b_scale),.(id1,id2)]
all_data_wide[Nmut==1,b_fitness_global_b0 := function_binding_dG2F(0,s_ddG_global1,best_global_par$b_dGwt,best_global_par$s_dGwt,
                                                                   best_global_par$b_bgr,best_global_par$b_scale),.(id1,id2)]
all_data_wide[Nmut==2,b_fitness_global_b0 := function_binding_dG2F(0,s_ddG_global1+s_ddG_global2,best_global_par$b_dGwt,best_global_par$s_dGwt,
                                                                   best_global_par$b_bgr,best_global_par$b_scale),.(id1,id2)]


#parameter/fitness scatter
# ggpairs(DT[is.global==T],columns = c("stability","s_fitness_global","s_ddG_global","binding","b_fitness_global","b_ddG_global"),
# columnLabels = c("Fstab","pred. Fstab","ddGstab","Fbind","pred. Fbind","ddGbind"))
ggpairs(all_data_wide[(!is.na(s_ddG_global1) | !is.na(b_ddG_global1)) & Nmut == 1],
        columns = c("s_fitness","s_fitness_global","s_ddG_global1","b_fitness","b_fitness_global","b_fitness_global_b0","b_ddG_global1"),
        aes(color=res_type1,alpha=0.3))
# columnLabels = c("Fstab","pred. Fstab","ddGstab","Fbind","pred. Fbind"))
ggsave("figures/dG/method2_dG_fromglobalparameters_top100.pdf",width=10,height=7)

#global relationship in bindingPCA and stabilityPCA space
ggplot(rbind(all_data_wide[Nmut==1],all_data_wide[Nmut==2][sample(x = .N,200)])) +
  geom_point(aes(b_fitness,s_fitness,color=factor(Nmut))) +
  geom_point(aes(b_fitness_global_b0,s_fitness_global),color="green",size=2) +
  geom_segment(aes(x=b_fitness_global_b0,y=s_fitness_global,xend = b_fitness,yend=s_fitness),alpha=0.3) +
  labs(x="fitness bindingPCA",y="fitness stabilityPCA",color="used in fitting")
ggsave("figures/dG/method2_binding_stability_relationship_top100.pdf",width=6,height=6)

ggplot(all_data_wide,aes(s_dGwt+s_ddG_global1)) + 
  geom_histogram() + 
  geom_vline(xintercept = s_dGwt,color="red")

ggplot(all_data_wide[Nmut==1],aes(b_dGwt + b_ddG_global1,s_dGwt+s_ddG_global1,color=res_type1)) + 
  geom_point() + 
  geom_hline(yintercept = s_dGwt,color="red") +
  geom_vline(xintercept = b_dGwt,color="red")
ggsave("figures/dG/method2_dG_fvsb_global_relationship_top100.pdf",width=6,height=6)








