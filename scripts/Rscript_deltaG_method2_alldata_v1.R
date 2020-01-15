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

#old
# DT = fread("dataset/GRB2_singles_dataset.txt")
#new DiMSum data
f_F1 = fread("dataset/DiMSum_last_version/GRB2_epPCR_GPD_fitness_singles.txt")
b_F1 = fread("dataset/DiMSum_last_version/GRB2_epPCR_CYC_fitness_singles.txt")
f_F2 = fread("dataset/DiMSum_last_version/GRB2_epPCR_GPD_fitness_doubles.txt")
b_F2 = fread("dataset/DiMSum_last_version/GRB2_epPCR_CYC_fitness_doubles.txt")

all_data = rbind(f_F1[Mut!="*" & !is.na(fitness),.(id1 = paste0(Pos,Mut),id2=NA,fitness,sigma,mean_count,type = "f_F1")],
                 b_F1[Mut!="*" & !is.na(fitness),.(id1 = paste0(Pos,Mut),id2=NA,fitness,sigma,mean_count,type = "b_F1")],
                 f_F2[Mut1!="*" & Mut2!="*" & !is.na(fitness_uncorr),.(id1 = paste0(Pos1,Mut1),id2=paste0(Pos2,Mut2),fitness=fitness_uncorr,sigma=sigma_uncorr,mean_count,type = "f_F2")],
                 b_F2[Mut1!="*" & Mut2!="*" & !is.na(fitness_uncorr),.(id1 = paste0(Pos1,Mut1),id2=paste0(Pos2,Mut2),fitness=fitness_uncorr,sigma=sigma_uncorr,mean_count,type = "b_F2")])
all_data[grep("1",type),Nmut:=1]
all_data[grep("2",type),Nmut:=2]

ggplot(all_data,aes(mean_count,fitness,color=type)) +
  geom_point() +
  scale_x_log10()

read_threshold = 20
all_data[,.N,mean_count > read_threshold]

all_data = all_data[mean_count > read_threshold]

#rescale fitness such that min(F) > 0 and wildtype ~ 1
all_data[,fitness := 0.15*fitness + 1]
all_data[,sigma := 0.15*sigma]

ggplot(all_data,aes(mean_count,fitness,color=type)) +
  geom_point() +
  scale_x_log10()

#parameters
R= 1.98*10^(-3) # kcal/mol
Temp= 310.15


##### wide table and other features (distance to GAB2 and RSA)
all_data_wide = merge(all_data[grep("^f",type),.(id1,id2,f_fitness=fitness,f_sigma=sigma)],
                      all_data[grep("^b",type),.(id1,id2,b_fitness=fitness,b_sigma=sigma)],
                      by=c("id1","id2"),all=T)
all_data_wide[is.na(id2),Nmut := 1]
all_data_wide[!is.na(id2),Nmut := 2]

all_data_wide[nchar(id1)==2,Pos1 := as.numeric(strsplit(id1,"")[[1]][1]),id1]
all_data_wide[nchar(id2)==2,Pos2 := as.numeric(strsplit(id2,"")[[1]][1]),id2]
all_data_wide[nchar(id1)==3,Pos1 := as.numeric(paste0(strsplit(id1,"")[[1]][1:2],collapse="")),id1]
all_data_wide[nchar(id2)==3,Pos2 := as.numeric(paste0(strsplit(id2,"")[[1]][1:2],collapse="")),id2]

### load distances between GRB2 and ligand
GRB2_GAB2_distances = fread("dataset/PDB_contactmap_2vwf_AB.txt")
GRB2_GAB2_distances[, WTAAPos2 := paste0(WT_AA2, Pos2)]
mindist = GRB2_GAB2_distances[,.(HAmin_ligand = min(HAmin),GAB2_AA = WTAAPos2[which.min(HAmin)]),Pos1]
all_data_wide[,GAB2_mindist1 := mindist[Pos1 %in% unlist(.SD),HAmin_ligand],Pos1,.SDcols = "Pos1"]
all_data_wide[,GAB2_mindist2 := mindist[Pos1 %in% Pos2,HAmin_ligand],Pos2]

### load RSA values
RSA_chainA = fread("dataset/2vwf_A.rsa",skip=8,nrows = 56)
names(RSA_chainA) = c("restype","aa","chain","Pos","RSA_all_abs","RSA_all_rel","RSA_ts_abs","RSA_ts_rel","RSA_mc_abs","RSA_mc_rel","RSA_np_abs","RSA_np_rel","RSA_pol_abs","RSA_pos_rel")
all_data_wide[,RSA1 := RSA_chainA[Pos %in% Pos1,RSA_all_rel],Pos1]
all_data_wide[,RSA2 := RSA_chainA[Pos %in% Pos2,RSA_all_rel],Pos2]
threshold_RSA = 10
threshold_ligand_dist = 5 
all_data_wide[RSA1 <= threshold_RSA,res_type1 := "core"]
all_data_wide[RSA1 > threshold_RSA & GAB2_mindist1 < threshold_ligand_dist,res_type1 := "bind"]
all_data_wide[RSA1 > threshold_RSA & GAB2_mindist1 > threshold_ligand_dist,res_type1 := "surf"]
all_data_wide[,res_type1 := factor(res_type1,levels=c("core","surf","bind"))]
all_data_wide[RSA2 <= threshold_RSA,res_type2 := "core"]
all_data_wide[RSA2 > threshold_RSA & GAB2_mindist2 < threshold_ligand_dist,res_type2 := "bind"]
all_data_wide[RSA2 > threshold_RSA & GAB2_mindist2 > threshold_ligand_dist,res_type2 := "surf"]
all_data_wide[,res_type2 := factor(res_type2,levels=c("core","surf","bind"))]
all_data_wide[Nmut==2,double_res_type := paste0(res_type1,res_type2,collapse="_"),.(res_type1,res_type2)]




function_folding_dG2F = function(f_dG,f_dGwt,f_bgr,f_scale) {
  # ff = (((1-f_bgr)/(1+exp((f_dGwt+f_dG)/R/Temp)) + f_bgr)) / f_scale
  # ff = (((1-f_bgr)/(1+exp((f_dGwt+f_dG)/R/Temp)) + f_bgr)) / (f_bgr + (1-f_bgr)*f_scale)
  ff = (1-f_bgr)/f_scale/(1+exp((f_dGwt+f_dG)/R/Temp)) + f_bgr
}
function_binding_dG2F = function(b_dG,f_dG,b_dGwt,f_dGwt,b_bgr,b_scale) {
  # bf = (((1-b_bgr)/(1+exp((b_dGwt+b_dG)/R/Temp)*(1+exp((f_dGwt+f_dG)/R/Temp))) + b_bgr)) / b_scale
  # bf = (((1-b_bgr)/(1+exp((b_dGwt+b_dG)/R/Temp)*(1+exp((f_dGwt+f_dG)/R/Temp))) + b_bgr)) / (b_bgr + (1-b_bgr)*b_scale)
  bf = (1-b_bgr)/b_scale/(1+exp((b_dGwt+b_dG)/R/Temp)*(1+exp((f_dGwt+f_dG)/R/Temp))) + b_bgr
}
function_folding_F2dG = function(f_fitness,f_bgr,f_scale) {
  # f_dG = R*Temp*log((1-f_bgr)/(f_fitness*f_scale-f_bgr)-1)
  # f_dG = R*Temp*log((1-f_bgr)/(f_fitness*(f_bgr + (1-f_bgr)*f_scale)-f_bgr)-1)
  f_dG = R*Temp*log((1-f_bgr)/f_scale/(f_fitness-f_bgr)-1)
}
function_binding_F2dG = function(b_fitness,f_dG,b_bgr,b_scale) {
  # b_dG = R*Temp*(log((1-b_bgr)/(b_fitness*b_scale-b_bgr)-1) - log(1+exp(f_dG/R/Temp)))
  # b_dG = R*Temp*(log((1-b_bgr)/(b_fitness*(b_bgr + (1-b_bgr)*b_scale)-b_bgr)-1) - log(1+exp(f_dG/R/Temp)))
  b_dG = R*Temp*(log((1-b_bgr)/b_scale/(b_fitness-b_bgr)-1) - log(1+exp(f_dG/R/Temp)))
}
## test these functions
X=data.table(f_dg = seq(-4,4,0.1))
X[,ff:=function_folding_dG2F(f_dg,-0.57,0.2,0.5),f_dg]
# X[,bf:=function_binding_dG2F(0,f_dg,-2,-1.3,0.2,1),f_dg]
ggplot(X,aes(f_dg,ff)) +
  geom_point()
# 
# ggplot(X,aes(f_dg,bf)) +
#   geom_point() 
# 
# ggplot(X,aes(bf,ff)) +
#   geom_point() #+
# scale_y_continuous(limits=c(0,1)) +
#   scale_x_continuous(limits=c(0,1))


#fit global relationship between deltaGs and fitess from binding and stabilityPCA data assuming deltadeltaG binding = 0 for all used variants
function_fit_allvars_bothassays = function(parameters) {
  
  f_dG_par = parameters[1:id_L]
  b_dG_par = parameters[(id_L+1):(2*id_L)]
  if (length(parameters) > (2*id_L)) {
    b_scale = parameters[2*id_L+1]
    f_scale = parameters[2*id_L+2]
    b_bgr = parameters[2*id_L+3]
    f_bgr = parameters[2*id_L+4]
  }
  
  #determine wild-type dG from set of parameters; this is for wild-type fitness being 1 !!! change if using fitness simply from log-ratios of counts
  f_dGwt = function_folding_F2dG(1,f_bgr,f_scale)
  b_dGwt = function_binding_F2dG(1,f_dGwt,b_bgr,b_scale)
  
  #mix and match mutants and variants
  # f_dG_for_f = as.matrix(f_var_mut_idx) %*% matrix(f_dG_par)
  f_dG_for_f = c(f_dG_par[f_singles_key],f_dG_par[f_doubles_key1]+f_dG_par[f_doubles_key2])
  ff = function_folding_dG2F(f_dG_for_f,f_dGwt,f_bgr,f_scale)
  
  # f_dG_for_b = as.matrix(b_var_mut_idx) %*% matrix(f_dG_par)
  # b_dG_for_b = as.matrix(b_var_mut_idx) %*% matrix(b_dG_par)
  f_dG_for_b = c(f_dG_par[b_singles_key],f_dG_par[b_doubles_key1]+f_dG_par[b_doubles_key2])
  b_dG_for_b = c(b_dG_par[b_singles_key],b_dG_par[b_doubles_key1]+b_dG_par[b_doubles_key2])
  bf = function_binding_dG2F(b_dG_for_b,f_dG_for_b,b_dGwt,f_dGwt,b_bgr,b_scale)
  
  
  #some bf/ff values might be NA because they are below background growth, set fitness and error NA to efit_global_relationshipclude in deviation calcuation
  binding_fitness2 = matrix(binding_fitness)
  # binding_fitness2[is.na(bf)] = NA
  stability_fitness2 = matrix(stability_fitness)
  # stability_fitness2[is.na(ff)] = NA
  binding_error2 = matrix(binding_error)
  binding_error2[is.na(bf)] = NA
  binding_error2[is.infinite(bf)] = NA
  stability_error2 = matrix(stability_error)
  stability_error2[is.na(ff)] = NA
  stability_error2[is.infinite(ff)] = NA
  
  #mean square deviation; divide by weights to correct for variants being NA
  # MSD = sum((binding_fitness2 - bf)^2 / binding_error2^2 + (stability_fitness2 - ff)^2 / stability_error2^2,na.rm=T) / sum(binding_error2^-2 + stability_error2^-2,na.rm=T)
  MSD = sum((binding_fitness2 - bf)^2 / binding_error2^2,na.rm=T) / sum(binding_error2^-2,na.rm=T) +
    sum((stability_fitness2 - ff)^2 / stability_error2^2,na.rm=T) / sum(stability_error2^-2,na.rm=T)
  # lambda*sum(abs(b_dG),na.rm=T)
  return(MSD)
}

## set up indicies for mutant:variant matching
setkey(all_data,id1,id2)
id_unique = unique(c(all_data_wide[!is.na(f_fitness) & !is.na(b_fitness),id1],all_data_wide[!is.na(f_fitness) & !is.na(b_fitness) & !is.na(id2),id2]))
id_unique = id_unique[!is.na(id_unique)]
# id_unique = names(sort(table(c(all_data$id1,all_data$id2)),decreasing = T)[1:100])
id_L = length(id_unique)

all_data[,id1_key := which(id_unique == id1),id1]
all_data[,id2_key := which(id_unique == id2),id2]

all_data[,has_keys := !is.na(id1_key) & (grepl("1",type) | !is.na(id2_key))]
ggplot(all_data,aes(fitness,color=type,linetype=has_keys)) +
  geom_density()
all_data[,.N,.(type,has_keys)]

f_singles = all_data[has_keys == T & type=="f_F1"]
f_singles_key = f_singles[,id1_key]
f_doubles = all_data[has_keys == T & type=="f_F2"]
f_doubles_key1 = f_doubles[,id1_key]
f_doubles_key2 = f_doubles[,id2_key]

#check whether doubles actualy have weight to change behaviour of singles
# f_singles[,doubles_sigma := f_doubles[id1 %in% .SD | id2 %in% .SD,sum(sigma^-2,na.rm=T)],id1,.SDcols="id1"]
# f_singles[,doubles_N := f_doubles[id1 %in% .SD | id2 %in% .SD,.N],id1,.SDcols="id1"]
# ggplot(f_singles,aes(sigma^-2,doubles_sigma,color=doubles_N)) +
# geom_point() +
# scale_x_log10() +
# scale_y_log10() + geom_abline()
# yes

b_singles = all_data[has_keys == T & type=="b_F1"]
b_singles_key = b_singles[,id1_key]
b_doubles = all_data[has_keys == T & type=="b_F2"]
b_doubles_key1 = b_doubles[,id1_key]
b_doubles_key2 = b_doubles[,id2_key]

# b_DT = all_data[grep("^b",type)]
# b_DT[,rank := 1:.N]
# b_L = nrow(b_DT)
# b_var_mut_idx = sparseMatrix(i = c(b_DT[,rank],b_DT[!is.na(id2_key),rank]),
#                              j = c(b_DT[,id1_key],b_DT[!is.na(id2_key),id2_key]),dims=c(b_L,id_L))

#fitness & error
# stability_fitness = f_DT[,fitness]
# stability_error = f_DT[,sigma]
# binding_fitness = b_DT[,fitness]
# binding_error = b_DT[,sigma]
stability_fitness = c(f_singles[,fitness],f_doubles[,fitness])
stability_error = c(f_singles[,sigma],f_doubles[,sigma])
binding_fitness = c(b_singles[,fitness],b_doubles[,fitness])
binding_error = c(b_singles[,sigma],b_doubles[,sigma])

# use regularization
# lambda = 1e-8

###
# fit parameters
registerDoMC(cores=Ncores)
global_models = foreach(m=1:25) %dopar% {
  #sample global parameters then guess good dGs
  # x = c(0.7 + 0.2*runif(1),0.7+0.2*runif(1),0.3+0.3*runif(2))
  # f_dG_guess = function_folding_F2dG(f_singles[,.(id1_key,fitness)][order(id1_key)][,fitness],
  #                                x[4],x[2])
  # f_dG_guess[is.na(f_dG_guess)] = mean(f_dG_guess[!is.na(f_dG_guess)])
  # b_dG_guess = function_binding_F2dG(b_singles[,.(id1_key,fitness)][order(id1_key)][,fitness],f_dG_guess,x[3],x[1])
  # b_dG_guess[is.na(b_dG_guess)] = mean(b_dG_guess[!is.na(b_dG_guess)])
  
  global_model = optim(par = c(rep(0,id_L),
                               rep(0,id_L),
                               0.9+0.1*runif(1),
                               0.9+0.1*runif(1),
                               0.1+0.3*runif(2)),
                       # model = optim(par = c(f_dG_guess,
                       #                       b_dG_guess,x),
                       fn = function_fit_allvars_bothassays,
                       method = "L-BFGS-B",
                       lower = c(rep(-10,id_L),rep(-10,id_L),1e-2,1e-2,1e-2,1e-2),
                       upper = c(rep(10,id_L),rep(10,id_L),1-1e-2,1-1e-2,1-1e-2,1-1e-2),
                       control = list(trace=3))
  # control = list())
  return(global_model)
}

#save
save(global_models,file = "dataset/dG/method2_Xvars_Xmodel.Rdata")

# evaluate parameters versus value of objective function
parameters = t(sapply(X=1:length(global_models),FUN = function(X){global_models[[X]]$par}))
objective = sapply(X=1:length(global_models),FUN = function(X){global_models[[X]]$value})
dt = data.table(log10(objective),parameters[,(2*id_L+1) : (2*id_L+4)])
names(dt) = c("objective","b_scale","f_scale","b_bgr","f_bgr")
ggpairs(dt)
ggsave("figures/dG/method2_parameter_dependencies_top100.pdf")
# best model
best_global_model = global_models[[which(rank(objective)==1)]]
best_global_par = data.table(b_scale = best_global_model$par[2*id_L+1],
                             f_scale = best_global_model$par[2*id_L+2],
                             b_bgr = best_global_model$par[2*id_L+3],
                             f_bgr = best_global_model$par[2*id_L+4])

#calcualte again dG wildtype values from final parameters
f_dGwt = function_folding_F2dG(1,best_global_par$f_bgr,best_global_par$f_scale)
b_dGwt = function_binding_F2dG(1,f_dGwt,best_global_par$b_bgr,best_global_par$b_scale)
best_global_par[,b_dGwt := b_dGwt]
best_global_par[,f_dGwt := f_dGwt]
print(best_global_par)
write.table(x = best_global_par,file = "dataset/dG/method2_top100vars_bestpar.txt",row.names=F,col.names=T,quote=F)

#transfer results to all_data table
ddG_global_pars = data.table(id = id_unique,f_ddG = best_global_model$par[1:id_L],b_ddG = best_global_model$par[(id_L+1):(2*id_L)])
all_data_wide[,f_ddG_global1 := ddG_global_pars[id == id1,f_ddG],id1]
all_data_wide[,f_ddG_global2 := ddG_global_pars[id == id2,f_ddG],id2]
all_data_wide[,b_ddG_global1 := ddG_global_pars[id == id1,b_ddG],id1]
all_data_wide[,b_ddG_global2 := ddG_global_pars[id == id2,b_ddG],id2]

all_data_wide[Nmut==1,f_fitness_global := function_folding_dG2F(f_ddG_global1,best_global_par$f_dGwt,
                                                                best_global_par$f_bgr,best_global_par$f_scale),.(id1,id2)]
all_data_wide[Nmut==2,f_fitness_global := function_folding_dG2F(f_ddG_global1+f_ddG_global2,best_global_par$f_dGwt,
                                                                best_global_par$f_bgr,best_global_par$f_scale),.(id1,id2)]
all_data_wide[Nmut==1,b_fitness_global := function_binding_dG2F(b_ddG_global1,f_ddG_global1,best_global_par$b_dGwt,best_global_par$f_dGwt,
                                                                best_global_par$b_bgr,best_global_par$b_scale),.(id1,id2)]
all_data_wide[Nmut==2,b_fitness_global := function_binding_dG2F(b_ddG_global1+b_ddG_global2,f_ddG_global1+f_ddG_global2,
                                                                best_global_par$b_dGwt,best_global_par$f_dGwt,best_global_par$b_bgr,best_global_par$b_scale),.(id1,id2)]
all_data_wide[Nmut==1,b_fitness_global_b0 := function_binding_dG2F(0,f_ddG_global1,best_global_par$b_dGwt,best_global_par$f_dGwt,
                                                                   best_global_par$b_bgr,best_global_par$b_scale),.(id1,id2)]
all_data_wide[Nmut==2,b_fitness_global_b0 := function_binding_dG2F(0,f_ddG_global1+f_ddG_global2,best_global_par$b_dGwt,best_global_par$f_dGwt,
                                                                   best_global_par$b_bgr,best_global_par$b_scale),.(id1,id2)]


#parameter/fitness scatter
# ggpairs(DT[is.global==T],columns = c("stability","f_fitness_global","f_ddG_global","binding","b_fitness_global","b_ddG_global"),
# columnLabels = c("Fstab","pred. Fstab","ddGstab","Fbind","pred. Fbind","ddGbind"))
ggpairs(all_data_wide[(!is.na(f_ddG_global1) | !is.na(b_ddG_global1)) & Nmut == 1],
        columns = c("f_fitness","f_fitness_global","f_ddG_global1","b_fitness","b_fitness_global","b_fitness_global_b0","b_ddG_global1"),
        aes(color=res_type1,alpha=0.3))
# columnLabels = c("Fstab","pred. Fstab","ddGstab","Fbind","pred. Fbind"))
ggsave("figures/dG/method2_dG_fromglobalparameters_top100.pdf",width=10,height=7)

#global relationship in bindingPCA and stabilityPCA space
ggplot(rbind(all_data_wide[Nmut==1],all_data_wide[Nmut==2][sample(x = .N,200)])) +
  geom_point(aes(b_fitness,f_fitness,color=factor(Nmut))) +
  geom_point(aes(b_fitness_global_b0,f_fitness_global),color="green",size=2) +
  geom_segment(aes(x=b_fitness_global_b0,y=f_fitness_global,xend = b_fitness,yend=f_fitness),alpha=0.3) +
  labs(x="fitness bindingPCA",y="fitness stabilityPCA",color="used in fitting")
ggsave("figures/dG/method2_binding_stability_relationship_top100.pdf",width=6,height=6)

ggplot(all_data_wide,aes(f_dGwt+f_ddG_global1)) + 
  geom_histogram() + 
  geom_vline(xintercept = f_dGwt,color="red")

ggplot(all_data_wide[Nmut==1],aes(b_dGwt + b_ddG_global1,f_dGwt+f_ddG_global1,color=res_type1)) + 
  geom_point() + 
  geom_hline(yintercept = f_dGwt,color="red") +
  geom_vline(xintercept = b_dGwt,color="red")
ggsave("figures/dG/method2_dG_fvsb_global_relationship_top100.pdf",width=6,height=6)



##########################################################
### now fit b_ddG and f_ddG for all remaining variants ###
## set up indicies for mutant:variant matching
id_unique = unique(c(all_data_wide[!is.na(f_fitness) & !is.na(b_fitness),id1],all_data_wide[!is.na(f_fitness) & !is.na(b_fitness) & !is.na(id2),id2]))
id_unique = id_unique[!is.na(id_unique)]
id_L = length(id_unique)

all_data[,id1_key := which(id_unique == id1),id1]
all_data[,id2_key := which(id_unique == id2),id2]

f_singles = all_data[type=="f_F1"]
f_singles_key = f_singles[,id1_key]
f_doubles = all_data[type=="f_F2"]
f_doubles_key1 = f_doubles[,id1_key]
f_doubles_key2 = f_doubles[,id2_key]

b_singles = all_data[type=="b_F1"]
b_singles_key = b_singles[,id1_key]
b_doubles = all_data[type=="b_F2"]
b_doubles_key1 = b_doubles[,id1_key]
b_doubles_key2 = b_doubles[,id2_key]

stability_fitness = c(f_singles[,fitness],f_doubles[,fitness])
stability_error = c(f_singles[,sigma],f_doubles[,sigma])
binding_fitness = c(b_singles[,fitness],b_doubles[,fitness])
binding_error = c(b_singles[,sigma],b_doubles[,sigma])

###
# fit parameters
# registerDoMC(cores=Ncores)
# local_models = foreach(m=1:25) %dopar% {
b_scale =  best_global_par$b_scale
f_scale =  best_global_par$f_scale
b_bgr = best_global_par$b_bgr
f_bgr = best_global_par$f_bgr
local_model = optim(par = c(rep(0,id_L),
                            rep(0,id_L)),
                    fn = function_fit_allvars_bothassays,
                    method = "L-BFGS-B",
                    lower = c(rep(-10,id_L),rep(-10,id_L)),
                    upper = c(rep(10,id_L),rep(10,id_L)),
                    control = list(trace=3))
# control = list())
# return(local_model)
# }

#transfer results to all_data table
ddG_local_pars = data.table(id = id_unique,f_ddG = local_model$par[1:id_L],b_ddG = local_model$par[(id_L+1):(2*id_L)])
all_data_wide[,f_ddG_local1 := ddG_local_pars[id == id1,f_ddG],id1]
all_data_wide[,f_ddG_local2 := ddG_local_pars[id == id2,f_ddG],id2]
all_data_wide[,b_ddG_local1 := ddG_local_pars[id == id1,b_ddG],id1]
all_data_wide[,b_ddG_local2 := ddG_local_pars[id == id2,b_ddG],id2]

#compare local versus global parameters
#should be the same
ddG_pars_all = merge(ddG_local_pars[,.(id,f_ddG_local=f_ddG,b_ddG_local=b_ddG)],
                     ddG_global_pars[,.(id,f_ddG_global=f_ddG,b_ddG_global=b_ddG)])
ggpairs(ddG_pars_all,columns = grep("ddG",names(ddG_pars_all),value=T))


all_data_wide[Nmut==1,f_fitness_local := function_folding_dG2F(f_ddG_local1,best_global_par$f_dGwt,
                                                               best_global_par$f_bgr,best_global_par$f_scale),.(id1,id2)]
all_data_wide[Nmut==2,f_fitness_local := function_folding_dG2F(f_ddG_local1+f_ddG_local2,best_global_par$f_dGwt,
                                                               best_global_par$f_bgr,best_global_par$f_scale),.(id1,id2)]
all_data_wide[Nmut==1,b_fitness_local := function_binding_dG2F(b_ddG_local1,f_ddG_local1,best_global_par$b_dGwt,best_global_par$f_dGwt,
                                                               best_global_par$b_bgr,best_global_par$b_scale),.(id1,id2)]
all_data_wide[Nmut==2,b_fitness_local := function_binding_dG2F(b_ddG_local1+b_ddG_local2,f_ddG_local1+f_ddG_local2,
                                                               best_global_par$b_dGwt,best_global_par$f_dGwt,best_global_par$b_bgr,best_global_par$b_scale),.(id1,id2)]
all_data_wide[Nmut==1,b_fitness_local_b0 := function_binding_dG2F(0,f_ddG_local1,best_global_par$b_dGwt,best_global_par$f_dGwt,
                                                                  best_global_par$b_bgr,best_global_par$b_scale),.(id1,id2)]
all_data_wide[Nmut==2,b_fitness_local_b0 := function_binding_dG2F(0,f_ddG_local1+f_ddG_local2,best_global_par$b_dGwt,best_global_par$f_dGwt,
                                                                  best_global_par$b_bgr,best_global_par$b_scale),.(id1,id2)]

write.table(x = all_data_wide,file = "dataset/dG/method2_top100vars_bestmodel.txt",quote=F,row.names=F,col.names=T)


ggplot(all_data_wide[Nmut==1]) + 
  geom_point(aes(b_dGwt + b_ddG_local1,f_dGwt+f_ddG_local1,color=res_type1,size=!is.na(f_ddG_global1))) + 
  geom_hline(yintercept = f_dGwt,color="red") +
  geom_vline(xintercept = b_dGwt,color="red")
ggsave("figures/dG/method2_dG_fvsb_local_relationship_top100.pdf",width=6,height=6)

#pleiotropy
DT = copy(all_data_wide[Nmut==1])
DT[,pleiotropy := "none"]
DT[b_ddG_local1 > 0.5 & f_ddG_local1 > 0.5,pleiotropy := "b_ddG > 0.5 & f_ddG > 0.5"]
DT[b_ddG_local1 > 0.5 & f_ddG_local1 > 1,pleiotropy := "b_ddG > 0.5 & f_ddG > 1"]
p1=ggplot(DT) + 
  # geom_density2d(aes(b_ddG_local1,f_ddG_local1,color=res_type1)) + 
  geom_point(aes(b_ddG_local1,f_ddG_local1,color=res_type1,shape=pleiotropy))

p2=ggplot(DT) +
  geom_point(aes(b_fitness,f_fitness,color=pleiotropy)) +
  geom_line(aes(b_fitness_global_b0,f_fitness_global),color="red",size=2) +
  geom_segment(inherit.aes = F,data=DT[pleiotropy!="none"],aes(x=b_fitness_global_b0,y=f_fitness_global,xend = b_fitness,yend=f_fitness),alpha=0.3)
p=grid.arrange(p1,p2,nrow=2)
ggsave(plot=p,file="figures/dG/method2_ddG_pleiotropy.pdf",width=6,height=7)



#   
# 
# function_fit_dGfromglobal = function(parameters) {
#   bf = function_binding_dG2F(parameters[1:L],parameters[(L+1):(L*2)],
#                              best_global_par$b_dGwt,
#                              best_global_par$f_dGwt,
#                              best_global_par$b_bgr,
#                              # best_global_par$b_s1,
#                              best_global_par$b_scale)
#   ff = function_folding_dG2F(parameters[(L+1):(L*2)],
#                              best_global_par$f_dGwt,
#                              best_global_par$f_bgr,
#                              # best_global_par$f_s1,
#                              best_global_par$f_scale)
#   binding_fitness2 = binding_fitness
#   binding_fitness2[is.na(bf)] = NA
#   stability_fitness2 = stability_fitness
#   stability_fitness2[is.na(ff)] = NA
#   binding_error2 = binding_error
#   binding_error2[is.na(bf)] = NA
#   stability_error2 = stability_error
#   stability_error2[is.na(ff)] = NA
#   
#   #mean square deviation; divide by weights to correct for variants being NA
#   MSD = sum((binding_fitness2 - bf)^2 / binding_error2^2 + (stability_fitness2 - ff)^2 / stability_error2^2,na.rm=T) / sum(binding_error2^-2 + stability_error2^-2,na.rm=T)
#   return(MSD)
# }
# 
# # DT2 = copy(GRB2_singles)
# binding_fitness = DT[,binding]
# stability_fitness = DT[,stability]
# binding_error = DT[,sigma_binding]
# stability_error = DT[,sigma_stab]
# L = nrow(DT)
# 
# fit_dGfromglobal = nlminb(start = c(rep(0,2*L)),
#                           # start = c(rnorm(n = L+2),0.3,0.3,0,1,0,1),
#                           objective = function_fit_dGfromglobal,
#                           control = list(iter.max = 500))
# print(fit_dGfromglobal)
# 
# DT[,b_ddG_local := fit_dGfromglobal$par[1:L]]
# DT[,f_ddG_local := fit_dGfromglobal$par[(L+1):(L*2)]]
# DT[,f_fitness_local := function_folding_dG2F(f_ddG_local,f_dGwt,
#                                              best_global_par$f_bgr,best_global_par$f_scale),f_ddG_local]
# DT[,b_fitness_local := function_binding_dG2F(0,f_ddG_local,b_dGwt,f_dGwt,
#                                              best_global_par$b_bgr,best_global_par$b_scale),f_ddG_local]

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





