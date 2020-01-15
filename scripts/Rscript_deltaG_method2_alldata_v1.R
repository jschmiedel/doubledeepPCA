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

##### load new DiMSum data
s_F1 = fread("dataset/DiMSum_GRB2/GRB2_epPCR_GPD_fitness_singles.txt")
b_F1 = fread("dataset/DiMSum_GRB2/GRB2_epPCR_CYC_fitness_singles.txt")
s_F2 = fread("dataset/DiMSum_GRB2/GRB2_epPCR_GPD_fitness_doubles.txt")
b_F2 = fread("dataset/DiMSum_GRB2/GRB2_epPCR_CYC_fitness_doubles.txt")

#combine
all_data = rbind(s_F1[Mut!="*" & !is.na(fitness) & Pos < 57,.(id1 = paste0(Pos,Mut),id2=NA,fitness,sigma,mean_count,type = "s_F1")],
                 b_F1[Mut!="*" & !is.na(fitness) & Pos < 57,.(id1 = paste0(Pos,Mut),id2=NA,fitness,sigma,mean_count,type = "b_F1")],
                 s_F2[Mut1!="*" & Mut2!="*" & !is.na(fitness_uncorr) & Pos1 < 57  & Pos2 < 57,.(id1 = paste0(Pos1,Mut1),id2=paste0(Pos2,Mut2),fitness=fitness_uncorr,sigma=sigma_uncorr,mean_count,type = "s_F2")],
                 b_F2[Mut1!="*" & Mut2!="*" & !is.na(fitness_uncorr) & Pos1 < 57  & Pos2 < 57,.(id1 = paste0(Pos1,Mut1),id2=paste0(Pos2,Mut2),fitness=fitness_uncorr,sigma=sigma_uncorr,mean_count,type = "b_F2")])
all_data[grep("1",type),Nmut:=1]
all_data[grep("2",type),Nmut:=2]

#apply read threshold
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

##### wide table (fitness and error values for both assays in same row)
all_data_wide = merge(all_data[grep("^s",type),.(id1,id2,s_fitness=fitness,s_sigma=sigma)],
                      all_data[grep("^b",type),.(id1,id2,b_fitness=fitness,b_sigma=sigma)],
                      by=c("id1","id2"),all=T)
all_data_wide[is.na(id2),Nmut := 1]
all_data_wide[!is.na(id2),Nmut := 2]

all_data_wide[nchar(id1)==2,Pos1 := as.numeric(strsplit(id1,"")[[1]][1]),id1]
all_data_wide[nchar(id2)==2,Pos2 := as.numeric(strsplit(id2,"")[[1]][1]),id2]
all_data_wide[nchar(id1)==3,Pos1 := as.numeric(paste0(strsplit(id1,"")[[1]][1:2],collapse="")),id1]
all_data_wide[nchar(id2)==3,Pos2 := as.numeric(paste0(strsplit(id2,"")[[1]][1:2],collapse="")),id2]

all_data_wide[,id := ifelse(Nmut==1,id1,paste0(id1,"_",id2)),.(id1,id2)]

##### load distances between GRB2 and ligand
GRB2_GAB2_distances = fread("dataset/PDB_contactmap_2vwf_AB.txt")
GRB2_GAB2_distances[, WTAAPos2 := paste0(WT_AA2, Pos2)]
mindist = GRB2_GAB2_distances[,.(HAmin_ligand = min(HAmin),GAB2_AA = WTAAPos2[which.min(HAmin)]),Pos1]
all_data_wide[,GAB2_mindist1 := mindist[Pos1 %in% unlist(.SD),HAmin_ligand],Pos1,.SDcols = "Pos1"]
all_data_wide[,GAB2_mindist2 := mindist[Pos1 %in% Pos2,HAmin_ligand],Pos2]

##### load RSA values
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
all_data_wide[Nmut==2,double_res_type := paste0(res_type1,"_",res_type2),.(res_type1,res_type2)]




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








