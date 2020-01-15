######## calculate free energies of folding and binding from stabilityPCA and deepPCA data ########
#### version: 1.0
#### created: 2019/11/27 by J??rn
#### last modified: 2019/11/27 by J??rn

### method 1: use only single mutant data from both assays

require(data.table)
require(ggplot2)
require(GGally)
require(foreach)
require(doMC)
theme_set(theme_bw(base_size=9))


#julia
setwd("deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/") #This is the only path common to any computer
#joern
setwd("/Users/jschmidel/GoogleDrive/deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/")
## Joern new Macbook
# setwd("/Users/jschmiedel/GoogleDrive/deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/")

#old
# DT = fread("dataset/GRB2_singles_dataset.txt")
#new DiMSum data
DT = fread("dataset/DiMSum_last_version/GRB2_singles_dataset.txt")
#new like old
DT[,binding := 0.24*binding + 1]
DT[,stability := 0.24*stability + 1]
DT[,sigma_binding := 0.24*sigma_binding]
DT[,sigma_stab := 0.24*sigma_stab]


#parameters
R= 1.98*10^(-3) # kcal/mol
Temp= 310.15

#fitness values b_fitness and f_fitness
#deltaG values b_dG and f_dG of each variant
#deltaG values of the wildtype b_dGwt and f_dGwt
#background growth b_bgr and f_bgr
#fitness to state probability transformation b_s1,b_s2,f_s1 and f_s2


function_folding_dG2F = function(f_dG,f_dGwt,f_bgr,f_s2) {
  # ff = (((1-f_bgr)/(1+exp((f_dGwt+f_dG)/R/Temp)) + f_bgr)) / f_s2
  ff = (((1-f_bgr)/(1+exp((f_dGwt+f_dG)/R/Temp)) + f_bgr)) / (f_bgr + (1-f_bgr)*f_s2)
}
function_binding_dG2F = function(b_dG,f_dG,b_dGwt,f_dGwt,b_bgr,b_s2) {
  # bf = (((1-b_bgr)/(1+exp((b_dGwt+b_dG)/R/Temp)*(1+exp((f_dGwt+f_dG)/R/Temp))) + b_bgr)) / b_s2
  bf = (((1-b_bgr)/(1+exp((b_dGwt+b_dG)/R/Temp)*(1+exp((f_dGwt+f_dG)/R/Temp))) + b_bgr)) / (b_bgr + (1-b_bgr)*b_s2)
}
function_folding_F2dG = function(f_fitness,f_bgr,f_s2) {
  # f_dG = R*Temp*log((1-f_bgr)/(f_fitness*f_s2-f_bgr)-1)
  f_dG = R*Temp*log((1-f_bgr)/(f_fitness*(f_bgr + (1-f_bgr)*f_s2)-f_bgr)-1)
}
function_binding_F2dG = function(b_fitness,f_dG,b_bgr,b_s2) {
  # b_dG = R*Temp*(log((1-b_bgr)/(b_fitness*b_s2-b_bgr)-1) - log(1+exp(f_dG/R/Temp)))
  b_dG = R*Temp*(log((1-b_bgr)/(b_fitness*(b_bgr + (1-b_bgr)*b_s2)-b_bgr)-1) - log(1+exp(f_dG/R/Temp)))
}
## test these functions
# X=data.table(f_dg = seq(-10,10,0.1))
# X[,ff:=function_folding_dG2F(f_dg,-1.3,0.2,0,1),f_dg]
# X[,bf:=function_binding_dG2F(0,f_dg,0,-1.3,0.2,0,1),f_dg]
# ggplot(X,aes(f_dg,ff)) + 
# geom_point() + 
# scale_y_continuous(limits=c(0,1))

# ggplot(X,aes(f_dg,bf)) + 
#   geom_point() + 
#   scale_y_continuous(limits=c(0,1))

# ggplot(X,aes(bf,ff)) + 
# geom_point() + 
# scale_y_continuous(limits=c(0,1)) + 
# scale_x_continuous(limits=c(0,1))


#fit global relationship between deltaGs and fitess from binding and stabilityPCA data assuming deltadeltaG binding = 0 for all used variants
function_fit_global_relationship = function(parameters) {
  
  f_dG = parameters[1:L]
  b_bgr = parameters[L+1]
  f_bgr = parameters[L+2]
  # b_s1 = parameters[L+3]
  b_s2 = parameters[L+3]
  # f_s1 = parameters[L+5]
  f_s2 = parameters[L+4]
  b_dG = parameters[(L+5):(2*L+4)]
  
  #determine wild-type dG from set of parameters; this is for wild-type fitness being 1 !!! change if using fitness simply from log-ratios of counts
  f_dGwt = function_folding_F2dG(1,f_bgr,f_s2)
  b_dGwt = function_binding_F2dG(1,f_dGwt,b_bgr,b_s2)
  #for new data
  # f_dGwt = function_folding_F2dG(0,f_bgr,f_s1,f_s2)
  # b_dGwt = function_binding_F2dG(0,f_dGwt,b_bgr,b_s1,b_s2)
  
  
  bf = function_binding_dG2F(b_dG=0,f_dG,b_dGwt,f_dGwt,b_bgr,b_s2)
  # bf = function_binding_dG2F(b_dG,f_dG,b_dGwt,f_dGwt,b_bgr,b_s2)
  ff = function_folding_dG2F(f_dG,f_dGwt,f_bgr,f_s2)
  
  #some bf/ff values might be NA because they are below background growth, set fitness and error NA to efit_global_relationshipclude in deviation calcuation
  binding_fitness2 = binding_fitness
  # binding_fitness2[is.na(bf)] = NA
  stability_fitness2 = stability_fitness
  # stability_fitness2[is.na(ff)] = NA
  binding_error2 = binding_error
  binding_error2[is.na(bf)] = NA
  binding_error2[is.infinite(bf)] = NA
  stability_error2 = stability_error
  stability_error2[is.na(ff)] = NA
  stability_error2[is.infinite(ff)] = NA
  
  #mean square deviation; divide by weights to correct for variants being NA
  MSD = sum((binding_fitness2 - bf)^2 / binding_error2^2 + (stability_fitness2 - ff)^2 / stability_error2^2,na.rm=T) / sum(binding_error2^-2 + stability_error2^-2,na.rm=T)
  # MSD = sum((binding_fitness2 - bf)^2 / binding_error2^2,na.rm=T) / sum(binding_error2^-2,na.rm=T) +
    # sum((stability_fitness2 - ff)^2 / stability_error2^2,na.rm=T) / sum(stability_error2^-2,na.rm=T) +
    # lambda*sum(abs(b_dG),na.rm=T)
  return(MSD)
}

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
# lambda = 1e-8
#old
# fit_global_relationship = nlminb(start = c(rep(0,L),0.3,0.3,0,0.9,0,0.9),
#                                  objective = function_fit_global_relationship,
#                                  lower = c(rep(-Inf,L),0.2,0.2,0,0,0,0),
#                                  upper = c(rep(Inf,L),0.4,0.4,1,1,1,1),
#                                  control = list(iter.max = 200))
#new DiMSum data
models = list()
objective = c()
for (i in 1:10) {
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
  fit_global_relationship = optim(par = c(rep(0,L),
                                          0.1+0.3*runif(2),
                                          0.9+0.1*runif(1),
                                          0.9+0.1*runif(1)),
                                  fn = function_fit_global_relationship,
                                  method = "L-BFGS-B",
                                  lower = c(rep(-10,L),1e-5,1e-5,1e-5,1e-5),
                                  upper = c(rep(10,L),1-1e-5,1-1e-5,1-1e-5,1-1e-5),
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
  objective[i] = fit_global_relationship$value
  models[[i]] = fit_global_relationship
  
}
hist(log10(objective))

all_models = t(sapply(X=1:length(models),FUN = function(X){models[[X]]$par}))
dt = data.table(log10(objective),all_models[,(L+1) : (L+4)])
names(dt) = c("objective","b_bgr","f_bgr","b_s2","f_s2")
ggpairs(dt)

fit_global_relationship = models[[which(rank(objective)==1)]]
table_global_parameters = data.table(b_bgr = fit_global_relationship$par[L+1],
                                     f_bgr = fit_global_relationship$par[L+2],
                                     # b_s1 = fit_global_relationship$par[L+3],
                                     b_s2 = fit_global_relationship$par[L+3],
                                     # f_s1 = fit_global_relationship$par[L+5],
                                     f_s2 = fit_global_relationship$par[L+4])

#calcualte again dG wildtype values from final parameters
#old
f_dGwt = function_folding_F2dG(1,table_global_parameters$f_bgr,table_global_parameters$f_s2)
b_dGwt = function_binding_F2dG(1,f_dGwt,table_global_parameters$b_bgr,table_global_parameters$b_s2)
#new DiMSum data
# f_dGwt = function_folding_F2dG(0,table_global_parameters$f_bgr,table_global_parameters$f_s1,table_global_parameters$f_s2)
# b_dGwt = function_binding_F2dG(0,f_dGwt,table_global_parameters$b_bgr,table_global_parameters$b_s1,table_global_parameters$b_s2)

table_global_parameters[,b_dGwt := b_dGwt]
table_global_parameters[,f_dGwt := f_dGwt]
print(table_global_parameters)

#plot results
DT[is.global==T,f_ddG_global := fit_global_relationship$par[1:L]]
# DT[is.global==T,b_ddG_global := fit_global_relationship$par[(L+5):(2*L+4)]]
DT[is.global==T,f_fitness_global := function_folding_dG2F(f_ddG_global,f_dGwt,
                                                          table_global_parameters$f_bgr,table_global_parameters$f_s2),f_ddG_global]
DT[is.global==T,b_fitness_global := function_binding_dG2F(0,f_ddG_global,b_dGwt,f_dGwt,
                                                          table_global_parameters$b_bgr,table_global_parameters$b_s2),f_ddG_global]

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

ggplot(DT,aes(f_dGwt+f_ddG_global)) + geom_histogram() + geom_vline(xintercept = f_dGwt,color="red")


## now fit b_ddG and f_ddG
function_fit_dGfromglobal = function(parameters) {
  bf = function_binding_dG2F(parameters[1:L],parameters[(L+1):(L*2)],
                             table_global_parameters$b_dGwt,
                             table_global_parameters$f_dGwt,
                             table_global_parameters$b_bgr,
                             # table_global_parameters$b_s1,
                             table_global_parameters$b_s2)
  ff = function_folding_dG2F(parameters[(L+1):(L*2)],
                             table_global_parameters$f_dGwt,
                             table_global_parameters$f_bgr,
                             # table_global_parameters$f_s1,
                             table_global_parameters$f_s2)
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
                                             table_global_parameters$f_bgr,table_global_parameters$f_s2),f_ddG_local]
DT[,b_fitness_local := function_binding_dG2F(0,f_ddG_local,b_dGwt,f_dGwt,
                                             table_global_parameters$b_bgr,table_global_parameters$b_s2),f_ddG_local]

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





