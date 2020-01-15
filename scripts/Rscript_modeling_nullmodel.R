######## 2-state and 3-state thermodynamics models ########
# Compare the effects of mutations in folding when the phenotype measured is the concentration of the mutated protein 
# or whether is the protein coplex amount.
# Explore different ranges of parameters.
# Find the best fit to our data to use as null model.

#### version: 1.0
#### last modified: 2019/11/11
#### created by: Julia
#### last modified by: Julia

require(ggplot2)
require(data.table)
library(RColorBrewer)
require(princurve)
theme_set(theme_classic(base_size=9))

#### Modeling - Function returns the amount of protein 1 folded (2-state) or the amount of bound complex (3-state, folded and bound)
# P1_total = total amount of protein 1 (folded + unfolded), the one is mutated.
# P2_total = total amount of protein 1 (folded + unfolded). Set to 1, since protein 2 is not mutated. The ratio btw 1 and 2 is what matters.
# P1_deltaG_wt = DeltaG of the WT protein 1. Range of -3~ -1 kcal/mol is reasonable.
# deltG_bind_wt = DeltaG of binding between protein 1 and 2. Range from -5 and -1. Use -2 as default.
# P1_deltadeltaG_folding, range from -1 to 5, since most mutations will increase the deltaG of folding.

myppi<- function(P1_total, P2_total, P1_deltaG_wt, P2_deltaG_wt, P1_deltadeltaG_folding, deltadeltaG_binding, deltG_bind_wt){ 
  # initialize parameters
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  deltG_bind_wt = deltG_bind_wt
  
  deltaG1= P1_deltaG_wt + P1_deltadeltaG_folding
  deltaG2= P2_deltaG_wt 
  
  deltaGb= deltG_bind_wt + deltadeltaG_binding
  
  k1 = exp(-deltaG1/(R*Temp))
  k2= exp(-deltaG2/(R*Temp))
  k3= exp(-deltaGb/(R*Temp))
  
  # Calculating the fraction of P1_folded in a two-state folding model
  sol_2s = k1*(P1_total)/(1 + k1)
  
  # Calculating the fraction of P1_folded when taking into acount the binding to P2 (three-state model)
  a1= k3/((1+k1^-1)* (1+k2^-1))
  b1= -(P1_total+P2_total)*k3/((1+k1^-1)* (1+k2^-1))-1
  c1= k3*P1_total*P2_total/((1+k1^-1)* (1+k2^-1))
  
  comp2= (-b1-(b1^2-4*a1*c1)^0.5)/(2*a1)
  comp1= (-b1+(b1^2-4*a1*c1)^0.5)/(2*a1) 
  
  if ( comp2>0 && comp2< min(P1_total, P2_total) ) {
    sol_3s=comp2
  } else if ( comp1>0 && comp1<min(P1_total, P2_total) ){
    sol_3s=comp1
  } else {
    sol_3s=NA
  }
  
  dt = data.table(amount_P1_folded = sol_2s, amount_complex = sol_3s, P1_deltadeltaG_folding = P1_deltadeltaG_folding)
  dt$amount_P1_folded_rel2wt = dt$amount_P1_folded/dt$amount_P1_folded[dt$P1_deltadeltaG_folding == 0]
  dt$amount_complex_rel2wt = dt$amount_complex/dt$amount_complex[dt$P1_deltadeltaG_folding == 0]
  return(dt)
}



#### Run simulations varying parameters
# default parameters
# deltadeltaG of folding always varies from -10 to 1 kcal/mol (most mutations have a negative effect on folding rather than stabilizing)
ddG_f_range = seq(-1, 10, 0.05)
# amount of P1 (relative to P1 which is 1)
P1_total = 1
# amount of P2 
P2_total = 1
# deltaG of P1
P1_deltaG_wt = -2
#  deltaG of P2
P2_deltaG_wt = -2
# deltaG of binding
deltG_bind_wt = -1



### Varying the P1 total amount relative to P2
P1_total_range = seq(0.8, 1.2, 0.1)

dt1 <- do.call("rbind",lapply(P1_total_range, function(x){
  dt = myppi(P1_total = x, P2_total = P2_total, P1_deltaG_wt = P1_deltaG_wt, P2_deltaG_wt = P2_deltaG_wt, deltG_bind_wt = deltG_bind_wt, 
             P1_deltadeltaG_folding = ddG_f_range, deltadeltaG_binding = rep(0, length(ddG_f_range)) )
  dt$P1_total = x
  dt
}))

p1 <- ggplot(dt1, aes(x=amount_complex_rel2wt, y=amount_P1_folded_rel2wt, color=P1_total, group=P1_total)) + geom_abline(linetype=2) +
  geom_point() + geom_line() + 
  scale_color_gradient2("[Total P1]", low = "tomato2", high ="steelblue2", mid = "khaki2", midpoint = median(P1_total_range)) +
  xlab("[Protein complex] relative to WT") + ylab("[Protein 1] relative to WT") +
  theme(legend.position = c(0.2, 0.75))
p1  


### Varying the P1 WT deltaG
P1_deltaG_wt_range = seq(-3, -1, 0.5)

dt2 <- do.call("rbind",lapply(P1_deltaG_wt_range, function(x){
  dt = myppi(P1_total = P1_total, P2_total = P2_total, P1_deltaG_wt = x, P2_deltaG_wt = P2_deltaG_wt, deltG_bind_wt = deltG_bind_wt, 
             P1_deltadeltaG_folding = ddG_f_range, deltadeltaG_binding = rep(0, length(ddG_f_range)) )
  dt$P1_deltaG_WT = x
  dt
}))

p2 <- ggplot(dt2, aes(x=amount_complex_rel2wt, y=amount_P1_folded_rel2wt, color=P1_deltaG_WT, group=P1_deltaG_WT)) + geom_abline(linetype=2) +
  geom_point() + geom_line() + 
  scale_color_gradient2("deltaG P1 WT", low = "tomato2", high ="steelblue2", mid = "khaki2", midpoint = median(P1_deltaG_wt_range)) +
  xlab("[Protein complex] relative to WT") + ylab("[Protein 1] relative to WT") +
  theme(legend.position = c(0.2, 0.75))
p2 


### Varying the P2 WT deltaG
P2_deltaG_wt_range = seq(-3, -1, 0.5)

dt3 <- do.call("rbind",lapply(P2_deltaG_wt_range, function(x){
  dt = myppi(P1_total = P1_total, P2_total = P2_total, P1_deltaG_wt = P1_deltaG_wt, P2_deltaG_wt = x, deltG_bind_wt = deltG_bind_wt, 
             P1_deltadeltaG_folding = ddG_f_range, deltadeltaG_binding = rep(0, length(ddG_f_range)))
  dt$P2_deltaG_WT = x
  dt
}))

p3 <- ggplot(dt3, aes(x=amount_complex_rel2wt, y=amount_P1_folded_rel2wt, color=P2_deltaG_WT, group=P2_deltaG_WT)) + geom_abline(linetype=2) +
  geom_point() + geom_line() + 
  scale_color_gradient2("deltaG P1 WT", low = "tomato2", high ="steelblue2", mid = "khaki2", midpoint = median(P2_deltaG_wt_range)) +
  xlab("[Protein complex] relative to WT") + ylab("[Protein 1] relative to WT") +
  theme(legend.position = c(0.2, 0.75))
p3 

### Varying the deltaG of binding
deltG_bind_wt_range = seq(-3, -0.25, 0.25)

dt4 <- do.call("rbind",lapply(deltG_bind_wt_range, function(x){
  dt = myppi(P1_total = P1_total, P2_total = P2_total, P1_deltaG_wt = P1_deltaG_wt, P2_deltaG_wt = P2_deltaG_wt, deltG_bind_wt = x, 
             P1_deltadeltaG_folding = ddG_f_range, deltadeltaG_binding = rep(0, length(ddG_f_range)))
  dt$deltG_bind_wt = x
  dt
}))

p4 <- ggplot(dt4, aes(x=amount_complex_rel2wt, y=amount_P1_folded_rel2wt, color=deltG_bind_wt, group=deltG_bind_wt)) + geom_abline(linetype=2) +
  geom_point() + geom_line() + 
  scale_color_gradient2("WT deltG binding", low = "tomato2", high ="steelblue2", mid = "khaki2", midpoint = median(deltG_bind_wt_range)) +
  xlab("[Protein complex] relative to WT") + ylab("[Protein 1] relative to WT") +
  theme(legend.position = c(0.2, 0.75))
p4  



#### Compare to our data from experiments
### load GRB2 single mutation data
# use  source data with read count threshold 20 and non-logged fitness values
singles_binding = fread("dataset/GRB2_CYC_singles_readT20.txt")
singles_stab = fread("dataset/GRB2_GPD_singles_readT20.txt")

# merge both datasets (binding and stability assays)
singles = merge(singles_binding[Mut!='*',.(Pos,WT_AA,Mut, sigma_binding = sigma,binding = fitness)],singles_stab[Mut!='*',.(Pos,WT_AA,Mut,sigma_stab= sigma,stability= fitness)],by=c("Pos","WT_AA","Mut"),all=T)
singles[,diff_binding_stability := (binding-stability)]
singles[,diff_stability_binding := (stability-binding)]
singles[,sigma_diff := sqrt((sigma_binding^2) + (sigma_stab^2))]
singles[,Mut := factor(Mut,levels=strsplit('RHKDESTNQCGPAVILMFYW','')[[1]])]

#few singles with infinite binding values, kick them out
singles[is.infinite(binding),binding := NA]
singles = singles[!is.na(binding)]
singles = singles[!is.na(stability)] 

ttest <- function(av, se, df = 2, mu = 1) {
  tstat <- (av-mu)/se
  # Using T-dist
  pval = 2*pt(abs(tstat), df, lower=FALSE)
  return(pval)
}
singles[,pval_bind := p.adjust(ttest(binding, sigma_binding), method = "fdr")]
singles[,pval_stab := p.adjust(ttest(stability, sigma_stab), method = "fdr")]
singles[,pval_diff := p.adjust(ttest(diff_binding_stability, sigma_diff, mu = 0), method = "fdr")]
singles[,affect_binding_only := pval_bind < 0.05 & pval_stab > 0.05 & diff_binding_stability < 0]

# remove those points that we know that only alter binding from the two assays
singles_subset = singles[affect_binding_only ==F, .(binding, stability)]
X = as.matrix(singles_subset)

ggplot(singles_subset, aes(x=binding, y=stability)) + geom_point() + geom_abline(linetype=2)

#### Vary the 4 free parameters to find the best fit
# iterate over the ranges of all the parameters defined above
# param_scan = setNames(data.table(do.call("rbind", lapply(P1_total_range, function(p1t){
#   do.call("rbind", lapply(P1_deltaG_wt_range, function(p1d){
#     do.call("rbind", lapply(P2_deltaG_wt_range, function(p2d){
#       do.call("rbind", lapply(deltG_bind_wt_range, function(dgb){
#         dt = myppi(P1_total = p1t, P2_total = P2_total, P1_deltaG_wt = p1d, P2_deltaG_wt = p2d, deltG_bind_wt = dgb, 
#                    P1_deltadeltaG_folding = ddG_f_range, deltadeltaG_binding = rep(0, length(ddG_f_range)))
#         fit_pc <- principal_curve(as.matrix(dt[, .(amount_complex_rel2wt, amount_P1_folded_rel2wt)]))
#         ord = fit_pc[[1]][order(fit_pc[[1]][,1], decreasing = TRUE), ]
#         YY = project_to_curve(X, ord, stretch = 0)
#         c(p1t, p1d, p2d, dgb, YY$dist)
#       }))
#     }))
#   }))
# }))), c("P1_total","P1_deltaG_wt", "P2_deltaG_wt", "deltG_bind_wt", "SSD"))
# param_scan[order(param_scan$SSD),]

# set a technical threshold were fitness lower to that is not detectable by the assay
thrs_fitness_stab = 0.3
thrs_fitness_bind = 0.25

param_scan_thrs = setNames(data.table(do.call("rbind", lapply(P1_total_range, function(p1t){
  do.call("rbind", lapply(P1_deltaG_wt_range, function(p1d){
    do.call("rbind", lapply(P2_deltaG_wt_range, function(p2d){
      do.call("rbind", lapply(deltG_bind_wt_range, function(dgb){
        dt = myppi(P1_total = p1t, P2_total = P2_total, P1_deltaG_wt = p1d, P2_deltaG_wt = p2d, deltG_bind_wt = dgb, 
                   P1_deltadeltaG_folding = ddG_f_range, deltadeltaG_binding = rep(0, length(ddG_f_range)))
        dt[amount_complex_rel2wt < thrs_fitness_bind, amount_complex_rel2wt:= thrs_fitness_bind]
        dt[amount_P1_folded_rel2wt < thrs_fitness_stab, amount_P1_folded_rel2wt:= thrs_fitness_stab]
        fit_pc <- principal_curve(as.matrix(dt[, .(amount_complex_rel2wt, amount_P1_folded_rel2wt)]))
        ord = fit_pc[[1]][order(fit_pc[[1]][,1], decreasing = TRUE), ]
        YY = project_to_curve(X, ord, stretch = 0)
        c(p1t, p1d, p2d, dgb, YY$dist)
      }))
    }))
  }))
}))), c("P1_total","P1_deltaG_wt", "P2_deltaG_wt", "deltG_bind_wt", "SSD"))

ordered_param_scan <- param_scan_thrs[order(param_scan_thrs$SSD),]


## Fit the best model with the lowest SSD value
dt_final = myppi(P1_total = as.numeric(ordered_param_scan[1,1]), P2_total = P2_total, P1_deltaG_wt = as.numeric(ordered_param_scan[1,2]), 
                 P2_deltaG_wt = as.numeric(ordered_param_scan[1,3]), deltG_bind_wt = as.numeric(ordered_param_scan[1,4]), 
                 P1_deltadeltaG_folding = ddG_f_range, deltadeltaG_binding = rep(0, length(ddG_f_range)))
dt_final[amount_complex_rel2wt < thrs_fitness_bind, amount_complex_rel2wt := thrs_fitness_bind]
dt_final[amount_P1_folded_rel2wt < thrs_fitness_stab, amount_P1_folded_rel2wt := thrs_fitness_stab]

p5 <- ggplot(singles, aes(x=binding, y=stability, color=affect_binding_only)) + geom_point() +
  geom_line(data=dt_final, aes(x=amount_complex_rel2wt, y=amount_P1_folded_rel2wt), color="red", size=2) 
p5


## For each data point calculate the distance to the projection on the null model expectation
# Propagate the error acordingly and do a t-test to see which points differ from the null model.
fit_pc_final <- as.matrix(dt_final[, .(amount_complex_rel2wt, amount_P1_folded_rel2wt)])
YY_final = project_to_curve(X, fit_pc_final)



plot(fit_pc_final, xlim=c(0,1.2), ylim=c(0,1.2)); points(X); points(YY_final$s, col="red"); whiskers(X, YY_final$s); 
hist(YY_final$dist_ind, breaks = 50)






