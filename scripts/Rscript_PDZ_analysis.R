### analyse PDZ data

## source functions
filelist = list.files('functions/')
sapply(paste0('functions/',filelist),source,.GlobalEnv)

# create directory structure
dir.create("processed_data/", showWarnings = FALSE)
dir.create("results/", showWarnings = FALSE)
dir.create("results/preprocessing/", showWarnings = FALSE)

require(data.table)
require(ggplot2)
require(GGally)
require(cowplot)
require(gridExtra)
theme_set(theme_classic(base_size=9))

## Colors
col_purple = "#9161A8"
col_blue =  "#0066CC"
col_orange = "#F7941E"
col_red = "#EF4136"


#### load data
singles = fread(file = "processed_data/PDZ_singles_alldata.txt")
# exp fitness values
singles[,s_fitness_exp := exp(s_fitness)]
singles[,b_fitness_exp := exp(b_fitness)]
singles[,diff_binding_stability := b_fitness_exp-s_fitness_exp]

singles = singles[!STOP & !STOP_readthrough]

#### repeat GRB2 analysis
##################################################################
############################ Figure 2 ############################
##################################################################

# detrimental fitness threshold
threshold_fitness = 0.75
fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}
### panel A: fitness distributions and enrichment of core/surface/ligand binding positions in lethal mutations
panel2A_1 = ggplot(singles,aes(s_fitness_exp)) +
  geom_density(adjust=.5, fill="grey95")  +
  geom_vline(xintercept = threshold_fitness,lty=2) +
  coord_cartesian(xlim = range(c(singles$b_fitness_exp,singles$s_fitness_exp),na.rm=T)) +
  scale_x_continuous(breaks = seq(0,1,0.25),expand=c(0,0)) +
  scale_y_continuous(expand = c(0,0.1)) +
  geom_text(inherit.aes = F,data=data.table(x=threshold_fitness,y=1.5,label="detrimental"),aes(x,y,label=label),hjust=1.1) +
  geom_segment(aes(x = threshold_fitness, y = 1.4, xend = threshold_fitness-0.2, yend = 1.4, colour = "segment"),arrow =  arrow(length = unit(0.03, "npc")),color="black") +
  labs(x = "fitness stabilityPCA") +
  ggtitle("stabilityPCA")
panel2A_2 = ggplot(singles,aes(b_fitness_exp)) +
  geom_density(adjust=.5, fill="grey95")  +
  geom_vline(xintercept = threshold_fitness,lty=2) +
  coord_cartesian(xlim = range(c(singles$b_fitness_exp,singles$s_fitness_exp),na.rm=T)) +
  scale_x_continuous(breaks = seq(0,1,0.25),expand=c(0,0)) +
  scale_y_continuous(expand = c(0,0.1),labels = fmt_dcimals(1)) +
  geom_text(inherit.aes = F,data=data.table(x=threshold_fitness,y=2,label="detrimental"),aes(x,y,label=label),hjust=1.1) +
  geom_segment(aes(x = threshold_fitness, y = 1.8, xend = threshold_fitness-0.2, yend = 1.8, colour = "segment"),arrow =  arrow(length = unit(0.03, "npc")),color="black") +
  labs(x = "fitness bindingPCA") +
  ggtitle("bindingPCA")
panel2A_3 = ggplot(singles[,.(lethal = sum(s_fitness_exp < threshold_fitness,na.rm=T)/sum(!is.na(s_fitness_exp)) / singles[,sum(s_fitness_exp < threshold_fitness,na.rm=T)/sum(!is.na(s_fitness_exp))]),type],
                   aes(x=type,y=lethal,fill=type)) +
  geom_bar(stat="identity", width = 0.8, fill="grey95", col="black") +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0.25,0.5,1,1.5,2,3)) +
  theme(legend.position='none') +
  coord_cartesian(ylim = c(0.4,2.25)) +
  labs(x="",y="enrichment detrimental mut.",fill="") +
  ggtitle("")
panel2A_4 = ggplot(singles[,.(lethal = sum(b_fitness_exp < threshold_fitness,na.rm=T)/sum(!is.na(b_fitness_exp)) / singles[,sum(b_fitness_exp < threshold_fitness,na.rm=T)/sum(!is.na(b_fitness_exp))]),type],
                   aes(x=type,y=lethal,fill=type)) +
  geom_bar(stat="identity", width = 0.8, fill="grey95", col="black") +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0.25,0.5,1,1.5,2,3)) +
  theme(legend.position=c(0.5,0.15), legend.direction = "horizontal") +
  coord_cartesian(ylim = c(0.4,2.25)) +
  labs(x="",y="enrichment detrimental mut.",fill="") +
  ggtitle("")

p2A = grid.arrange(grobs = list(panel2A_1,panel2A_2,panel2A_3,panel2A_4),nrow=2,
                   layout_matrix = rbind(c(1,1,1,1,1,3,3),
                                         c(2,2,2,2,2,4,4)))
### panel B: comparison stability and binding fitness colouring by distance to the ligand residues

p2B = ggplot(singles,
             aes(b_fitness_exp, s_fitness_exp,color=log10(HAmin_ligand))) +
  geom_point(size=2) +
  geom_abline(linetype=2) +
  scale_color_gradient2("min dist to ligand [??]", low = col_orange, high = col_purple, mid = "grey90", midpoint = 0.9, breaks = log10(c(5,10,15,20)), labels = c(5,10,15,20)) +
  scale_x_continuous(breaks = seq(0,1.25,0.25),expand=c(0,0.0), limits = c(0.3, 1.25)) +
  scale_y_continuous(breaks = seq(0,1.25,0.25),expand=c(0,0.0),  limits = c(0.3, 1.25)) +
  labs(color="min.dist[A]",x="fitness bindingPCA",y="fitness stabilityPCA") +
  theme(legend.direction = "horizontal", legend.position = c(0.8,0.075))
p2B

### figure 2 assembled
fig2 <- grid.arrange(p2A, p2B, layout_matrix = rbind(c(1,1,1,2,2),
                                                     c(1,1,1,2,2)))
ggsave("results/PDZ_fig2_fitness_distributions.pdf", fig2, width = 9, height = 3.5)


##################################################################
############################ Figure 3 ############################
##################################################################

### panel A: heatmaps of mutation effects
setkey(singles,Mut,Pos)

singles_plot = copy(singles)
singles_plot[,fitness_diff := sign(diff_binding_stability) * min(abs(diff_binding_stability),0.5),.(Pos,Mut)]

panel3A_1 = ggplot(singles,aes(Pos,Mut,fill=s_fitness_exp)) +
  geom_raster() +
  geom_raster(inherit.aes = F, data=unique(singles[,.(Pos,WT_AA)]),aes(Pos,WT_AA),fill="black") +
  geom_raster(inherit.aes = F, data=singles[,.(mean(s_fitness_exp,na.rm=T)),Pos],aes(Pos,y=-0.5,fill=V1)) +
  geom_rect(inherit.aes = F,data=unique(singles[,.(Pos,close=HAmin_ligand <= threshold_ligand_dist)])[order(Pos)][,.(idx = rleid(close),close,Pos)][close==T][,.(minP = min(Pos),maxP = max(Pos)),idx],
            aes(xmin = minP-0.5,xmax=maxP+0.5,ymin=-1,ymax = 20.5),color="black",fill="white",lty=2,alpha=0) +
  scale_x_continuous(expand=c(0,0.02), breaks = seq(5,85,5)) +
  coord_cartesian(ylim=c(-0.5,20)) +
  scale_fill_gradient2(low=col_red,mid="grey80",high=col_blue,na.value="white",midpoint = 1,limits=range(c(singles$b_fitness_exp,singles$s_fitness_exp),na.rm=T),breaks = seq(0,1.25,0.25)) +
  labs(x="",y="aa mutation",fill="stability") +
  ggtitle("stabilityPCA")

panel3A_2 = ggplot(singles,aes(Pos,Mut,fill=b_fitness_exp)) +
  geom_raster() +
  geom_raster(inherit.aes = F, data=unique(singles[,.(Pos,WT_AA)]),aes(Pos,WT_AA),fill="black") +
  geom_raster(inherit.aes = F, data=singles[is.finite(b_fitness_exp),.(mean(b_fitness_exp,na.rm=T)),Pos],aes(Pos,y=-0.5,fill=V1)) +
  geom_rect(inherit.aes = F,data=unique(singles[,.(Pos,close=HAmin_ligand <= threshold_ligand_dist)])[order(Pos)][,.(idx = rleid(close),close,Pos)][close==T][,.(minP = min(Pos),maxP = max(Pos)),idx],
            aes(xmin = minP-0.5,xmax=maxP+0.5,ymin=-1,ymax = 20.5),color="black",fill="white",lty=2,alpha=0) +
  scale_x_continuous(expand=c(0,0.02), breaks = seq(5,85,5)) +
  coord_cartesian(ylim=c(-0.5,20)) +
  scale_fill_gradient2(low=col_red,mid="grey80",high=col_blue,na.value="white",midpoint = 1,limits=range(c(singles$b_fitness_exp,singles$s_fitness_exp),na.rm=T),breaks = seq(0,1.25,0.25)) +
  labs(x="",y="aa mutation",fill="binding") +
  ggtitle("bindingPCA")

panel3A_3 = ggplot(singles_plot,aes(Pos,Mut,fill=fitness_diff)) +
  geom_raster() +
  geom_raster(inherit.aes = F, data=unique(singles[,.(Pos,WT_AA)]),aes(Pos,WT_AA),fill="black") +
  geom_raster(inherit.aes = F, data=singles[is.finite(diff_binding_stability),.(mean(diff_binding_stability,na.rm=T)),Pos],aes(Pos,y=-0.5,fill=V1)) +
  geom_rect(inherit.aes = F,data=unique(singles[,.(Pos,close=HAmin_ligand <= threshold_ligand_dist)])[order(Pos)][,.(idx = rleid(close),close,Pos)][close==T][,.(minP = min(Pos),maxP = max(Pos)),idx],
            aes(xmin = minP-0.5,xmax=maxP+0.5,ymin=-1,ymax = 20.5),color="black",fill="white",lty=2,alpha=0) +
  # geom_text(inherit.aes = F,data=singles[HAmin_ligand < threshold_ligand_dist,.(Pos,GAB2_AA)],aes(x=Pos,y=-2,label=GAB2_AA), size = 2) +
  scale_x_continuous(expand=c(0,0.02), breaks = seq(5,85,5)) +
  coord_cartesian(ylim=c(-2.5,20)) +
  scale_fill_gradient2(low=col_red,mid="grey80",high=col_blue,na.value="white",midpoint = 0,limits=c(-0.5,0.5),breaks = seq(-1,1,0.5)) + # fix this setting upper bound
  labs(x="position",y="aa mutation",fill="bin-sta") +
  ggtitle("Difference")


p3A = grid.arrange(grobs = list(panel3A_1,panel3A_2,panel3A_3),nrow=3)
ggsave(plot=p3A,file = "results/PDZ_fig3_panelA_heatmaps.pdf",height=unit(9,"cm"),width=unit(7,"cm"))


### >>move this to analysis script
ggplot(singles,aes(b_fitness_exp,s_fitness_exp,color=log10(HAmin_ligand))) +
  geom_point(size=2) +
  scale_color_gradient2("min dist to ligand [??]", low = col_orange, high = col_purple, mid = "grey90", midpoint = 0.9, breaks = log10(c(5,10,15,20)), labels = c(5,10,15,20)) +
  geom_abline(intercept = c(-0.4,0.4)) +
  labs(x="fitness b_fitnessPCA",y="fitness s_fitnessPCA")
ggsave("results/PDZ_s_fitness_b_fitness_fitness_scatter.pdf",width=7,height=6)

# set of mutations that increase b_fitnessPCA fitness but has negative s_fitnessPCA fitness
singles[s_fitness_exp-b_fitness_exp < -0.3 & b_fitness_exp > -0.2,.(Pos,WT_AA,Mut,s_fitness_exp,s_sigma,b_fitness_exp,b_sigma,HAmin_ligand,type)][order(HAmin_ligand)]
#proximal sites:
#   all mutations position 12
#distant sites: 
#   D22W
#   G23W
#   G25A/P
#   T75P

# potentially allosteric sites
singles[s_fitness_exp-b_fitness_exp > 0.4,.(Pos,WT_AA,Mut,s_fitness_exp,s_sigma,b_fitness_exp,b_sigma,HAmin_ligand,type,RSA_unbound,RSA_bound)][order(HAmin_ligand)]
# positions 26, 57,35,38,80 (all >9.7A from CRIPT)

#average over positions
singles_avgPos=singles[,.(s_fitness = median(s_fitness,na.rm=T),s_fitness_sd = sd(s_fitness,na.rm=T),
                          s_fitness_q10 = quantile(s_fitness,0.25,na.rm=T),s_fitness_q90 = quantile(s_fitness,0.75,na.rm=T),
                          b_fitness = median(b_fitness,na.rm=T),b_fitness_sd = sd(b_fitness,na.rm=T),
                          b_fitness_q10 = quantile(b_fitness,0.25,na.rm=T),b_fitness_q90 = quantile(b_fitness,0.75,na.rm=T),
                          type=unique(type),HAmin_ligand = unique(HAmin_ligand),N= .N,WT_AA=unique(WT_AA)),.(Pos)]
ggplot(singles_avgPos,aes(b_fitness,s_fitness,color=type)) +
  geom_point() +
  geom_errorbar(aes(ymin = s_fitness_q10,ymax=s_fitness_q90)) +
  geom_errorbarh(aes(xmin = b_fitness_q10,xmax=b_fitness_q90)) +
  geom_abline() +
  geom_vline(xintercept = 0,lty=2) +
  geom_hline(yintercept = 0,lty=2)
singles_avgPos[s_fitness > 0 | b_fitness > 0]
#positions where mutations increase stability are first and last positions in the domain
#positions where mutations increase only binding are 11/12 and 67