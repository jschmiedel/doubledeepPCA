PDZ_EC = fread("/Users/jschmidel/Dropbox (Personal)/Science/CRG/DMS2struct/DMS2struct_analysis/PDZ/processed_data/PWI_EC.txt")

dataset_dir = "/Users/jschmidel/Dropbox (Personal)/Science/CRG/DMS2struct/DMS2struct_analysis/PDZ/"
protein_seq = scan(paste0(dataset_dir,"dataset/PDZ.fasta"),what="character")[2]
protein_seq_split = strsplit(protein_seq,"")[[1]]
PDZ_EC[,WT_AA1 := protein_seq_split[Pos1],Pos1]
PDZ_EC[,WT_AA2 := protein_seq_split[Pos2],Pos2]

PDZ_singles = fread("processed_data/PDZ_singles_alldata.txt")
unique(PDZ_singles[b_fitness - s_fitness > 0.4 & HAmin_ligand > 10 & !STOP_readthrough,Pos])

ggplot(PDZ_EC,aes(Pos1,Pos2,fill=pLL)) + 
  geom_raster() +
  scale_fill_gradient2(midpoint=0) +
  geom_rect(inherit.aes = F,data=unique(PDZ_singles[,.(Pos,close=HAmin_ligand <= threshold_ligand_dist)])[order(Pos)][,.(idx = rleid(close),close,Pos)][close==T][,.(minP = min(Pos),maxP = max(Pos)),idx],
          aes(xmin = minP-0.5,xmax=maxP+0.5,ymin= 0,ymax = 84),color="red",fill="white",lty=2,alpha=0) +
  geom_rect(inherit.aes = F,data=unique(PDZ_singles[STOP_readthrough == F,.(maxF = max(b_fitness-s_fitness,na.rm=T),HAmin_ligand=unique(HAmin_ligand)),Pos][,.(Pos,dist_improv=maxF > 0.4 & HAmin_ligand > 10)])[order(Pos)][,.(idx = rleid(dist_improv),dist_improv,Pos)][dist_improv==T][,.(minP = min(Pos),maxP = max(Pos)),idx],
            aes(ymin = minP-0.5,ymax=maxP+0.5,xmin= 0,xmax = 84),color="green",fill="white",lty=2,alpha=0) +
  geom_rect(inherit.aes = F,data=unique(PDZ_singles[STOP_readthrough == F,.(minF = min(b_fitness-s_fitness,na.rm=T),HAmin_ligand=unique(HAmin_ligand)),Pos][,.(Pos,dist_improv=minF < -0.4 & HAmin_ligand > 10)])[order(Pos)][,.(idx = rleid(dist_improv),dist_improv,Pos)][dist_improv==T][,.(minP = min(Pos),maxP = max(Pos)),idx],
            aes(ymin = minP-0.5,ymax=maxP+0.5,xmin= 0,xmax = 84),color="purple",fill="white",lty=2,alpha=0)
  

binding_site = unique(PDZ_singles[HAmin_ligand <= 5,Pos])
distal_improving = unique(PDZ_singles[STOP_readthrough == F & b_fitness-s_fitness > 0.4 & HAmin_ligand > 10 & Mut != "*",Pos])
distal_worsening = unique(PDZ_singles[STOP_readthrough == F & b_fitness-s_fitness < -0.4 & HAmin_ligand > 10 & Mut != "*",Pos])

X = PDZ_EC[Pos1 %in% binding_site | Pos2 %in% binding_site]
X[,min_dist := min(PDZ_singles[Pos == Pos1,HAmin_ligand],PDZ_singles[Pos == Pos2,HAmin_ligand]),.(Pos1,Pos2)]
X[,max_dist := max(PDZ_singles[Pos == Pos1,HAmin_ligand],PDZ_singles[Pos == Pos2,HAmin_ligand]),.(Pos1,Pos2)]
X[,Pos1_dist := PDZ_singles[Pos == Pos1,HAmin_ligand],Pos1]
X[,Pos2_dist := PDZ_singles[Pos == Pos2,HAmin_ligand],Pos2]
X[,type := "neutral"]
X[Pos1 %in% distal_improving | Pos2 %in% distal_improving,type := "better"]
X[Pos1 %in% distal_worsening | Pos2 %in% distal_worsening,type := "worse"]

ggplot(X[abs(Pos1-Pos2)>5 ],aes(max_dist,pLL,color=type)) + 
  geom_point() +
  geom_smooth(method="loess") +
  geom_hline(yintercept = 0.5)
X[abs(Pos1-Pos2)>5 & type == "better" & max_dist > 10 & pLL > 0.75]
#pair 53-70 has strong pLL
X[(Pos1-Pos2)< -5 & type == "better" & max_dist > 10 & pLL > 0.45][order(pLL)]
#positions 22,25,34 and 53(+54) all connected to binding sites positions via high pLL values

ggplot(X[abs(Pos1-Pos2)>5],aes(max_dist,nmi,color=type)) + 
  geom_point() +
  geom_smooth(method="loess")
X[abs(Pos1-Pos2)>5 & type == "worse" & max_dist > 10 & nmi > 0.1]
#pair 4/5 versus 12 has strong MI

X[(Pos1-Pos2)< -5 & type == "better" & max_dist > 10 & nmi > 0.025]
# positions 22,23&25 show up
PDZ_singles[Pos==2]

X[(Pos1-Pos2)< -5 & max_dist > 10 & pLL > 0.45][order(Pos1)]


Y = PDZ_singles[STOP==F & STOP_readthrough ==F,.(q25 = quantile(b_fitness-s_fitness,0.1,na.rm=T),q75 = quantile(b_fitness-s_fitness,0.9,na.rm=T),HAmin_ligand=unique(HAmin_ligand)),Pos]
ggplot(Y,aes(HAmin_ligand,q25,color=Pos)) + geom_point() + geom_smooth()

Z = merge(PDZ_singles[STOP!="*" & STOP_readthrough == F,.(Pos,WT_AA,Mut,s_fitness,b_fitness,diff_bs = b_fitness - s_fitness,type,HAmin_ligand)],
          # X[!(Pos1 %in% binding_site) & abs(Pos1-Pos2) > 2,.(max_pLL = max(pLL),sum_pLL = sum(pLL)),Pos1],by.x="Pos",by.y="Pos1")
          X[,.(max_pLL = max(nmi),sum_pLL = sum(nmi)),Pos1],by.x="Pos",by.y="Pos1")
Z[,diff_bs_mean := mean(diff_bs,na.rm=T),Pos]
Z[,diff_bs_sd := sd(diff_bs,na.rm=T),Pos]
Z[,diff_bs_mean_abs := mean(abs(diff_bs),na.rm=T),Pos]
ggplot(Z) +
  geom_point(aes(max_pLL,diff_bs)) +
  geom_pointrange(aes(x=max_pLL,y=diff_bs_mean,ymin=diff_bs_mean-diff_bs_sd,ymax=diff_bs_mean+diff_bs_sd),color="red") +
  geom_label(inherit.aes=F,data=unique(Z[,.(max_pLL,diff_bs_mean,Pos)]),aes(max_pLL,diff_bs_mean,label=Pos))

ggplot(Z) +
  geom_point(aes(sum_pLL,diff_bs)) +
  geom_pointrange(aes(x=sum_pLL,y=diff_bs_mean,ymin=diff_bs_mean-diff_bs_sd,ymax=diff_bs_mean+diff_bs_sd),color="red") +
  geom_label(inherit.aes=F,data=unique(Z[,.(sum_pLL,diff_bs_mean,Pos)]),aes(sum_pLL,diff_bs_mean,label=Pos))



contactmap = fread("dataset/PDZ/PDB_contactmap_pdb1be9_A.txt")
C = contactmap[Pos1 %in% Z[HAmin_ligand > 0,Pos] & Pos2 %in% Z[HAmin_ligand > 5,Pos] & HAmin < 5 & abs(Pos1 - Pos2) > 0]
C[,max_pLL1 := Z[Pos %in% unlist(.SD),unique(sum_pLL)],Pos1,.SDcols = c("Pos1")]
C[,max_pLL2 := Z[Pos %in% unlist(.SD),unique(sum_pLL)],Pos2,.SDcols = c("Pos2")]
C[,diff_bs_mean_abs1 := Z[Pos %in% unlist(.SD),unique(diff_bs_mean_abs)],Pos1,.SDcols = c("Pos1")]
C[,diff_bs_mean_abs2 := Z[Pos %in% unlist(.SD),unique(diff_bs_mean_abs)],Pos2,.SDcols = c("Pos2")]

ggplot(Z) +
  geom_segment(inherit.aes=F,data=C[!(Pos1 %in% binding_site) & !(Pos2 %in% binding_site)],aes(x=max_pLL1,xend=max_pLL2,y=diff_bs_mean_abs1,yend=diff_bs_mean_abs2)) +
  geom_label(inherit.aes=F,data=unique(Z[HAmin_ligand > 5,.(max_pLL=sum_pLL,diff_bs_mean_abs,diff_bs_mean,Pos,HAmin_ligand)]),
             aes(max_pLL,diff_bs_mean_abs,
                 label=Pos,
                 fill=diff_bs_mean,
                 color=Pos %in% unique(C[Pos1 %in% binding_site & HAmin < 5, Pos2]))) +
  scale_fill_gradient2(midpoint = 0) +
  scale_color_manual(values=c("black","darkgrey")) +
  labs(color="close to BS")



  # geom_smooth()
  # geom_boxplot()
  
