### Infered GRB2 deltaG of binging against in vitro measurements ###
# Use in vitro data from Malagrino et al. Scientific Reports, 2019

require(data.table)
library(ggplot2)
require(gdata)
require(gridExtra)

setwd("~/Google Drive/PhD/Projects/deepStructure/005-Manuscripts/001-ProteinLigandInterfaces_v3/scripts/")

### Load data 
# deltaGs inferred from Joer's models
dG <- fread("../dataset/dG/method2_top100vars_bestmodel.txt")
dG_singles <- dG[Nmut == 1, ]

R= 1.98*10^(-3) # kcal/mol
Temp= 310.15

# Malagrino data
Malagrino <- read.xls("../../../002-ExternalData/Malagrino2019/Malagrino2019_Kds_deltaGs_GAB2wt.xlsx", sheet = 1, header = TRUE)[, 1:7]
Malagrino = Malagrino[!(Malagrino$id == "WT"),]
Malagrino$id1 <- do.call("c",lapply(as.character(Malagrino$id), FUN = function(x) {substr(x, 2, nchar(x))}))
Malagrino$deltaG_malagrino = log(Malagrino$kd)*R*Temp

DS <- merge(dG_singles, Malagrino, by ="id1")


r1 <- round(cor(DS$b_ddG_local1, DS$deltaG_malagrino),2)
pval1 <- round(cor.test(DS$b_ddG_local1, DS$deltaG_malagrino)$p.value, 2)
p1 <- ggplot(DS, aes(x=b_ddG_local1, y=deltaG_malagrino)) + theme_classic() +
  geom_point() +
  geom_text(aes(label=id, x=b_ddG_local1 +0.01 , y= deltaG_malagrino + 0.02)) +
  geom_abline(intercept = 0, slope = -1, linetype=2) +
  xlab("inferred deltaG binding (local)") + ylab("Malagrino et al. deltaG") +
  ggtitle(paste0("R = ", r1, ", pval = ", pval1))
p1


r2 <- round(cor(DS$b_fitness, DS$deltaG_malagrino),2)
pval2 <- round(cor.test(DS$b_fitness, DS$deltaG_malagrino)$p.value, 2)
p2 <- ggplot(DS, aes(x=b_fitness, y=deltaG_malagrino)) + theme_classic() +
  geom_point() +
  geom_text(aes(label=id, x=b_fitness +0.001 , y= deltaG_malagrino + 0.02)) +
  xlab("fitness (bindingPCA)") + ylab("Malagrino et al. deltaG") +
  ggtitle(paste0("R = ", r2, ", pval = ", pval2))
p2

p <- grid.arrange(p1,p2, nrow=1)
ggsave("../figures/dG/method1_newdata/20191218_ComparisonInferredDeltaGs_Malagrino2019.pdf", p, width = 10, height = 4.5)
