##### preprocess data for dG estimation
function_dG_prepare_dataset = function(
  name = "GRB2_epPCR",
  DMS_file_list = c("dataset/DiMSum_GRB2/GRB2_epPCR_GPD_fitness_singles.txt",
                    "dataset/DiMSum_GRB2/GRB2_epPCR_GPD_fitness_doubles.txt",
                    "dataset/DiMSum_GRB2/GRB2_epPCR_CYC_fitness_singles.txt",
                    "dataset/DiMSum_GRB2/GRB2_epPCR_CYC_fitness_doubles.txt"),
  PDB_interaction_file = "dataset/PDB_contactmap_2vwf_AB.txt",
  RSA_file = "dataset/2vwf_A.rsa",
  read_threshold = 20,
  execute=TRUE
) {

  if (!execute) {return()}

  require(data.table)
  require(ggplot2)
  ##### load new DiMSum data
  s_F1 = fread(DMS_file_list[1])
  s_F2 = fread(DMS_file_list[2])
  b_F1 = fread(DMS_file_list[3])
  b_F2 = fread(DMS_file_list[4])
  
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
    geom_vline(xintercept = read_threshold) +
    scale_x_log10()
  ggsave(paste0("results/dG/",name,"_dG_prepare_dataset_readthresholding.pdf"))
  all_data[,.N,mean_count > read_threshold]
  all_data = all_data[mean_count > read_threshold]
  
  #rescale fitness such that min(F) = 0.1 0 and wildtype = 1
  minF = all_data[,abs(min(fitness,na.rm=T))]
  all_data[,fitness := fitness/(minF/0.9) + 1]
  all_data[,sigma := sigma/(minF/0.9)]
  
  # ggplot(all_data,aes(mean_count,fitness,color=type)) +
  #   geom_point() +
  #   scale_x_log10()
  
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
  interaction_distances = fread(PDB_interaction_file)
  interaction_distances[, WTAAPos2 := paste0(WT_AA2, Pos2)]
  mindist = interaction_distances[,.(HAmin_ligand = min(HAmin),GAB2_AA = WTAAPos2[which.min(HAmin)]),Pos1]
  all_data_wide[,interaction_mindist1 := mindist[Pos1 %in% unlist(.SD),HAmin_ligand],Pos1,.SDcols = "Pos1"]
  all_data_wide[,interaction_mindist2 := mindist[Pos1 %in% Pos2,HAmin_ligand],Pos2]
  
  ##### load RSA values
  RSA_chainA = fread(RSA_file,skip=8,nrows = 56)
  names(RSA_chainA) = c("restype","aa","chain","Pos","RSA_all_abs","RSA_all_rel","RSA_ts_abs","RSA_ts_rel","RSA_mc_abs","RSA_mc_rel","RSA_np_abs","RSA_np_rel","RSA_pol_abs","RSA_pos_rel")
  all_data_wide[,RSA1 := RSA_chainA[Pos %in% Pos1,RSA_all_rel],Pos1]
  all_data_wide[,RSA2 := RSA_chainA[Pos %in% Pos2,RSA_all_rel],Pos2]
  threshold_RSA = 10
  threshold_ligand_dist = 5 
  all_data_wide[RSA1 <= threshold_RSA,res_type1 := "core"]
  all_data_wide[RSA1 > threshold_RSA & interaction_mindist1 < threshold_ligand_dist,res_type1 := "bind"]
  all_data_wide[RSA1 > threshold_RSA & interaction_mindist1 > threshold_ligand_dist,res_type1 := "surf"]
  all_data_wide[,res_type1 := factor(res_type1,levels=c("core","surf","bind"))]
  all_data_wide[RSA2 <= threshold_RSA,res_type2 := "core"]
  all_data_wide[RSA2 > threshold_RSA & interaction_mindist2 < threshold_ligand_dist,res_type2 := "bind"]
  all_data_wide[RSA2 > threshold_RSA & interaction_mindist2 > threshold_ligand_dist,res_type2 := "surf"]
  all_data_wide[,res_type2 := factor(res_type2,levels=c("core","surf","bind"))]
  all_data_wide[Nmut==2,double_res_type := paste0(res_type1,"_",res_type2),.(res_type1,res_type2)]
  
  #### restrict to variants with stability and fitness measured
  all_data_wide2 = all_data_wide[!is.na(s_fitness) & !is.na(b_fitness),.(id,id1,id2,Pos1,Pos2,Nmut,
                                                                         s_fitness,s_sigma,b_fitness,b_sigma,
                                                                         interaction_mindist1,interaction_mindist2,
                                                                         RSA1,RSA2,res_type1,res_type2,double_res_type)]
  
  ## stratify doubles into 10 groups for cross validation
  all_data_wide2[Nmut==2,tenfold_Xval := ceiling(runif(.N)*10)]
  
  write.table(all_data_wide2,file = paste0("processed_data/",name,"_dG_dataset.txt"),row.names=F,quote=F)
}