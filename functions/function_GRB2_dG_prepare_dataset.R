##### preprocess data for dG estimation
function_GRB2_dG_prepare_dataset <- function(
                                             dataset_name = "GRB2",
                                             DMS_file_list = c(
                                                "processed_data/GRB2_wildtype_alldata.txt",
                                                "processed_data/GRB2_singles_alldata.txt",
                                                "processed_data/GRB2_doubles_alldata.txt"
                                              ),
                                             PDB_interaction_file,
                                             read_threshold = 20,
                                             execute = TRUE) {
  if (!execute) {
    return()
  }

  require(data.table)
  require(ggplot2)
  require(gridExtra)

  ##### load new DiMSum data
  wildtype <- fread(DMS_file_list[1])
  singles <- fread(DMS_file_list[2])
  doubles <- fread(DMS_file_list[3])


  #compare s_fitness from NM and PE
  X = rbind(singles[,.(s_fitness, s_fitness_EP, s_fitness_NM, type = "single")],
    doubles[,.(s_fitness, s_fitness_EP, s_fitness_NM, type = "double")])
  ggplot(melt(X,id.vars = "type"), aes(value, color = variable, linetype = type)) +
    geom_density() 
  ggsave("results/preprocessing/GRB2_sfitness_distribution.pdf")

  #
  doubles[,s_fitness1 := singles[Pos %in% Pos1 & Mut %in% Mut1, s_fitness],.(Pos1, Mut1)]
  doubles[,s_fitness2 := singles[Pos %in% Pos2 & Mut %in% Mut2, s_fitness],.(Pos2, Mut2)]
  ggplot(doubles, aes(s_fitness1 * s_fitness2, s_fitness)) +
   geom_hex() +
   geom_smooth() +
   geom_abline(color = "red") +
   # scale_x_log10() +
   scale_fill_distiller(direction = 1, trans = "log10")
  ggsave("results/preprocessing/GRB2_sfitness_doubleshex.pdf")

  doubles[,s_fitness1_EP := singles[Pos %in% Pos1 & Mut %in% Mut1, s_fitness_EP],.(Pos1, Mut1)]
  doubles[,s_fitness2_EP := singles[Pos %in% Pos2 & Mut %in% Mut2, s_fitness_EP],.(Pos2, Mut2)]
  ggplot(doubles, aes(s_fitness1_EP * s_fitness2_EP, s_fitness_EP)) +
   geom_hex() +
   geom_smooth() +
   geom_abline(color = "red") +
   # scale_x_log10() +
   scale_fill_distiller(direction = 1, trans = "log10")
  ggsave("results/preprocessing/GRB2_sfitness_doubleshex_EP.pdf")

  doubles[,s_fitness1_NM := singles[Pos %in% Pos1 & Mut %in% Mut1, s_fitness_NM],.(Pos1, Mut1)]
  doubles[,s_fitness2_NM := singles[Pos %in% Pos2 & Mut %in% Mut2, s_fitness_NM],.(Pos2, Mut2)]
  ggplot(doubles, aes(s_fitness1_NM * s_fitness2_NM, s_fitness_NM)) +
   geom_hex() +
   geom_smooth() +
   geom_abline(color = "red") +
   # scale_x_log10() +
   scale_fill_distiller(direction = 1, trans = "log10")
  ggsave("results/preprocessing/GRB2_sfitness_doubleshex_NM.pdf")

  doubles[,b_fitness1 := singles[Pos %in% Pos1 & Mut %in% Mut1, b_fitness],.(Pos1, Mut1)]
  doubles[,b_fitness2 := singles[Pos %in% Pos2 & Mut %in% Mut2, b_fitness],.(Pos2, Mut2)]
  ggplot(doubles, aes(b_fitness1 * b_fitness2, b_fitness)) +
   geom_hex() +
   geom_smooth() +
   geom_abline(color = "red") +
   # scale_x_log10() +
   scale_fill_distiller(direction = 1, trans = "log10")
  ggsave("results/preprocessing/GRB2_bfitness_doubleshex.pdf")


  # combine
  all_data <- rbind(
    wildtype[Mut != "*" & WT_AA != "*", .(
      Nmut = 0,
      Pos1 = Pos,
      Pos2 = NA,
      id1 = paste0(Pos, Mut),
      id2 = NA,
      s_fitness = s_fitness ,
      s_sigma,
      b_fitness = b_fitness ,
      b_sigma
    )],
    singles[Mut != "*" & WT_AA != "*", .(
      Nmut = 1,
      Pos1 = Pos,
      Pos2 = NA,
      id1 = paste0(Pos, Mut),
      id2 = NA,
      s_fitness = s_fitness ,
      s_sigma,
      b_fitness = b_fitness ,
      b_sigma
    )],
    doubles[Mut1 != "*" & WT_AA1 != "*" & Mut2 != "*" & WT_AA2 != "*", .(
      Nmut = 2,
      Pos1,
      Pos2,
      id1 = paste0(Pos1, Mut1),
      id2 = paste0(Pos2, Mut2),
      s_fitness = s_fitness ,
      s_sigma,
      b_fitness = b_fitness ,
      b_sigma
    )]
  )

  all_data[, id := ifelse(Nmut == 2, paste0(id1, "_", id2), id1), .(id1, id2)]

  p1 <- ggplot(all_data[Nmut >= 1], aes(s_fitness, color = factor(Nmut))) +
    geom_density(adjust = 0.25) +
    geom_vline(xintercept = c(0.15, 0.99))
  p2 <- ggplot(all_data[Nmut >= 1], aes(b_fitness, color = factor(Nmut))) +
    geom_density(adjust = 0.25) +
    geom_vline(xintercept = c(0.275, 0.99))
  p <- grid.arrange(p1, p2, nrow = 2)
  ggsave(
    plot = p, "results/preprocessing/GRB2_dG_fitness_dist.pdf",
    width = 6,
    height = 8
  )


  ### how many variants with fitness values?
  print(all_data[
    , .N,
    .(Nmut,
      s_fitness_values = !is.na(s_fitness),
      b_fitness_values = !is.na(b_fitness)
    )
  ])

  

  #define test set
  #10% of double mutants observed in both assays
  set.seed(1603)
  all_data[,.N,.(Nmut,!is.na(s_fitness) & !is.na(b_fitness))]
  all_data[,test_set := Nmut == 2 & !is.na(s_fitness) & !is.na(b_fitness) & runif(1) <= 0.1,id]
  all_data[,.N,test_set]
  
  
  ### how often are variants present in double mutants?
  ## count number of doubles in training data
  tmp = rbind(all_data[Nmut == 2 & !is.na(s_fitness) & test_set == F, .(id = id1)],
    all_data[Nmut == 2 & !is.na(s_fitness) & test_set == F, .(id = id2)])[,.N,id]
  all_data[Nmut == 1 & !is.na(s_fitness), Nsid := tmp[id %in% id1, N], id1] 

  tmp = rbind(all_data[Nmut == 2 & !is.na(b_fitness) & test_set == F, .(id = id1)],
    all_data[Nmut == 2 & !is.na(b_fitness) & test_set == F, .(id = id2)])[,.N,id]
  all_data[Nmut == 1 & !is.na(b_fitness), Nbid := tmp[id %in% id1, N], id1] 

  tmp = rbind(all_data[Nmut == 2 & !is.na(b_fitness) & !is.na(s_fitness) & test_set == F, .(id = id1)],
    all_data[Nmut == 2 & !is.na(b_fitness) & !is.na(s_fitness) & test_set == F, .(id = id2)])[,.N,id]
  all_data[Nmut == 1 & !is.na(b_fitness) & !is.na(s_fitness), Nid := tmp[id %in% id1, N], id1] 

  all_data[Nmut == 1, quantile(Nsid, na.rm = T)]
  all_data[Nmut == 1, quantile(Nbid, na.rm = T)]
  all_data[Nmut == 1, quantile(Nid, na.rm = T)]


  # save
  write.table(all_data,
    file = paste0("processed_data/", dataset_name, "_dG_dataset.txt"),
    row.names = F,
    quote = F
  )


  
  # #### restrict to variants with stability and fitness measured
  # all_data2 <- all_data[!is.na(s_fitness) & !is.na(b_fitness), .(
  #   id, id1, id2, Pos1, Pos2, Nmut,
  #   s_fitness, s_sigma, b_fitness, b_sigma
  # )]

  # print(all_data2[
  #   , .N,
  #   .(Nmut,
  #     s_fitness_values = !is.na(s_fitness),
  #     b_fitness_values = !is.na(b_fitness)
  #   )
  # ])

  # write.table(all_data2,
  #   file = paste0("processed_data/", dataset_name, "_dG_dataset.txt"),
  #   row.names = F,
  #   quote = F
  # )

  #################################################
  ### prepare data for Otwinowski Julia scripts ###
  #################################################

  ## all variants covered in bindingPCA assays
  Odata_allvars <- rbind(
    wildtype[, .(
      ham = 0,
      Mut = "0X",
      DNA = total_reads_input,
      SelAll = total_reads_output
    )],
    singles[!is.na(b_fitness) & !STOP, .(
      ham = 1,
      Mut = paste0(WT_AA, Pos, Mut),
      DNA = total_reads_input,
      SelAll = total_reads_output
    )],
    doubles[!is.na(b_fitness) & !STOP, .(
      ham = 2,
      Mut = paste0(WT_AA1, Pos1, Mut1, "-", WT_AA2, Pos2, Mut2),
      DNA = total_reads_input,
      SelAll = total_reads_output
    )]
  )

  write.table(Odata_allvars,
    file = paste0("processed_data/", dataset_name, "_dG_dataset_Otwinowski.txt"),
    row.names = F,
    quote = F
  )

  ## only variants covered in both assays
  # Odata <- rbind(
  #   wildtype[, .(
  #     ham = 0,
  #     Mut = "0X",
  #     DNA = total_reads_input,
  #     SelAll = total_reads_output
  #   )],
  #   singles[!is.na(s_fitness) & !is.na(b_fitness) & !STOP, .(
  #     ham = 1,
  #     Mut = paste0(WT_AA, Pos, Mut),
  #     DNA = total_reads_input,
  #     SelAll = total_reads_output
  #   )],
  #   doubles[!is.na(s_fitness) & !is.na(b_fitness) & !STOP, .(
  #     ham = 2,
  #     Mut = paste0(WT_AA1, Pos1, Mut1, "-", WT_AA2, Pos2, Mut2),
  #     DNA = total_reads_input,
  #     SelAll = total_reads_output
  #   )]
  # )

  # write.table(Odata,
  #   file = paste0("processed_data/", dataset_name, "_dG_dataset_Otwinowski.txt"),
  #   row.names = F,
  #   quote = F
  # )
}
