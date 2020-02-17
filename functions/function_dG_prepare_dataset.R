##### preprocess data for dG estimation
function_dG_prepare_dataset <- function(
                    dataset_name = "GRB2",
                    DMS_file_list,
                    PDB_interaction_file = "dataset/PDB_contactmap_2vwf_AB.txt",
                    RSA_file = "dataset/2vwf_A.rsa",
                    read_threshold = 20,
                    execute = TRUE) {
  if (!execute) {
    return()
  }

  require(data.table)
  require(ggplot2)
  ##### load new DiMSum data
  wildtype <- fread(DMS_file_list[1])
  singles <- fread(DMS_file_list[2])
  doubles <- fread(DMS_file_list[3])

  # combine
  all_data <- rbind(
    singles[Mut != "*" & WT_AA != "*", .(
      Nmut = 1,
      Pos1 = Pos,
      Pos2 = NA,
      id1 = paste0(Pos, Mut),
      id2 = NA,
      s_fitness = s_fitness + 1,
      s_sigma,
      b_fitness = b_fitness + 1,
      b_sigma
    )],
    doubles[Mut1 != "*" & WT_AA1 != "*" & Mut2 != "*" & WT_AA2 != "*", .(
      Nmut = 2,
      Pos1,
      Pos2,
      id1 = paste0(Pos1, Mut1),
      id2 = paste0(Pos2, Mut2),
      s_fitness = s_fitness + 1,
      s_sigma,
      b_fitness = b_fitness + 1,
      b_sigma
    )])

  all_data[, id := ifelse(Nmut == 1, id1, paste0(id1, "_", id2)), .(id1, id2)]

  ### how many variants with fitness values?
  print(all_data[, .N,
                  .(Nmut,
                    s_fitness_values = !is.na(s_fitness),
                    b_fitness_values = !is.na(b_fitness))])

  ### how often are variants present in double mutants?
  X <- rbind(all_data[Nmut == 2 & !is.na(s_fitness) & !is.na(b_fitness),
                     .(id = id1)],
    all_data[Nmut == 2 & !is.na(s_fitness) & !is.na(b_fitness), .(id = id2)])
  X[, .N, id][, quantile(N)]

  #save
  write.table(all_data, 
    file = paste0("processed_data/", dataset_name, "_dG_dataset_allvars.txt"),
    row.names = F,
    quote = F)


  #### restrict to variants with stability and fitness measured
  all_data2 <- all_data[!is.na(s_fitness) & !is.na(b_fitness), .(
    id, id1, id2, Pos1, Pos2, Nmut,
    s_fitness, s_sigma, b_fitness, b_sigma
  )]

  print(all_data2[, .N,
                  .(Nmut,
                    s_fitness_values = !is.na(s_fitness),
                    b_fitness_values = !is.na(b_fitness))])

  write.table(all_data2, 
    file = paste0("processed_data/", dataset_name, "_dG_dataset.txt"),
    row.names = F,
    quote = F)

  #################################################
  ### prepare data for Otwinowski Julia scripts ###
  #################################################

  ## all variants covered in bindingPCA assays
  Odata_allvars = rbind(wildtype[, .(
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
    )])

  write.table(Odata_allvars, 
    file = paste0("processed_data/", dataset_name, "_dG_dataset_allvars_Otwinowski.txt"),
    row.names = F,
    quote = F)

  ## only variants covered in both assays
  Odata = rbind(wildtype[, .(
      ham = 0,
      Mut = "0X",
      DNA = total_reads_input,
      SelAll = total_reads_output
    )],
    singles[!is.na(s_fitness) & !is.na(b_fitness) & !STOP, .(
      ham = 1,
      Mut = paste0(WT_AA,Pos, Mut),
      DNA = total_reads_input,
      SelAll = total_reads_output
    )],
    doubles[!is.na(s_fitness) & !is.na(b_fitness) & !STOP, .(
      ham = 2,
      Mut = paste0(WT_AA1, Pos1, Mut1, "-", WT_AA2, Pos2, Mut2),
      DNA = total_reads_input,
      SelAll = total_reads_output
    )])

  write.table(Odata, 
    file = paste0("processed_data/", dataset_name, "_dG_dataset_Otwinowski.txt"),
    row.names = F,
    quote = F)
}
