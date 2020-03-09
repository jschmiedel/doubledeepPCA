function_PDZ_dG_prepare_dataset = function(
  dataset_name,
  DMS_file_list
) {

  require(data.table)
  require(ggplot2)
  require(GGally)
  require(gridExtra)

  ##### load data
  wildtype_ddPCA <- fread(DMS_file_list[1])
  singles_ddPCA <- fread(DMS_file_list[2])
  McLaughlin <- fread(DMS_file_list[3])

  McLaughlin[, Pos := as.integer(gsub("[A-Z]","",mutant)) - 310, mutant]
  McLaughlin[, WT_AA := strsplit(mutant, "")[[1]][1], mutant]
  McLaughlin[, Mut := strsplit(mutant, "")[[1]][5], mutant]

  # combine
  all_data <- rbind(wildtype_ddPCA[,.(
      Nmut = 0,
      Pos,
      id1 = paste0(Pos, Mut),
      s_fitness = s_fitness + 1,
      s_sigma,
      b1_fitness = b_fitness + 1,
      b1_sigma = b_sigma,
      b2_fitness = 1,
      b2_sigma = 0.01,
      b3_fitness = 1,
      b3_sigma = 0.01
    )],
    merge(
      singles_ddPCA[Mut != "*" & WT_AA != "*" & !is.na(s_fitness) & !is.na(b_fitness), .(
        Nmut = 1,
        Pos,
        id1 = paste0(Pos, Mut),
        s_fitness = s_fitness + 1,
        s_sigma,
        b1_fitness = b_fitness + 1,
        b1_sigma = b_sigma
      )],
      McLaughlin[Mut != "*" & WT_AA != "*", .(
        Nmut = 1,
        Pos,
        id1 = paste0(Pos, Mut),
        b2_fitness = exp(CRIPT),
        b2_sigma = 0.05, ## guess, range is 1.5x larger than ddPCA
        b3_fitness = exp(Tm2F),
        b3_sigma = 0.05 ## guess, range is similar to ddPCA
      )]))

  p = ggpairs(all_data,
    columns = grep("fitness", names(all_data), value = T),
    aes(alpha = 0.01))
  ggsave(plot = p, "results/preprocessing/PDZ_dG_data.pdf")


  p1 = ggplot(all_data, aes(s_fitness)) +
    geom_density() +
    geom_vline(xintercept = c(0.05, 0.95))

  p2 = ggplot(all_data, aes(b1_fitness)) +
    geom_density() +
    geom_vline(xintercept = c(0.1, 0.95))

  p3 = ggplot(all_data, aes(b2_fitness)) +
    geom_density() +
    geom_vline(xintercept = c(0.25, 1.025))

  p4 = ggplot(all_data, aes(b3_fitness)) +
    geom_density() +
    geom_vline(xintercept = c(0.425, 0.95))


  P = grid.arrange(p1, p2, p3, p4, nrow = 2)
  ggsave(plot = P, "results/preprocessing/PDZ_fitness_distributions.pdf")

  #save
  write.table(all_data, 
    file = paste0("processed_data/", dataset_name, "_dG_dataset.txt"),
    row.names = F,
    quote = F)
}