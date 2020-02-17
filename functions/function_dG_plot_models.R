function_dG_plot_models <- function(
                      models_name,
                      dataset_name,
                      model_dir = "processed_data/",
                      print_dir = "results/dG/",
                      which_model = "best"
) {
  require(data.table)
  require(ggplot2)
  require(GGally)
  theme_set(theme_bw(base_size = 9))

  # load Rdata file with models
  load(file = file.path(model_dir, paste0(models_name, ".Rdata")))

  ## extract parameters and objective from models
  parameters <- t(sapply(X = 1:length(models), FUN = function(X) {
    models[[X]]$par
  }))
  # global parameters
  global_par_idx <- grep("[sb]_", colnames(parameters))
  global_par <- data.table(parameters[, global_par_idx])
  # which/how many variants?
  Nvar <- (ncol(parameters) - length(global_par_idx)) / 2
  id_var <- names(models[[1]]$par[1:Nvar])
  # optim objective
  objective <- sapply(X = 1:length(models), FUN = function(X) {
    models[[X]]$value
  })

  # add wild-type dG values if missing
  if (!("s_dGwt" %in% names(global_par))) {
    global_par[, s_dGwt := function_folding_F2dG(
      s_fitness = ifelse(names(global_par) %in% s_Fwt, global_par$s_Fwt, 1), 
      s_bgr = s_bgr, 
      s_scale = s_scale
    )]
  }
  if (!("b_dGwt" %in% names(global_par))) {
    global_par[, b_dGwt := function_binding_F2dG(
      b_fitness = ifelse(names(global_par) %in% b_Fwt, global_par$b_Fwt, 1),  
      s_dGwt = s_dGwt, 
      b_bgr = b_bgr, 
      b_scale = b_scale
    )]
  }


  #### first plot relationship between global parameters and objective
  dt <- data.table(objective = log10(objective), global_par)
  p <- ggpairs(dt)
  ggsave(
    plot = p, 
    filename = file.path(print_dir, paste0(models_name, "_globalpar.pdf")), 
    width = 10, 
    height = 10
  )

  #### next plot dG relationship for select model
  if (which_model == "best") {
    sel_model <- which.min(objective)
  } else { # if which_model is integer
    sel_model <- which_model
  }
  global_par_sel <- global_par[sel_model, ]

  table_dg = data.table(
    variant = names(parameters[sel_model,1:Nvar]),
    s_ddG = parameters[sel_model,1:Nvar],
    b_ddG = parameters[sel_model,(Nvar + 1):(2 * Nvar)]
  )


  ## plot global relationship between stability and binding versus fitness of variants
  ## only single mutants
  # all_data <- fread(paste0("processed_data/", dataset_name, "_dG_dataset.txt"))[id %in% id_var]
  # setkey(all_data, id)
  # all_data[names(table_s_dG), s_ddG := table_s_dG]
  # all_data[names(table_b_dG), b_ddG := table_b_dG]
  ## all variants
  all_data <- fread(paste0("processed_data/", dataset_name, "_dG_dataset.txt"))
  all_data[, s_ddG1 := table_dg[variant == id1,s_ddG],id1]
  all_data[, b_ddG1 := table_dg[variant == id1,b_ddG],id1]
  all_data[, s_ddG2 := table_dg[variant == id2,s_ddG],id2]
  all_data[, b_ddG2 := table_dg[variant == id2,b_ddG],id2]
 

  # calculate stabilityPCA fitness from dG values
  all_data[, s_fitness_pred := function_folding_dG2F(
      s_dG = sum(c(s_ddG1, s_ddG2), na.rm = T),
      s_dGwt = global_par_sel$s_dGwt,
      s_bgr = ifelse("s_scale" %in% names(global_par), 
        global_par_sel$s_bgr, 0), 
      s_scale = ifelse("s_scale" %in% names(global_par), 
        global_par_sel$s_scale, 1)
    ), 
    .(s_ddG1, s_ddG2)
  ]
  
  # calculate bindingPCA fitness if b_ddG = 0; i.e. global relationship governed by s_ddG
  all_data[, 
    b_fitness_pred_bdG0 := function_binding_dG2F(
      b_dG = 0,
      s_dG = sum(c(s_ddG1, s_ddG2), na.rm = T),
      b_dGwt = global_par_sel$b_dGwt,
      s_dGwt = global_par_sel$s_dGwt,
      b_bgr = global_par_sel$b_bgr,
      b_scale = global_par_sel$b_scale
    ), 
    .(s_ddG1, s_ddG2)
  ]
  # calculate bindingPCA fitness with estimated b_ddG values
  all_data[, 
    b_fitness_pred := function_binding_dG2F(
      b_dG = sum(c(b_ddG1, b_ddG2), na.rm = T),
      s_dG = sum(c(s_ddG1, s_ddG2), na.rm = T),
      b_dGwt = global_par_sel$b_dGwt,
      s_dGwt = global_par_sel$s_dGwt,
      b_bgr = global_par_sel$b_bgr,
      b_scale = global_par_sel$b_scale
    ),
    .(b_ddG1, s_ddG1, b_ddG2, s_ddG2)
  ]


  ## relationship between parameters and fitness
  # load structural property file to color residues by type
  variant_structuralproperties = fread(paste0("processed_data/",dataset_name,"_variant_structuralproperties.txt"))
  all_data[,res_type1 := variant_structuralproperties[Pos == Pos1,type],Pos1]

  ## singles
  p <- ggpairs(
    all_data[is.na(id2), cbind(.SD,
        s_ddG = s_ddG1,
        b_ddG = b_ddG1,
        res_type1
      ),
      .SDcols = c(
        "s_fitness", 
        "s_fitness_pred", 
        "b_fitness", 
        "b_fitness_pred"
      )
    ],
    columns = c(
      "s_ddG", 
      "b_ddG", 
      "s_fitness", 
      "s_fitness_pred", 
      "b_fitness", 
      "b_fitness_pred"
    ),
    aes(color = res_type1, alpha = 0.3), title = which_model
  )
  ggsave(
    plot = p, 
    file = file.path(print_dir, paste0(models_name, "_parameters_fitness.pdf")), 
    width = 10, 
    height = 10
  )

  ## all doubles
  p <- ggpairs(
    all_data[!is.na(id2),cbind(.SD,
        s_ddG = s_ddG1 + s_ddG2,
        b_ddG = b_ddG1 + b_ddG2
      ),
      .SDcols = c(
        "s_fitness", 
        "s_fitness_pred", 
        "b_fitness", 
        "b_fitness_pred"
      )
    ],
    columns = c(
      "s_ddG", 
      "b_ddG", 
      "s_fitness", 
      "s_fitness_pred", 
      "b_fitness", 
      "b_fitness_pred"
    ),
    aes(alpha = 0.1), title = which_model
  )
  ggsave(
    plot = p, 
    file = file.path(print_dir, paste0(models_name, "_parameters_fitness_alldoubles.pdf")), 
    width = 10, 
    height = 10
  )


  ## global relationship in bindingPCA and stabilityPCA space
  if (sum(grep("s_scale", names(global_par))) > 0) {
    p <- ggplot(all_data[is.na(id2)]) +
      geom_point(aes(
        b_fitness, 
        s_fitness
      )) +
      geom_point(aes(
        b_fitness_pred_bdG0, 
        s_fitness_pred), 
        color = "red",
        size = 2
      ) +
      geom_segment(aes(
        x = b_fitness_pred_bdG0, 
        y = s_fitness_pred, 
        xend = b_fitness, 
        yend = s_fitness), 
        alpha = 0.3
      ) +
      labs(
        x = "fitness bindingPCA", 
        y = "fitness stabilityPCA", 
        color = "used in fitting", 
        title = which_model
      )
    ggsave(
      plot = p, 
      file = file.path(print_dir, paste0(models_name, "_fitness_globalrelationship.pdf")), 
      width = 6, 
      height = 6)
  }
}