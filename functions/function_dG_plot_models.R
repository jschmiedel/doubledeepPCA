function_dG_plot_models <- function(
                      models_name,
                      dataset_name,
                      model_dir = "processed_data/",
                      print_dir = "results/dG/",
                      data_file = ".Rdata",
                      which_model = "best"
) {
  require(data.table)
  require(ggplot2)
  require(GGally)
  theme_set(theme_bw(base_size = 9))

  if (data_file == ".Rdata") {
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
    dt <- data.table(objective = log10(objective), global_par)

    #### next plot dG relationship for select model
    if (which_model == "best") {
      sel_model <- which.min(objective)
    } else { # if which_model is integer
      sel_model <- which_model
    }
    global_par_sel <- global_par[sel_model, ]

    table_dg <- data.table(
      variant = names(parameters[sel_model,1:Nvar]),
      s_ddG = parameters[sel_model,1:Nvar],
      b_ddG = parameters[sel_model,(Nvar + 1):(2 * Nvar)]
    )
  } else if (data_file == "data.table") {
    # dt <- fread(file.path(model_dir, paste0(models_name, ".txt")))
    models_dt <- fread(file.path(model_dir, paste0(models_name, ".txt")))
    global_par <- models_dt[, .SD, .SDcols = grep("^[sb]",names(models_dt))]
    objective <- models_dt[, objective]
    dt <- cbind(objective = log10(objective),global_par)

    if (which_model == "best") {
      sel_model <- which.min(objective)
    } else { # if which_model is integer
      sel_model <- which_model
    }
    global_par_sel <- global_par[sel_model, ]

    s_ddg <- data.table(variant = sapply(X = grep("s_ddg", names(models_dt)),
      FUN = function(X) {strsplit(names(models_dt)[X], "_")[[1]][1]}),
      s_ddg = models_dt[sel_model, unlist(.SD), .SDcols = grep("s_ddg", names(models_dt))])
    b_ddg <- data.table(variant = sapply(X = grep("b_ddg", names(models_dt)),
      FUN = function(X) {strsplit(names(models_dt)[X], "_")[[1]][1]}),
      b_ddg = models_dt[sel_model, unlist(.SD), .SDcols = grep("b_ddg", names(models_dt))])
    
    table_dg <- merge(s_ddg, b_ddg, all = T)
  }

  #### first plot relationship between global parameters and objective
  
  p <- ggpairs(dt)
  ggsave(
    plot = p, 
    filename = file.path(print_dir, paste0(models_name, "_globalpar.pdf")), 
    width = 10, 
    height = 10
  )

  


  ## plot global relationship between stability and binding versus fitness of variants
  ## only single mutants
  # all_data <- fread(paste0("processed_data/", dataset_name, "_dG_dataset.txt"))[id %in% id_var]
  # setkey(all_data, id)
  # all_data[names(table_s_dG), s_ddG := table_s_dG]
  # all_data[names(table_b_dG), b_ddG := table_b_dG]
  ## all variants
  all_data <- fread(paste0("processed_data/", dataset_name, "_dG_dataset.txt"))
  all_data[, s_ddg1 := table_dg[variant == id1,s_ddg],id1]
  all_data[, b_ddg1 := table_dg[variant == id1,b_ddg],id1]
  all_data[, s_ddg2 := table_dg[variant == id2,s_ddg],id2]
  all_data[, b_ddg2 := table_dg[variant == id2,b_ddg],id2]
 

  # calculate stabilityPCA fitness from dG values
  all_data[, s_fitness_pred := function_folding_dG2F(
      s_ddg = sum(c(s_ddg1, s_ddg2), na.rm = T),
      s_dgwt = global_par_sel$s_dgwt,
      s_f0 = ifelse("s_f0" %in% names(global_par), 
        global_par_sel$s_f0, 0), 
      s_fwt = ifelse("s_fwt" %in% names(global_par), 
        global_par_sel$s_fwt, 1)
    ), 
    .(s_ddg1, s_ddg2)
  ]
  
  # calculate bindingPCA fitness if b_ddG = 0; i.e. global relationship governed by s_ddG
  all_data[, 
    b_fitness_pred_bdg0 := function_binding_dG2F(
      b_ddg = 0,
      s_ddg = sum(c(s_ddg1, s_ddg2), na.rm = T),
      b_dgwt = global_par_sel$b_dgwt,
      s_dgwt = global_par_sel$s_dgwt,
      b_f0 = global_par_sel$b_f0,
      b_fwt = global_par_sel$b_fwt
    ), 
    .(s_ddg1, s_ddg2)
  ]
  # calculate bindingPCA fitness with estimated b_ddG values
  all_data[, 
    b_fitness_pred := function_binding_dG2F(
      b_ddg = sum(c(b_ddg1, b_ddg2), na.rm = T),
      s_ddg = sum(c(s_ddg1, s_ddg2), na.rm = T),
      b_dgwt = global_par_sel$b_dgwt,
      s_dgwt = global_par_sel$s_dgwt,
      b_f0 = global_par_sel$b_f0,
      b_fwt = global_par_sel$b_fwt
    ),
    .(b_ddg1, s_ddg1, b_ddg2, s_ddg2)
  ]


  ## relationship between parameters and fitness
  # load structural property file to color residues by type
  variant_structuralproperties = fread(paste0("processed_data/",dataset_name,"_variant_structuralproperties.txt"))
  all_data[,res_type1 := variant_structuralproperties[Pos == Pos1,type],Pos1]

  ## singles
  p <- ggpairs(
    all_data[is.na(id2), cbind(.SD,
        s_ddg = s_ddg1,
        b_ddg = b_ddg1,
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
      "s_ddg", 
      "b_ddg", 
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
        s_ddg = s_ddg1 + s_ddg2,
        b_ddg = b_ddg1 + b_ddg2
      ),
      .SDcols = c(
        "s_fitness", 
        "s_fitness_pred", 
        "b_fitness", 
        "b_fitness_pred"
      )
    ],
    columns = c(
      "s_ddg", 
      "b_ddg", 
      "s_fitness", 
      "s_fitness_pred", 
      "b_fitness", 
      "b_fitness_pred"
    ),
    aes(alpha = 0.1), title = which_model,
    lower = list(continuous = 'density')
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
        b_fitness_pred_bdg0, 
        s_fitness_pred), 
        color = "red",
        size = 2
      ) +
      geom_segment(aes(
        x = b_fitness_pred_bdg0, 
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