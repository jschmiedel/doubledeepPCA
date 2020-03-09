## analyse GRB2 dG estimation models

require(data.table)
require(ggplot2)
require(GGally)

files <- list.files("processed_data/tmp/")

files <- files[grep("GRB2", files)]

files_otwin <- files[grep("Otwin", files)]
files <- files[!grepl("Otwin", files)]

files1 <- files[!grepl("allvars", files)]
files2 <- files[grepl("allvars", files)]

for (i in seq_along(files1)) {
    output <- fread(file.path("processed_data/tmp/", files1[i]))
    tmp <- strsplit(files1[i], "_")[[1]]
    output[, method := as.integer(gsub("method", "", tmp[1]))]
    output[, iteration := as.integer(gsub("\\.txt", "", tmp[length(tmp)]))]
    output[, dataset := "dual_observed"]
    if (i == 1) {
        DT <- output
    } else {
        DT <- rbind(DT,
            output,
            fill = TRUE
        )
    }
}
for (i in seq_along(files2)) {
    output <- fread(file.path("processed_data/tmp/", files2[i]))
    tmp <- strsplit(files2[i], "_")[[1]]
    output[, method := as.integer(gsub("method", "", tmp[1]))]
    output[, iteration := as.integer(gsub("\\.txt", "", tmp[length(tmp)]))]
    output[, dataset := "all_vars"]
    DT <- rbind(DT,
        output,
        fill = TRUE
    )
}

DT[, modelnr := .N:1, .(method, dataset, iteration)]

# how many models?
DT[, .N, .(method, dataset, iteration)][, mean(N), .(method, dataset)]
# for method 2 all_vars, there wasn't enough time to finish simulations

DT[, objective := log10(objective)]
DT = DT[is.finite(objective)]

# rank models
setkey(DT, objective)
DT[, rank := 1:.N, .(method, iteration, dataset)]

# plot for one process how models improved
ggpairs(DT[ method == 1 & dataset == "dual_observed" & iteration == 9,
    .SD,
    .SDcols = grep("^[osb]", names(DT), value = T)
    ])
ggsave("results/dG/GRB2_dG_method1_it9.pdf")

# compare iterations for same method/dataset
n <- grep("^[osb]", names(DT), value = T)
# for (i in seq(1, 3)) {
for (i in 1) {
    # for (d in c(1, 2)) {
    d = 1
        if (i == 3) {
            n0 <- n[!grepl("^s_f", n)]
        } else {
            n0 <- n
        }
        p <- ggpairs(DT[method == i & dataset == ifelse(d == 1, "dual_observed", "all_vars") & rank == 1,
            .SD,
            .SDcols = n0
            ])#[between(
                # objective,
            # mean(objective) - sd(objective),
                # -Inf,
                # mean(objective) + sd(objective)
            # )])
        ggsave(
            plot = p,
            filename = paste0("results/dG/GRB2", 
                ifelse(d == 1, "", "_allvars"),"_method", i, "_allit.pdf")
        )
for (j in 1:10) {
        best_model = DT[method == i & dataset == ifelse(d == 1, "dual_observed", "all_vars") & iteration == j][order(objective)][1]
        ddG_lines = grep("ddG", names(best_model), value = T)
        variant = sapply(X = seq_along(ddG_lines), 
            FUN = function(X){strsplit(ddG_lines[X], "_")[[1]][1]})
        type = sapply(X = seq_along(ddG_lines), 
            FUN = function(X){paste0(strsplit(ddG_lines[X], "_")[[1]][2:3],
                collapse = "_")})
        table_dg = data.table(
            value = best_model[,unlist(.SD),
            .SDcols = grep("ddG", names(best_model), value = T)], 
            variant,
            type)

        global_par_sel = best_model[,.SD,.SDcols = grep("^[sb]",names(best_model))]

        all_data <- fread(paste0("processed_data/GRB2_dG_dataset.txt"))
        all_data[, s_ddG1 := table_dg[variant == id1 & type == "s_ddG", value], id1]
        all_data[, b_ddG1 := table_dg[variant == id1 & type == "b_ddG", value], id1]
        all_data[, s_ddG2 := table_dg[variant == id2 & type == "s_ddG", value], id2]
        all_data[, b_ddG2 := table_dg[variant == id2 & type == "b_ddG", value], id2]


          # calculate stabilityPCA fitness from dG values
        all_data[, s_fitness_pred := function_folding_dG2F(
                s_ddg = sum(c(s_ddG1, s_ddG2), na.rm = T),
                s_dgwt = global_par_sel$s_dgwt,
                s_f0 = global_par_sel$s_f0, 
                s_fwt = global_par_sel$s_fwt
            ), 
            .(s_ddG1, s_ddG2)
        ]

          # calculate bindingPCA fitness if b_ddG = 0; i.e. global relationship governed by s_ddG
        all_data[, 
            b_fitness_pred_bdG0 := function_binding_dG2F(
                b_ddg = 0,
                s_ddg = sum(c(s_ddG1, s_ddG2), na.rm = T),
                b_dgwt = global_par_sel$b_dgwt,
                s_dgwt = global_par_sel$s_dgwt,
                b_f0 = global_par_sel$b_f0,
                b_fwt = global_par_sel$b_fwt
            ), 
            .(s_ddG1, s_ddG2)
        ]
          # calculate bindingPCA fitness with estimated b_ddG values
        all_data[, 
            b_fitness_pred := function_binding_dG2F(
                b_ddg = sum(c(b_ddG1, b_ddG2), na.rm = T),
                s_ddg = sum(c(s_ddG1, s_ddG2), na.rm = T),
                b_dgwt = global_par_sel$b_dgwt,
                s_dgwt = global_par_sel$s_dgwt,
                b_f0 = global_par_sel$b_f0,
                b_fwt = global_par_sel$b_fwt
            ),
            .(b_ddG1, s_ddG1, b_ddG2, s_ddG2)
        ]


          ## relationship between parameters and fitness
          # load structural property file to color residues by type
        variant_structuralproperties = 
            fread("processed_data/GRB2_variant_structuralproperties.txt")
        all_data[,res_type1 := 
            variant_structuralproperties[Pos == Pos1, type], Pos1]

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
            aes(color = res_type1, alpha = 0.3)
        )
        ggsave(
            plot = p, 
            file = paste0("results/dG/GRB2", 
                ifelse(d==1, "", "_allvars"), "_method", 
                i, "_iteration", j, "_parameters_fitness.pdf"), 
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
            aes(alpha = 0.1)
        )
        ggsave(
            plot = p, 
            file = paste0("results/dG/GRB2", 
                ifelse(d==1, "", "_allvars"), "_method", 
                i, "_iteration", j, "_parameters_fitness_doubles.pdf"), 
            width = 10, 
            height = 10
        )


          ## global relationship in bindingPCA and stabilityPCA space
        if (sum(grep("s_fwt", names(global_par_sel))) > 0) {
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
                color = "used in fitting"
            )
            ggsave(
                plot = p, 
                file = paste0("results/dG/GRB2", 
                    ifelse(d==1, "", "_allvars"), "_method", 
                    i, "_iteration", j, "_fitness_globalrelationship.pdf"),
                width = 6, 
                height = 6)
        }
    }
}

#### exclude weird iterations
DT[method == 2 & dataset == "dual_observed" & rank == 1,
.SD,
.SDcols = grep("^[iosb]", names(DT), value = T)
]
# first two iterations, run separately as tests are clearly off
# that's true for all "dual_observed" dataset methods 1-3
DT = DT[!(dataset == "dual_observed" & iteration %in% c(1,2))]

DT[method == 1 & dataset == "all_vars" & rank == 1,
.SD,
.SDcols = grep("^[iosb]", names(DT), value = T)
]
# the four best iterations seem very weird, having -10 wildtype dG binding and -13 b_f0
DT = DT[!(method == 1 & dataset == "all_vars" & iteration %in% c(70,4,97,71))]

DT[method == 2 & dataset == "all_vars" & rank == 1,
.SD,
.SDcols = grep("^[iosb]", names(DT), value = T)
]
# all fine

DT[method == 4 & dataset == "all_vars" & rank == 1,
.SD,
.SDcols = grep("^[iosb]", names(DT), value = T)
]

# [between(
#     objective,
#     mean(objective) - sd(objective),
#     mean(objective) + sd(objective)
# )]

DT[ method == 1 & dataset == "dual_observed" & iteration == 61,
.SD,
.SDcols = grep("^[iosb]", names(DT), value = T)
]

