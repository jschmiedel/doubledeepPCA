## analyse PDZ dG estimation models

require(data.table)
require(ggplot2)
require(GGally)
require(gridExtra)
theme_set(theme_bw(base_size = 9))

files <- list.files("processed_data/tmp/")

files <- files[grep("PDZ", files)]

for (i in seq_along(files)) {
    output <- fread(file.path("processed_data/tmp/", files[i]))
    tmp <- strsplit(files[i], "_")[[1]]

    output[, iteration := as.integer(gsub("\\.txt", "", tmp[length(tmp)]))]

    if (i == 1) {
        DT <- output
        } else {
        DT <- rbind(DT,
            output,
            fill = TRUE
        )
    }
}

DT[, modelnr := .N:1, .(iteration)]

# how many models?
DT[, .N, .(iteration)][, mean(N)]

DT[, objective := log10(objective)]
DT = DT[is.finite(objective)]
setkey(DT, objective)
DT[, rank := 1:.N, .(iteration)]

DT[, .SD, .SDcols = grep("^[ormibs]",names(DT))]

# save data.table with models
write.table(DT, 
    file = "processed_data/dG_PDZ_10models.txt",
    quote = F, row.names = F)

# plot for one process how models improved
ggpairs(DT[iteration == 8,
    .SD,
    .SDcols = grep("^[osb]", names(DT), value = T)
])
ggsave("results/dG/PDZ_dG_it8.pdf")


# compare iterations for same method/dataset
n <- grep("^[osb]", names(DT), value = T)
p <- ggpairs(DT[rank == 1,
        .SD,
        .SDcols = n
        ][between(
            objective,
            -Inf,
            mean(objective) + sd(objective)
        )],
        aes(alpha = 0.01))
ggsave(
    plot = p,
    filename = paste0("results/dG/PDZ_allmodels.pdf"),
    width = 12,
    height = 12
)

best_model = DT[order(objective)][1]
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

all_data <- fread(paste0("processed_data/PDZ_dG_dataset.txt"))
all_data[, s_ddG := table_dg[variant == id1 & type == "s_ddG", value], id1]
all_data[, b12_ddG := table_dg[variant == id1 & type == "b12_ddG", value], id1]
all_data[, b3_ddG := table_dg[variant == id1 & type == "b3_ddG", value], id1]


  # calculate stabilityPCA fitness from dG values
all_data[, s_fitness_pred := function_folding_dG2F(
        s_ddg = s_ddG,
        s_dgwt = global_par_sel$s1_dgwt,
        s_f0 = global_par_sel$s_f0, 
        s_fwt = global_par_sel$s_fwt
    ), 
    .(s_ddG)
]

###### b1
all_data[, 
    b1_fitness_pred_bdG0 := function_binding_dG2F(
        b_ddg = 0,
        s_ddg = s_ddG,
        b_dgwt = global_par_sel$b1_dgwt,
        s_dgwt = global_par_sel$s1_dgwt,
        b_f0 = global_par_sel$b1_f0,
        b_fwt = global_par_sel$b1_fwt
    ), 
    .(s_ddG)
]
all_data[, 
    b1_fitness_pred := function_binding_dG2F(
        b_ddg = b12_ddG,
        s_ddg = s_ddG,
        b_dgwt = global_par_sel$b1_dgwt,
        s_dgwt = global_par_sel$s1_dgwt,
        b_f0 = global_par_sel$b1_f0,
        b_fwt = global_par_sel$b1_fwt
    ),
    .(b12_ddG, s_ddG)
]

###### b2
all_data[, 
    b2_fitness_pred_bdG0 := function_binding_dG2F(
        b_ddg = 0,
        s_ddg = s_ddG,
        b_dgwt = global_par_sel$b2_dgwt,
        s_dgwt = global_par_sel$s23_dgwt,
        b_f0 = global_par_sel$b2_f0,
        b_fwt = global_par_sel$b2_fwt
    ), 
    .(s_ddG)
]
all_data[, 
    b2_fitness_pred := function_binding_dG2F(
        b_ddg = b12_ddG,
        s_ddg = s_ddG,
        b_dgwt = global_par_sel$b2_dgwt,
        s_dgwt = global_par_sel$s23_dgwt,
        b_f0 = global_par_sel$b2_f0,
        b_fwt = global_par_sel$b2_fwt
    ),
    .(b12_ddG, s_ddG)
]

###### b3
all_data[, 
    b3_fitness_pred_bdG0 := function_binding_dG2F(
        b_ddg = 0,
        s_ddg = s_ddG,
        b_dgwt = global_par_sel$b3_dgwt,
        s_dgwt = global_par_sel$s23_dgwt,
        b_f0 = global_par_sel$b3_f0,
        b_fwt = global_par_sel$b3_fwt
    ), 
    .(s_ddG)
]
all_data[, 
    b3_fitness_pred := function_binding_dG2F(
        b_ddg = b3_ddG,
        s_ddg = s_ddG,
        b_dgwt = global_par_sel$b3_dgwt,
        s_dgwt = global_par_sel$s23_dgwt,
        b_f0 = global_par_sel$b3_f0,
        b_fwt = global_par_sel$b3_fwt
    ),
    .(b3_ddG, s_ddG)
]

  ## relationship between parameters and fitness
  # load structural property file to color residues by type
PDZstruct <- fread("processed_data/PDZ_singles_alldata.txt")
PDZstruct[,ID := paste0(Pos, Mut, collapse = ""),.(Pos,Mut)]
all_data[, type := PDZstruct[ID == id1,type], id1]

  ## singles
p <- ggpairs(
    all_data[sample(.N, 250)],
    columns = c(
        "s_ddG", 
        "b12_ddG", 
        "b3_ddG", 
        "s_fitness", 
        "s_fitness_pred", 
        "b1_fitness", 
        "b1_fitness_pred",
        "b1_fitness_pred_bdG0", 
        "b2_fitness",
        "b2_fitness_pred",
        "b2_fitness_pred_bdG0", 
        "b3_fitness",
        "b3_fitness_pred",
        "b3_fitness_pred_bdG0"
    ),
    # aes(alpha = 0.01, color = type)
    aes(alpha = 0.01, color = b12_ddG > 1 & s_ddG < 1)
)
ggsave(
    plot = p, 
    file = paste0("results/dG/PDZ_parameters_fitness.pdf"), 
    width = 20, 
    height = 20
)

# only ddG values
p <- ggpairs(
    all_data,
    columns = c(
        "s_ddG", 
        "b12_ddG", 
        "b3_ddG"
    ),
    aes(alpha = 0.1, color = type)
)
ggsave(
    plot = p, 
    file = paste0("results/dG/PDZ_ddG_value_scatter.pdf"), 
    width = 6, 
    height = 6
)


## global relationship in bindingPCA and stabilityPCA space
p1 <- ggplot(all_data) +
    geom_point(aes(
            b1_fitness, 
            s_fitness, 
            color = b12_ddG,
        ), alpha = 1) +
    geom_point(aes(
            b1_fitness_pred_bdG0, 
            s_fitness_pred), 
        color = "red",
        size = 2
    ) +
    geom_segment(aes(
            x = b1_fitness_pred_bdG0, 
            y = s_fitness_pred, 
            xend = b1_fitness, 
            yend = s_fitness
        ), 
        alpha = 0.1
    ) +
    scale_color_gradient2(midpoint = 0, mid = "grey") +
    labs(
        x = "fitness bindingPCA", 
        y = "fitness stabilityPCA"
    )
p2 <- ggplot(all_data) +
    geom_point(aes(
            b1_fitness, 
            s_fitness, 
            color = s_ddG,
        ), alpha = 1) +
    geom_point(aes(
            b1_fitness_pred_bdG0, 
            s_fitness_pred), 
        color = "red",
        size = 2
    ) +
    geom_segment(aes(
            x = b1_fitness_pred_bdG0, 
            y = s_fitness_pred, 
            xend = b1_fitness, 
            yend = s_fitness
        ), 
        alpha = 0.1
    ) +
    scale_color_gradient2(midpoint = 0, mid = "grey") +
    labs(
        x = "fitness bindingPCA", 
        y = "fitness stabilityPCA"
    )

# ggsave(
#     plot = p, 
#     file = paste0("results/dG/PDZ_globalrel_ddPCA.pdf"),
#     width = 8, 
#     height = 4)

p3 <- ggplot(all_data) +
    geom_point(aes(
            b2_fitness, 
            s_fitness,
            color = b12_ddG,
        ), alpha = 1) +
    geom_point(aes(
            b2_fitness_pred_bdG0, 
            s_fitness_pred), 
        color = "red",
        size = 2
    ) +
    geom_segment(aes(
            x = b2_fitness_pred_bdG0, 
            y = s_fitness_pred, 
            xend = b2_fitness, 
            yend = s_fitness), 
        alpha = 0.1
    ) +
    scale_color_gradient2(midpoint = 0, mid = "grey") +
    labs(
        x = "fitness CRIPT binding", 
        y = "fitness stabilityPCA"
    )
p4 <- ggplot(all_data) +
    geom_point(aes(
            b2_fitness, 
            s_fitness,
            color = s_ddG,
        ), alpha = 1) +
    geom_point(aes(
            b2_fitness_pred_bdG0, 
            s_fitness_pred), 
        color = "red",
        size = 2
    ) +
    geom_segment(aes(
            x = b2_fitness_pred_bdG0, 
            y = s_fitness_pred, 
            xend = b2_fitness, 
            yend = s_fitness), 
        alpha = 0.1
    ) +
    scale_color_gradient2(midpoint = 0, mid = "grey") +
    labs(
        x = "fitness CRIPT binding", 
        y = "fitness stabilityPCA"
    )

# p <- grid.arrange(p1, p2, p3, p4,
#     nrow = 2)
# ggsave(
#     plot = p, 
#     file = paste0("results/dG/PDZ_globalrel_CRIPT.pdf"),
#     width = 8, 
#     height = 8)

p5 <- ggplot(all_data) +
    geom_point(aes(
            b3_fitness, 
            s_fitness,
            color = b3_ddG,
        ), alpha = 1) +
    geom_point(aes(
            b3_fitness_pred_bdG0, 
            s_fitness_pred), 
        color = "red",
        size = 2
    ) +
    geom_segment(aes(
            x = b3_fitness_pred_bdG0, 
            y = s_fitness_pred, 
            xend = b3_fitness, 
            yend = s_fitness), 
        alpha = 0.1
    ) +
    scale_color_gradient2(midpoint = 0, mid = "grey") +
    labs(
        x = "fitness Tm2F binding", 
        y = "fitness stabilityPCA"
    )

p6 <- ggplot(all_data) +
    geom_point(aes(
            b3_fitness, 
            s_fitness,
            color = s_ddG,
        ), alpha = 1) +
    geom_point(aes(
            b3_fitness_pred_bdG0, 
            s_fitness_pred), 
        color = "red",
        size = 2
    ) +
    geom_segment(aes(
            x = b3_fitness_pred_bdG0, 
            y = s_fitness_pred, 
            xend = b3_fitness, 
            yend = s_fitness), 
        alpha = 0.1
    ) +
    scale_color_gradient2(midpoint = 0, mid = "grey") +
    labs(
        x = "fitness Tm2F binding", 
        y = "fitness stabilityPCA"
    )

p <- grid.arrange(p1, p2, p3, p4, p5, p6,
    nrow = 3)
ggsave(
    plot = p, 
    file = paste0("results/dG/PDZ_globalrel_stability_binding.pdf"),
    width = 10, 
    height = 12)


p <- ggplot(all_data) +
    geom_point(aes(
            b2_fitness, 
            b1_fitness,
            color = b12_ddG,
        ), alpha = 1) +
    # geom_point(aes(
    #         b3_fitness_pred_bdG0, 
    #         s_fitness_pred), 
    #     color = "red",
    #     size = 2
    # ) +
    # geom_segment(aes(
    #         x = b3_fitness_pred_bdG0, 
    #         y = s_fitness_pred, 
    #         xend = b3_fitness, 
    #         yend = s_fitness), 
    #     alpha = 0.1
    # ) +
    scale_color_gradient2(midpoint = 0, mid = "grey") +
    labs(
        x = "fitness CRIPT binding", 
        y = "fitness bindingPCA"
    )
ggsave(
    plot = p, 
    file = paste0("results/dG/PDZ_globalrel_bindingPCA_CRIPT.pdf"),
    width = 6, 
    height = 6)

#relationship between b12_ddG and b3_ddG
p <- ggplot(all_data,
    aes(b12_ddG,
        b3_ddG,
        color = type
    )) +
    geom_point() +
    geom_abline()

ggsave(plot = p, 
    file = "results/dG/PDZ_b12_b3_ddG.pdf")
    



