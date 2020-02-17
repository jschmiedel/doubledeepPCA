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

# how many models?
DT[, .N, .(method, dataset, iteration)][, mean(N), .(method, dataset)]
# for method 2 all_vars, there wasn't enough time to finish simulations

DT[, objective := log10(objective)]
DT = DT[is.finite(objective)]

# rank models
setkey(DT, objective)
DT[, rank := 1:.N, .(method, iteration, dataset)]

# plot for one process how models improved
ggpairs(DT[ method == 2 & dataset == "dual_observed" & iteration == 10,
    .SD,
    .SDcols = grep("^[osb]", names(DT), value = T)
])
ggsave("results/dG/example_method2_it10.pdf")

# compare iterations for same method/dataset
n <- grep("^[osb]", names(DT), value = T)
for (i in c(1, 2, 3)) {
    for (d in c(1, 2)) {
        if (i == 3) {
            n0 <- n[!grepl("^s_f", n)]
        } else {
            n0 <- n
        }
        p <- ggpairs(DT[method == i & dataset == ifelse(d == 1, "dual_observed", "all_vars") & rank == 1,
            .SD,
            .SDcols = n0
        ][between(
            objective,
            # mean(objective) - sd(objective),
            -Inf,
            mean(objective) + sd(objective)
        )])
        ggsave(
            plot = p,
            filename = paste0("results/dG/example_method", i, "_", ifelse(d == 1, "", "allvars_"), "allit.pdf")
        )
    }
}

#### exclude weird iterations
DT[method == 1 & dataset == "dual_observed" & rank == 1,
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

DT[method == 3 & dataset == "all_vars" & rank == 1,
    .SD,
    .SDcols = grep("^[iosb]", names(DT), value = T)
]

# [between(
#     objective,
#     mean(objective) - sd(objective),
#     mean(objective) + sd(objective)
# )]

DT[ method == 3 & dataset == "dual_observed" & iteration == 56,
    .SD,
    .SDcols = grep("^[iosb]", names(DT), value = T)
]

