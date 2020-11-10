# get (minimal) pairwise distances of multi-mutation variants
get_pairdistance <- function(
  vd,
  dist = "scHAmin"
) {
  cm <- fread(file.path(dataset_folder, "data/contactmap_cis.txt"))
  vd[, tmp := gsub("[A-Z]", "", aa_subs), aa_subs]
  mult <- unique(vd[, grep("_", tmp, value = T)])
  for (m in mult) {
    as <- as.integer(strsplit(m, "_")[[1]])
    idx <- utils::combn(length(as), 2)
    mindist <- c()
    for (i in 1:ncol(idx)) {
      mindist[i] <- cm[Pos1 == as[idx[1, i]] & Pos2 == as[idx[2, i]],
        unlist(.SD),
        .SDcols = dist]
    }
    vd[tmp == m, min_pairdist := min(mindist)]
  }
  vd[, tmp := NULL]
  return(vd)
}
