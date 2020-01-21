#ggpairs helper functions

ggpairs_log10_diagonal <- function(data, mapping, ...) {
  ggally_densityDiag(data, mapping, ...) + scale_x_log10()
}
ggpairs_log10_points <- function(data, mapping, ...) {
  ggally_points(data, mapping, alpha=0.5, ...) + scale_x_log10() + scale_y_log10() + geom_abline(color="red")
}
ggpairs_log10_density <- function(data, mapping, ...) {
  ggally_density(data, mapping, ...) + scale_x_log10() + scale_y_log10() + geom_abline(color="red")
}
ggpairs_density <- function(data, mapping, ...) {
  ggally_density(data, mapping, ...) + geom_abline(color="red")
}
ggpairs_log10_hex <- function(data, mapping, ...) {
  geom_hex(data, mapping, ...) + scale_x_log10() + scale_y_log10() + geom_abline(color="red") +
    scale_fill_gradient(low="grey",high="dodgerblue3",trans="log10")
}
ggpairs_log10_cor <- function(data, mapping, ...) {
  # preprocess the data for the correlation calculations
  data[[quo_name(mapping$x)]] <- log10(data[[quo_name(mapping$x)]])
  data[[quo_name(mapping$y)]] <- log10(data[[quo_name(mapping$y)]])
  
  ggally_cor(data, mapping, ...) + # grids will not match. hide them
    theme(
      panel.grid.major = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA)
    )
}