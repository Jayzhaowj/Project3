root_dir <- "/Users/johnn/Documents/Research/Project3/hierarchical/eeg/"
plot_dir <- paste0(root_dir, "plots/for_paper")

load(paste0(root_dir, "eeg_results.RData"))

### load some parameters
n_I <- 9
w <- seq(0.001, 0.499, by = 0.001)
label <- c("F3", "C3", "P3", "Fz", "Cz", "Pz", "F4", "C4", "P4")
P <- 15
n_t <- 3600

draw_density_hier_eeg <- function(w, index, P, n_t, s, ...){
  constant1 <- 84.375/3600
  x_coord <- seq(constant1*P, 84.375 - constant1*P, length.out = n_t - 2 * P)
  constant2 <- (256/6)
  y_coord <- constant2 * w
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                                   "yellow", "#FF7F00", "red", "#7F0000"))
  filled.contour(x_coord, y_coord, s[(P+1):(n_t-P), ], xlab = 'time',
                 ylab = 'frequency',
                 color.palette = jet.colors, ...)
  cat('Graph has been drawn!\n')
}


#####################################
## draw spectral density plots
#####################################
#### compute spectral density #####
#### for each channel
library(Rfast)
#max <- max(unlist(lapply(s, function(x) nth(x, 2, descending=TRUE))))
#zlim <- c(range(s, s_sample, na.rm=TRUE)[1], max)
zlim <- c(range(s, s_quantile, s_mean_quantile))

for(index in 1:n_I){
  #zlim <- c(range(s[index, , ], s_quantile), max)
  png(filename = paste0(plot_dir, '/est_mean_ch_', label[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_eeg(w = w, index = index, P = P,
                        n_t = n_t, s = s[index, , ],
                        main = bquote("log spectral density: "*.(label[index])),
                        zlim = zlim)
  dev.off()

  png(filename = paste0(plot_dir, '/est_lb_ch_', label[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_eeg(w = w, index = index, P = P,
                        n_t = n_t, s = s_quantile[[index]][1, , ],
                        main = bquote("log spectral density: "*.(label[index])),
                        zlim = zlim)
  dev.off()


  png(filename = paste0(plot_dir, '/est_ub_ch_', label[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_eeg(w = w, index = index, P = P,
                        n_t = n_t, s = s_quantile[[index]][2, , ],
                        main = bquote("log spectral density: "*.(label[index])),
                        zlim = zlim)
  dev.off()
  cat("\n The channel of", label[index], "has been drawn. \n")

}


#### for mean of baseline spectral density
index <- 1
png(filename = paste0(plot_dir, '/est_', index, 'mean.png'))
par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
draw_density_hier_eeg(w = w, index = index, P = P,
                      main = "Estimated baseline of log SD",
                      n_t = n_t, s = s_mean[1, , ], zlim = zlim)
dev.off()

##### 95% credible interval of baseline spectral density
index <- 1
png(filename = paste0(plot_dir, '/est_lb_', index, 'mean.png'))
par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
draw_density_hier_eeg(w = w, index = index, P = P,
                      main = "Estimated baseline of log SD",
                      n_t = n_t, s = s_mean_quantile[1, , ], zlim = zlim)
dev.off()


index <- 1
png(filename = paste0(plot_dir, '/est_ub_', index, 'mean.png'))
par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
draw_density_hier_eeg(w = w, index = index, P = P,
                      main = "Estimated baseline of log SD",
                      n_t = n_t, s = s_mean_quantile[2, , ], zlim = zlim)
dev.off()



png(filename = paste0(plot_dir, "/ll_fwd.png"))
par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
plot(1:15, ll_fwd, pch = 4, xlab = "model order", ylab = "log (likelihood)", type = "l")
points(1:15, ll_fwd, pch = 1)
dev.off()
