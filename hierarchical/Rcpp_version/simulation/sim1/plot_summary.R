dev.off()
graphics.off()
index <- 1
png(paste0(plot_dir, "sim1_sc", index, ".png"), width = 2000, height=900)
draw_density_ci <- function(x_coord, y_coord, s, s_quantile, index){
  cex.axis <- 3
  cex.lab <- 4
  cex.main <- 4.5
  line <- 4.5
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  ww <- (3 + mar.orig[2L]) * par("csi") * 2.54
  
  dev.new()
  #layout(matrix(4:1, ncol = 4L), widths = c(2.1, 2, 2, lcm(ww-1.5)))
  cex.axis <- 2.2
  cex.lab <- 4
  cex.main <- 4.5
  line <- 4.5
  layout(matrix(c(4:2, 1, 1, 1), nrow = 2, byrow = TRUE), widths = c(1, 1, 1), heights = c(1,lcm(ww)))
  layout.show(n=4)
  par(las = 1)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  levels <- pretty(zlim, 20)
  levels_range <- c(range(levels)[1] - .5, range(levels)[2]+1)
  
  plot.new()
  #plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = 'i',
  #            yaxs = 'i')
  #rect(0, levels[-length(levels)], 1, levels[-1L], col=jet.colors(length(levels)-1))
  #axis(4, cex.axis = cex.axis)
  plot.window(ylim = c(0, 3), xlim = levels_range, xaxs = 'i',
              yaxs = 'i')
  
  rect(levels[-length(levels)]+.5, 0, levels[-1L]+.5, 3, col=jet.colors(length(levels)-1))
  axis(1, at = levels+.5, labels = levels, cex.axis = cex.axis, font.axis=2)
  
  
  xlim <- range(x_coord)
  ylim <- range(y_coord)
  
  mar <- mar.orig
  mar[4L] <- 1.8
  mar[c(1L, 3L)] <- 7.1
  mar[2L] <- 9.1
  par(mar = mar)
  ## upper bound
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s_quantile[[index]][2, (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  title(main = bquote("log SD: "*hat(s)[.(index)]), xlab = "time", ylab = "frequency", line = line,
        cex.lab = cex.lab, cex.main = cex.main)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, cex.axis=cex.axis, font.axis=2)
  box()
  
  ## mean
  mar <- mar.orig
  mar[4L] <- 1.8
  mar[c(1L, 3L)] <- 7.1
  mar[2L] <- 9.1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s[index, (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  title(main = bquote("log SD: "*hat(s)[.(index)]), xlab = "time", ylab = "frequency", line = line,
        cex.lab = cex.lab, cex.main = cex.main)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, cex.axis=cex.axis, font.axis=2)
  box()
  
  ## lower bound
  mar <- mar.orig
  mar[4L] <- 1.8
  mar[c(1L, 3L)] <- 7.1
  mar[2L] <- 9.1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s_quantile[[index]][1, (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  title(main = bquote("log SD: "*hat(s)[.(index)]), xlab = "time", ylab = "frequency", line = line,
        cex.lab = cex.lab, cex.main = cex.main)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, cex.axis=cex.axis, font.axis=2)
  box()
}
dev.off()




graphics.off()

draw_density_te <- function(x_coord, y_coord, s, s_true, index, main1, main2){
  cex.axis <- 3
  cex.lab <- 4
  cex.main <- 4.5
  line <- 4.5
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  ww <- (3 + mar.orig[2L]) * par("csi") * 2.54
  dev.new()
  #layout(matrix(4:1, ncol = 4L), widths = c(2.1, 2, 2, lcm(ww-1.5)))
  cex.axis <- 2.2
  cex.lab <- 4
  cex.main <- 4.5
  line <- 4.5
  layout(matrix(c(3:2, 1, 1), nrow = 2, byrow = TRUE), heights = c(1,lcm(ww+1)))
  layout.show(n=3)
  par(las = 1)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  levels <- pretty(zlim, 20)
  levels_range <- c(range(levels)[1] - .5, range(levels)[2]+1)
  
  plot.new()
  #plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = 'i',
  #            yaxs = 'i')
  #rect(0, levels[-length(levels)], 1, levels[-1L], col=jet.colors(length(levels)-1))
  #axis(4, cex.axis = cex.axis)
  plot.window(ylim = c(0, 3), xlim = levels_range, xaxs = 'i',
              yaxs = 'i')
  
  rect(levels[-length(levels)]+.6, 0, levels[-1L]+.6, 3, col=jet.colors(length(levels)-1))
  axis(1, at = levels+.6, labels = levels, cex.axis = cex.axis, font.axis=2)
  
  
  xlim <- range(x_coord)
  ylim <- range(y_coord)
  
  
  
  ## mean
  mar <- mar.orig
  mar[4L] <- 1.8
  mar[c(1L, 3L)] <- 7.1
  mar[2L] <- 9.1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s[index, (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  title(main = main1, xlab = "time", ylab = "frequency", line = line,
        cex.lab = cex.lab, cex.main = cex.main)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, cex.axis=cex.axis, font.axis=2)
  box()
  
  ## true
  mar <- mar.orig
  mar[4L] <- 1.8
  mar[c(1L, 3L)] <- 7.1
  mar[2L] <- 9.1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s_true[index, (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  title(main = main2, xlab = "time", ylab = "frequency", line = line,
        cex.lab = cex.lab, cex.main = cex.main)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, cex.axis=cex.axis, font.axis=2)
  box()
}

setwd("/Users/johnn/Documents/Research/Project3/hierarchical/Rcpp_version/simulation/sim1/")
load("simulation_data.RData")
load("simulation_results.RData")
library(PARCOR)



set.seed(1234)
n_t <- 1024
n_I <- 5
P <- 2
#wt <- 0.01
vt <- 0.5
et <- .8
sim <- y.sim
yt <- sim$yt
true_ar <- sim$phit
true_ar_mean <- array(cbind((0.1*(1:n_t)/n_t + 0.85)*cos(2*pi/sim$lambdat_mu), -0.9^2), dim = c(1, n_t, 2))


## span of frequence
w <- seq(0.001, 0.499, by = 0.001)

### True spectral density
s_true <- cp_sd_uni(phi=true_ar, sigma2 = rep(et, n_t), w=w)
s_mean_true <- cp_sd_uni(phi=true_ar_mean, sigma2=rep(et, n_t), w=w)

## span of frequence
w <- seq(0.001, 0.499, by = 0.001)

zlim <- range(s, s_quantile, s_mean_quantile, s_true, s_mean_true)
P <- 5

x_coord <- seq(0, n_t-2*P-1, by = 1)
y_coord <- w

plot_dir <- "/Users/johnn/Documents/Research/Project3/hierarchical/for_papers/sim1/"
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                                 "yellow", "#FF7F00", "red", "#7F0000"))



zlim <- range(s, s_mean, s_mean_quantile, s_true, s_mean_true)

draw_density_ci(x_coord = x_coord, y_coord = y_coord, s = s, s_quantile = s_quantile, 
                index = 1)

draw_density_te(x_coord = x_coord, y_coord = y_coord, s = s, s_true = s_true, 
                main1 = bquote("log SD: "*hat(s)[.(index)]),
                main2 = bquote("log SD: "*s[.(index)]),
                index = 1)

draw_density_te(x_coord = x_coord, y_coord = y_coord, s = s, s_true = s_true, 
                main1 = bquote("log SD: "*hat(s)[.(index)]),
                main2 = bquote("log SD: "*s[.(index)]),
                index = 4)

draw_density_te(x_coord = x_coord, y_coord = y_coord, s = s, s_true = s_true, 
                main1 = bquote("log SD: "*hat(s)[.(index)]),
                main2 = bquote("log SD: "*s[.(index)]),
                index = 5)


draw_density_te(x_coord = x_coord, y_coord = y_coord, s = s_mean, s_true = s_mean_true, 
                main1 = bquote("log SD: "*hat(bar(s))),
                main2 = bquote("log SD: "*bar(s)),
                index = 1)




