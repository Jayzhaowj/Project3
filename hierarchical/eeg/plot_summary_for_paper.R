draw_density_ci <- function(x_coord, y_coord, s, s_quantile, index, main = "baseline log SD"){
  cex.axis <- 2.8
  cex.lab <- 4
  cex.main <- 4.5
  line <- 4
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  ww <- (3 + mar.orig[2L]) * par("csi") * 2.54
  
  dev.new()
  #layout(matrix(4:1, ncol = 4L), widths = c(2.1, 2, 2, lcm(ww-1.5)))

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
  .filled.contour(x_coord, y_coord, s_quantile[2, (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  title(main = main, xlab = "time", ylab = "frequency", line = line,
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
  title(main = main, xlab = "time", ylab = "frequency", line = line,
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
  .filled.contour(x_coord, y_coord, s_quantile[1, (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  title(main = main, xlab = "time", ylab = "frequency", line = line,
        cex.lab = cex.lab, cex.main = cex.main)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, cex.axis=cex.axis, font.axis=2)
  box()
}
dev.off()

graphics.off()
draw_density_te <- function(x_coord, y_coord, s, index1, index2, index3, main){
  cex.axis <- 2.8
  cex.lab <- 4
  cex.main <- 4.5
  line <- 4
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  ww <- (3 + mar.orig[2L]) * par("csi") * 2.54
  dev.new()
  #layout(matrix(4:1, ncol = 4L), widths = c(2.1, 2, 2, lcm(ww-1.5)))

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
  
  rect(levels[-length(levels)]+.6, 0, levels[-1L]+.6, 3, col=jet.colors(length(levels)-1))
  axis(1, at = levels+.6, labels = levels, cex.axis = cex.axis, font.axis=2)
  
  
  xlim <- range(x_coord)
  ylim <- range(y_coord)
  
  
  
  ## F4
  mar <- mar.orig
  mar[4L] <- 1.8
  mar[c(1L, 3L)] <- 7.1
  mar[2L] <- 9.1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s[index3, (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  title(main = main[index3], xlab = "time", ylab = "frequency", line = line,
        cex.lab = cex.lab, cex.main = cex.main)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, cex.axis=cex.axis, font.axis=2)
  box()
  
  ## Pz
  mar <- mar.orig
  mar[4L] <- 1.8
  mar[c(1L, 3L)] <- 7.1
  mar[2L] <- 9.1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s[index2, (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  title(main = main[index2], xlab = "time", ylab = "frequency", line = line,
        cex.lab = cex.lab, cex.main = cex.main)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, cex.axis=cex.axis, font.axis=2)
  box()
  
  ## Cz
  mar <- mar.orig
  mar[4L] <- 1.8
  mar[c(1L, 3L)] <- 7.1
  mar[2L] <- 9.1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s[index1, (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  title(main = main[index1], xlab = "time", ylab = "frequency", line = line,
        cex.lab = cex.lab, cex.main = cex.main)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, cex.axis=cex.axis, font.axis=2)
  box()
}

root_dir <- "/Users/johnn/Documents/Research/Project3/hierarchical/eeg/"
plot_dir <- paste0(root_dir, "plots/for_paper")

load(paste0(root_dir, "eeg_results.RData"))

### load some parameters
n_I <- 9
w <- seq(0.001, 0.499, by = 0.001)
label <- c("F3", "C3", "P3", "Fz", "Cz", "Pz", "F4", "C4", "P4")
P <- 15
n_t <- 3600


constant1 <- 84.375/3600
x_coord <- seq(constant1*P, 84.375 - constant1*P, length.out = n_t - 2 * P)
constant2 <- (256/6)
y_coord <- constant2 * w
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                                 "yellow", "#FF7F00", "red", "#7F0000"))

zlim <- range(s, s_mean_quantile, s_mean)
## draw Channel Cz, Pz, and F4
draw_density_te(x_coord = x_coord, y_coord = y_coord, s = s, index1 = 5, index2 = 6, index3 = 7, main = label)

## draw estimated mean of eeg channel
draw_density_ci(x_coord = x_coord, y_coord = y_coord, s = s_mean, 
                s_quantile = s_mean_quantile, index = 1)

## draw Channel Cz
draw_density_ci(x_coord = x_coord, y_coord = y_coord, s = s, 
                s_quantile = s_quantile[[5]], index = 5, main = label[5])

## draw Channel Pz
draw_density_ci(x_coord = x_coord, y_coord = y_coord, s = s, 
                s_quantile = s_quantile[[6]], index = 6, main = label[6])

## draw Channel F4
draw_density_ci(x_coord = x_coord, y_coord = y_coord, s = s, 
                s_quantile = s_quantile[[7]], index = 7, main = label[7])
