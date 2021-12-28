draw_density_ci <- function(x_coord, y_coord, s, s_quantile, index, main,
                            xlim, ylim, zlim){

  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  ww <- (3 + mar.orig[2L]) * par("csi") * 2.54
  loc <- c(0, 1)
  dev.new()
  #layout(matrix(4:1, ncol = 4L), widths = c(2.1, 2, 2, lcm(ww-1.5)))
  cex.axis <- 2.8
  cex.lab <- 4
  cex.main <- 4.5
  line <- 4
  layout(matrix(c(4:2, 1, 1, 1), nrow = 2, byrow = TRUE), widths = c(1, 1, 1), heights = c(1,lcm(ww)))
  layout.show(n=4)
  par(las = 1)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  levels <- pretty(zlim, 20)
  levels_range <- c(range(levels)[1] - .1, range(levels)[2]+.1)
  
  plot.new()
  #plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = 'i',
  #            yaxs = 'i')
  #rect(0, levels[-length(levels)], 1, levels[-1L], col=jet.colors(length(levels)-1))
  #axis(4, cex.axis = cex.axis)
  plot.window(ylim = c(0, 3), xlim = levels_range, xaxs = 'i',
              yaxs = 'i')
  
  rect(levels[-length(levels)]+.05, 0, levels[-1L]+.05, 3, col=jet.colors(length(levels)-1))
  axis(1, at = levels+.05, labels = levels, cex.axis = cex.axis, font.axis=2)
  

  mar <- mar.orig
  mar[4L] <- 1.8
  mar[c(1L, 3L)] <- 7.1
  mar[2L] <- 9.1
  par(mar = mar)
  ## upper bound
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s_quantile[[index]][2, (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  title(main = main, xlab = "time (s)", ylab = "frequency (Hz)", line = line,
        cex.lab = cex.lab, cex.main = cex.main)
  abline(v=loc, col = 'black', lty = 2)
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
  .filled.contour(x_coord, y_coord, s[(P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  title(main = main, xlab = "time (s)", ylab = "frequency (Hz)", line = line,
        cex.lab = cex.lab, cex.main = cex.main)
  abline(v=loc, col = 'black', lty = 2)
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
  title(main = main, xlab = "time (s)", ylab = "frequency (Hz)", line = line,
        cex.lab = cex.lab, cex.main = cex.main)
  abline(v=loc, col = 'black', lty = 2)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, cex.axis=cex.axis, font.axis=2)
  box()
}






  



root_dir <- "/Users/johnn/Documents/Research/Project3/hierarchical/walking_EEG/"
######## subject id #######
subject_id <- 25
load(file = paste0(root_dir, "/data/S", subject_id, ".RData"))


####### name of cluster area #######
cond_type <- "stand_pull"
cluster_area <- EEG_data[[cond_type]]$names
times <- EEG_data[[cond_type]]$times
raw_data <- EEG_data[[cond_type]]$data
data <- raw_data

# ### do the first difference
# data_diff <- apply(raw_data, c(1,3), diff)
# data <- aperm(data_diff, c(2,1,3))
# diff_times <- as.numeric(times[1, 1:n_t+1])

n_t <- dim(data)[2]
P <- 10


######### load results ############
s_all <- list(NA)
cond_type <- "stand_pull"
load(file = paste0(root_dir, "/results/", cond_type, "_S", subject_id, ".RData"))
s_all[[eval(cond_type)]] <- s_mean

cond_type <- "walk_pull"
load(file = paste0(root_dir, "/results/", cond_type, "_S", subject_id, ".RData"))
s_all[[eval(cond_type)]] <- s_mean

cond_type <- "stand_rotate"
load(file = paste0(root_dir, "/results/", cond_type, "_S", subject_id, ".RData"))
s_all[[eval(cond_type)]] <- s_mean

cond_type <- "walk_rotate"
load(file = paste0(root_dir, "/results/", cond_type, "_S", subject_id, ".RData"))
s_all[[eval(cond_type)]] <- s_mean

######### plot directory ########
plot_dir <- paste0(root_dir, "/results/plots/subject_", subject_id, "/abbvie/")


######### span of frequency #########
w <- seq(0.001, 0.499, by = 0.001)

zlim <- c(-0.8, 0.8)
ylim <- c(1, 50)
xlim <- c(-.5, 1.5)

x_coord <- times[(P+1):(n_t-P)]/1000
y_coord <- w*128
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                                 "yellow", "#FF7F00", "red", "#7F0000"))



#########################################################
#### stand pull and walk pull
########################################################
png(paste0(plot_dir, "stand_walk_pull.png"), width = 1200, height = 1000)
mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
on.exit(par(par.orig))
ww <- (3 + mar.orig[2L]) * par("csi") * 2.54
loc <- c(0, 1)
#layout(matrix(4:1, ncol = 4L), widths = c(2.1, 2, 2, lcm(ww-1.5)))
cex.axis <- 1.55
cex.lab <- 3
cex.main <- 3
line <- 4

layout(cbind(matrix(c(1:8), nrow = 4, byrow = TRUE), 
             matrix(9:16, nrow = 4, byrow = TRUE)),
       widths = c(1.2, 1, 1, 1.2), heights = c(1.2,1,1,1))
layout.show(n=16)
# par(las = 1)
# mar <- mar.orig
# mar[1L] <- 1
# mar[4L] <- mar[2L]
# mar[2L] <- 1
# par(mar = mar, mai=c(0.5,0.5,0.5,0))
levels <- pretty(zlim, 20)
levels_range <- c(range(levels)[1] - .1, range(levels)[2]+.1)
# 
# plot.new()
# #plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = 'i',
# #            yaxs = 'i')
# #rect(0, levels[-length(levels)], 1, levels[-1L], col=jet.colors(length(levels)-1))
# #axis(4, cex.axis = cex.axis)
# plot.window(ylim = c(0, 3), xlim = levels_range, xaxs = 'i',
#             yaxs = 'i')
# 
# rect(levels[-length(levels)]+.05, 0, levels[-1L]+.05, 3, col=jet.colors(length(levels)-1))
# axis(1, at = levels+.05, labels = levels, cex.axis = cex.axis, font.axis=2)
# 


## upper bound
index <- 1
for(i in 1:(length(cluster_area))){
  # ## stand pull
  if(index == 1){
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, 1.8, .8, 0))
  }else if(index > 1 & index < 5){
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, 1.8, 0.1, 0))
  }else if(index == 5){
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, .4, .8, .6))
  }else{
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, .4, 0.1,.6))
  }
  s_mean_diff <- t(apply(s_all$stand_pull[[i]][1, , ], 1, 
                         function(x) x - s_all$stand_pull[[i]][1, which(times==-500), ]))
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s_mean_diff[ (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  abline(v=loc, col = 'black', lty = 2)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, at = c(1, 4, 8, 13, 30, 50), 
       labels = as.character(c(1, 4, 8, 13, 30, 50)),
       las = 2,
       cex.axis=cex.axis, font.axis=2)
  box()
  if(index == 1 | index == 5){
    title(main = "stand pull", cex.main = cex.main)
  }
  if(index >=1 && index < 5){
    title(ylab = cluster_area[i], cex.lab = cex.lab, line = line)
  }
  
  ### walk pull
  if(index == 1){
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, .6, .8, .4))
  }else if(index > 1 & index < 5){
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, .6, 0.1, .4))
  }else if(index == 5){
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, 0, .8, 1.8))
  }else{
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, 0, 0.1, 1.8))
  }
  
  s_mean_diff <- t(apply(s_all$walk_pull[[i]][1, , ], 1, 
                         function(x) x - s_all$walk_pull[[i]][1, which(times==-500), ]))
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s_mean_diff[ (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  abline(v=loc, col = 'black', lty = 2)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, at = c(1, 4, 8, 13, 30, 50), 
       labels = as.character(c(1, 4, 8, 13, 30, 50)),
       las = 2,
       cex.axis=cex.axis, font.axis=2)
  box()
  if(index == 1 | index == 5){
    title(main = "walk pull", cex.main = cex.main)
  }
  if(index > 4){
    mtext(cluster_area[i], cex = 2, line = 3, side = 4)
  }
  index <- index + 1
}
# title(xlab = "time (s)", ylab = "frequency (Hz)", line = line,
#       cex.lab = cex.lab, cex.main = cex.main)
# 
dev.off()

#########################################################
#### stand rotate and walk rotate
########################################################
png(paste0(plot_dir, "stand_walk_rotate.png"), width = 1200, height = 1000)
mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
on.exit(par(par.orig))
ww <- (3 + mar.orig[2L]) * par("csi") * 2.54
loc <- c(0, 0.5)
#layout(matrix(4:1, ncol = 4L), widths = c(2.1, 2, 2, lcm(ww-1.5)))
cex.axis <- 1.55
cex.lab <- 3
cex.main <- 3
line <- 4

layout(cbind(matrix(c(1:8), nrow = 4, byrow = TRUE), 
             matrix(9:16, nrow = 4, byrow = TRUE)),
       widths = c(1.2, 1, 1, 1.2), heights = c(1.2,1,1,1))
layout.show(n=16)
# par(las = 1)
# mar <- mar.orig
# mar[1L] <- 1
# mar[4L] <- mar[2L]
# mar[2L] <- 1
# par(mar = mar, mai=c(0.5,0.5,0.5,0))
levels <- pretty(zlim, 20)
levels_range <- c(range(levels)[1] - .1, range(levels)[2]+.1)
# 
# plot.new()
# #plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = 'i',
# #            yaxs = 'i')
# #rect(0, levels[-length(levels)], 1, levels[-1L], col=jet.colors(length(levels)-1))
# #axis(4, cex.axis = cex.axis)
# plot.window(ylim = c(0, 3), xlim = levels_range, xaxs = 'i',
#             yaxs = 'i')
# 
# rect(levels[-length(levels)]+.05, 0, levels[-1L]+.05, 3, col=jet.colors(length(levels)-1))
# axis(1, at = levels+.05, labels = levels, cex.axis = cex.axis, font.axis=2)
# 


## upper bound
index <- 1
for(i in 1:(length(cluster_area))){
  # ## stand pull
  if(index == 1){
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, 1.8, .8, 0))
  }else if(index > 1 & index < 5){
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, 1.8, 0.1, 0))
  }else if(index == 5){
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, .4, .8, .6))
  }else{
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, .4, 0.1,.6))
  }
  s_mean_diff <- t(apply(s_all$stand_rotate[[i]][1, , ], 1, 
                         function(x) x - s_all$stand_rotate[[i]][1, which(times==-500), ]))
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s_mean_diff[ (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  abline(v=loc, col = 'black', lty = 2)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, at = c(1, 4, 8, 13, 30, 50), 
       labels = as.character(c(1, 4, 8, 13, 30, 50)),
       las = 2,
       cex.axis=cex.axis, font.axis=2)
  box()
  if(index == 1 | index == 5){
    title(main = "stand rotate", cex.main = cex.main)
  }
  if(index >=1 && index < 5){
    title(ylab = cluster_area[i], cex.lab = cex.lab, line = line)
  }
  
  ### walk pull
  if(index == 1){
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, .6, .8, .4))
  }else if(index > 1 & index < 5){
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, .6, 0.1, .4))
  }else if(index == 5){
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, 0, .8, 1.8))
  }else{
    mar <- mar.orig
    # mar[4L] <- 1.8
    # mar[c(1L, 3L)] <- 3.1
    # mar[2L] <- 9.1
    par(mar = mar, mai = c(0.3, 0, 0.1, 1.8))
  }
  
  s_mean_diff <- t(apply(s_all$walk_rotate[[i]][1, , ], 1, 
                         function(x) x - s_all$walk_rotate[[i]][1, which(times==-500), ]))
  plot.new()
  plot.window(xlim, ylim, "", xaxs="i", yaxs="i", asp = NA)
  .filled.contour(x_coord, y_coord, s_mean_diff[ (P+1):(n_t-P), ], as.double(levels), col = jet.colors(length(levels)-1))
  abline(v=loc, col = 'black', lty = 2)
  Axis(x_coord, side = 1, cex.axis=cex.axis, font.axis=2)
  Axis(y_coord, side = 2, at = c(1, 4, 8, 13, 30, 50), 
       labels = as.character(c(1, 4, 8, 13, 30, 50)),
       las = 2,
       cex.axis=cex.axis, font.axis=2)
  box()
  if(index == 1 | index == 5){
    title(main = "walk rotate", cex.main = cex.main)
  }
  if(index > 4){
    mtext(cluster_area[i], cex = 2, line = 3, side = 4)
  }
  index <- index + 1
}
# title(xlab = "time (s)", ylab = "frequency (Hz)", line = line,
#       cex.lab = cex.lab, cex.main = cex.main)
# 
dev.off()


#### draw scale
filled.contour(s_mean_diff, zlim = zlim, color.palette = jet.colors)
## draw density of supplementary motor
index <- 7
s_mean_diff <- t(apply(s_mean[[index]][1, , ], 1, 
                       function(x) x - s_mean[[index]][1, which(times==-500), ]))
draw_density_ci(x_coord = x_coord, y_coord = y_coord, s = s_mean_diff, 
                s_quantile = s_mean_diff_quantile, index = 7, main = cluster_area[7],
                xlim = xlim, ylim = ylim, zlim = zlim)
draw_density_te(x_coord = x_coord, y_coord = y_coord, s = s_mean, 
                index1 = 7, index2 = 3, index3 = 5, 
                main = cluster_area, xlim = xlim, ylim = ylim, zlim = zlim)


## left sensorimotor 
index <- 3
s_mean_diff <- t(apply(s_mean[[index]][1, , ], 1, 
                       function(x) x - s_mean[[index]][1, which(times==-500), ]))
draw_density_ci(x_coord = x_coord, y_coord = y_coord, s = s_mean_diff, 
                s_quantile = s_mean_diff_quantile, index = index, main = cluster_area[index],
                xlim = xlim, ylim = ylim, zlim = zlim)

## right sensorimotor 
index <- 5
s_mean_diff <- t(apply(s_mean[[index]][1, , ], 1, 
                       function(x) x - s_mean[[index]][1, which(times==-500), ]))
draw_density_ci(x_coord = x_coord, y_coord = y_coord, s = s_mean_diff, 
                s_quantile = s_mean_diff_quantile, index = index, main = cluster_area[index],
                xlim = xlim, ylim = ylim, zlim = zlim)
