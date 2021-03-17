####################################################
###### draw spectral density for walking EEG #######
####################################################
draw_density_hier_wEEG <- function(w, times, P, n_t, s, ...){
  constant_hz <- 128
  y_coord <- w*128
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                                   "yellow", "#FF7F00", "red", "#7F0000"))
  filled.contour(times[(P+1):(n_t-P)]/1000, y_coord, s[(P+1):(n_t-P), ], 
                 xlab = 'time (s)',
                 ylab = 'frequency (Hz)',
                 color.palette = jet.colors, ...)
  cat('Graph has been drawn!\n')
}





##################################
###### load data #########
library(PARCOR)
library(snowfall)
######## load dataset #########
root_dir <- getwd()

######## subject id #######
subject_id <- 6
load(file = paste0(root_dir, "/data/S", subject_id, ".RData"))


####### name of cluster area #######
cond_type <- "walk_rotate"
cluster_area <- EEG_data[[cond_type]]$names
times <- EEG_data[[cond_type]]$times
raw_data <- EEG_data[[cond_type]]$data

### do the first difference
data_diff <- apply(raw_data, c(1,3), diff)
data <- aperm(data_diff, c(2,1,3))

n_t <- dim(data)[2]
P <- 10
diff_times <- as.numeric(times[1, 1:n_t+1])

######### load results ############
load(file = paste0(root_dir, "/results/", cond_type, "_S", subject_id, ".RData"))

######### plot directory ########
plot_dir <- paste0(root_dir, "/results/plots/subject_", subject_id, "/", cond_type)

######### span of frequency #########
w <- seq(0.001, 0.499, by = 0.001)

######### range of zlim ##########
zlim <- c(-2, 10)

for(index in 1:dim(data)[1]){
  png(filename = paste0(plot_dir, '/est_mean_ch_', cluster_area[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_wEEG(w = w, P = P, times = times,
                         n_t = n_t, s = s_mean[[index]][1, , ],
                         main = bquote("log spectral density: "*.(cluster_area[index])),
                         zlim = zlim)
  dev.off()
  
  png(filename = paste0(plot_dir, '/est_lb_ch_', cluster_area[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_wEEG(w = w, P = P, times = times,
                        n_t = n_t, s = s_mean_quantile[[index]][1, , ],
                        main = bquote("log spectral density: "*.(cluster_area[index])),
                        zlim = zlim)
  dev.off()
  
  
  png(filename = paste0(plot_dir, '/est_ub_ch_', cluster_area[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_wEEG(w = w,  P = P, times=times,
                        n_t = n_t, s = s_mean_quantile[[index]][2, , ],
                        main = bquote("log spectral density: "*.(cluster_area[index])),
                        zlim = zlim)
  dev.off()
}

###########################################################################################

####### name of cluster area #######
cond_type <- "walk_pull"
cluster_area <- EEG_data[[cond_type]]$names
times <- EEG_data[[cond_type]]$times
raw_data <- EEG_data[[cond_type]]$data

### do the first difference
data_diff <- apply(raw_data, c(1,3), diff)
data <- aperm(data_diff, c(2,1,3))

n_t <- dim(data)[2]
P <- 10
diff_times <- as.numeric(times[1, 1:n_t+1])

######### load results ############
load(file = paste0(root_dir, "/results/", cond_type, "_S", subject_id, ".RData"))

######### plot directory ########
plot_dir <- paste0(root_dir, "/results/plots/subject_", subject_id, "/", cond_type)


########


for(index in 1:dim(data)[1]){
  png(filename = paste0(plot_dir, '/est_mean_ch_', cluster_area[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_wEEG(w = w, P = P, times = times,
                         n_t = n_t, s = s_mean[[index]][1, , ],
                         main = bquote("log spectral density: "*.(cluster_area[index])),
                         zlim = zlim)
  dev.off()
  
  png(filename = paste0(plot_dir, '/est_lb_ch_', cluster_area[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_wEEG(w = w, P = P, times = times,
                         n_t = n_t, s = s_mean_quantile[[index]][1, , ],
                         main = bquote("log spectral density: "*.(cluster_area[index])),
                         zlim = zlim)
  dev.off()
  
  
  png(filename = paste0(plot_dir, '/est_ub_ch_', cluster_area[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_wEEG(w = w,  P = P, times=times,
                         n_t = n_t, s = s_mean_quantile[[index]][2, , ],
                         main = bquote("log spectral density: "*.(cluster_area[index])),
                         zlim = zlim)
  dev.off()
}

###########################################################################################

####### name of cluster area #######
cond_type <- "stand_rotate"
cluster_area <- EEG_data[[cond_type]]$names
times <- EEG_data[[cond_type]]$times
raw_data <- EEG_data[[cond_type]]$data

### do the first difference
data_diff <- apply(raw_data, c(1,3), diff)
data <- aperm(data_diff, c(2,1,3))

n_t <- dim(data)[2]
P <- 10
diff_times <- as.numeric(times[1, 1:n_t+1])

######### load results ############
load(file = paste0(root_dir, "/results/", cond_type, "_S", subject_id, ".RData"))

######### plot directory ########
plot_dir <- paste0(root_dir, "/results/plots/subject_", subject_id, "/", cond_type)


########


for(index in 1:dim(data)[1]){
  png(filename = paste0(plot_dir, '/est_mean_ch_', cluster_area[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_wEEG(w = w, P = P, times = times,
                         n_t = n_t, s = s_mean[[index]][1, , ],
                         main = bquote("log spectral density: "*.(cluster_area[index])),
                         zlim = zlim)
  dev.off()
  
  png(filename = paste0(plot_dir, '/est_lb_ch_', cluster_area[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_wEEG(w = w, P = P, times = times,
                         n_t = n_t, s = s_mean_quantile[[index]][1, , ],
                         main = bquote("log spectral density: "*.(cluster_area[index])),
                         zlim = zlim)
  dev.off()
  
  
  png(filename = paste0(plot_dir, '/est_ub_ch_', cluster_area[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_wEEG(w = w,  P = P, times=times,
                         n_t = n_t, s = s_mean_quantile[[index]][2, , ],
                         main = bquote("log spectral density: "*.(cluster_area[index])),
                         zlim = zlim)
  dev.off()
}




###########################################################################################

####### name of cluster area #######
cond_type <- "stand_pull"
cluster_area <- EEG_data[[cond_type]]$names
times <- EEG_data[[cond_type]]$times
raw_data <- EEG_data[[cond_type]]$data

### do the first difference
data_diff <- apply(raw_data, c(1,3), diff)
data <- aperm(data_diff, c(2,1,3))

n_t <- dim(data)[2]
P <- 10
diff_times <- as.numeric(times[1, 1:n_t+1])

######### load results ############
load(file = paste0(root_dir, "/results/", cond_type, "_S", subject_id, ".RData"))

######### plot directory ########
plot_dir <- paste0(root_dir, "/results/plots/subject_", subject_id, "/", cond_type)


########


for(index in 1:dim(data)[1]){
  png(filename = paste0(plot_dir, '/est_mean_ch_', cluster_area[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_wEEG(w = w, P = P, times = times,
                         n_t = n_t, s = s_mean[[index]][1, , ],
                         main = bquote("log spectral density: "*.(cluster_area[index])),
                         zlim = zlim)
  dev.off()
  
  png(filename = paste0(plot_dir, '/est_lb_ch_', cluster_area[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_wEEG(w = w, P = P, times = times,
                         n_t = n_t, s = s_mean_quantile[[index]][1, , ],
                         main = bquote("log spectral density: "*.(cluster_area[index])),
                         zlim = zlim)
  dev.off()
  
  
  png(filename = paste0(plot_dir, '/est_ub_ch_', cluster_area[index], '.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier_wEEG(w = w,  P = P, times=times,
                         n_t = n_t, s = s_mean_quantile[[index]][2, , ],
                         main = bquote("log spectral density: "*.(cluster_area[index])),
                         zlim = zlim)
  dev.off()
}
