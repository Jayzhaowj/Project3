library(PARCOR)
library(snowfall)
library(pracma)
######## load dataset #########
setwd("/Users/johnn/Documents/Research/Project3/hierarchical/walking_EEG/")
root_dir <- getwd()
plot_dir <- paste0(getwd(), "/results/plots/subject_25/abbvie/")
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                                 "yellow", "#FF7F00", "red", "#7F0000"))



col_band_sd <- function(phi, SIGMA, P_max, w_band, n_I, index){
  n_t <- dim(phi)[2]
  spectral_density <- matrix(NA, nrow = n_I, ncol = n_t-2*P_max)
  coherence <- matrix(NA, nrow = dim(index)[2], ncol = n_t-2*P_max)
  partial_coherence <- matrix(NA, nrow = dim(index)[2], ncol = n_t-2*P_max)
  tmp <- cp_sd(phi = phi, SIGMA=SIGMA, w = w_band)
  sd <- simplify2array(tmp$sd)
  f_spec <- apply(sd, c(1, 2, 4), function(x) trapz(w_band, x))
  f_dens <- array(apply(f_spec, 3, abs), dim = c(n_I, n_I, n_t))
  g_spec <- apply(f_spec, 3, function(x) pracma::inv(x))
  g_spec <- array(g_spec, dim = c(n_I, n_I, n_t))
  g_dens <- array(apply(g_spec, 3, abs), dim = c(n_I, n_I, n_t))
  for(i in 1:dim(index)[2]){
    if(i == 1){
      spectral_density[i, ] <- log(f_dens[index[1, i], index[1, i], (P_max+1):(n_t-P_max)])
      spectral_density[i+1, ] <- log(f_dens[index[2, i], index[2, i], (P_max+1):(n_t-P_max)])
    }else if (i < n_I-1){
      spectral_density[index[2, i], ] <- log(f_dens[index[2, i], index[2, i], (P_max+1):(n_t-P_max)])
    }
    coherence[i, ] <- (f_dens[index[1, i], index[2, i], (P_max+1):(n_t-P_max)])^2/(f_dens[index[1,i], index[1,i], (P_max+1):(n_t-P_max)]*f_dens[index[2, i], index[2, i], (P_max+1):(n_t-P_max)])
    partial_coherence[i, ] <- (g_dens[index[1, i], index[2, i], (P_max+1):(n_t-P_max)])^2/(g_dens[index[1,i], index[1,i], (P_max+1):(n_t-P_max)]*g_dens[index[2, i], index[2, i], (P_max+1):(n_t-P_max)])
  }
  return(list("sd" = spectral_density, "coh" = coherence, "pcoh" = partial_coherence))
}

col_band_dcoh <- function(phi, SIGMA, P_max, w_band, n_I, index){
  n_t <- dim(phi)[2]
  DTF <- matrix(NA, nrow = dim(index)[2], ncol = n_t-2*P_max)
  PDC <- matrix(NA, nrow = dim(index)[2], ncol = n_t-2*P_max)
  tmp <- cp_sd_dcoh(phi = phi, SIGMA=SIGMA, w = w_band)
  tmp_PHI <- apply(simplify2array(tmp$PHI), c(1,2,4), function(x) trapz(w_band, x))
  tmp_PHI_inv <- apply(simplify2array(tmp$PHI_inv), c(1,2,4), function(x) trapz(w_band, x))
  pdc <- array(dim = dim(tmp_PHI))
  dtf <- array(dim = dim(tmp_PHI_inv))
  
  ### PDC
  for(i in 1:dim(tmp_PHI)[3]){
    PHI_norm <- sqrt(diag(t(tmp_PHI[, , i]) %*% tmp_PHI[, , i]))
    pdc[, , i] <- t(apply(tmp_PHI[, , i], 1, function(x) x/PHI_norm))
  }
  
  ### DTF 
  for(i in 1:dim(tmp_PHI)[3]){
    PHI_inv_norm <- sqrt(diag(tmp_PHI[, , i] %*% t(tmp_PHI[, , i])))
    dtf[, , i] <- apply(tmp_PHI[, , i], 2, function(x) x/PHI_norm)
  }
  for(i in 1:dim(index)[2]){
    DTF[i, ] <- abs(dtf[index_all[1, i], index_all[2, i], (P_max+1):(n_t-P_max)])
    PDC[i, ] <- abs(pdc[index_all[1, i], index_all[2, i], (P_max+1):(n_t-P_max)])
  }
  return(list("PDC" = PDC, "DTF" = DTF))
}


######## subject id #######
subject_id <- 25
load(file = paste0(root_dir, "/data/S", subject_id, ".RData"))


####### name of cluster area #######
cond_type <- "stand_pull"
cluster_area <- EEG_data[[cond_type]]$names
times <- EEG_data[[cond_type]]$times
raw_data <- EEG_data[[cond_type]]$data
data_mean <- apply(raw_data, 1:2, mean)
data_median <- apply(raw_data, 1:2, median)
data_sd <- apply(raw_data, 1:2, sd)
is_difference <- "FALSE"

####### check the plots ###############
par(mfrow = c(2, 2))
plot(data_mean[1, ], type = 'l', main = "mean", ylab = "values", xlab = "time")
plot(data_median[1, ], type = 'l', main = "median", ylab = "values", xlab = "time")
plot(data_sd[1, ], type = 'l', main = "sd", ylab = "values", xlab = "time")
plot(raw_data[1, , 2], type = 'l', main = "epoch2", ylab = "values", xlab = "time")

####### fit model ###########
n_t <- dim(data_median)[2]
K <- dim(data_median)[1]
S_0 <- 10*diag(K)
P <- 15

#### set up discount factor ####
grid_seq <- seq(0.97, 1, 0.001)
delta <- array(dim = c(length(grid_seq), K^2, P))
#### discount factor for PARCOR model
for(i in 1:P){
  delta[, , i] <- matrix(rep(grid_seq, times = K), ncol = K)
}

result <- run_parcor_parallel(data_mean, delta = delta, P = P,
                              DIC = TRUE, S_0 = S_0, 
                              uncertainty = FALSE)
P_opt <- which.min(result$DIC)
ar_coef <- result$ar_coef[[P_opt]]$forward
w <- seq(0.001, 0.499, by = 0.001)
constant_hz <- 128
x_coord <- times[(P+1):(n_t-P)]/1000
y_coord <- w*constant_hz

####### delta band ########
index <- combn(K, 2)
w_band <- w[1 < constant_hz*w & constant_hz*w < 4]
band_delta <- col_band_sd(ar_coef, array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), P, w_band, K, index)

index_all <- t(expand.grid(1:K, 1:K))
band_delta_dtf <- col_band_dcoh(ar_coef, array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), P, w_band, K, index_all)

####### theta band #########
index <- combn(K, 2)
w_band <- w[4 < constant_hz*w & constant_hz*w < 8]
band_theta <- col_band_sd(ar_coef, array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), P, w_band, K, index)

index_all <- t(expand.grid(1:K, 1:K))
band_theta_dtf <- col_band_dcoh(ar_coef, array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), P, w_band, K, index_all)

####### alpha band #########
w_band <- w[8 < constant_hz*w & constant_hz*w < 13]
band_alpha <- col_band_sd(ar_coef, array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), P, w_band, K, index)
band_alpha_dtf <- col_band_dcoh(ar_coef, array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), P, w_band, K, index_all)

####### draw plots ########
####### delta band ########
##### coherence
type <- "Delta"
region <- "SMA"
cluster_area_abb <- c("LO", "RO", "LS", "AC", "RS", "PP", "SMA", "AP")
region_ind <- which(cluster_area == "supplementary motor")
suppl <- c(6, 12, 17, 21, 24, 26, 28)
col <- jet.colors(length(suppl))
png(paste0(plot_dir, cond_type, "/coh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
par(mfrow = c(1, 1), cex = 1.5)
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of coh between ", region, " and rest of clusters under stand pull") 
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_delta$coh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_delta$coh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)],
        col = col[i])
}
legend(x = 1.2, y = 1.05, legend = cluster_area_abb[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####### partial coherence
png(paste0(plot_dir, cond_type, "/pcoh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)

main <- paste0(type, " band of pcoh between ", region, " and rest of clusters under stand pull")
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_delta$pcoh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[i],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_delta$pcoh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 1.2, y = 1.05, legend = cluster_area_abb[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####################################################
####### theta band ########
##### coherence
type <- "Theta"
region <- "SMA"
region_ind <- which(cluster_area == "supplementary motor")
png(paste0(plot_dir, cond_type, "/coh_", region, "_", type, "_band.png"), 
    width = 985, height = 680)
suppl <- c(6, 12, 17, 21, 24, 26, 28)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of coh between ", region, " and rest of clusters under stand pull") 
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_theta$coh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_theta$coh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 1.2, y = 1.05, legend = cluster_area_abb[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####### partial coherence
png(paste0(plot_dir, cond_type, "/pcoh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of pcoh between ", region, " and rest of clusters under stand pull")
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_theta$pcoh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_theta$pcoh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 1.2, y = 1.05, legend = cluster_area_abb[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()

##### partial direct coherence
type <- "Theta"
region <- "SMA"
region_ind <- which(cluster_area == "supplementary motor")
png(paste0(plot_dir, cond_type, "/pdc_from_", region, "_", type, "_band.png"), 
    width = 985, height = 680)
suppl <- which(index_all[2, ] == region_ind)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of PDC from ", region, " to rest of clusters under stand pull") 
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_theta_dtf$PDC[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_theta_dtf$PDC[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 1.2, y = 1.05, legend = cluster_area_abb[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.1)
abline(v = loc, col = 'black', lty = 2)
dev.off()

###################################################
type <- "Theta"
region <- "SMA"
region_ind <- which(cluster_area == "supplementary motor")
png(paste0(plot_dir, cond_type, "/pdc_to_", region, "_", type, "_band.png"), 
    width = 985, height = 680)
suppl <- which(index_all[1, ] == region_ind)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of PDC from rest of clusters to ", region, " under stand pull") 
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_theta_dtf$PDC[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", col = col[1],
     ylim = c(range(band_theta_dtf$PDC)[1], range(band_theta_dtf$PDC)[2]+0.02), main = main)
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_theta_dtf$PDC[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend("topright", legend = cluster_area_abb[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.1)
abline(v = loc, col = 'black', lty = 2)
dev.off()

####################################################
##### PDC 
#############################
w <- seq(0.001, 0.499, 0.001)
sd <- cp_sd(phi = ar_coef, SIGMA = array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), w = w)
for(i in 1:dim(index_all)[2]){
  tmp_sd <- get_sd(sd$PDC, index_all[1, i], index_all[2, i], 4)[(P+1):(n_t - P),  ]
  png(paste0(plot_dir, cond_type, "/pdc_from_", cluster_area_abb[index_all[2, i]], "_to_", cluster_area_abb[index_all[1, i]], ".png"))
  par(par(mar=c(5,6,4,1)+.1), cex.lab = 2, cex.axis = 1.5, cex.main = 2.5)
  filled.contour(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
                 y_coord, tmp_sd[which(x_coord==-0.5):which(x_coord==1.5), ],
                 plot.title = {title(xlab="time (s)", ylab = "frequency (Hz)", main = bquote("PDC: "*.(cluster_area_abb[index_all[2, i]])%->%.(cluster_area_abb[index_all[1, i]])))
                   abline(v=loc, col = 'black', lty = 2)},
                 color.palette = jet.colors, zlim = c(0, 1))
  dev.off()
  cat("iterations: ", i, "/", dim(index_all)[2], "\r")
}


####################################################
##### DTF
#############################
w <- seq(0.001, 0.499, 0.001)
sd <- cp_sd(phi = ar_coef, SIGMA = array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), w = w)
for(i in 1:dim(index_all)[2]){
  tmp_sd <- get_sd(sd$DTF, index_all[1, i], index_all[2, i], 4)[(P+1):(n_t - P),  ]
  png(paste0(plot_dir, cond_type, "/dtf_from_", cluster_area_abb[index_all[2, i]], "_to_", cluster_area_abb[index_all[1, i]], ".png"))
  par(par(mar=c(5,6,4,1)+.1), cex.lab = 2, cex.axis = 1.5, cex.main = 2.5)
  filled.contour(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
                 y_coord, tmp_sd[which(x_coord==-0.5):which(x_coord==1.5), ],
                 plot.title = {title(xlab="time (s)", ylab = "frequency (Hz)", main = bquote("DTF: "*.(cluster_area_abb[index_all[2, i]])%->%.(cluster_area_abb[index_all[1, i]])))
                   abline(v=loc, col = 'black', lty = 2)},
                 color.palette = jet.colors, zlim = c(0, 1))
  dev.off()
  cat("iterations: ", i, "/", dim(index_all)[2], "\r")
}




####################################################

####### alpha band ########
##### coherence
type <- "Alpha"
region <- "LS"
region_ind <- which(cluster_area == "left sensorimotor")
png(paste0(plot_dir, cond_type, "/coh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
suppl <- c(2, 8, 14, 15, 16, 17, 18)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of coh between ", region, " and rest of clusters")
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_alpha$coh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_alpha$coh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 0.85, y = 1.0, legend = cluster_area[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####### partial coherence
png(paste0(plot_dir, cond_type, "/pcoh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of pcoh between ", region, " and rest of clusters")
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_theta$pcoh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_theta$pcoh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 0.85, y = 1.0, legend = cluster_area[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####################################################

##### coherence
type <- "Alpha"
region <- "RS"
region_ind <- which(cluster_area == "right sensorimotor")
png(paste0(plot_dir, cond_type, "/coh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
suppl <- c(4, 10, 15, 19, 23, 24, 25)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of coh between ", region, " and rest of clusters")
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_alpha$coh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_alpha$coh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 0.85, y = 1.0, legend = cluster_area[-7], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####### partial coherence
png(paste0(plot_dir, cond_type, "/pcoh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of pcoh between ", region, " and rest of clusters")
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_alpha$pcoh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_alpha$pcoh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 0.85, y = 1.0, legend = cluster_area[-7], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()




####### name of cluster area #######
cond_type <- "walk_pull"
cluster_area <- EEG_data[[cond_type]]$names
times <- EEG_data[[cond_type]]$times
raw_data <- EEG_data[[cond_type]]$data
data_mean <- apply(raw_data, 1:2, mean)
data_median <- apply(raw_data, 1:2, median)
data_sd <- apply(raw_data, 1:2, sd)
is_difference <- "FALSE"

####### check the plots ###############
par(mfrow = c(2, 2))
plot(data_mean[1, ], type = 'l', main = "mean", ylab = "values", xlab = "time")
plot(data_median[1, ], type = 'l', main = "median", ylab = "values", xlab = "time")
plot(data_sd[1, ], type = 'l', main = "sd", ylab = "values", xlab = "time")
plot(raw_data[1, , 5], type = 'l', main = "epoch2", ylab = "values", xlab = "time")

####### fit model ###########
n_t <- dim(data_median)[2]
K <- dim(data_median)[1]
S_0 <- 10*diag(K)
P <- 15

#### set up discount factor ####
grid_seq <- seq(0.97, 1, 0.001)
delta <- array(dim = c(length(grid_seq), K^2, P))
#### discount factor for PARCOR model
for(i in 1:P){
  delta[, , i] <- matrix(rep(grid_seq, times = K), ncol = K)
}

result <- run_parcor_parallel(data_mean, delta = delta, P = P,
                              DIC = TRUE, S_0 = S_0, 
                              uncertainty = FALSE)
P_opt <- which.min(result$DIC)
ar_coef <- result$ar_coef[[P_opt]]$forward
w <- seq(0.001, 0.499, by = 0.001)
constant_hz <- 128
x_coord <- times[(P+1):(n_t-P)]/1000
y_coord <- w*constant_hz
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                                 "yellow", "#FF7F00", "red", "#7F0000"))
####### delta band ########
index <- combn(K, 2)
w_band <- w[1 < constant_hz*w & constant_hz*w < 4]
band_delta <- col_band_sd(ar_coef, array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), P, w_band, K, index)

####### theta band #########
index <- combn(K, 2)
w_band <- w[4 < constant_hz*w & constant_hz*w < 8]
band_theta <- col_band_sd(ar_coef, 
                          array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), 
                          P, w_band, K, index)

index_all <- t(expand.grid(1:K, 1:K))
band_theta_dtf <- col_band_dcoh(ar_coef, array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), P, w_band, K, index_all)

####### alpha band #########
w_band <- w[8 < constant_hz*w & constant_hz*w < 13]
band_alpha <- col_band_sd(ar_coef, array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), P, w_band, K, index)
band_alpha_dtf <- col_band_dcoh(ar_coef, array(result$St_fwd[[P_opt]], dim = c(K, K, n_t)), P, w_band, K, index_all)


####### draw plots ########
####### delta band ########
##### coherence
type <- "Delta"
region <- "SMA"
region_ind <- which(cluster_area == "left sensorimotor")
suppl <- c(6, 12, 17, 21, 24, 26, 28)
col <- jet.colors(length(suppl))
png(paste0(plot_dir, cond_type, "/coh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of coh between ", region, " and rest of clusters") 
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_delta$coh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_delta$coh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 0.85, y = 1.0, legend = cluster_area[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####### partial coherence
png(paste0(plot_dir, cond_type, "/pcoh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)

main <- paste0(type, " band of pcoh between ", region, " and rest of clusters")
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_delta$pcoh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[i],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_delta$pcoh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 0.85, y = 1.0, legend = cluster_area[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####################################################
####### theta band ########
##### coherence
type <- "Theta"
region <- "SMA"
region_ind <- which(cluster_area == "supplementary motor")
png(paste0(plot_dir, cond_type, "/coh_", region, "_", type, "_band.png"), 
    width = 985, height = 680)
suppl <- c(6, 12, 17, 21, 24, 26, 28)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of coh between ", region, " and rest of clusters") 
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_theta$coh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_theta$coh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 0.85, y = 1.0, legend = cluster_area[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####### partial coherence
png(paste0(plot_dir, cond_type, "/pcoh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of pcoh between ", region, " and rest of clusters")
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_theta$pcoh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_theta$pcoh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 0.85, y = 1.0, legend = cluster_area[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####################################################

####### alpha band ########
##### coherence
type <- "Alpha"
region <- "LS"
region_ind <- which(cluster_area == "left sensorimotor")
png(paste0(plot_dir, cond_type, "/coh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
suppl <- c(2, 8, 14, 15, 16, 17, 18)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of coh between ", region, " and rest of clusters")
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_alpha$coh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_alpha$coh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 0.85, y = 1.0, legend = cluster_area[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####### partial coherence
png(paste0(plot_dir, cond_type, "/pcoh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of pcoh between ", region, " and rest of clusters")
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_theta$pcoh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_theta$pcoh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 0.85, y = 1.0, legend = cluster_area[-region_ind], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####################################################

##### coherence
type <- "Alpha"
region <- "RS"
region_ind <- which(cluster_area == "right sensorimotor")
png(paste0(plot_dir, cond_type, "/coh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
suppl <- c(4, 10, 15, 19, 23, 24, 25)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of coh between ", region, " and rest of clusters")
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_alpha$coh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_alpha$coh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 0.85, y = 1.0, legend = cluster_area[-7], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
####### partial coherence
png(paste0(plot_dir, cond_type, "/pcoh_", region, "_", type, "_band.png"),
    width = 985, height = 680)
par(mfrow = c(1, 1), cex = 1.5)
plot.new()
xlim <- c(-0.5, 1.5)
ylim <- c(0, 1)
loc <- c(0, 1)
main <- paste0(type, " band of pcoh between ", region, " and rest of clusters")
plot(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
     band_alpha$pcoh[suppl[1], which(x_coord==-0.5):which(x_coord==1.5)], 
     type = 'l', xlab = "time", ylab = "value", yaxt = 'n', col = col[1],
     ylim = c(0, 1), main = main)
axis(side = 2, at = c(0, 0.5, 1))
for(i in 2:length(suppl)){
  lines(x_coord[which(x_coord==-0.5):which(x_coord==1.5)], 
        band_alpha$pcoh[suppl[i], which(x_coord==-0.5):which(x_coord==1.5)], col = col[i])
}
legend(x = 0.85, y = 1.0, legend = cluster_area[-7], 
       col = col, lty = rep(1, length(suppl)), bty = 'n', cex = 1.2)
abline(v = loc, col = 'black', lty = 2)
dev.off()
