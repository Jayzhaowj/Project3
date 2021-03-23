# dir1 <- '/soe/wjzhao/project/Project3/hierarchical/Rcpp_version/'
# plot_dir <- "/soe/wjzhao/project/Project3/hierarchical/eeg/plots/"
# plot_dir <- "/Users/johnn/Documents/Research/Project3/hierarchical/eeg/plots/new/"
# setwd(dir1)
# source(paste0(dir1, 'hier_PARCOR_cpp.R'))
# source(paste0(dir1, "draw_density_RcppVer.R"))
library(PARCOR)
library(snowfall)




################# EEG #####################################
########       Read data       ########

### set up data directory
#mydir <- "/Users/jay/Documents/Statistics/Papers/Research/log/11:30:17/data/"
mydir <- "/soe/wjzhao/project/project1/EEG/data/"
#mydir <- "C:/Users/johnn/Desktop/EEG/data/"

label <- c("F3", "C3", "P3", "Fz", "Cz", "Pz", "F4", "C4", "P4")
dlmd7_data <- scan(paste0(mydir, 'dlmd7.dat'))  ## F3
dlmd8_data <- scan(paste0(mydir, 'dlmd8.dat'))  ## C3
dlmd9_data <- scan(paste0(mydir, 'dlmd9.dat'))  ## P3
dlmd11_data <- scan(paste0(mydir, 'dlmd11.dat')) ## Fz
dlmd12_data <- scan(paste0(mydir, 'dlmd12.dat')) ## Cz
dlmd13_data <- scan(paste0(mydir, 'dlmd13.dat')) ## Pz
dlmd15_data <- scan(paste0(mydir, 'dlmd15.dat')) ## F4
dlmd16_data <- scan(paste0(mydir, 'dlmd16.dat')) ## C4
dlmd17_data <- scan(paste0(mydir, 'dlmd17.dat')) ## P4


dlmd7_data <- dlmd7_data - mean(dlmd7_data)
dlmd8_data <- dlmd8_data - mean(dlmd8_data)
dlmd9_data <- dlmd9_data - mean(dlmd9_data)
dlmd11_data <- dlmd11_data - mean(dlmd11_data)
dlmd12_data <- dlmd12_data - mean(dlmd12_data)
dlmd13_data <- dlmd13_data - mean(dlmd13_data)
dlmd15_data <- dlmd15_data - mean(dlmd15_data)
dlmd16_data <- dlmd16_data - mean(dlmd16_data)
dlmd17_data <- dlmd17_data - mean(dlmd17_data)


##### generate dataset ######
data_all <- t(rbind(dlmd7_data, dlmd8_data, dlmd9_data,
                  dlmd11_data, dlmd12_data, dlmd13_data,
                  dlmd15_data, dlmd16_data, dlmd17_data))

sample_size <- 1000
n_t <- 3600
n_I <- 9
P <- 15

delta <- seq(0.99, 0.999, by = 0.001)
delta_matrix <- as.matrix(expand.grid(delta, delta))
result_parcor <- hparcor(yt = data_all, P = P,
                         delta = delta_matrix, sample_size = sample_size,
                         chains = 1, DIC = TRUE, uncertainty = TRUE)

##########################################
##### optimal model order ############
##########################################
P_opt <- which.min(result_parcor$DIC_fwd)
cat("Optimal model order: ", P_opt)

##########################################
##### estimation of sigma2   ############
##########################################
est_sigma2 <- rep(result_parcor$sigma2t_fwd[n_t-P, P_opt], n_t)

###################################
#### discount factor
###################################
print(result_parcor$best_delta_fwd)

## compute ar coefficients
coef_parcor <- run_dl(phi_fwd = result_parcor$phi_fwd,
                             phi_bwd = result_parcor$phi_bwd)
coef <- coef_parcor[[P_opt]]$forward

coef_parcor_mean <- run_dl(phi_fwd = result_parcor$mu_fwd,
                                  phi_bwd = result_parcor$mu_bwd)

coef_mean <- coef_parcor_mean[[P_opt]]$forward

##############################################
##### 95% credible interval
##############################################
coef_sample <- lapply(1:sample_size, function(x) compute_TVAR_hier(phi_fwd = result_parcor$phi_fwd_sample[x, , ,],
                                                              phi_bwd = result_parcor$phi_bwd_sample[x, , ,],
                                                              P_opt = P_opt))
coef_sample <- simplify2array(coef_sample)
coef_mean_sample <- lapply(1:sample_size, function(x) compute_TVAR_hier(phi_fwd = result_parcor$mu_fwd_sample[x, , ,],
                                                                   phi_bwd = result_parcor$mu_bwd_sample[x, , ,],
                                                                   P_opt = P_opt))
coef_mean_sample <- simplify2array(coef_mean_sample)
coef_mean_sample <- coef_mean_sample[1, , , , drop = FALSE]

## compute spectral density
## span of frequence
w <- seq(0.001, 0.499, by = 0.001)
sfInit(parallel = TRUE, cpus=10, type="SOCK")
sfLibrary(PARCOR)
sfExport("w", "coef_sample", "est_sigma2")
s_sample <- sfLapply(1:(sample_size), function(x) cp_sd_uni(w=w,
                                                            phi = coef_sample[, , , x],
                                                            sigma2 = est_sigma2))
sfStop()

s_sample <- simplify2array(s_sample)
####################
sfInit(parallel = TRUE, cpus=10, type="SOCK")
sfLibrary(PARCOR)
sfExport("w", "coef_mean_sample", "est_sigma2", "n_t", "P_opt")
s_mean_sample <- sfLapply(1:(sample_size), function(x) cp_sd_uni(w=w,
                                                                 phi=array(coef_mean_sample[1, , , x],
                                                                           dim = c(1, n_t, P_opt)),
                                                                 sigma2=est_sigma2))
sfStop()
s_mean_sample <- simplify2array(s_mean_sample)
### spectral density of each time series
s <- cp_sd_uni(w = w,
                phi = coef,
                sigma2 = est_sigma2)

### average of spectral density
s_mean <- cp_sd_uni(w=w,
                     phi= coef_mean,
                     sigma2 = est_sigma2)

######################################################
### compute 95% credible interval of spectral density
######################################################
s_quantile <- rep(list(NA), n_I)
for(i in 1:n_I){
  s_quantile[[i]] <- apply(s_sample[i, , , ], 1:2, quantile, c(0.025, 0.975))
  cat("\n iterations: ", i, "/", n_I)
}

s_mean_quantile <- apply(s_mean_sample[1, , , ], 1:2, quantile, c(0.025, 0.975))





### draw ar coefficients
# png(paste0(plot_dir, "ar_coef.png"), width = 1500, height = 1400)
# layout(matrix(1:4, 2, 2))
# layout.show(n = 4)
# par(cex.lab = 3, cex.axis = 2.5, cex.main = 3)
# for(i in 1:2){
#   for(j in 1:P){
#     plot(coef[i, (P+1):(n_t-P), j], xlab = "time", ylab = "value",
#          type = 'l', main = bquote(phi[.(i)*.(j)]),
#          ylim = range(coef, true_ar, na.rm = TRUE))
#     lines(true_ar[i, (P+1):(n_t-P), j], type = 'l', col = 'red')
#   }
# }
# dev.off()




