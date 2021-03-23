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
data <- raw_data
is_difference <- "FALSE"

# ### do the first difference
# data_diff <- apply(raw_data, c(1,3), diff)
# data <- aperm(data_diff, c(2,1,3))
# is_difference <- "TRUE"

####### potential model order ########
P <- 10

####### sample size ##########
sample_size <- 100
####### construct discount factor #########
delta <- seq(0.9, 0.99, by = 0.01)
delta_matrix <- as.matrix(expand.grid(delta, delta))

####### set up parameters for number of time series and number of time points ######
n_I <- dim(data)[1]
n_t <- dim(data)[2]

####### fit hierarchical dynamic linear model ########
sfInit(parallel = TRUE, cpus=n_I, type = "SOCK")
sfLibrary(PARCOR)
sfExport("data", "P", "delta_matrix", "sample_size")
ptm <- proc.time()
result <- sfLapply(1:n_I, function(x) hparcor(yt=data[x, , ], P=P, delta = delta_matrix, 
                                              sample_size=sample_size, chains=5, DIC = TRUE, uncertainty = TRUE))
consumed_time <- proc.time()-ptm
sfStop()
########################################################

#################################
##### retrieve results ##########
#################################
P_opt <- numeric(n_I)
est_sigma2 <- matrix(NA, nrow = n_I, ncol = n_t)
coef_mean <- rep(list(NA), n_I)
coef_mean_sample <- rep(list(NA), n_I)
for(i in 1:n_I){
  P_opt[i] <- which.min(result[[i]]$DIC_fwd)
  est_sigma2[i, ] <- rep(result[[i]]$sigma2t_fwd[n_t-P, P_opt[i]], n_t)
  #### transform estimated mean of PARCOR coefficients to estimated mean of TVAR coefficients
  coef_parcor_mean <- run_dl(phi_fwd = result[[i]]$mu_fwd[1, , , drop=FALSE],
                             phi_bwd = result[[i]]$mu_bwd[1, , , drop=FALSE])
  coef_mean[[i]] <- coef_parcor_mean[[P_opt[i]]]$forward
  
  
  #### transform sample of PARCOR coefficients to estimated mean of TVAR coefficients
  mu_fwd_sample <- result[[i]]$mu_fwd_sample[, 1, , ]
  mu_bwd_sample <- result[[i]]$mu_bwd_sample[, 1, , ]
  tmp <- lapply(1:sample_size, function(x) compute_TVAR_hier(phi_fwd = mu_fwd_sample[x, , ,drop=FALSE],
                                                             phi_bwd = mu_bwd_sample[x, , ,drop=FALSE],
                                                             P_opt = P_opt[i]))
  coef_mean_sample[[i]] <- simplify2array(tmp)
}

###########################################################

####### compute spectral density #########
####### set up span of frequency
w <- seq(0.001, 0.499, by = 0.001)
s_mean <- rep(list(NA), n_I)
s_mean_sample <- rep(list(NA), n_I)
s_mean_quantile <- rep(list(NA), n_I)
sfInit(parallel=TRUE, cpus=10, type="SOCK")
sfLibrary(PARCOR)
sfExport("w", "coef_mean_sample", "est_sigma2", "n_t", "P_opt")
for(i in 1:n_I){
  sfExport("i")
  tmp <- sfLapply(1:sample_size, function(x) cp_sd_uni(w=w,
                                                       phi=array(coef_mean_sample[[i]][1, , ,x], 
                                                                 dim = c(1, n_t, P_opt[i])),
                                                       sigma2=est_sigma2[i,]))
  s_mean_sample[[i]] <- simplify2array(tmp)[1, , , ]
  cat("\n Iterations of computing spectral density: ", i, "/", n_I, "\n")
  s_mean[[i]] <- cp_sd_uni(w=w,
                           phi=coef_mean[[i]],
                           sigma2=est_sigma2[i,])
  s_mean_quantile[[i]] <- apply(s_mean_sample[[i]], 1:2, quantile, c(0.025, 0.975))
}
sfStop()


sink(file=paste0(root_dir, "/results/", cond_type, "_S", subject_id, ".txt"))
cat("\n The first order difference:", is_difference, "\n")
cat("\n Optimal model order:", P_opt, "\n")
cat("\n Computation time: ", consumed_time, "\n")
cat("\n The discount factor range:", range(delta), "\n")
for(i in 1:length(result)){
  cat("\n Selected forward discount factor of area", cluster_area[i], ":", result[[i]]$best_delta_fwd, "\n")
  cat("\n Selected backward discount factor of area", cluster_area[i], ":", result[[i]]$best_delta_bwd, "\n")
}
sink()

save("s_mean", "s_mean_quantile", "is_difference", file = paste0(root_dir, "/results/", cond_type, "_S", subject_id, ".RData"))