# dir1 <- '/soe/wjzhao/project/Project3/hierarchical/Rcpp_version/'
# setwd(dir1)
plot_dir <- "/soe/wjzhao/project/Project3/hierarchical/plots/Rcpp/"
# source(paste0(dir1, 'hier_PARCOR_cpp.R'))
# source(paste0(dir1, "draw_density_RcppVer.R"))


library(matrixcalc)
library(PARCOR)
library(snowfall)

########################################
##### simulation 4
########################################
##################
## draw time series
##################
draw_ts <- function(yt, ind, ...){
  par(cex = 1.5)
  plot(yt, type = 'l', xlab = 'time', ylab = 'value',
       main = bquote("ts"[.(ind)]), ...)
}




gen.sim <- function(n_t, n_I, P=2, vt=0.1, et){
  yt <- matrix(nrow = n_t, ncol = n_I)
  mu <- matrix(nrow = n_t, ncol = P)
  #betat <- array(NA, dim = c(n_I, n_t, P))
  phit <- array(NA, dim=c(n_I, n_t, P))
  lambdat <- matrix(NA, nrow = n_I, ncol = n_t)
  delta_i <- c(0, 0, 0, 1, 5)
  phit[, , 2] <- -0.9^2
  rt <- 0.1*(1:n_t)/n_t + 0.85
  y_init <- matrix(0, nrow = P, ncol = n_I)
  yt_tmp <- rbind(y_init, yt)
  ## generate the common part lambda_mu
  lambdat_mu <- 15*(1:n_t)/n_t + 5

  ## generate the ar coefficients phit
  for(ts in 1:n_I){
    for(i in 1:n_t){
      if(i == 1){
        lambdat[ts, i] <- delta_i[ts] + rnorm(1, mean = 0, sd = vt)
      }else{
        lambdat[ts, i] <- lambdat_mu[i-1] + delta_i[ts] + rnorm(n = 1, mean = 0, sd = vt)
      }
    }
  }

  for(ts in 1:n_I){
    phit[ts, , 1] <- rt*cos(2*pi/lambdat[ts, ])

  }

  ## genrate the data sets
  for(ts in 1:n_I){
    for(i in 1:n_t+P){
      yt_tmp[i, ts] <- t(as.matrix(yt_tmp[i-(1:P), ts]))%*%as.matrix(phit[ts, i-P, ]) + rnorm(n=1, sd = et)
    }
  }
  return(list(yt = yt_tmp[(P+1):(n_t+P), ],
              phit=phit,
              lambdat=lambdat,
              lambdat_mu = lambdat_mu))
}

########################
#### generate simulation data
########################

set.seed(1234)
n_t <- 1024
n_I <- 5
P <- 2
#wt <- 0.01
vt <- 0.5
et <- .8
sim <- gen.sim(n_t = n_t, n_I = n_I, P = P, vt = vt, et = et)
yt <- sim$yt
true_ar <- sim$phit
true_ar_mean <- array(cbind((0.1*(1:n_t)/n_t + 0.85)*cos(2*pi/sim$lambdat_mu), -0.9^2), dim = c(1, n_t, 2))

#################################
#### draw time series and
#################################
###
layout(matrix(1:(n_I+1), nrow = 2, ncol = 3, byrow = TRUE))
for(i in 1:n_I){
  draw_ts(yt[, i], ind = i)
}

#########################
###### fit model #######
#########################
F2t <- matrix(1, nrow = n_I, ncol = n_I)
F2t[1:4, 2:5] <- diag(4)
F2t[5, ] <- -1*F2t[5, ]
F2t[5,1] <- 1

delta <- seq(0.99, 0.999, 0.002)
delta_matrix <- as.matrix(expand.grid(delta, delta))
# delta_matrix <- cbind(replicate(2, delta_matrix_tmp[, 1]),
#                       delta_matrix_tmp)
# result_parcor <- hier_parcor(yt = yt, delta = delta_matrix,
#                              P = 5,
#                              F2 = F2t)

sample_size <- 1000
ptm <- proc.time()
result_parcor <- hparcor(yt = yt, delta = delta_matrix,
                         P = 5, F2 = F2t, sample_size = sample_size,
                         chains = 10, DIC = TRUE, uncertainty = TRUE)
print(proc.time() - ptm)
###########################
##### Optimal model order
###########################
P_opt <- which.min(result_parcor$DIC_fwd)
cat("Optimal model order: ", P_opt)
print(result_parcor$DIC_fwd)
print(result_parcor$best_pred_dens_fwd)
###########################
##### selected discount factor
###########################
print(result_parcor$best_delta_fwd)


###########################
#### estimated innovation variance
###########################
sigma2 <- rep(result_parcor$sigma2t_fwd[n_t-P, P_opt], n_t)

## compute ar coefficients
coef_parcor <- run_dl(phi_fwd = result_parcor$phi_fwd,
                      phi_bwd = result_parcor$phi_bwd)


coef <- coef_parcor[[P_opt]]$forward

## compute ar coefficients
coef_parcor_mean <- run_dl(phi_fwd = result_parcor$mu_fwd,
                           phi_bwd = result_parcor$mu_bwd)

coef_mean <- coef_parcor_mean[[P_opt]]$forward


##
coef_sample <- lapply(1:sample_size, function(x) compute_TVAR_hier(phi_fwd = result_parcor$phi_fwd_sample[x, , ,],
                                                                   phi_bwd = result_parcor$phi_bwd_sample[x, , ,],
                                                                   P_opt = P_opt))
coef_sample <- simplify2array(coef_sample)
coef_mean_sample <- lapply(1:sample_size, function(x) compute_TVAR_hier(phi_fwd = result_parcor$mu_fwd_sample[x, , ,],
                                                                        phi_bwd = result_parcor$mu_bwd_sample[x, , ,],
                                                                        P_opt = P_opt))
coef_mean_sample <- simplify2array(coef_mean_sample)


coef_quantile <- apply(coef_sample, 1:3, quantile, c(0.025, 0.5, 0.975))


## compute spectral density
## span of frequence
w <- seq(0.001, 0.499, by = 0.001)

### spectral density of each time series
s <- cp_sd_uni(w = w,
                phi = coef,
                sigma2 = sigma2)

### average of spectral density
s_mean <- cp_sd_uni(w=w,
                    phi= coef_mean,
                    sigma2 = sigma2)

### True spectral density
s_true <- cp_sd_uni(phi=true_ar, sigma2 = rep(et, n_t), w=w)
s_mean_true <- cp_sd_uni(phi=true_ar_mean, sigma2=rep(et,n_t), w=w)


####################################
#### sample spectral density
####################################
sfInit(parallel = TRUE, cpus=10, type="SOCK")
sfLibrary(PARCOR)
sfExport("w", "coef_sample", "sigma2")
s_sample <- sfLapply(1:(sample_size), function(x) cp_sd_uni(w=w,
                                                            phi = coef_sample[, , , x],
                                                            sigma2 = sigma2))
sfStop()

s_sample <- simplify2array(s_sample)

sfInit(parallel = TRUE, cpus=10, type="SOCK")
sfLibrary(PARCOR)
sfExport("w", "coef_mean_sample", "sigma2", "n_t", "P_opt")
s_mean_sample <- sfLapply(1:(sample_size), function(x) cp_sd_uni(w=w,
                                                            phi = array(coef_mean_sample[1, , , x],
                                                                        dim = c(1, n_t, P_opt)),
                                                            sigma2 = sigma2))
sfStop()

s_mean_sample <- simplify2array(s_mean_sample)



### draw ar coefficients

sim_index <- "simulation_4"
index <- 1:(n_t-2*P)
png(paste0(plot_dir, sim_index, "/scale/ar_coef.png"), width = 2500, height = 2000)
layout(matrix(1:10, 5, 2, byrow = TRUE))
layout.show(n = 10)
par(cex.lab = 2.5, cex.axis = 2.5, cex.main = 3)
for(i in 1:n_I){
  for(j in 1:P_opt){
    plot(coef[i, (P+1):(n_t-P), j], xlab = "time", ylab = "value",
         type = 'l', main = bquote(phi[.(i)*.(j)]),
         ylim = range(coef[i, (P+1):(n_t-P), j],
                      coef_mean[1, (P+1):(n_t-P), j],
                      true_ar[i, (P+1):(n_t-P), j],
                      coef_quantile[, i, (P+1):(n_t-P), j],
                      na.rm = TRUE))

    polygon(c(index,
              rev(index)),
            c(coef_quantile[3, i, (P+1):(n_t-P), j],
              rev(coef_quantile[1, i, (P+1):(n_t-P), j])),
            col="skyblue", border = NA)
    lines(coef[i, (P+1):(n_t-P), j], col = 'black', lty = 1)
    lines(true_ar[i, (P+1):(n_t-P), j], type = 'l', col = 'red', lty = 3)
    lines(coef_mean[1, (P+1):(n_t-P), j], type = 'l',
          col = 'blue', lty = 2)

    if(j == 1 | i==1 | i==3){
      loc <- "bottomright"
    }else{
      loc <- "topright"
    }
    legend(loc, legend=c("true", "estimated", "mean"),
           lty = c(3,1,2), col = c("red", "black", "blue"),
           cex = 3)
  }
}
dev.off()


png(paste0(plot_dir, sim_index, "/scale/BLF_scree.png"), width = 1500, height = 1400)
layout(matrix(1, 1, 1))
par(cex = 2)
plot(result_parcor$best_pred_dens_fwd, type = 'l', xlab = 'P', ylab = 'log-likelihood')
dev.off()

# png(paste0(plot_dir, sim_index, "/scale/ar_coef_mean.png"), width = 1500, height = 1400)
# layout(matrix(1:4, 2, 2))
# layout.show(n = 4)
# par(cex.lab = 3, cex.axis = 2.5, cex.main = 3)
# for(i in 1:2){
#   for(j in 1:P){
#     plot(coef_mean[1, (P+1):(n_t-P), j], xlab = "time", ylab = "value",
#          type = 'l', main = bquote(phi[.(i)*.(j)]),
#          ylim = range(coef, true_ar, na.rm = TRUE))
#     lines(true_ar[i, (P+1):(n_t-P), j], type = 'l', col = 'red')
#   }
# }
# dev.off()


######################################################
### compute 95% credible interval of spectral density
######################################################
s_quantile <- rep(list(NA), n_I)
for(i in 1:n_I){
  s_quantile[[i]] <- apply(s_sample[i, , , ], 1:2, quantile, c(0.025, 0.975))
}

s_mean_quantile <- apply(s_mean_sample[1, , , ], 1:2, quantile, c(0.025, 0.975))

#####################################
## draw spectral density plots
#####################################
#### compute spectral density #####
#### for time series for all 5 time series
zlim <- range(s, s_quantile, s_mean_quantile, s_true, s_mean_true)
for(i in 1:n_I){
  index <- i
  # zlim <- range(s[index, (P+1):(n_t-P),],
  #               s_quantile[, (P+1):(n_t-P),],
  #               s_true[index, (P+1):(n_t-P), ])
  png(filename = paste0(plot_dir, sim_index, '/scale/est_', index, 'st.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier(w = w, index = index, P = P,
                    n_t = n_t, s = s[index, , ], zlim = zlim)
  dev.off()


  png(filename = paste0(plot_dir, sim_index, '/scale/true_', index, 'st.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier(w = w, index = index, P = P,
                    n_t = n_t, s = s_true[index, , ], zlim = zlim)
  dev.off()


  png(filename = paste0(plot_dir, sim_index, '/scale/ql_', index, 'st.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier(w = w, index = index, P = P,
                    n_t = n_t, s = s_quantile[1, , ], zlim = zlim)
  dev.off()


  png(filename = paste0(plot_dir, sim_index, '/scale/qu_', index, 'st.png'))
  par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  draw_density_hier(w = w, index = index, P = P,
               n_t = n_t, s = s_quantile[2, , ], zlim = zlim)
  dev.off()
}




#### for baseline of all five time series
index <- 1
# zlim <- range(s_mean[1, ,], s_mean_true[1, , ])
png(filename = paste0(plot_dir, sim_index, '/scale/est_', index, 'mean.png'))
par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
draw_density_hier(w = w, index = index, P = P,
                  n_t = n_t, s = s_mean[1, , ], zlim = zlim)
dev.off()


png(filename = paste0(plot_dir, sim_index, '/scale/true_', index, 'mean.png'))
par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
draw_density_hier(w = w, index = index, P = P,
                  n_t = n_t, s = s_mean_true[1, , ], zlim = zlim)
dev.off()




#### for baseline of all five time series
index <- 1
# zlim <- range(s_mean[1, ,], s_mean_true[1, , ])
png(filename = paste0(plot_dir, sim_index, '/scale/est_lb_', index, 'mean.png'))
par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
draw_density_hier(w = w, index = index, P = P,
                  n_t = n_t, s = s_mean_quantile[1, , ], zlim = zlim)
dev.off()

#### for baseline of all five time series
index <- 1
# zlim <- range(s_mean[1, ,], s_mean_true[1, , ])
png(filename = paste0(plot_dir, sim_index, '/scale/est_ub_', index, 'mean.png'))
par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
draw_density_hier(w = w, index = index, P = P,
                  n_t = n_t, s = s_mean_quantile[2, , ], zlim = zlim)
dev.off()
