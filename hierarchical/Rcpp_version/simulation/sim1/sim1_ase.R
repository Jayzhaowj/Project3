plot_dir <- "/Users/johnn/Documents/Research/Project3/hierarchical/plots/Rcpp/simulation_4/scale/"


library(matrixcalc)
library(PARCOR)

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
vt <- 0.5
et <- 0.8
n_sim <- 50
sim <- rep(list(NA), n_sim)
yt <- rep(list(NA), n_sim)
true_ar <- rep(list(NA), n_sim)
true_ar_mean <- rep(list(NA), n_sim)
for(i in 1:n_sim){
  sim[[i]] <- gen.sim(n_t=n_t, n_I=n_I, P=P, vt=vt, et=et)
  yt[[i]] <- sim[[i]]$yt
  true_ar[[i]] <- sim[[i]]$phit
  true_ar_mean[[i]] <- (0.1*(1:n_t)/n_t + 0.85)*cos(2*pi/sim[[i]]$lambdat_mu)
}

#######
##
######
mse <- matrix(NA, nrow = n_I, ncol = n_sim)
## compute spectral density
## span of frequence
w <- seq(0.001, 0.499, by = 0.001)


#########################
###### fit model #######
#########################

delta <- seq(0.99, 0.999, 0.001)
delta_matrix <- as.matrix(expand.grid(delta, delta))
ll <- matrix(NA, nrow = n_I, ncol = n_sim)
P_opt <- numeric(n_sim)

sample_size <- 1000
ptm <- proc.time()

for(i in 1:n_sim){
  result_parcor <- hparcor(yt = yt[[i]], delta = delta_matrix,
                           P = 5, sample_size = sample_size,
                           chains = 10, DIC = TRUE, uncertainty = FALSE)
  ll[, i] <- result_parcor$best_pred_dens_fwd
  P_opt[i] <- which.min(result_parcor$DIC_fwd)
  sigma2 <- rep(result_parcor$sigma2t_fwd[n_t-P, P_opt[i]], n_t)
  coef_parcor <- run_dl(phi_fwd = result_parcor$phi_fwd,
                        phi_bwd = result_parcor$phi_bwd)


  coef <- coef_parcor[[P_opt[i]]]$forward

  ## compute ar coefficients
  coef_parcor_mean <- run_dl(phi_fwd = result_parcor$mu_fwd,
                             phi_bwd = result_parcor$mu_bwd)

  coef_mean <- coef_parcor_mean[[P_opt[i]]]$forward

  ### spectral density of each time series
  s <- cp_sd_uni(w = w,
                 phi = coef,
                 sigma2 = sigma2)
  ### True spectral density
  s_true <- cp_sd_uni(phi=true_ar[[i]], sigma2 = rep(et, n_t), w=w)
  mse[, i] <- apply((s-s_true)^2, 1, mean)
  cat("\n iterations of simulation: ", i, "/", n_sim, "\n")
}

print(proc.time() - ptm)


##### draw log-likelihood
png(paste0(plot_dir, 'll_sim1.png'))
par(par(mar=c(5,6,4,1)+.1), cex.lab = 2, cex.axis = 1.5, cex.main = 2.5)
plot(ll[, 1], type = 'l', xlab = 'model order',
     ylab = 'log(likelihood)', ylim = range(ll),
     main = "Simulation 1", cex = 2)
#     main = expression(paste(phi["1, 1, 2, t"], " = -0.8")), cex = 2)
for(i in 2L:50L){
  lines(ll[, i], type = 'l', col = i)
}
dev.off()


