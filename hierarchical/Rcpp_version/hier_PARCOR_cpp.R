dir <- '/Users/johnn/Documents/Research/Project3/hierarchical/Rcpp_version/'
setwd(dir)
library(Rcpp)
library(matrixcalc)
library(mvnfast)
sourceCpp("ffbs.cpp")
sourceCpp("PARCOR_to_AR.cpp")

# compute_pred_dens <- function(yt, ft, Qt, mt, Ct, nt, n_t, type, m, P){
#   ## dimension of yt: n_t * n_I
#   ## dimension of F1: n_t * n_I
#   ## dimension of fnt: n_I * 1
#   ## dimension of mnt: n_I * n_t
#   if(type == 1){
#     lbound <- P + 1
#     ubound <- n_t
#     sign <- 1
#   }else{
#     lbound <- 1
#     ubound <- n_t - P
#     sign <- -1
#   }
#   ll <- 0
#   for(i in lbound:ubound){
#     # fnt <- F1[i - sign*m, drop = FALSE] * mnt[, i]
#     # Qnt <- diag(F1[i - sign*m, ]) %*% Cnt[, , i] %*% diag(F1[i - sign*m, ])
#     # if(i == lbound){
#     #   ll <- tryCatch({ll + dmvt(X = yt[i, ], mu = ft[, i], 
#     #                             sigma = Qt[, , i], df = 1, log = TRUE)}, error = function(e) {
#     #                               cat("type: ", type, "\n")
#     #                               cat("iteration: ", i, "\n")
#     #                               print(Qt)
#     #                               print(F1[i-sign*m, ])
#     #                               # print(sigma2t[i])
#     #                               # print(sigma2t[ubound])
#     #                               print(Ct[, , i])})
#     # }else{
#       # ll <- tryCatch({ll + dmvt(X = yt[i, ], mu = ft[, i], 
#       #                           sigma = Qt[, , i], df = nt[i-1], log = TRUE)}, error = function(e) {
#       #                             cat("type: ", type, "\n")
#       #                             cat("iteration: ", i, "\n")
#       #                             print(Qt)
#       #                             print(F1[i-sign*m, ])
#       #                             # print(sigma2t[i])
#       #                             # print(sigma2t[ubound])
#       #                             print(Ct[, , i])})
#     # }
#     
#     ll <- tryCatch({ll + dmvn(X = yt[i, ], mu = ft[, i], sigma = Qt[, , i], log = TRUE)})
#   }
#   return(ll)
# }
# 
# 
# sample_ll <- function(yt, F1, F2, ft, Qt, mt, Ct, nt, n_t, type, m, 
#                       sample_size = 1000L, P){
#   n_I <- dim(yt)[2]
#   # obtain sample
#   if(type == 1){
#     lbound <- P + 1
#     ubound <- n_t
#     sign <- 1
#   }else{
#     lbound <- 1
#     ubound <- n_t - P
#     sign <- -1
#   }
#   ll <- numeric(sample_size)
#   for(i in lbound:ubound){
#     # if(i == lbound){
#     #   mt_sample <- rmvt(n = sample_size, mu = mt[, i], sigma = Ct[, , i], df = 1)
#     #   ft_sample <- t(F1[i-sign*m, drop=FALSE] * t(mt_sample))
#     #   ll_sample <- tryCatch({dmvt(X = ft_sample, mu = as.numeric(yt[i, ]),
#     #                               sigma = Qt[, , i], df = 1, log = TRUE)}, error = function(e) {
#     #                                 cat("type: ", type, "\n")
#     #                                 cat("iteration: ", i, "\n")
#     #                                 print(Qt)
#     #                                 print(F1[i-sign*m, ])
#     #                                 # print(sigma2t[i])
#     #                                 # print(sigma2t[ubound])
#     #                                 print(Ct[, , i])})
#     # }else{
#       # if(!is.positive.definite(Ct[, , i-1])){
#       #   Ct[, , i-1] <- 0.5*Ct[, , i-1] + 0.5*t(Ct[, , i-1])
#       # }
#       # mt_sample <- tryCatch(rmvt(n = sample_size, mu = mt[, i], 
#       #                       sigma = Ct[, , i], df = nt[i-1]), 
#       #                       error = function(e) {print(Ct[, , i])})
#       # ft_sample <- t(F1[i-sign*m, drop=FALSE] * t(mt_sample))
#       # ll_sample <- tryCatch({dmvt(X = ft_sample, mu = as.numeric(yt[i, ]), 
#       #                             sigma = Qt[, , i], df = nt[i-1], log = TRUE)}, error = function(e) {
#       #                               cat("type: ", type, "\n")
#       #                               cat("iteration: ", i, "\n")
#       #                               print(Qt)
#       #                               print(F1[i-sign*m, ])
#       #                               # print(sigma2t[i])
#       #                               # print(sigma2t[ubound])
#       #                               print(Ct[, , i])})
#     # }
#       mt_sample <- tryCatch(rmvn(n = sample_size, mu = mt[,i], sigma = Ct[, ,i]),
#                             error = function(e){cat("wrong in sampling \n")
#                                                 print(Ct[, , i])})
#       at_sample <- F2%*%t(mt_sample)
#       ft_sample <- diag(F1[i - sign*m, ])%*%at_sample
#       ll_sample <- tryCatch({dmvn(X = ft_sample, mu = as.numeric(yt[i, ]),
#                                   sigma = Qt[, , i], log = TRUE)}, error = function(e){
#                                     cat("Wrong in computing density \n")
#                                     print(Qt[, , i])})
#     
#     ll <- ll + ll_sample
#     cat("iterations of sampling: ", i, "/", ubound - lbound+1, "\r")
#   }
#   return(mean(ll))
# }
# 
# 
# compute_pDIC <- function(yt, F1, ft, Qt, mt, Ct, nt, n_t, type, m, 
#                          sample_size, DIC = TRUE, P){
#   ll <- compute_pred_dens(yt = yt, ft = ft, Qt = Qt, mt = mt, Ct = Ct,
#                           nt = nt, n_t = n_t, type = type, m = m, P=P)
#   if(DIC){
#     pDIC_sample <- sample_ll(yt = yt, F1 = F1, ft = ft, Qt = Qt, mt = mt,
#                              Ct = Ct, nt = nt, n_t = n_t, 
#                              type = type, m = m, sample_size = sample_size, P = P)
#   }else{
#     pDIC_sample <- 0
#   }
# 
#   
#   pDIC <- 2*(ll - pDIC_sample)
#   return(list(ll = ll, pDIC = pDIC))
# }

# hier_parcor <- function(yt, delta, P, F2, sample_size = 1000L, DIC = TRUE){
#   ## dimension of yt: n_t * n_I
#   ## dimension of delta: num_delta * 2
#   n_t <- dim(yt)[1]
#   I <- dim(yt)[2]
#   resid_fwd <- array(0, dim = c(n_t, I, P+1))
#   resid_bwd <- array(0, dim = c(n_t, I, P+1))
#   resid_fwd[, , 1] <- yt
#   resid_bwd[, , 1] <- yt
#   phi_fwd <- array(0, dim = c(I, n_t, P))
#   phi_bwd <- array(0, dim = c(I, n_t, P))
#   mu_fwd <- array(0, dim = c(I, n_t, P))
#   mu_bwd <- array(0, dim = c(I, n_t, P))
#   sigma2t_fwd <- matrix(0, nrow = n_t, ncol = P)
#   sigma2t_bwd <- matrix(0, nrow = n_t, ncol = P)
#   num_delta <- dim(delta)[1]
#   best_delta_fwd <- matrix(0, nrow = P, ncol = 2)
#   best_delta_bwd <- matrix(0, nrow = P, ncol = 2)
#   best_pred_dens_fwd <- numeric(P)
#   best_pred_dens_bwd <- numeric(P)
#   pDIC_fwd <- numeric(P)
#   pDIC_bwd <- numeric(P)
#   for(j in 1:P){
#     for(i in 1:num_delta){
#       tmp_fwd <- forward_filter_backward_smooth(yt = t(resid_fwd[, , j]), 
#                                                 F1 = t(resid_bwd[, , j]), F2 = F2, 
#                                                 delta1 = delta[i, 1], delta2 = delta[i, 2], 
#                                                 n_t = n_t, I = I, m = j, type = 1, P = P)
#       tmp_bwd <- forward_filter_backward_smooth(yt = t(resid_bwd[, , j]), 
#                                                 F1 = t(resid_fwd[, , j]), F2 = F2, 
#                                                 delta1 = delta[i, 1], delta2 = delta[i, 2], 
#                                                 n_t = n_t, I = I, m = j, type = 0, P = P)
#       tmp_DIC_fwd <- compute_pDIC(yt = resid_fwd[, , j], F1 = resid_bwd[, , j], 
#                                   ft = tmp_fwd$ft, Qt = tmp_fwd$Qt, mt = tmp_fwd$akt,
#                                   Ct = tmp_fwd$Rkt, nt = tmp_fwd$nt,
#                                   n_t = n_t, type = 1, m = j, DIC = DIC,
#                                   #sigma2t = tmp_fwd$sigma2t, 
#                                   sample_size = sample_size, P = P)
#       tmp_DIC_bwd <- compute_pDIC(yt = resid_bwd[, , j], F1 = resid_fwd[, , j], 
#                                   ft = tmp_bwd$ft, Qt = tmp_bwd$Qt, mt = tmp_bwd$akt,
#                                   Ct = tmp_bwd$Rkt, nt = tmp_bwd$nt, 
#                                   n_t = n_t, type = 0, m = j, DIC = DIC,
#                                   #sigma2t = tmp_bwd$sigma2t,
#                                   sample_size = sample_size, P = P)
#       if(i == 1){
#         pDIC_fwd[j] <- tmp_DIC_fwd$pDIC
#         pDIC_bwd[j] <- tmp_DIC_bwd$pDIC
#         best_pred_dens_fwd[j] <- tmp_DIC_fwd$ll
#         best_pred_dens_bwd[j] <- tmp_DIC_bwd$ll
#         
#         best_fwd <- tmp_fwd
#         best_bwd <- tmp_bwd
#         best_delta_fwd[j, ] <- delta[i, ]
#         best_delta_bwd[j, ] <- delta[i, ]
#         #best_pred_dens_fwd[j] <- tmp_DIC_fwd$ll
#         #pDIC_fwd[j] <- tmp_DIC_fwd$pDIC
#         
#         #best_pred_dens_bwd[j] <- tmp_DIC_bwd$ll
#         #pDIC_bwd[j] <- tmp_DIC_bwd$pDIC
#       }else{
#         if(best_pred_dens_fwd[j] < tmp_DIC_fwd$ll){
#           best_fwd <- tmp_fwd
#           best_delta_fwd[j, ] <- delta[i, ]
#           best_pred_dens_fwd[j] <- tmp_DIC_fwd$ll
#           pDIC_fwd[j] <- tmp_DIC_fwd$pDIC
#         }
#         if(best_pred_dens_bwd[j] < tmp_DIC_bwd$ll){
#           best_bwd <- tmp_bwd
#           best_delta_bwd[j, ] <- delta[i, ]
#           best_pred_dens_bwd[j] <- tmp_DIC_bwd$ll
#           pDIC_bwd[j] <- tmp_DIC_bwd$pDIC
#         }
#       }
#     }
#     
#     #browser()
#     resid_fwd[, , j+1] <- t(best_fwd$residuals) 
#     resid_bwd[, , j+1] <- t(best_bwd$residuals)
#     phi_fwd[, , j] <- best_fwd$mnt
#     phi_bwd[, , j] <- best_bwd$mnt
#     mu_fwd[, , j] <- best_fwd$mnkt
#     mu_bwd[, , j] <- best_bwd$mnkt
#     sigma2t_fwd[, j] <- best_fwd$sigma2t
#     sigma2t_bwd[, j] <- best_bwd$sigma2t
#     DIC_fwd <- 2*(cumsum(pDIC_fwd) - best_pred_dens_fwd)
#     DIC_bwd <- 2*(cumsum(pDIC_bwd) - best_pred_dens_bwd)
#     cat('\niter: ', j)
#   }
#   return(list(resid_fwd = resid_fwd, resid_bwd = resid_bwd,
#               phi_fwd = phi_fwd, phi_bwd = phi_bwd, 
#               mu_fwd = mu_fwd, mu_bwd = mu_bwd,
#               sigma2t_fwd = sigma2t_fwd, sigma2t_bwd = sigma2t_bwd,
#               best_delta_fwd = best_delta_fwd,
#               best_delta_bwd = best_delta_bwd,
#               best_pred_dens_fwd = best_pred_dens_fwd,
#               best_pred_dens_bwd = best_pred_dens_bwd,
#               pDIC_fwd = pDIC_fwd,
#               pDIC_bwd = pDIC_bwd,
#               DIC_fwd = DIC_fwd,
#               DIC_bwd = DIC_bwd))
# }



hparcor <- function(yt, delta, P, F2, sample_size = 1000L, chains, DIC = TRUE, uncertainty = TRUE){
  ## dimension of yt: n_t * n_I
  ## dimension of delta: num_delta * 2
  n_t <- dim(yt)[1]
  I <- dim(yt)[2]
  resid_fwd <- array(0, dim = c(n_t, I, P+1))
  resid_bwd <- array(0, dim = c(n_t, I, P+1))
  resid_fwd[, , 1] <- yt
  resid_bwd[, , 1] <- yt
  phi_fwd <- array(0, dim = c(I, n_t, P))
  phi_bwd <- array(0, dim = c(I, n_t, P))
  mu_fwd <- array(0, dim = c(I, n_t, P))
  mu_bwd <- array(0, dim = c(I, n_t, P))
  sigma2t_fwd <- matrix(0, nrow = n_t, ncol = P)
  sigma2t_bwd <- matrix(0, nrow = n_t, ncol = P)
  num_delta <- dim(delta)[1]
  best_delta_fwd <- matrix(0, nrow = P, ncol = 2)
  best_delta_bwd <- matrix(0, nrow = P, ncol = 2)
  best_pred_dens_fwd <- numeric(P)
  best_pred_dens_bwd <- numeric(P)
  pDIC_fwd <- numeric(P)
  pDIC_bwd <- numeric(P)
  
  phi_fwd_sample <- array(NA, dim = c(sample_size, I, n_t, P))
  phi_bwd_sample <- array(NA, dim = c(sample_size, I, n_t, P))
  mu_fwd_sample <- array(NA, dim = c(sample_size, I, n_t, P))
  mu_bwd_sample <- array(NA, dim = c(sample_size, I, n_t, P))
  for(j in 1:P){
    best_fwd <- ffbs_DIC(yt = t(resid_fwd[, , j]), 
                         F1 = t(resid_bwd[, , j]), F2 = F2, 
                         n_t = n_t, I = I, m = j, type = 1, P = P,
                         delta = delta, DIC = DIC, sample_size = sample_size,
                         chains = chains, uncertainty=TRUE)
    best_bwd <- ffbs_DIC(yt = t(resid_bwd[, , j]), 
                         F1 = t(resid_fwd[, , j]), F2 = F2, 
                         n_t = n_t, I = I, m = j, type = 0, P = P,
                         delta = delta, DIC = DIC, sample_size = sample_size,
                         chains = chains, uncertainty=TRUE)
    
    phi_fwd_sample[, , , j] <- best_fwd$mnt_sample
    phi_bwd_sample[, , , j] <- best_bwd$mnt_sample
    mu_fwd_sample[, , , j] <- best_fwd$mnkt_sample
    mu_bwd_sample[, , , j] <- best_bwd$mnkt_sample
    
    pDIC_fwd[j] <- best_fwd$pDIC
    pDIC_bwd[j] <- best_bwd$pDIC
    best_delta_fwd[j, ] <- best_fwd$delta
    best_delta_bwd[j, ] <- best_bwd$delta
    best_pred_dens_fwd[j] <- best_fwd$ll
    best_pred_dens_bwd[j] <- best_bwd$ll
    
    ### obtain (j+1) residuals
    resid_fwd[, , j+1] <- t(best_fwd$residuals) 
    resid_bwd[, , j+1] <- t(best_bwd$residuals)
    phi_fwd[, , j] <- best_fwd$mnt
    phi_bwd[, , j] <- best_bwd$mnt
    mu_fwd[, , j] <- best_fwd$mnkt
    mu_bwd[, , j] <- best_bwd$mnkt
    sigma2t_fwd[, j] <- best_fwd$sigma2t
    sigma2t_bwd[, j] <- best_bwd$sigma2t
    DIC_fwd <- 2*(cumsum(pDIC_fwd) - best_pred_dens_fwd)
    DIC_bwd <- 2*(cumsum(pDIC_bwd) - best_pred_dens_bwd)
    cat('\niter: ', j)
    
    
  }
  return(list(resid_fwd = resid_fwd, resid_bwd = resid_bwd,
              phi_fwd = phi_fwd, phi_bwd = phi_bwd, 
              mu_fwd = mu_fwd, mu_bwd = mu_bwd,
              sigma2t_fwd = sigma2t_fwd, sigma2t_bwd = sigma2t_bwd,
              best_delta_fwd = best_delta_fwd,
              best_delta_bwd = best_delta_bwd,
              best_pred_dens_fwd = best_pred_dens_fwd,
              best_pred_dens_bwd = best_pred_dens_bwd,
              pDIC_fwd = pDIC_fwd,
              pDIC_bwd = pDIC_bwd,
              DIC_fwd = DIC_fwd,
              DIC_bwd = DIC_bwd,
              phi_fwd_sample = phi_fwd_sample,
              phi_bwd_sample = phi_bwd_sample,
              mu_fwd_sample = mu_fwd_sample,
              mu_bwd_sample = mu_bwd_sample))
}


compute_TVAR <- function(phi_fwd, phi_bwd, P_opt){
  result <- PAR_to_AR_fun(phi_fwd = phi_fwd, 
                          phi_bwd = phi_bwd)
  return(result[[P_opt]]$forward)
}

