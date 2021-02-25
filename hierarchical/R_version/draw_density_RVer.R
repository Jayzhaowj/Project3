## calculate the spectral density, output log
cal.tfr <- function(phi, sigma2, w, n_t, P){
  exp_part <- matrix(exp(-2i*pi*w*(1:P)), ncol = 1)
  s <- sapply(1:n_t, function(x) sigma2[x]/Mod(1 - t(phi[, x, drop = FALSE]) %*% exp_part)^2)
  return(log(s))
}

## draw spectral density graph
draw.density <- function(phi, sigma2, interval = 0.001, ...){
  w <- seq(0.0001, 0.4999, by = interval)
  P <- dim(phi)[1]
  n_t <- dim(phi)[2]
  s <- matrix(NA, nrow = n_t, ncol = length(w))
  for(i in 1:length(w)){
    s[, i] <- cal.tfr(phi, sigma2, w[i], n_t, P)
  }
  x_coord <- seq(0, 1, length.out = n_t)
  y_coord <- w
  cat('Calculation has been completed! \n')
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", 
                                   "yellow", "#FF7F00", "red", "#7F0000"))
  filled.contour(x_coord, y_coord, s, xlab = 'time', 
                 ylab = 'frequency', color.palette = jet.colors, ...)
  cat('Graph has been drawn!\n')
  return(s)
}

  ## compute coefficient of the forward and backward TVAR model from PARCOR
#ad.est <- function(phi_f, phi_b, P){
#  akm <- rep(list(NA), P)
#  dkm <- rep(list(NA), P)
#  akmm <- rep(list(NA), P)
#  dkmm <- rep(list(NA), P)
#  for(i in 1:P){
#    akmm[[i]] <- rep(list(NA), i)
#    dkmm[[i]] <- rep(list(NA), i)
#  }
#  akm <- phi_f 
#  dkm <- phi_b
#  akmm[[1]][[1]] <- phi_f[[1]]
#  dkmm[[1]][[1]] <- phi_b[[1]]
#  for(j in 2:P){
#    akmm[[j]][[j]] <- phi_f[[j]]
#    dkmm[[j]][[j]] <- phi_b[[j]]
#    for(i in 1:(j-1)){
#     akmm[[j]][[i]] <- akmm[[j-1]][[i]] - phi_f[[j]]*dkmm[[j-1]][[j-i]]
#      dkmm[[j]][[i]] <- dkmm[[j-1]][[i]] - phi_b[[j]]*akmm[[j-1]][[j-i]]
#    }
#  }
#  return(list(akmm = akmm, dkmm = dkmm))
#}

#compute.coef <- function(coef, I, P, n_t){
#  n_tt <- (1+P):(n_t-P)
#  coef_new <- rep(list(NA), I)
  ## i for the number of process, j for the lag
#  for(i in 1:I){
#    tmp <- matrix(nrow = P, ncol = length(n_tt))
#    for(j in 1:P){
#      tmp[j, ] <- coef$akmm[[P]][[j]][i, (1+P):(n_t-P)]
#    }
#    coef_new[[i]] <- tmp
#  }
#  return(coef_new)
#}


