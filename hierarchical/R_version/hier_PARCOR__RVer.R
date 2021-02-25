update.seq.fun <- function(F1, F2, delta1, delta2, yt, n_t, I, m, type){
  mk_0 <- rep(0, I)
  Ck_0 <- diag(I)
  n_0 <- 1
  d_0 <- 1
  nt <- numeric(n_t)
  dt <- numeric(n_t)
  at <- matrix(nrow = I, ncol = n_t)
  Rt <- rep(list(NA), n_t)
  mt <- matrix(nrow = I, ncol = n_t)
  Ct <- rep(list(NA), n_t)
  St <- rep(list(NA), n_t)
  Skt <- rep(list(NA), n_t)
  akt <- matrix(nrow = I, ncol = n_t)
  Rkt <- rep(list(NA), n_t)
  mkt <- matrix(nrow = I, ncol = n_t)
  Ckt <- rep(list(NA), n_t)
  ft <- matrix(nrow = I, ncol = n_t)
  Qt <- rep(list(NA), n_t)
  V2t <- rep(list(NA), n_t)
  if(type == 'forward'){
    lbound <- m + 1
    ubound <- n_t
    sign <- 1
  }else{
    lbound <- 1
    ubound <- n_t - m
    sign <- -1
  }
  #Ckt[[lbound-1]] <- Ck_0
  #mkt[, lbound-1] <- mk_0
  
  for(i in lbound:ubound){
    ## prior distribution update
    F1t <- diag(F1[, i - sign*m])
    if(i == lbound){
      akt[, i] <- mk_0
      #Rkt[[i]] <- Ck_0 + Wt[[i]]
      Rkt[[i]] <- Ck_0/delta1
    }else{
      akt[, i] <- mkt[, i-1, drop = FALSE]
      #Rkt[[i]] <- Ckt[[i-1]] + Wt[[i]]
      Rkt[[i]] <- Ckt[[i-1]]/delta1
    }
    at[, i] <- F2%*%akt[, i, drop = FALSE]
    #Rt[[i]] <- F2%*%Rkt[[i]]%*%t(F2) + V2t[[i]]
    Rt[[i]] <- F2%*%Rkt[[i]]%*%t(F2)/delta2
    V2t[[i]] <- (1-delta2)/delta2 * F2%*%Rkt[[i]]%*%t(F2)
    
    ## predictive distribution update
    ft[, i] <- F1t%*%at[, i, drop = FALSE]
    #Qt[[i]] <- F1t%*%Rt[[i]]%*%t(F1t) + V1t[[i]]
    
    if(i == lbound){
      Qt[[i]] <- F1t%*%Rt[[i]]%*%t(F1t) + d_0/n_0 * diag(I)
    }else{
      Qt[[i]] <- F1t%*%Rt[[i]]%*%t(F1t) + dt[i-1]/nt[i-1]*diag(I)
    }
    
    
    ## posterior distribution 
    # evolution equation
    Skt[[i]] <- Rkt[[i]]%*%t(F1t%*%F2)
    inv_Qt <- chol2inv(chol(Qt[[i]]))
    if(i == lbound){
      nt[i] <- n_0
      dt[i] <- d_0
    }else{
      nt[i] <- nt[i-1] + I
      dt[i] <- dt[i-1] + t(yt[, i, drop = FALSE] - ft[, i, drop = FALSE]) %*%inv_Qt%*%(yt[, i, drop = FALSE] - ft[, i, drop = FALSE])
    }
    
    
    mkt[, i] <- akt[, i, drop = FALSE] + Skt[[i]]%*%inv_Qt%*%(yt[, i, drop = FALSE] - ft[, i, drop = FALSE])
    Ckt[[i]] <- Rkt[[i]] - Skt[[i]]%*%inv_Qt%*%t(Skt[[i]])
    
    # structural equation
    St[[i]] <- Rt[[i]]%*%t(F1t)
    mt[, i] <- at[, i, drop = FALSE] + St[[i]]%*%inv_Qt%*%(yt[, i, drop = FALSE] - ft[, i, drop = FALSE])
    Ct[[i]] <- Rt[[i]] - St[[i]]%*%inv_Qt%*%t(St[[i]])
  }
  return(list(mt = mt, mkt = mkt, Ct = Ct, Ckt = Ckt, 
              Rt = Rt, Rkt = Rkt, at = at, akt = akt,
              Qt = Qt, yt = yt, ft = ft, V1t=dt[ubound]/nt[ubound], 
              V2t = V2t))
}

smooth.fun <- function(mt, mkt, Ct, Ckt, at, akt, I, Rt, Rkt, n_t, m, F1, F2, V1t, V2t, type){
  mnt <- matrix(nrow = I, ncol = n_t)
  Cnt <- rep(list(NA), n_t)
  mnkt <- matrix(nrow = I, ncol = n_t)
  Cnkt <- rep(list(NA), n_t)
  Ant <- rep(list(NA), n_t)
  Ankt <- rep(list(NA), n_t)
  
  if(type == 'forward'){
    lbound <- m + 1
    ubound <- n_t
    sign <- 1
  }else{
    lbound <- 1
    ubound <- n_t - m
    sign <- -1
  }

  mnt[, ubound] <- mt[, ubound]
  Cnt[[ubound]] <- Ct[[ubound]]
  mnkt[, ubound] <- mkt[, ubound]
  Cnkt[[ubound]] <- Ckt[[ubound]]
  
  for(i in (ubound-1):(lbound)){
    F1t <- diag(F1[, i-sign*m])
    #V02t_s <- V1t[[i]] + F1t%*%V2t[[i]]%*%t(F1t)
    V02t_s <- V1t*diag(I) + F1t%*%V2t[[i]]%*%t(F1t)
    inv_V02t_s <- chol2inv(chol(V02t_s))
    Ant[[i]] <- F2%*%Ckt[[i]]%*%t((diag(I) - V2t[[i]]%*%t(F1t)%*%inv_V02t_s%*%F1t)%*%F2)
    Ankt[[i]] <- Ckt[[i]]
    inv_Rtp1 <- chol2inv(chol(Rt[[i+1]]))
    inv_Rktp1 <- chol2inv(chol(Rkt[[i+1]]))
    mnt[, i] <- mt[, i, drop = FALSE] + t(Ant[[i]])%*%inv_Rtp1%*%(mnt[, i+1, drop = FALSE] - at[, i+1, drop = FALSE])
    Cnt[[i]] <- Ct[[i]] - t(Ant[[i]])%*%inv_Rtp1%*%(Rt[[i+1]] - Cnt[[i+1]])%*%inv_Rtp1%*%Ant[[i]]
    mnkt[, i] <- mkt[, i, drop = FALSE] + t(Ankt[[i]])%*%inv_Rktp1%*%(mnkt[, i+1, drop = FALSE] - akt[, i+1, drop = FALSE])
    Cnkt[[i]] <- Ckt[[i]] - t(Ankt[[i]])%*%inv_Rktp1%*%(Rkt[[i+1]] - Cnkt[[i+1]])%*%inv_Rktp1%*%Ankt[[i]]
  }
  return(list(mnt = mnt, mnkt = mnkt, Cnt = Cnt, Cnkt = Cnkt, Ant = Ant, Ankt = Ankt))
}

update.fun <- function(f, par_coef, b, m, type, n_t){
  index <- m
  if(type == 'forward'){
    lbound <- m+1
    ubound <- n_t
    f_new <- f
    for(i in lbound:ubound){
      f_new[, i] <- f[, i] - b[, i-m]*par_coef[, i]
    }
    return(f_new)
  }
  if(type == 'backward'){
    lbound <- 1
    ubound <- n_t - m
    b_new <- b
    for(i in (index+1):(n_t-index)){
      b_new[, i] <- b[, i] - f[, i+m]*par_coef[, i]
    }
    return(b_new)
  }
}


parcor.fun <- function(x.sim, delta1, delta2, P, F2){
  I <- dim(x.sim)[1]
  n_t <- dim(x.sim)[2]
  f_ini <- rep(list(NA), P+1)
  b_ini <- rep(list(NA), P+1)
  f_ini[[1]] <- x.sim
  b_ini[[1]] <- x.sim
  phi_f <- rep(list(NA), P)
  mu_f <- rep(list(NA), P)
  phi_b <- rep(list(NA), P)
  mu_b <- rep(list(NA), P)
  tmp_f <- rep(list(NA), P)
  tmp2_f <- rep(list(NA), P)
  tmp_b <- rep(list(NA), P)
  tmp2_b <- rep(list(NA), P)
  sigma2 <- numeric(P)
  for(j in 1:P){
    cat('\niter:', j)
    tmp_f[[j]] <- update.seq.fun(F1 = b_ini[[j]], F2 = F2, 
                                 delta1 = delta1, 
                                 delta2 = delta2, yt = f_ini[[j]], 
                                 n_t, I,
                                 m=j, type = 'forward')
    sigma2[j] <- tmp_f[[j]]$V1t
    tmp2_f[[j]] <- smooth.fun(tmp_f[[j]]$mt, tmp_f[[j]]$mkt, tmp_f[[j]]$Ct,
                              tmp_f[[j]]$Ckt, tmp_f[[j]]$at, tmp_f[[j]]$akt, 
                              I, tmp_f[[j]]$Rt, tmp_f[[j]]$Rkt, n_t, m = j,
                              b_ini[[j]], F2, V1t = tmp_f[[j]]$V1t,
                              V2t = tmp_f[[j]]$V2t, type = 'forward')
    phi_f[[j]] <- tmp2_f[[j]]$mnt
    mu_f[[j]] <- tmp2_f[[j]]$mnkt
    f_ini[[j+1]] <- update.fun(f = f_ini[[j]], par_coef = phi_f[[j]], 
                               b = b_ini[[j]], m = j, type = 'forward', n_t)
    #cat('\n', f_ini[[j+1]][, 1:10])
    tmp_b[[j]] <- update.seq.fun(F1 = f_ini[[j]], F2 = F2, 
                                 delta1 = delta1, 
                                 delta2 = delta2, yt = b_ini[[j]], 
                                 n_t, I, m=j, type = 'backward')
    tmp2_b[[j]] <- smooth.fun(tmp_b[[j]]$mt, tmp_b[[j]]$mkt, tmp_b[[j]]$Ct,
                              tmp_b[[j]]$Ckt, tmp_b[[j]]$at, tmp_b[[j]]$akt, 
                              I, tmp_b[[j]]$Rt, tmp_b[[j]]$Rkt, n_t, m = j, 
                              f_ini[[j]], F2, V1t = tmp_b[[j]]$V1t,
                              V2t = tmp_b[[j]]$V2t, type = 'backward')
    phi_b[[j]] <- tmp2_b[[j]]$mnt
    mu_b[[j]] <- tmp2_b[[j]]$mnkt
    b_ini[[j+1]] <- update.fun(f = f_ini[[j]], par_coef = phi_b[[j]], 
                               b = b_ini[[j]], m = j, type = 'backward', n_t)
  }
  return(list(tmp_f = tmp_f, tmp2_f = tmp2_f, phi_f = phi_f, 
              f_ini = f_ini, tmp_b = tmp_b, tmp2_b = tmp2_b, 
              phi_b = phi_b, b_ini = b_ini, mu_f = mu_f, mu_b = mu_b,
              sigma2 = sigma2))
}

## compute coefficient of the forward and backward TVAR model from PARCOR
ad.est <- function(phi_f, phi_b, P){
  akmm <- rep(list(NA), P)
  dkmm <- rep(list(NA), P)
  for(i in 1:P){
    akmm[[i]] <- rep(list(NA), i)
    dkmm[[i]] <- rep(list(NA), i)
  }
  akmm[[1]][[1]] <- phi_f[[1]]
  dkmm[[1]][[1]] <- phi_b[[1]]
  for(j in 2:P){
    akmm[[j]][[j]] <- phi_f[[j]]
    dkmm[[j]][[j]] <- phi_b[[j]]
    for(i in 1:(j-1)){
      akmm[[j]][[i]] <- akmm[[j-1]][[i]] - phi_f[[j]]*dkmm[[j-1]][[j-i]]
      dkmm[[j]][[i]] <- dkmm[[j-1]][[i]] - phi_b[[j]]*akmm[[j-1]][[j-i]]
    }
  }
  return(list(akmm = akmm, dkmm = dkmm))
}

compute.coef <- function(coef, I, P, n_t){
  n_tt <- (1+P):(n_t-P)
  coef_new <- rep(list(NA), I)
  ## i for the number of process, j for the lag
  for(i in 1:I){
    tmp <- matrix(nrow = P, ncol = length(n_tt))
    for(j in 1:P){
      tmp[j, ] <- coef$akmm[[P]][[j]][i, (1+P):(n_t-P)]
    }
    coef_new[[i]] <- tmp
  }
  return(coef_new)
}