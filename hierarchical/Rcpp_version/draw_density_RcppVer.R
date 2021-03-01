cal.tfr <- function(phi, sigma2, w, n_t, P){
  ## phi is ith time series coefficients
  ## return the log spectral density
  exp_part <- exp(-2i * pi * w * (1:P))
  #browser()
  return(log(sigma2) - log((Mod(1 - apply(t(phi) * exp_part, 2, sum, na.rm = TRUE)))^2))
}


compute_sd <- function(phi, sigma2, w){
  P <- dim(phi)[3]
  n_t <- dim(phi)[2]
  n_I <- dim(phi)[1]
  s <- array(NA, dim = c(n_I, n_t, length(w)))

  for(i in 1:n_I){
    for(j in 1:length(w)){
      s[i, , j] <- cal.tfr(phi[i, , ], sigma2, w[j], n_t, P)
    }
    #browser()
  }
  cat("Calculation has been completed! \n")
  return(s)
}



draw.density <- function(w, index, P, n_t, s, ...){
  x_coord <- seq(0, 1, length.out = n_t-2*P)
  y_coord <- w
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                                   "yellow", "#FF7F00", "red", "#7F0000"))
  filled.contour(x_coord, y_coord, s[index, (P+1):(n_t-P), ], xlab = 'time',
                 ylab = 'frequency', main = bquote("log sc: "*"f"[.(index)]),
                 color.palette = jet.colors, ...)
  cat('Graph has been drawn!\n')
}
