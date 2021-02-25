#####################################
## Camparison between west and Yang's
#####################################

#####################################
## Simulation data 
#####################################
## x_t = a_t x_t-1 -0.81x_t-2 + epsilon_t
## a_t = 0.8(1-0.5cos(pit/1024))
## epsilon ~ N(0, 1)

## West


dir1 <- '/Users/jay/Documents/Statistics/Papers/Research/Advancement/codes/'
source(paste0(dir1, 'west_TVAR.R'))
#load(paste0(dir1, 'west_sim1_pl.RData'))
set.seed(1)
x1 = 1
x2 = 2
t = 1:1024
x.sim = matrix(nrow = 1, ncol = 1024)
x.sim[1, 1] = x1
x.sim[1, 2] = x2
at1 = 0.8*(1 - 0.5 * cos(pi*t/1024))
for(i in 3:1024){
  x.sim[1, i] = (at1[i]*x.sim[1, i-1] - 0.81*x.sim[1, i-2] + rnorm(1, 0, 1))
}


### true spectral density
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/TRUE_TVAR(2)_spectral_density.png')
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/advancement/revised/images/TRUE_TVAR2_spectral_density.png')
par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
draw.density(rbind(at1, -0.81), rep(1, 1024), main = 'True spectrum', cex.lab = 2, zlim = c(-3, 5), cex.axis = 1.5)
dev.off()

png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/advancement/revised/images/TVAR(2)_ts.png')
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/advancement/revised/images/TVAR(2)_ts.png')
par(cex.lab = 1.5, cex.axis = 1.5)
plot(ts(as.numeric(x.sim)), ylab = 'values')
dev.off()
########################
## initialization 
########################
m_0 <- matrix(0, 2, 1)
C_0 <- 100*diag(2)
k_0 <- 3 
d_0 <- 1
#delta = 1
#beta = 0.935

n_t <- 1024
m <- matrix(0, 2, n_t)
m[, 2] <- m_0
C <- rep(list(0), n_t)
C[[2]] <- C_0
X <- matrix(0, 2, n_t)
X[, 1] <- c(1, 2)
d <- numeric(n_t)
d[2] <- d_0
k <- numeric(n_t)
k[2] <- k_0
## beta from 0.8 to 1
################################
beta1 <- c(0.8, seq(0.9, 1, 0.01))
result1 <- rep(list(0), length(beta1))
ll1 <- numeric(length(beta1))
for(i in 1:length(beta1)){
  result1[[i]] <- update.seq(m, x.sim, d, k, C, P = 2, beta1[i], delta = 1)
  ll1[i] <- result1[[i]]$ll
}
beta1[which(ll1 == max(ll1))]

## beta from 0.97 to 1 by 0.001
## delta from 0.8 to 1
##################################
beta2 <- seq(0.96, 1, 0.001)
delta <- seq(0.8, 1, 0.01)
par_grid <- expand.grid(beta2, delta)
result2 <- rep(list(0), dim(par_grid)[1])
ll2 <- numeric(dim(par_grid)[1])
ptm <- proc.time()
for(i in 1:dim(par_grid)[1]){
  result2[[i]] = update.seq(m, x.sim, d, k, C, P = 2, par_grid[i, 1], delta = par_grid[i, 2])
  ll2[i] = result2[[i]]$ll
  cat('\n:', i)
}
par_grid[which(ll2 == max(ll2)), ]

## Smoothing
result_retro <- retrofilter(result2[[which(ll2 == max(ll2))]], 
                           beta = par_grid[which(ll2 == max(ll2)), 1], 
                           delta = par_grid[which(ll2 == max(ll2)), 2], P = 2, n_t = n_t)
time1 <- proc.time() - ptm

## Plot estimated spectral density
##################################
phi.post = result_retro$m[, 2:1024]
sigma.post = (result_retro$d/(result_retro$k))[2:1024]

## Draw spectral density
##################################
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/estimated_TVAR(2)_by_west_pl.png')
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/advancement/revised/images/estimated_TVAR2_by_west_pl.png')
par(cex.lab = 2, cex.axis = 1.5, cex.main = 2)
draw.density(phi.post, sigma.post, main = 'TVAR(2) model', cex.lab = 1.5, zlim = c(-3, 5))
dev.off()

#save.image(paste0(dir1, 'west_sim1.RData'))
load(paste0(dir1, 'west_sim1.RData'))
## large variance of prior
save.image(paste0(dir1, 'west_sim1_pl.RData'))

load(paste0(dir1, 'west_sim1_pl.RData'))
## Yang's method

source(paste0(dir1,'yang_PARCOR.R'))
n_t <- 1024
P <- 2
mu_0 <- rep(0, P)
c_0 <- rep(100, P)
#c_0 <- rep(0.2^2, P)
v_0 <- rep(3, P)
kappa_0 <- rep(1, P)


gamma <- seq(0.96, 1, by = 0.001)
delta <- seq(0.8, 1, by = 0.01)
par_grid <- expand.grid(gamma, delta)
ptm <- proc.time()

result <- parcor.fun(par_grid, x.sim = x.sim[1, ],
                    mu_0, c_0, v_0, kappa_0, P)
time2 <- proc.time() - ptm

phi_f <- result$phi_f
phi_b <- result$phi_b
sigma_f <- result$sigma_f

akm <- matrix(phi_f[1, ], nrow = 1)
dkm <- matrix(phi_b[1, ], nrow = 1)
fb_coef <- rep(list(0), P)
fb_coef[[1]] <- list(akm = akm, dkm = dkm)
for(i in 2:P){
  fb_coef[[i]] <- ad.est(akmm = fb_coef[[i-1]]$akm, dkmm = fb_coef[[i-1]]$dkm, phi_f[i, ], phi_b[i, ], i)
}

phi.post1 <- fb_coef[[P]]$akm[, (1+P):(n_t-P)]
sigma.post1 <- sigma_f[P, (1+P):(n_t-P)]

## Draw spectral density
##################################
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/estimated_TVAR(2)_by_yang_pl.png')
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/advancement/revised/images/estimated_TVAR2_by_yang_pl.png')
par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
draw.density(phi.post1, sigma.post1, main = 'PARCOR model', zlim = c(-3, 5))
dev.off()


save.image(paste0(dir1, 'yang_sim1_pl.RData'))
#save.image(paste0(dir1, 'yang_sim1.RData'))
load(paste0(dir1, 'yang_sim1_pl.RData'))




## Piecewise Stationary AR Process
## simulation from t = 1 to 1024
## if 1<= t <= 512 x_t = 0.9x_(t-1) + epsilon_t;
## if 513 <= t <= 768, x_t = 1.69x_(t-1) - 0.81 x_(t-2) + epsilon_t;
## if 769 <= t <= 1024, x_t = 1.32x_(t-1) - 0.81 x_(t-2) + epsilon_t; 
##########################################################################
set.seed(1)
x_ini <- 0
n_t <- 1024
x.sim <- length(n_t)
x.sim[1] <- 0.9*x_ini + rnorm(1)
for(i in 2:512){
  x.sim[i] <- 0.9*x.sim[i-1] + rnorm(1)
}
for(i in 513:768){
  x.sim[i] <- 1.69*x.sim[i-1] - 0.81*x.sim[i-2] + rnorm(1)
}
for(i in 769:1024){
  x.sim[i] <- 1.32*x.sim[i-1] - 0.81*x.sim[i-2] + rnorm(1)
}

plot(seq(0, 1, length.out = 1024), x.sim, type = 'l', xlab = 'time', ylab = 'value', main = 'piecewise AR process')


## true piecewise AR process
n_t <- 1024
phi_true_PAR <- matrix(ncol = n_t, nrow = 2)
phi_true_PAR[1, 1:512] <- 0.9
phi_true_PAR[2, 1:512] <- 0
phi_true_PAR[1, 513:768] <- 1.69
phi_true_PAR[2, 513:1024] <- -0.81
phi_true_PAR[1, 769:1024] <- 1.32
sigma_true_PAR = rep(1, n_t)

## draw true piecewise AR process spectral density
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/estimated_piecewise_true.png')
draw.density(phi_true_PAR, sigma_true_PAR, main = NULL, cex.lab = 1.5, zlim = c(-3, 6))
dev.off()


########################
## initialization 
########################
m_0 <- matrix(0, 2, 1)
#C_0 = 0.2^2*diag(2)
C_0 <- 100*diag(2)
k_0 <- 3 
d_0 <- 1
delta <- 1


n_t <- 1024
m <- matrix(0, 2, n_t)
m[, 2] <- m_0
C <- rep(list(0), n_t)
C[[2]] <- C_0
X <- matrix(0, 2, n_t)
X[, 1] <- c(1, 2)
d <- numeric(n_t)
d[2] <- d_0
k <- numeric(n_t)
k[2] <- k_0

## fit 
beta2 <- seq(0.96, 1, 0.001)
delta <- seq(0.8, 1, 0.01)
par_grid <- expand.grid(beta2, delta)
result2 <- rep(list(0), dim(par_grid)[1])
ll2 <- numeric(dim(par_grid)[1])
for(i in 1:dim(par_grid)[1]){
  result2[[i]] <- update.seq(m, x.sim, d, k, C, P = 2, par_grid[i, 1], delta = par_grid[i, 2])
  ll2[i] <- result2[[i]]$ll
  cat('\n:', i)
}
par_grid[which(ll2 == max(ll2)), ]

## Smoothing
result_retro = retrofilter(result2[[which(ll2 == max(ll2))]], 
                           beta = par_grid[which(ll2 == max(ll2)), 1], 
                           delta = par_grid[which(ll2 == max(ll2)), 2], P = 2, n_t = n_t)

## Plot estimated spectral density
##################################
phi.post = result_retro$m[, 2:1024]
sigma.post = (result_retro$d/(result_retro$k))[2:1024]
## Draw spectral density
##################################
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/estimated_piecewise_by_west_pl.png')
draw.density(phi.post, sigma.post, main = NULL, cex.lab = 1.5, zlim = c(-3, 6))
dev.off()


save.image(paste0(dir1, 'west_sim2_pl.RData'))
#save.image(paste0(dir1, 'west_sim2.RData'))

## yang's method
source(paste0(dir1,'yang_PARCOR.R'))

n_t <- 1024
P <- 2
mu_0 <- rep(0, P)
c_0 <- rep(100, P)
#c_0 = rep(0.2^2, P)
v_0 <- rep(3, P)
kappa_0 <- rep(1, P)


gamma <- seq(0.99, 1, by = 0.001)
delta <- seq(0.8, 1, by = 0.1)
par_grid <- expand.grid(gamma, delta)
result <- parcor.fun(par_grid, x.sim = x.sim,
                    mu_0, c_0, v_0, kappa_0, P)

phi_f <- result$phi_f
phi_b <- result$phi_b
sigma_f <- result$sigma_f

akm <- matrix(phi_f[1, ], nrow = 1)
dkm <- matrix(phi_b[1, ], nrow = 1)
fb_coef <- rep(list(0), P)
fb_coef[[1]] <- list(akm = akm, dkm = dkm)
for(i in 2:P){
  fb_coef[[i]] <- ad.est(akmm = fb_coef[[i-1]]$akm, dkmm = fb_coef[[i-1]]$dkm, phi_f[i, ], phi_b[i, ], i)
}

phi.post1 <- fb_coef[[P]]$akm[, (1+P):(n_t-P)]
sigma.post1 <- sigma_f[P, (1+P):(n_t-P)]

## Draw spectral density
##################################
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/estimated_piecewise_by_yang_pl.png')
draw.density(phi.post1, sigma.post1, main = NULL, cex.lab = 1.5, zlim = c(-3, 6))
dev.off()

save.image(paste0(dir1, 'yang_sim2_pl.RData'))
#save.image(paste0(dir1, 'yang_sim2.RData'))


## TVAR(2)
set.seed(1)
x1 = 1
x2 = 2
t = 1:1024
x.sim = matrix(nrow = 2, ncol = 1024)
x.sim[1, 1] = x1
x.sim[1, 2] = x2
x.sim[2, 1] = x1
x.sim[2, 2] = x2
at1 = 0.8*(1 - 0.5 * cos(pi*t/1024))
at2 = 0.9*(1 - 0.1 * cos(pi*t/1024))
for(i in 3:1024){
  x.sim[1, i] = (at1[i]*x.sim[1, i-1] - 0.81*x.sim[1, i-2] + rnorm(1, 0, 1))
  x.sim[2, i] = at2[i]*x.sim[2, i-1] - 0.9*x.sim[2, i-2] + rnorm(1, 0, 1)
}

png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/images/TVAR(2)_ts1.png')
plot(x.sim[1, ], type = 'l', xlab = 'time', ylab = expression(x[t]^1), cex.lab = 1.5, mgp = c(2, 1, 0))
dev.off()
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/images/TVAR(2)_ts2.png')
plot(x.sim[2, ], type = 'l', xlab = 'time', ylab = expression(x[t]^2), cex.lab = 1.5, mgp = c(2, 1, 0))
dev.off()

plot(at1, type = 'l', xlab = 'time', ylab = expression(a[1]), ylim = c(0.2, 1.3))
plot(at2, type = 'l', xlab = 'time', ylab = expression(a[2]), ylim = c(0.2, 1.3))

## hierarchical partrial correlation coefficients
source(paste0(dir1, 'hier_PARCOR.R'))
n_t <- 1024
V1t <- diag(2)
mk_0 <- c(0, 0) 
Ck_0 <- 10*diag(2)
V2t <- 0.0001*diag(2)
Wt <- rep(list(0.00001*diag(2)), n_t)
F2t <- rbind(c(1, 1), c(1, -1))
Gt <- rep(list(diag(2)), n_t)
result_parcor <- parcor.fun(x.sim, Wt, V1t, V2t, mk_0, Ck_0, 2, F2t, Gt)
coef_parcor <- ad.est(result_parcor$phi_f, result_parcor$phi_b, P = 2)
coef <- compute.coef(coef_parcor, I=2, P = 2, n_t)
coef_mu_tmp <- ad.est(result_parcor$mu_f, result_parcor$mu_b, P=2)
coef_mu <- compute.coef(coef_mu_tmp, I = 2, P = 2, n_t)

save.image(paste0(dir1, 'hier_sim1.RData'))
load(paste0(dir1, 'hier_sim1.RData'))
#png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/true_1st.png')
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/advancement/revised/images/true_1st.png')
par(cex.lab = 1.5, cex.axis = 1.5)
draw.density(rbind(at1, -0.81), rep(1, 1024), main = NULL, zlim = c(-3, 6))
dev.off()
#png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/est_1st.png')
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/advancement/revised/images/est_1st.png')
par(cex.lab = 1.5, cex.axis = 1.5)
draw.density(coef[[1]], rep(1, 1024), main = NULL, zlim = c(-3, 6), cex.lab = 1.5)
dev.off()

#png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/true_2st.png')
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/advancement/revised/images/true_2st.png')
par(cex.lab = 1.5, cex.axis = 1.5)
draw.density(rbind(at2, -0.9), rep(1, 1024), main = NULL, zlim = c(-3, 6), cex.lab = 1.5)
dev.off()

#png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/est_2st.png')
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/advancement/revised/images/est_2st.png')
par(cex.lab = 1.5, cex.axis = 1.5)
draw.density(coef[[2]], rep(1, 1024), main = NULL, zlim = c(-3, 6), cex.lab = 1.5)
dev.off()

#png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/mean_simulation.png')
png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/advancement/revised/images/mean_simulation.png')
par(cex.lab = 1.5, cex.axis = 1.5)
draw.density(coef_mu[[1]], rep(1, 1024), main = NULL, zlim = c(-3, 6), cex.lab = 1.5)
dev.off()



### EEG data analysis

mydir <- "/Users/jay/Documents/Statistics/Papers/Research/11:30/data/"
dlmd7_data <- scan(paste0(mydir, 'dlmd7.dat'))
dlmd8_data <- scan(paste0(mydir, 'dlmd8.dat'))
dlmd9_data <- scan(paste0(mydir, 'dlmd9.dat'))
dlmd11_data <- scan(paste0(mydir, 'dlmd11.dat'))
dlmd12_data <- scan(paste0(mydir, 'dlmd12.dat'))
dlmd13_data <- scan(paste0(mydir, 'dlmd13.dat'))
dlmd15_data <- scan(paste0(mydir, 'dlmd15.dat'))
dlmd16_data <- scan(paste0(mydir, 'dlmd16.dat'))
dlmd17_data <- scan(paste0(mydir, 'dlmd17.dat'))

dlmd7_data <- dlmd7_data - mean(dlmd7_data)
dlmd8_data <- dlmd8_data - mean(dlmd8_data)
dlmd9_data <- dlmd9_data - mean(dlmd9_data)
dlmd11_data <- dlmd11_data - mean(dlmd11_data)
dlmd12_data <- dlmd12_data - mean(dlmd12_data)
dlmd13_data <- dlmd13_data - mean(dlmd13_data)
dlmd15_data <- dlmd15_data - mean(dlmd15_data)
dlmd16_data <- dlmd16_data - mean(dlmd16_data)
dlmd17_data <- dlmd17_data - mean(dlmd17_data)


data_all <- rbind(dlmd7_data, dlmd8_data, dlmd9_data,
                  dlmd11_data, dlmd12_data, dlmd13_data,
                  dlmd15_data, dlmd16_data, dlmd17_data)

png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/presentation/advancement/revised/images/ts_csz.png')
par(cex.lab = 2, cex.axis = 2, cex.main = 2)
plot(dlmd12_data, type = 'l', yaxt = 'n', xlab = 'time (secs)', main = 'Channel Cz', ylab = '', xaxt = 'n')
axis(1, c(0, 1200, 2400, 3600), labels = c(0, 27.9, 55.8, 83.72))
dev.off()
n_t <- 3600
I <- 9
var_EEG <- c(var(dlmd7_data), var(dlmd8_data), var(dlmd9_data),
             var(dlmd11_data), var(dlmd12_data), var(dlmd13_data),
             var(dlmd15_data), var(dlmd16_data), var(dlmd17_data))
V1t <- diag(c(var(dlmd7_data), var(dlmd8_data), var(dlmd9_data),
              var(dlmd11_data), var(dlmd12_data), var(dlmd13_data),
              var(dlmd15_data), var(dlmd16_data), var(dlmd17_data)))
mk_0 <- rep(0, I) 
Ck_0 <- 10*diag(I)
V2t <- 0.0001*diag(I)
Wt <- rep(list(0.00001*diag(I)), n_t)
F2t_m <- data.frame(a = gl(I, 1))
F2t <- as.matrix(model.matrix(~a, F2t_m, contrasts = list(a = "contr.sum")))
Gt <- rep(list(diag(I)), n_t)
result_parcor <- parcor.fun(data_all, Wt, V1t, V2t, mk_0, Ck_0, P = 12, F2t, Gt)
coef_parcor <- ad.est(result_parcor$phi_f, result_parcor$phi_b, P = 12)
coef <- compute.coef(coef_parcor, I = 9, P = 12, n_t)
coef_parcor_mean <- ad.est(result_parcor$mu_f, result_parcor$mu_b, P = 12)
coef_mean <- compute.coef(coef_parcor_mean, I = 9, P = 12, n_t)

#save.image(paste0(dir1, 'hier_eeg1.RData'))
load(paste0(dir1, 'hier_eeg1.RData'))
draw.density.eeg <- function(phi, sigma2, interval = 0.001, ...){
  w <- seq(0.0001, 0.4999, by = interval)
  P <- dim(phi)[1]
  n_t <- dim(phi)[2]
  s <- matrix(nrow = n_t, ncol = length(w))
  for(i in 1:length(w)){
    s[, i] <- cal.tfr(phi, sigma2, w[i], n_t, P)
  }
  constant1 <- 83.72/3600
  x_coord <- seq(constant1*P, 83.72 - constant1*P, length.out = n_t)
  constant2 <- 3600/(256/5)
  y_coord <- constant2/(1/w)
  
  cat('Calculation has been completed! \n')
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", 
                                   "yellow", "#FF7F00", "red", "#7F0000"))
  filled.contour(x_coord, y_coord, s, xlab = 'time (s)', 
                 ylab = 'frequency (Hz)', color.palette = jet.colors, ...)
  cat('Graph has been drawn!\n')
}

name2 <- c('F3', 'C3', 'P3', 'Fz', 'Cz', 'Pz', 'F4', 'C4', 'P4')

for(i in 1:9){
  #name <- paste0('Sd', name2[i], 'by_hier.pdf')
  #pdf(paste0('/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/', name), height = 8, width = 13)
  name <- paste0('Sd', name2[i], 'by_hier.png')
  png(filename = paste0('/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/', name))
  par(cex.lab = 1.5, cex.axis = 2)
  draw.density.eeg(phi = coef[[1]], sigma2 = rep(var_EEG[i], 3600), cex.lab = 1.5, cex.main = 2,
               main = paste0('Channel ', name2[i]), zlim = c(0, 16))
  dev.off()
}

png(filename = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/mean_PARCOR_EEG.png')
par(cex.lab = 1.5, cex.axis = 2)
draw.density.eeg(phi = coef_mean[[1]], sigma2 = rep(mean(var_EEG), 3600), cex.lab = 1.5, zlim = c(0, 16))
dev.off()

### multivariate version

dir1 <- '/Users/jay/Documents/Statistics/Papers/Research/Advancement/codes/'
source(paste0(dir1, 'multi_PARCOR.R'))

### example 1

set.seed(1)
x1 <- 1
x2 <- 2
t <- 1:1024
n_t <- 1024
x.sim <- matrix(nrow = 2, ncol = 1024)
x.sim[1, 1] <- x1
x.sim[1, 2] <- x2
x.sim[2, 1] <- x1
x.sim[2, 2] <- x2

rt1 <- seq(0.85, 0.95, length.out = n_t)
rt2 <- seq(0.95, 0.85, length.out = n_t)
lt1 <- seq(5, 20, length.out = n_t)
lt2 <- seq(15, 5, length.out = n_t)

at11 <- rt1*cos(2*pi/lt1)
at22 <- rt2*cos(2*pi/lt2)
at12 <- rt1^2
at21 <- rt2^2

At1 <- rep(list(NA), n_t)
for(i in 3:1024){
  At1[[i]] <- matrix(c(at11[i], 0, 0, at22[i]), nrow = 2, ncol = 2)
  At2 <- matrix(c(-at12[i], 0, 0, -at21[i]), 2, 2)
  x.sim[, i] <- At1[[i]]%*%x.sim[, i - 1, drop = FALSE] +  At2%*%x.sim[, i - 2, drop = FALSE] + matrix(rnorm(2, 0, 1), nrow = 2, ncol = 1)
  #x.sim[, i] <- At1[[i]]%*%x.sim[, i - 1, drop = FALSE] + matrix(rnorm(2, 0, 1), nrow = 2, ncol = 1)
}
#par(mfrow = c(1, 1))
#plot(x.sim[1, ], type = 'l')
#plot(x.sim[2, ], type = 'l')


n_t <- 1024
V1t <- diag(2)
#V1t <- diag(apply(x.sim, 1, var))
mk_0 <- c(0, 0, 0, 0) 
Ck_0 <- 100*diag(4)
Wt <- rep(list(0.00001*diag(4)), n_t)
Gt <- rep(list(diag(4)), n_t)
result_parcor <- parcor.fun(x.sim, Wt, V1t, mk_0, Ck_0, 2, Gt)
coef_parcor <- ad.est(result_parcor$phi_f, result_parcor$phi_b, P = 2, 2, n_t)
result_coef <- compute.coef(coef_parcor$akmm, coef_parcor$dkmm, n_t, 2, 4)


true_phi <- rep(list(NA), 2)
true_phi[[1]] <- rbind(at11, 0, 0, at22)
true_phi[[2]] <- rbind(-at12, rep(0, 1024), rep(0, 1024), -at21)
draw.density(true_phi, diag(2), target_dir = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/',
             type = 'true1_', zlim = c(-3, 6))

draw.density(result_coef, V1t, 
             target_dir = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/', type = 'est1_', 
             zlim = c(-3, 6))
#################


### example 2

set.seed(1)
x1 <- 1
x2 <- 2
t <- 1:1024
n_t <- 1024
x.sim <- matrix(nrow = 2, ncol = 1024)
x.sim[1, 1] <- x1
x.sim[1, 2] <- x2
x.sim[2, 1] <- x1
x.sim[2, 2] <- x2

rt1 <- seq(0.85, 0.95, length.out = n_t)
rt2 <- seq(0.95, 0.85, length.out = n_t)
lt1 <- seq(5, 20, length.out = n_t)
lt2 <- seq(15, 5, length.out = n_t)

at11 <- rt1*cos(2*pi/lt1)
at22 <- rt2*cos(2*pi/lt2)
at12 <- rt1^2
at21 <- rt2^2

At1 <- rep(list(NA), n_t)
for(i in 3:1024){
  At1[[i]] <- matrix(c(at11[i], 0, -.5, at22[i]), nrow = 2, ncol = 2)
  At2 <- matrix(c(-at12[i], 0, 0, -at21[i]), 2, 2)
  x.sim[, i] <- At1[[i]]%*%x.sim[, i - 1, drop = FALSE] +  At2%*%x.sim[, i - 2, drop = FALSE] + matrix(rnorm(2, 0, 1), nrow = 2, ncol = 1)
  #x.sim[, i] <- At1[[i]]%*%x.sim[, i - 1, drop = FALSE] + matrix(rnorm(2, 0, 1), nrow = 2, ncol = 1)
}
#par(mfrow = c(1, 1))
#plot(x.sim[1, ], type = 'l')
#plot(x.sim[2, ], type = 'l')


n_t <- 1024
V1t <- diag(2)
#V1t <- diag(apply(x.sim, 1, var))
mk_0 <- c(0, 0, 0, 0) 
Ck_0 <- 100*diag(4)
Wt <- rep(list(0.00001*diag(4)), n_t)
Gt <- rep(list(diag(4)), n_t)
result_parcor <- parcor.fun(x.sim, Wt, V1t, mk_0, Ck_0, 2, Gt)
coef_parcor <- ad.est(result_parcor$phi_f, result_parcor$phi_b, P = 2, 2, n_t)
result_coef <- compute.coef(coef_parcor$akmm, coef_parcor$dkmm, n_t, 2, 4)


true_phi <- rep(list(NA), 2)
true_phi[[1]] <- rbind(at11, 0, -.5, at22)
true_phi[[2]] <- rbind(-at12, rep(0, 1024), rep(0, 1024), -at21)

draw.density(true_phi, diag(2), target_dir = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/', 
             type = 'true2_', zlim = c(-3, 6))


draw.density(result_coef, V1t, target_dir = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/', 
             type = 'est2_', zlim = c(-2, 6))
############

### example 3

set.seed(1)
x1 <- 1
x2 <- 2
t <- 1:1024
n_t <- 1024
x.sim <- matrix(nrow = 2, ncol = 1024)
x.sim[1, 1] <- x1
x.sim[1, 2] <- x2
x.sim[2, 1] <- x1
x.sim[2, 2] <- x2

rt1 <- seq(0.85, 0.95, length.out = n_t)
rt2 <- seq(0.95, 0.85, length.out = n_t)
lt1 <- seq(5, 20, length.out = n_t)
lt2 <- seq(15, 5, length.out = n_t)

at11 <- rt1*cos(2*pi/lt1)
at22 <- rt2*cos(2*pi/lt2)
at12 <- rt1^2
at21 <- rt2^2

At1 <- rep(list(NA), n_t)
for(i in 3:1024){
  At1[[i]] <- matrix(c(at11[i], 0, -0.8, at22[i]), nrow = 2, ncol = 2)
  At2 <- matrix(c(-at12[i], 0, 0, -at21[i]), 2, 2)
  x.sim[, i] <- At1[[i]]%*%x.sim[, i - 1, drop = FALSE] +  At2%*%x.sim[, i - 2, drop = FALSE] + matrix(rnorm(2, 0, 1), nrow = 2, ncol = 1)
  #x.sim[, i] <- At1[[i]]%*%x.sim[, i - 1, drop = FALSE] + matrix(rnorm(2, 0, 1), nrow = 2, ncol = 1)
}
#par(mfrow = c(1, 1))
#plot(x.sim[1, ], type = 'l')
#plot(x.sim[2, ], type = 'l')


n_t <- 1024
V1t <- diag(2)
#V1t <- diag(apply(x.sim, 1, var))
mk_0 <- c(0, 0, 0, 0) 
Ck_0 <- 100*diag(4)
Wt <- rep(list(0.00001*diag(4)), n_t)
Gt <- rep(list(diag(4)), n_t)
result_parcor <- parcor.fun(x.sim, Wt, V1t, mk_0, Ck_0, 2, Gt)
coef_parcor <- ad.est(result_parcor$phi_f, result_parcor$phi_b, P = 2, 2, n_t)
result_coef <- compute.coef(coef_parcor$akmm, coef_parcor$dkmm, 1024, 2, 4)


true_phi <- rep(list(NA), 2)
true_phi[[1]] <- rbind(at11, 0, -0.8, at22)
true_phi[[2]] <- rbind(-at12, rep(0, 1024), rep(0, 1024), -at21)
draw.density(true_phi, diag(2), target_dir = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/', 
             type = 'true3_', zlim = c(-3, 7))

draw.density(result_coef, V1t, target_dir = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/', 
             type = 'est3_', zlim = c(-3, 7))

#########################################


## EEG channel 11 and 12
####################
mydir <- "/Users/jay/Documents/Statistics/Papers/Research/11:30/data/"
dlmd11_data <- scan(paste0(mydir, 'dlmd11.dat'))
dlmd12_data <- scan(paste0(mydir, 'dlmd12.dat'))
dlmd11_data <- dlmd11_data - mean(dlmd11_data)
dlmd12_data <- dlmd12_data - mean(dlmd12_data)

data_all <- rbind(dlmd11_data, dlmd12_data)

n_t <- 3600
I <- 2
var_EEG <- c(var(dlmd11_data), var(dlmd12_data))
V1t <- diag(c(var(dlmd11_data), var(dlmd12_data)))
#V1t <- diag(2)
mk_0 <- c(0, 0, 0, 0) 
Ck_0 <- 10*diag(4)
Wt <- rep(list(0.0001*diag(4)), n_t)
Gt <- rep(list(diag(4)), n_t)
result_parcor <- parcor.fun(data_all, Wt, V1t, mk_0, Ck_0, P = 5, Gt)
coef_parcor <- ad.est(result_parcor$phi_f, result_parcor$phi_b, P = 5, 2, n_t)
result_coef <- compute.coef(coef_parcor$akmm, coef_parcor$dkmm, n_t, P = 5, 4)


draw.density(result_coef, V1t, target_dir = '/Users/jay/Documents/Statistics/Papers/Research/Advancement/images/', 
             type = 'eeg1_multivariate_')
