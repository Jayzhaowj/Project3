dir1 <- '/Users/johnn/Documents/Research/Project3/hierarchical/'
plot_dir <- "/Users/johnn/Documents/Research/Project3/hierarchical/plots/"
source(paste0(dir1, 'R_version/hier_PARCOR_RVer.R'))
source(paste0(dir1, "R_version/draw_density_RVer.R"))
library(PARCOR)


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



n_t <- 1024
V1t <- rep(list(diag(2)), n_t)
#mk_0 <- c(0, 0) 
#Ck_0 <- 10*diag(2)
V2t <- rep(list(0.00001*diag(2)), n_t)
Wt <- rep(list(0.00001*diag(2)), n_t)
F2t <- rbind(c(1, 1), c(1, -1))
#Gt <- rep(list(diag(2)), n_t)
result_parcor <- parcor.fun(x.sim = x.sim, delta1 = 0.99, delta2 = 0.99,
                            P = 2, F2 = F2t)
coef_parcor <- ad.est(result_parcor$phi_f, result_parcor$phi_b, P = 2)
coef <- compute.coef(coef_parcor, I=2, P = 2, n_t)
coef_mu_tmp <- ad.est(result_parcor$mu_f, result_parcor$mu_b, P=2)
coef_mu <- compute.coef(coef_mu_tmp, I = 2, P = 2, n_t)

result_parcor$sigma2[2]

#####################################
## draw spectral density plots
#####################################
#### compute spectral density #####
png(filename = paste0(plot_dir, 'est_1st.png'))
par(cex.lab = 1.5, cex.axis = 1.5)
draw.density(coef[[1]], rep(result_parcor$sigma2[1], 1024), main = NULL, zlim = c(-8, 4), cex.lab = 1.5)
dev.off()



png(filename = paste0(plot_dir, 'est_2st.png'))
par(cex.lab = 1.5, cex.axis = 1.5)
draw.density(coef[[2]], rep(result_parcor$sigma2[2], 1024), main = NULL, zlim = c(-8, 4), cex.lab = 1.5)
dev.off()

png(filename = paste0(plot_dir, 'est_mean.png'))
par(cex.lab = 1.5, cex.axis = 1.5)
draw.density(coef_mu[[1]], rep(1, 1024), main = NULL, zlim = c(-8, 6), cex.lab = 1.5)
dev.off()

png(filename = paste0(plot_dir, 'true_1st.png'))
par(cex.lab = 1.5, cex.axis = 1.5)
draw.density(rbind(at1, -0.81), rep(1, 1024), main = NULL, zlim = c(-3, 8))
dev.off()



png(filename = paste0(plot_dir, 'true_2st.png'))
par(cex.lab = 1.5, cex.axis = 1.5)
draw.density(rbind(at2, -0.9), rep(1, 1024), main = NULL, zlim = c(-3, 8), cex.lab = 1.5)
dev.off()
