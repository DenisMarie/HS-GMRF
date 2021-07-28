# R source code

# script for simulating data with \Sigma_g = \Sigma_{g,half}, \rho = 0.9, and
# setting \beta_g=(5, 5/sqrt(10), ..., 5/sqrt(10))
library(mvtnorm)
library(coda)
library(invgamma)
library(bayesreg)
library(Matrix)
library(bayesreg)
library(igraph)
library(mvnfast)

source("UtilFunctions.R")
source("MainFunction.R")

n <- 100 # the number of observations
p <- 140 # the number of predictors
rho <- 0.9 # the correlation 
k <- 10 # the size of groups
nb.gp <- p/k # the number of groups
rho.seq <- c(rep(rho, nb.gp/2), rep(0, nb.gp/2)) # vector containing the correlation associated to each group (here the same rho)

# the associated graph
G <- adjacency(k = k, nb.gp = nb.gp, rho.seq = rho.seq)

# the associated contrast matrix
C <- CQ_const(G = G, signX = matrix(1, ncol = p, nrow = p))

# the matrix of predictors: matrix of dimension n x p
X <- matrix(0, nrow = n, ncol = p)
for (g in 0:(nb.gp-1)){
  X[,g*k+1] <- rnorm(n, 0, 1)
  for (i in 2:k)
    X[,g*k+(i)] <- rnorm(n, mean = rho.seq[g+1]*X[,g*k+(1)],
                         sd = sqrt(1-rho.seq[g+1]^2))
}

# the true coefficients
beta_tmp <- rep(c(5, rep(5/sqrt(10), k-1)), times = nb.gp)
beta_true <- rep(0,p)
ind_gp <- c(1,3,5,8,10)
ind_sel <- rep(1:nb.gp, each = k)%in%ind_gp
beta_true[ind_sel] <- beta_tmp[ind_sel]

# the residual variance
se2 <- sum(beta_true[ind_sel]^2)/length(ind_gp)

# the response variable: vector of dimension n
Y <- as.vector(X %*% beta_true + rnorm(n, 0, sqrt(se2)))


# to run HS-GMRF model on a simulated datasete
niter <- 1000
nburn <- 500
res.gmrf <- HS_GMRF(Y = Y, X = X, C = C, a = 1, b = 1,
                    niter = niter)

# to analyze the results
beta_hat <- colMeans(res.gmrf$beta_store[-(1:nburn),])
ci_gmrf <- t(apply(res.gmrf$beta_store[-(1:nburn),], 2, function(c) HPDinterval(as.mcmc(c), probs = 0.95)))
sel_gmrf <- which(apply(sign(ci_gmrf), 1, prod) == 1)
output(beta_true, sel_gmrf)
MSE(beta_true, beta_hat)

plot(beta_true, t = "l")
lines(beta_hat, col =2)
# to compare with HS
res.hs <- bayesreg(Y ~ X,
               data = data.frame(X, Y),
               n.samples =  niter, prior = "hs")
beta_hs <- rowMeans(res.hs$beta[, -(1:nburn)])
ci_hs <- t(apply(res.hs$beta[, -(1:nburn)], 1, function(c) HPDinterval(as.mcmc(c), probs = 0.95)))
sel_hs <- which(apply(sign(ci_hs), 1, prod) == 1)
output(beta_true, sel_hs)
MSE(beta_true, beta_hs)

lines(beta_hs, col =3)

# Fused Lasso  ------------------------------------------------------------

out <- fusedlasso1d(y = Y, X = X)
plot(out, style = "path")
plot(coef(out, lambda=sqrt(n*log(p)))$beta, t = "l")
lines(beta_true, col=2)
