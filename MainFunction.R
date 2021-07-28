# R source code for HS-GMRF model

# Arguments to input:
# - Y = vector of observations (n x 1)
# - X = matrix of predictors (n x p)
# - C = contrast matrix of dimension q x p
# - niter = number of iterations for the MCMC algorithm
# - a, b = hyperparameters of the inverse Gamma prior of the residual variance

# Output of function:
# a list containing the niter samples of beta, tau, lambda, se2, mu 


HS_GMRF <- function(Y = Y, X = X, C = C,
                        a = 1, b = 1, niter = 1000){
  p <- ncol(C)
  q <- nrow(C)
  n <- length(Y)
  XprimeX <- crossprod(X)
  # to store samples
  mu_store <- rep(0, niter)
  tau_store <- matrix(0, ncol = q, nrow = niter)
  nu_store <- matrix(0, ncol = q, nrow = niter)
  beta_store <- matrix(0, ncol = p, nrow = niter)
  lambda_store <- rep(0,niter)
  se2_store <- rep(0,niter)
  omega_store<- rep(0,niter)
  # intial values
  lambda_store[1] <- 1
  se2_store[1] <- 1
  beta_store[1,] <- rep(1,p)
  nu_store[1,] <- rep(1,q)
  tau_store[1,] <- rep(1,q)
  omega_store[1] <- 1
  mu_store[1] <- mean(Y) 
  phi <- rep(1,q)
  
  for (i in 2:niter){
    # update of tau and nu
    for (j in 1:q){
      tau_store[i,j] <- max(1E-5,rinvgamma(1,
                                            shape = 1,
                                            rate = phi[j]^2/(2*lambda_store[i-1])+1/nu_store[i-1,j]))
     
      nu_store[i,j] <- rinvgamma(1,
                                    shape = 1,
                                    rate = 1+1/tau_store[i,j])
      
    }
    Q_tmp <- crossprod(C,diag(1/tau_store[i,]))%*%C
    
    # update of beta
    Sigma_tmp <- solve(XprimeX/se2_store[i-1]+
                         Q_tmp/lambda_store[i-1])
    beta_store[i,] <- rmvn(1, 
                            mu = (1/se2_store[i-1])*Sigma_tmp%*%crossprod(X,Y-mu_store[i-1]), 
                            sigma = Sigma_tmp)
    
    phi <- C%*%c(beta_store[i,])
    
    # update of mu
    mu_store[i] <- rnorm(1, (t(rep(1,n))/n)%*%(Y-X%*%beta_store[i,]), sd = sqrt(se2_store[i-1]/n))

    # update of lambda and omega
    lambda_store[i] <- max(1E-5,rinvgamma(1,
                                           shape =(p+1)/2,
                                           rate = 1/omega_store[i-1]+crossprod(beta_store[i,],Q_tmp)%*%beta_store[i,]/2))
    omega_store[i] <- rinvgamma(1,
                                shape =1,
                                rate = 1/se2_store[i-1]+1/lambda_store[i])
    
    #update of se2
    se2_store[i] <- rinvgamma(1,
                              shape = (n+1)/2+a,
                              rate = 0.5*crossprod(Y-X%*%beta_store[i,]-mu_store[i]) +
                                1/omega_store[i]+b )
  }
  return(list(beta_store = beta_store,lambda_store = lambda_store,
              tau_store = tau_store, se2_store =se2_store , mu_store=mu_store))
}
