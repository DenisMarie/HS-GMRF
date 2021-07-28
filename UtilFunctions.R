# Functions called in the main code

# Function to create the graph/adjacency matrix

# Arguments to input:
# - k: the size of groups
# - nb.gp: the number of groups
# - rho.seq: vector of dimension nb.gp containing the correlation associated to each group

# Output of function:
# - G: matrix of dimension p x p representing the graph

adjacency <- function( k, nb.gp, rho.seq){
  adj <- matrix(0, ncol = k, nrow = k)
  adj[1,] <- 1; adj[,1] <- 1
  diag(adj) <- 1
  G <- kronecker(rho.seq*diag(nb.gp),adj) 
  G[G!=0] <- 1
  diag(G) <- 1
  return(G)
}

# Function to create the contrast matrix

# Arguments to input:
# - G = matix representing the graph
# - signX = matrix containing the sample correlation between covariates

# Output of function:
# - C = contrast matrix of dimension q x p

CQ_const <- function(G = G, signX = FALSE){
  p <- nrow(G)
  q <- sum(G[upper.tri(G)]!=0) # number of connected parameters
  ind <- apply(G, 1, function(x){which(x!=0)})
  C <- NULL
  for(i in 1:p)
  {
    ind_tmp <- ind[[i]]
    if (sum(ind_tmp!=i & ind_tmp>i) != 0){
      ind_tmp <- ind_tmp[ind_tmp!=i & ind_tmp>i]
      tmp <- matrix(0, ncol=p, nrow=length(ind_tmp))
      tmp[, i] <- 1
      for (j in 1:length(ind_tmp))
      {
        tmp[j, ind_tmp[j]] <- -1*signX[i, ind_tmp[j]]
      }
      C <- rbind(C, tmp)
    }else{
      if (length(ind_tmp) == 1){
        tmp <- rep(0,p)
        tmp[i] <- 1
        C <- rbind(C, tmp)
      }
    }
  }
  sub.graph <- clusters(graph.adjacency(G))$membership
  ind.gp <- as.numeric(names(table(sub.graph))[table(sub.graph)>1])
  ind.plus <-  match( ind.gp,sub.graph)
  zero <- matrix(0, ncol = p, nrow = length(ind.plus))
  if (length(ind.plus) == 0) zero <- c(1, rep(p-1)) else{
    for (i in 1:length(ind.plus)) zero[i,ind.plus[i]] <- 1
  }
  C <- rbind(C, zero)
  return(C)
}

# Function to compute Matthews correlation coefficient

# Arguments to input:
# - beta = vector of true coefficients
# - ind_sel = indices of the selected variables

# Output of function:
# - MCC: Matthews correlation coefficient

output <- function(beta, ind_sel){
  aux <- rep(0,length(beta))
  aux[ind_sel] <- 1
  TP <- sum(aux[beta != 0] > 0.5)
  FN <- sum(aux[beta != 0] < 0.5)
  FP <- sum(aux[beta == 0] > 0.5)
  TN <- sum(aux[beta == 0] < 0.5)
  MCC <- (TP*TN -FP*FN) / sqrt((TP+FP)* (TP+FN)*(TN +FP)*(TN +FN))
  return(MCC)
}

# Function to compute the mean squared error of the regression coefficients

# Arguments to input:
# - beta_true = vector of true coefficients
# - beta_hat = vector of estimated coefficients

# Output of function:
# - the mean squared error of the regression coefficients

MSE <- function(beta_true, beta_hat){
   return(mean((beta_true-beta_hat)^2))
}

