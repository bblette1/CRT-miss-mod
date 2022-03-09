# Check probit mixed model Gibbs sampler
rm(list = ls())
library(MCMCpack)
library(mvtnorm)

N <- 2000
p <- 3
n <- 50
nis <- 40
numiter <- 1000
burn <- numiter / 2
thin <- 1

betares <- array(NA, dim = c(100, p))
bsigres <- rep(NA, 100)

for (g in 1:100) {
  
  Beta <- array(NA, dim = c(nsim - burn, p))
  Sigmab <- rep(NA, nsim - burn)
  
  id <- rep(1:n, nis)
  X1 <- rnorm(N, 0, 1)
  X2 <- rbinom(N, 1, 0.5)
  X <- cbind(1, X1, X2)
  b <- rnorm(n, 0, 0.4)
  yprobit <- 0.5 + 1.2*X1 - 0.6*X2 + rep(b, nis)
  yprob <- pnorm(yprobit)
  y <- rbinom(N, 1, yprob)

  # Prior mean for beta
  beta0 <- c(0.25, 1.5, -0.5)
  # Prior precision for beta
  T0 <- diag(.01, 3)
  # Priors for Random Effects
  nu0 <- 3				      # DF for IWish prior on Sigmab  
  c0 <- 0.01			# Scale matrix for Wishart Prior on Sigmae
  
  # Initial values
  beta <- c(0.25, 1.5, -0.5)
  z <- rep(0, N)
  b <- rep(0, n) # Random effects
  taub <- 2 # Random effect precision
  
  # Posterior Var(beta) outside of loop
  vbeta <- solve(T0 + crossprod(X, X))
  
  # Algorithm
  for (h in 1:numiter) {
    
    # Draw latent variable z
    muz <- X %*% beta + rep(b, nis)
    z[y == 0] <- qnorm(runif(N, 0, pnorm(0, muz)), muz)[y == 0]
    z[y == 1] <- qnorm(runif(N, pnorm(0, muz), 1), muz)[y == 1]
    z[z == -Inf] <- -100
    
    # Update beta
    mbeta <- vbeta %*% (T0 %*% beta0 + crossprod(X, z - rep(b, nis)))
    beta <- c(rmvnorm(1, mbeta, vbeta))
    
    # Update b
    #taub12<-taub1/(1-rhob^2)  		                    # Prior precision of b1|b2
    #mb12<-rhob*sqrt(taub2/taub1)*b2 	                # Prior mean of b1|b2
    #vb1<-1/(nis+taub12)                              # Posterior var of b1|b2,y
    #mb1<-vb1*(taub12*mb12+tapply(z-X%*%beta,id,sum)) # Posterior mean of b1|b2,y
    #b1<-rnorm(n,mb1,sqrt(vb1))
    
    # Update b
    vb <- 1 / (c(taub) + nis)
    #vb <- 1 / nis
    #vb <- 1 / (c(taub) + tapply(z, id, sum))
    #mb <- vb*(c(taub) + tapply((z - X %*% beta), id, sum))
    mb <- vb*(tapply((z - X %*% beta), id, sum))
    b <- rnorm(n, mb, sqrt(vb))
    
    # Update taub
    Sigmab <- riwish(nu0 + n, c0 + crossprod(b))
    taub <- 1 / Sigmab
    
    # Output
    if (h> burn & h%%thin==0) {
      j<-(h-burn)/thin
      Beta[j,]<-beta
      Sigmab[j]<-1/taub
    }
    if (h%%1000==0) print(g)
  
  }
  
  betares[g, ] <- colMeans(Beta)
  bsigres[g] <- mean(Sigmab)
  
}

colMeans(betares) - c(0.5, 1.2, -0.6)
mean(bsigres)
