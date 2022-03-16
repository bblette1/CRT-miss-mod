# Test random effects BART code on its own
library(BART)
library(dbarts)
library(MixRF)

# Simulator function
simulator <- function(trial, ICC_out, ICC_mod, num_clusters) {
  
  # Data simulation
  size_clusters <- rpois(num_clusters, 100)
  
  # Create data frame
  df <- data.frame(cluster_ID = rep(1:num_clusters, size_clusters))
  df <- df %>% group_by(cluster_ID) %>%
    mutate(ind_ID = row_number(cluster_ID))
  
  # Treatment arm
  treated <- sample(1:num_clusters, num_clusters/2, replace = FALSE)
  df$A <- 1*(df$cluster_ID %in% treated)
  
  # Covariate
  n <- dim(df)[1]
  df$X <- rnorm(n, 0, 1)
  
  # Modifier
  alpha0_var <- pi^2 * ICC_mod / 3 / (1 - ICC_mod)
  alpha0 <- rnorm(num_clusters, 0, sqrt(alpha0_var))
  Mlogit <- 0.5 + rep(alpha0, size_clusters)
  Mprob <- exp(Mlogit) / (1 + exp(Mlogit))
  df$Mfull <- rbinom(n, 1, Mprob)
  
  # Outcome
  Y_resvar <- 3
  alpha1_var <- Y_resvar * ICC_out / (1 - ICC_out)
  alpha1 <- rnorm(num_clusters, 0, sqrt(alpha1_var))
  #df$Y <- 1 + 1.5*df$A + df$Mfull - 0.75*df$A*df$Mfull +
    #0.7*df$X +
    #0.7*df$X*df$A +
    #0.7*df$X*df$Mfull +
    #0.7*df$X*df$Mfull*df$A +
    #0.8*df$X*df$A - 0.4*df$X*df$Mfull +
    #rep(alpha1, size_clusters) + rnorm(n, 0, sqrt(Y_resvar))
  
  df$Y <- 1 + 1.5*df$A + df$Mfull - 0.75*df$A*df$Mfull +
    #0.7*df$X +
    #0.7*df$X*df$A +
    #0.7*df$X*df$Mfull +
    #0.7*df$X*df$Mfull*df$A +
    0.8*df$X*df$A - 0.4*df$X*df$Mfull + exp(0.1*df$A*df$X*df$Mfull)
    + cos(df$X) + sin(df$X*df$Mfull) +
    #rep(alpha1, size_clusters) +
    rnorm(n, 0, sqrt(Y_resvar))
  
  #testmod <- rbart_vi(Mfull ~ . - cluster_ID, n.trees = 75,
                      #df[, c(1, 3, 4, 5, 6)],
                      #group.by = cluster_ID, keepTrees = TRUE,
                      #n.chains = 2, power = 2)
  trainids <- sample(1:n, 1500)
  train <- df[trainids, ]
  testdat <- df[-trainids, ]
  
  #testmod <- rbart_vi(Y ~ . - cluster_ID, n.trees = 75,
                      #train[, c(1, 3:6)],
                      #group.by = cluster_ID, keepTrees = TRUE,
                      #power = 2, base = 0.95)
  #mean(testmod$tau^2)
  testmod <- bart(Y ~ ., train[, c(1, 3:6)], verbose = FALSE, keeptrees = T)
  
  preds <- predict(testmod, testdat)
  
  testmod2 <- geeglm(Y ~ A*X*Mfull, family = "gaussian", data = train,
                     id = cluster_ID, corstr = "exchangeable")
  
  preds2 <- predict(testmod2, testdat)
  
  #testmod3 <- MixRF(train$Y, cbind(train$A, train$X, train$Mfull),
                    #random = "(1 | cluster_ID)",
                    #train, MaxIterations = 50)
  
  #testmod3 <- gbmm(cbind(train$A, train$X, train$Mfull), train$Y,
                   #u.train = train$cluster_ID)
  
  out1 <- mean(testdat$Y - colMeans(preds))
  out2 <- mean(testdat$Y - preds2)
  
  return(c(out1, out2))
  
}

simres <- array(NA, dim = c(20, 2))
for (i in 1:20) {
  simres[i, ] <- simulator(i, 0.1, 0.1, 100)
}
colMeans(simres)







gbmm=function(
  x.train, y.train,
  x.test=matrix(0,0,0), type='wbart',
  u.train=NULL, B=NULL,
  ntype=as.integer(
    factor(type, levels=c('wbart', 'pbart'))),
  ##factor(type, levels=c('wbart', 'pbart', 'lbart'))),
  sparse=FALSE, theta=0, omega=1,
  a=0.5, b=1, augment=FALSE, rho=NULL,
  xinfo=matrix(0,0,0), usequants=FALSE,
  rm.const=TRUE,
  sigest=NA, sigdf=3, sigquant=0.90,
  k=2, power=2, base=0.95,
  ##sigmaf=NA,
  lambda=NA, tau.num=c(NA, 3, 6)[ntype],
  offset=NULL, ##w=rep(1, length(y.train)),
  ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
  ndpost=1000L, nskip=100L,
  keepevery=c(1L, 10L, 10L)[ntype],
  printevery=100L, transposed=FALSE,
  hostname=FALSE,
  mc.cores = 1L, nice = 19L, seed = 99L
)
{
  if(is.na(ntype))
    stop("type argument must be set to either 'wbart' or 'pbart'")
  ##stop("type argument must be set to either 'wbart', 'pbart' or 'lbart'")
  
  n = length(y.train)
  
  if(length(u.train)==0)
    stop("the random effects indices must be provided")
  
  c.train=integer(n) ## changing from R/arbitrary indexing to C/0
  r.train=unique(u.train)
  ## for(i in 1:n) {
  ##     print(c(i=i, u.train=u.train[i]))
  ##     c.train[i]=which(u.train[i]==r.train)-1
  ## }
  for(i in 1:n) c.train[i]=which(u.train[i]==r.train)-1
  r.train=unique(c.train)
  L=length(r.train)
  n.train=integer(L) ## n_i for each i
  for(i in 1:L) n.train[i]=length(which(r.train[i]==c.train))
  
  if(!transposed) {
    temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                           xinfo=xinfo, rm.const=rm.const)
    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    ## if(length(x.test)>0)
    ##     x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
    if(length(x.test)>0) {
      x.test = bartModelMatrix(x.test)
      x.test = t(x.test[ , temp$rm.const])
    }
    rm.const <- temp$rm.const
    grp <- temp$grp
    rm(temp)
  }
  else {
    rm.const <- NULL
    grp <- NULL
  }
  
  if(n!=ncol(x.train))
    stop('The length of y.train and the number of rows in x.train must be identical')
  
  p = nrow(x.train)
  np = ncol(x.test)
  if(length(rho)==0) rho=p
  if(length(rm.const)==0) rm.const <- 1:p
  if(length(grp)==0) grp <- 1:p
  
  check <- unique(sort(y.train))
  
  if(length(check)==2) {
    if(!all(check==0:1))
      stop('Binary y.train must be coded as 0 and 1')
    if(type=='wbart')
      stop("The outcome is binary so set type to 'pbart'")
    ##stop("The outcome is binary so set type to 'pbart' or 'lbart'")
  }
  
  ## check <- c('wbart', 'pbart', 'lbart')
  
  ## if(!(type %in% check))
  ##     stop("type argument must be set to either 'wbart', 'pbart' or 'lbart'")
  
  if(length(offset)==0) {
    offset=mean(y.train)
    if(type=='pbart') offset=qnorm(offset)
    ##else if(type=='lbart') offset=qlogis(offset)
  }
  
  u=double(L)
  u[1]=NaN
  
  if(type=='wbart') {
    ##y.train = y.train-offset
    
    if(!is.na(sigest) && !is.na(lambda) && lambda==0) {
      ##no op: sigma is fixed and known at given sigest value
    }
    else if(is.na(lambda)) {
      if(is.na(sigest)) {
        if(p < n) {
          temp = lme(y.train~., random=~1|factor(u.train),
                     data.frame(t(x.train),u.train,y.train))
          sigest = summary(temp)$sigma
          u = c(temp$coefficients$random[[1]])
          if(length(B)==0) B=2*sd(u)
        }
        else sigest = sd(y.train)
      }
      qchi = qchisq(1-sigquant, sigdf)
      lambda = (sigest^2)*qchi/sigdf #lambda parameter for sigma prior
    } else {
      sigest=sqrt(lambda)
    }
    
    if(is.na(tau.num)) {
      tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
    } else {
      tau=tau.num/(k*sqrt(ntree))
    }
  } else {
    lambda=1
    sigest=1
    tau=tau.num/(k*sqrt(ntree))
    ## tau=1-tau.interval
    
    ## if(type=='pbart')
    ##     tau=qnorm(1-0.5*tau)/(k*sqrt(ntree))
    ## else if(type=='lbart')
    ##     tau=qlogis(1-0.5*tau)/(k*sqrt(ntree))
  }
  
  ## hot deck missing imputation
  ## must be conducted here since it would
  ## cause trouble with multi-threading on the C++ side
  
  check=(np>0 && np==n)
  
  for(i in 1:n)
    for(j in 1:p) {
      if(check) check=((is.na(x.train[j, i]) && is.na(x.test[j, i])) ||
                         (!is.na(x.train[j, i]) && !is.na(x.test[j, i]) &&
                            x.train[j, i]==x.test[j, i]))
      
      while(is.na(x.train[j, i])) {
        h=sample.int(n, 1)
        x.train[j, i]=x.train[j, h]
      }
    }
  
  if(check) x.test=x.train
  else if(np>0) {
    for(i in 1:np)
      for(j in 1:p)
        while(is.na(x.test[j, i])) {
          h=sample.int(np, 1)
          x.test[j, i]=x.test[j, h]
        }
  }
  
  if(length(B)==0) B <- sigest
  
  if(.Platform$OS.type!='unix') hostname <- FALSE
  else if(hostname)
    hostname <- system('hostname', intern=TRUE)
  
  ptm <- proc.time()
  
  res = .Call("cgbmm",
              ntype, ##as.integer(factor(type, levels=check))-1,
              n,  #number of observations in training data
              p,  #dimension of x
              np, #number of observations in test data
              x.train,   #pxn training data x
              y.train,   #pxn training data x
              x.test,    #p*np test data x
              c.train,
              n.train,
              u,      ## random effects, if estimated
              L,
              B,
              ntree,
              numcut,
              ndpost*keepevery,
              nskip,
              keepevery,
              power,
              base,
              offset,
              tau,
              sigdf,
              lambda,
              sigest,
              ##w,
              sparse,
              theta,
              omega,
              grp,
              a,
              b,
              rho,
              augment,
              printevery,
              xinfo
  )
  
  res$proc.time <- proc.time()-ptm
  res$hostname <- hostname
  
  if(type=='wbart')
    res$yhat.train.mean <- apply(res$yhat.train, 2, mean)
  else {
    if(type=='pbart') res$prob.train = pnorm(res$yhat.train)
    ##else if(type=='lbart') res$prob.train = plogis(res$yhat.train)
    
    res$prob.train.mean <- apply(res$prob.train, 2, mean)
  }
  
  if(np>0) {
    if(type=='wbart')
      res$yhat.test.mean <- apply(res$yhat.test, 2, mean)
    else {
      if(type=='pbart') res$prob.test = pnorm(res$yhat.test)
      ##else if(type=='lbart') res$prob.test = plogis(res$yhat.test)
      
      res$prob.test.mean <- apply(res$prob.test, 2, mean)
    }
  }
  
  res$offset = offset
  names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
  dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
  dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
  res$varcount.mean <- apply(res$varcount, 2, mean)
  res$varprob.mean <- apply(res$varprob, 2, mean)
  res$rm.const <- rm.const
  res$sigest <- sigest
  res$B <- B
  res$u <- u
  attr(res, 'class') <- type
  return(res)
}
