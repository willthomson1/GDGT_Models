### This file contains functions to 
### - fit the GP regression model (fitGPR),
### - compute nearest neighbour distances (Dnear),
### - make predictions from a GPR model (predictGPR),
### - fit the forward model (fitFWD)
### - make predictions from the forward model (predictFWD)
### - extract full posterior predictive distributions (posteriorFWD)


### You will need to install Python and also the GPy library for Python 
### (run '!pip install GPy' in Python's command line)
### both are free to install. This code is based on Python 3.6.


### point to your Python distribution
Sys.setenv(RETICULATE_PYTHON = YOUR_PYTHON_PATH)
###

require(reticulate)

GPy <- import("GPy")
np <- import("numpy")

### import weights and nodes for 500-point Gauss-Hermite quadrature

weightsAndNodes <- read.csv("ghWeightsNodes.csv")[,2:3]

### This is a simple wrapper around the GPy functions for fitting GP regression models
### It takes the modern data and returns an environment containing information about the
### fitted GP regression model
### load_previous allows one to load the already-optimised GP model object (so is faster)

fitGPR <- function(modern.data,modern.temperatures, load_existing = TRUE){
  if(load_existing){
    mb <- np$load("mb.npy")
  } else{
  KK <- GPy$kern$RBF(input_dim = ncol(modern.data),ARD = TRUE)
  mb <- GPy$models$GPRegression(as.matrix(modern.data),as.matrix(modern.temperatures),KK)
  mb$optimize()
  }
  
  return(mb)
}

### This function computes nearest neighbour distances (weighted by the lengthscales of
### the kernel of a fitted GPR object (obtained via fitGPR)

Dnear <- function(newX,model){
  if(ncol(newX) != ncol(model$X)){
    stop("newX has the wrong number of columns")
  }
    K <- model$kern$K(as.matrix(model$X),as.matrix(newX))
    dists <- -log(K / as.numeric(model$kern$variance))
    return(apply(dists,2,min))
}

### This function predicts mean temperatures and standard deviations of predictions
### for the points in newX given the model object (obtained via fitGPR)

predictGPR <- function(newX,model){
  if(ncol(newX) != ncol(model$X)){
    stop("newX has the wrong number of columns")
  }
  
  pred <- mb$predict(newX)
  
  list(means = pred[[1]],sds = sqrt(pred[[2]]))
}

### This function fits the forward model. It takes the modern data and returns a model
### object for use with other functions. Specifying load_existing will load a previously
### trained model object:
### - load_existing = 1 loads an MOGP model based on GDGTs 0-3
### - load_existing = 2 loads an MOGP model based on GDGTs 0-5
### - load_existing = NULL fits the model using GPy (can take some time)

fitFWD <- function(modern.data,modern.temp,load_existing = c(1,2,NULL)){
  
  if(load_existing == 1){
    message('Loading MOGP model object based on GDGTs 0-3')
    mf <- np$load('mf4.npy')
  } else if(load_existing == 2){
    message('Loading MOGP model object based on all 6 GDGTs')
    mf <- np$load('mf6.npy')
  } else{
  
  message("Loading required package robCompositions")
  require(robCompositions)
  message("Imputing zeros")
  modern.data[modern.data == 0] <- NA
  modern.data <- impCoda(modern.data)$xImp
  message("ilr transforming the data")
  modern.data.ilr <- pivotCoord(modern.data)
  
  message("Setting up Multi-Output GP model")
  
  KK <- GPy$kern$Matern32(1)
  icm <- GPy$util$multioutput$ICM(input_dim = 1,
                                  num_outputs = ncol(modern.data.ilr),
                                  kernel = KK)
  temp.list <- lapply(1:ncol(modern.data.ilr),function(j) as.matrix(modern.temp))
  ilr.list <- lapply(1:ncol(modern.data.ilr),function(j) as.matrix(modern.data.ilr[,j]))
  mf <- GPy$models$GPCoregionalizedRegression(temp.list,ilr.list,kernel = icm)
  ### horrible hacky way to fix kernel variance parameter
  mf$constraints$add('fixed',c(0L,1L,2L))
  mf$constraints$remove('fixed',c(1L,2L))
  
  message("Optimising hyperparameters; this might take some time (tens of minutes) depending on your machine")
  mf$optimize()
  }
  
  return(mf)
  
}

### make predictions from the forward model
### inputs: newX :- a matrix or data.frame of GDGT values
###         model :- a model object obtained via fitFWD()
###         prior :- a 2-vector containing the mean and sd of the Gaussian prior on
###                  temperature. Defaults to (15,10)
###         PofXgivenT :- a list containing means, invcovs, dets of p(X|T_j) for each
###                       Gauss-Hermite node T_j
###         returnFullPosterior :- one of:
###                       - FALSE (default): only return means and variances
###                       - A vector of indices for which full posterior should be computed
###                       - TRUE: return full posterior for every new point          

predictFWD <- function(newX,
                       model,
                       prior = c(15,10),
                       PofXgivenT = NULL,
                       returnFullPosterior = FALSE,
                       transformed = F){
  
  dd <- max(model$Y_metadata[[1]]) + 1
  npred <- nrow(newX)
  
  if(ncol(newX) != (dd + 1)){
    stop("newX has the wrong number of columns")
  }
  
  if(returnFullPosterior){
    returnFullPosterior <- 1:npred
  }
  
  whichzerorows <- NULL
  
  if(!transformed){
    if(npred > 2*dd){
      message("Loading required package robCompositions")
      require(robCompositions)
      message("Imputing zeros")
      newX[newX == 0] <- NA
      newX <- impCoda(newX)$xImp
    } else{
      message("Not enough data points to impute zeros; removing rows containing zeros")
      whichzerorows <- which(apply(newX,1,function(x) any(x == 0)))
    }
    message("ilr transforming the data")
    newX <- as.matrix(pivotCoord(newX))
    
  }
  
  
  ## 500 node Gauss-Hermite quadrature (straightforward to use fastGHquad package to
  ## change this if desired)
  
  n_nodes <- 500
  
  xx <- sqrt(2) * prior[2] *weightsAndNodes$x + prior[1]
  
  if(!is.null(returnFullPosterior)){
    priorAtNodes <- dnorm(xx,prior[1],prior[2])
  }
  ww <- weightsAndNodes$w
  
  if(is.null(PofXgivenT)){
    warning("For speed on repeated runs, it is recommended to provide PofXgivenT, 
            which can be obtained via getPofXgivenT()")
  
  inds <- as.integer(0:(dd-1))
  noise_dict <- dict(list(output_index = matrix(inds,dd,1)))
  
  message("Computing p(X|T) at each quadrature node...")
  pb <- txtProgressBar(0,n_nodes)
  
  means <- matrix(NA,n_nodes,dd)
  invcovs <- array(NA,c(n_nodes,dd,dd))
  dets <- rep(NA,n_nodes)
  for(j in 1:n_nodes){
    X <- rep(xx[j],5)
    X <- cbind(X,inds)
    tmpp <- model$predict(X,Y_metadata=noise_dict,full_cov = TRUE)
    means[j,] <- tmpp[[1]]
    cholInvCov <- chol(tmpp[[2]])
    invtmp = chol2inv(cholInvCov)
    invcovs[j,,] = invtmp
    dets[j] = prod(diag(cholInvCov)^2)
    setTxtProgressBar(pb,j)
  }
  
  message("DONE")
  } else{
    means = PofXgivenT$means
    invcovs = PofXgivenT$invcovs
    dets = PofXgivenT$dets
  }
  
  posterior_means <- rep(NA,npred)
  posterior_vars <- rep(NA,npred)
  full_posteriors <- list()
  Zout <- rep(NA,npred)
  
  message("Computing p(T|X) for new data...")
  pb <- txtProgressBar(0,npred)
  
  for(i in which(!((1:npred)%in%whichzerorows))){
    ff <- rep(NA,n_nodes)
    xi <- newX[i,]
    for (j in 1:n_nodes){
      ### evaluate multivariate Gaussian density at i-th  ######################
      ### composition, at j-th temperature node, p(x_i|T_j) ####################
      ##########################################################################
      qf <- t(xi - means[j,])%*%invcovs[j,,]%*%(xi - means[j,])   ##############
      ff[j] <- exp(-0.5 * qf) / sqrt((2 * pi)^dd * dets[j])       ##############
      ##########################################################################
    }
      
    ## compute normalising factor, int p(t) dt, by Gauss-Hermite quadrature
    Z <- t(ww)%*%ff
      
    mu <- t(ww)%*%(ff*xx) / Z ## Gauss-Hermite quadrature again
    posterior_means[i] <- mu
    posterior_vars[i] <- t(ww)%*%(ff * (xx - rep(mu,n_nodes))^2) / Z
    Zout[i] <- Z
    if(i %in% returnFullPosterior){
      full_posteriors[[i]] <- data.frame(xx = xx,posterior = (ff * priorAtNodes) / rep(Z,n_nodes))
    }
    setTxtProgressBar(pb,i)
  }
  message("DONE")
  
  if(!is.null(whichzerorows)){
    message(paste("Predictions not made for points",whichzerorows,
                  "because they contained zero entries"))
  }
  
  return(list(mean = posterior_means,
              variance = posterior_vars,
              full_posteriors = full_posteriors,
              Z = Zout,
              transformedData = newX))
  
}

#### Function to obtain densities p(X|T) at the quadrature nodes 

getPofXgivenT <- function(model){
  dd <- max(model$Y_metadata[[1]]) + 1
  
  ## 500 node Gauss-Hermite quadrature (straightforward to use fastGHquad package to
  ## change this if desired)
  
  n_nodes <- 500
  
  xx <- sqrt(2) * prior[2] *weightsAndNodes$x + prior[1]
  ww <- weightsAndNodes$w
  
  inds <- as.integer(0:(dd-1))
  noise_dict <- dict(list(output_index = matrix(inds,dd,1)))
  
  message("Computing p(X|T) at each quadrature node...")
  pb <- txtProgressBar(0,n_nodes)
  
  means <- matrix(NA,n_nodes,dd)
  invcovs <- array(NA,c(n_nodes,dd,dd))
  dets <- rep(NA,n_nodes)
  for(j in 1:n_nodes){
    X <- rep(xx[j],5)
    X <- cbind(X,inds)
    tmpp <- model$predict(X,Y_metadata=noise_dict,full_cov = TRUE)
    means[j,] <- tmpp[[1]]
    cholInvCov <- chol(tmpp[[2]])
    invtmp = chol2inv(cholInvCov)
    invcovs[j,,] = invtmp
    dets[j] = prod(diag(cholInvCov)^2)
    setTxtProgressBar(pb,j)
  }
  
  return(list(means = means,invcovs = invcovs,dets = dets))

}


#### Compute the (unnormalised) posterior predictive density at the specified
#### points for the data points in newX.
#### Z is a vector of normalising constants (which can be obtained via predictFWD()).
#### transformed is a logical input indicating whether newX contains transformed data.

posteriorFWD <- function(newX,model,points = seq(-10,60,len = 200),
                                     prior = c(15,10),Z = NULL, transformed = FALSE){
  dd <- max(model$Y_metadata[[1]]) + 1
  
  npoints <- length(points)
  npred <- nrow(newX)
  
  priorAtPoints <- dnorm(points,prior[1],prior[2])
  
  inds <- as.integer(0:(dd-1))
  noise_dict <- dict(list(output_index = matrix(inds,dd,1)))
  whichzerorows <- NULL
  
  if(!transformed){
    if(npred > 2*dd){
      message("Loading required package robCompositions")
      require(robCompositions)
      message("Imputing zeros")
      newX[newX == 0] <- NA
      newX <- impCoda(newX)$xImp
      } else{
        whichzerorows <- which(apply(newX,1,function(x) any(x == 0)))
      }
      message("ilr transforming the data")
      newX <- as.matrix(pivotCoord(newX))
  
  }
  
  
  means <- matrix(NA,npoints,dd)
  invcovs <- array(NA,c(npoints,dd,dd))
  dets <- rep(NA,npoints)
  for(j in 1:npoints){
    X <- rep(points[j],5)
    X <- cbind(X,inds)
    tmpp <- model$predict(X,Y_metadata=noise_dict,full_cov = TRUE)
    means[j,] <- tmpp[[1]]
    cholInvCov <- chol(tmpp[[2]])
    invtmp = chol2inv(cholInvCov)
    invcovs[j,,] = invtmp
    dets[j] = prod(diag(cholInvCov)^2)
  }
  
  PPD <- matrix(NA,npoints,npred)
  
  for(i in which(!((1:npred)%in%whichzerorows))){
    ff <- rep(NA,npoints)
    xi <- newX[i,]
    for (j in 1:npoints){
      ### evaluate multivariate Gaussian density at i-th  ######################
      ### composition, at j-th temperature node, p(x_i|T_j) ####################
      ##########################################################################
      qf <- t(xi - means[j,])%*%invcovs[j,,]%*%(xi - means[j,])   ##############
      ff[j] <- exp(-0.5 * qf) / sqrt((2 * pi)^dd * dets[j])       ##############
      ##########################################################################
    }
    PPD[i,] <- ff*priorAtPoints / (ifelse(!is.null(Z),Z[i],1))
  }
  
  if(!is.null(whichzerorows)){
    message(paste("Predictions not made for points",whichzerorows,
                  "because they contained zero entries"))
  }
  return(PPD)
}
