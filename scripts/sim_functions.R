####----ABOUT---------------------------------------------------------------####
####   
####   Evaluating alternative methods for modeling trap         
####   efficiencies of outmigrating juvenile salmonids          
####
####   Authors:  Walden et al.                               
####
####   Content: R-Code for functions used in simulation      
####
####   This is the file with all functions needed to run the
####   simulation shown in Section X of the manuscript.
####

####
####----1. Function to generate beta-distributed data-----------------------####
####

## This uses names of variables declared in "scripts/sim_variables.R"
simulateDataset <- function(nobs, type=c("train", "test", "tune"), 
                            efficiency=c("Max efficiency 0.1", "Max efficiency 0.25"), 
                            seeds=NA) {
  
  # Get number of observations depending on dataset type
  thisnobs <- ifelse(type == "test", ntest, nobs)
  
  ## Create vector of parameter values
  Beta <- NA
  Beta0 <- ifelse(efficiency == "Max efficiency 0.1", Beta0_maxpt1, Beta0_maxpt25)
  Beta[1] <- Beta1
  Beta[2] <- Beta2
  Beta[3] <- Beta3
  Beta[4] <- Beta4
  Beta[5] <- Beta5
  Beta[6] <- Beta6
  Beta[7] <- 1          # needed for year random effect
  
  # While-loop checks that generated dataset y values are all less than max 
  # trap efficiency
  thiscondition <- TRUE
  while(thiscondition) {
 
    # If seeds aren't provided, then generate seeds, otherwise retrieve
    # seeds.
    if(any(is.na(seeds))) {
      seed1to6 <- round(runif(6, 1, 999999999),0)
      # While loop checks that there is at least one obs per random effect
      # group
      thiscond2 <- TRUE
      while (thiscond2) {
        seedyr <- round(runif(1, 1, 999999999),0)
        set.seed(seedyr)
        year <- sample.int(nyear, size=thisnobs, replace=TRUE)
        thiscond2 <- length(unique(year)) != nyear
      } # end while loop
    } else {
      seed1to6 <- as.numeric(seeds[1:6])
      seedyr <- as.numeric(seeds[7])
      set.seed(seedyr)
      year <- sample.int(nyear, size=thisnobs, replace=TRUE)
    } # end if-else loop
    
    # Create fixed-effects covariate vectors
    set.seed(seed1to6[2])
    Q <- runif(thisnobs, Qmin, Qmax)
    Q <- scale(Q)
    set.seed(seed1to6[3])
    C <- runif(thisnobs, degCmin, degCmax)
    C <- scale(C)
    set.seed(seed1to6[4])
    effort <- runif(thisnobs, effortmin, effortmax)
    effort <- scale(effort)
    set.seed(seed1to6[5])
    precip <- sample(c(0, 1), size=thisnobs, replace=TRUE)    # yes-no
    
    # Create random effects, one for each group
    set.seed(seed1to6[1])
    yearRE <- rnorm(nyear, 0, raneffSD)

    # Create matrix of covariates
    covariates<-matrix(c(Q, C, effort, precip, 
                         year), nrow=thisnobs, ncol=5)    

    # Create design matrix
    designMat<-cbind(covariates[,1], covariates[,2], covariates[,3], 
                     covariates[,4], 
                     covariates[,3] * covariates[,3],  
                     covariates[,1] * covariates[,2],
                     yearRE[covariates[,5]] + Beta0)
    
    # Calculate mean mu 
    lin.mu <- designMat%*%Beta
    mu <- inverselogit(lin.mu)

    # Calculate beta distribution parameters
    a <- mu / phiTrue
    b <- (1 - mu) / phiTrue
  
    # Realized y
    set.seed(seed1to6[6])
    y <- rbeta(thisnobs, a, b)
    
    # Check realized y to exit while loop
    thiscondition <- ifelse(efficiency == "Max efficiency 0.1", 
                            length(which(y > 0.1)) > 0, 
                            length(which(y > 0.25)) > 0)
  } # Exit while loop
  
  # Collate output
  seedcol <- as.numeric(c(seed1to6, seedyr, rep(NA, thisnobs-7)))
  results <- cbind(y, covariates, a, b, seedcol)
  return(results)
} # End function

####
####----2. Function to create, compile, and run Nimble model----------------####
####

run_MCMC_allcode <- function(model, data, monitors, inits, niter, nburnin) {
  
  myModel <- readBUGSmodel(
    model,
    data=data,
    buildDerivs=FALSE
  )
  
  CmyModel <- compileNimble(myModel)
  myMCMC <- buildMCMC(CmyModel, monitors=monitors)
  CmyMCMC <- compileNimble(myMCMC)
  
  results <- runMCMC(
    CmyMCMC,
    inits=inits,
    niter=niter,
    nburnin=nburnin,
    thin=10
  )
  
  return(results)
}

####
####----3. Function to parallelize Nimble MCMC chains-----------------------####
####

do <- function(cores, pop, fun, ...){
  require(doParallel)
  require(foreach)
  require(nimble)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  results <- foreach(pop, 
                     .packages=c("nimble")) %dopar% fun(...)
  stopCluster(cl)
  return(results)
}

####
####----4. Calculate f-statistic: proportion of mcmc chain with same -------####
####       sign as the mean
####

fstat <- function(x) {
  if(mean(x) < 0) {
    temp <- mean(x < 0, na.rm=TRUE)
    return(temp)
  }
  if(mean(x) > 0) {
    temp <- mean(x > 0, na.rm=TRUE)
    return(temp)
  }
  if(mean(x) == 0) {
    return(0)
  }
}

####
####----5. Calculate log-likelihood for beta-distributed variable-----------####
####       (Weinhold et al. 2020)
####

betaLogLik <- function(y, mu, phi) {
  lgamma(phi) - lgamma(mu*phi) - lgamma((1 - mu)*phi) + (mu*phi - 1)*log(y) + 
    ((1 - mu)*phi - 1)*log(1 - y)
}

####
####----6. Loss function for beta-distributed variable used in MLE----------####
####       (Weinhold et al. 2020)
####

betaLoss <- function(par, y, mu) {
  -sum(lgamma(par) - lgamma(mu*par) - lgamma((1 - mu)*par) + 
         (mu*par - 1)*log(y) + ((1 - mu)*par - 1)*log(1 - y))
}

####
####----7. Get predictions from ranger random forest object-----------------####
####----   when pdp or iml packages used

pred.fun.ranger <- function(object, newdata) { 
  results <- predict(object, newdata)$predictions
  return(results)
}

####
####----8. Get root mean squared error (RMSE) using beta distribution-------####
####       Standardized residuals defined by Ferrari and Cribari-Neto 2004  ####

rmse <- function(y_test, y_predicted, phi) {
  std.residuals <- (y_test - y_predicted) / 
    sqrt( (y_predicted * (1 - y_predicted)) / (1 + phi))
  out <- sqrt(mean(std.residuals^2))
  return(out)
}

####
####----9. Parallel bootstrapping of boosted beta regression----------------####
####                                                                        

run_bootstrap <- function(nboot, mstopmu, mstopphi, maindata, mainfmu, 
                          mainfphi, testdata, i, nyear) {
  inverselogit <- function(alpha) {
    exp(alpha) / ( exp(alpha) + 1 )
  }
  logit <- function(p) {
    log(p/(1-p))
  }
  results <- matrix(data=NA, nrow=300, ncol=1)
  thiscondition <- FALSE
  while(thiscondition == FALSE){
    bootingfunc <- function() {
      tryCatch(
        expr={
              thiscond <- TRUE
              while (thiscond) {
                temp <- sample(1:i, i, replace=TRUE) # get index of obs to keep
                bootdat <- maindata[temp, ] # get dataset
                thiscond <- length(unique(bootdat$year)) != nyear
              }
              results <- FDboostLSS(list(mu=mainfmu, phi=y ~ 1), 
                                      timeformula=NULL, 
                                      families=BetaLSS(), 
                                      data=bootdat, 
                                      control=boost_control(mstop=c(mu=mstopmu, 
                                                                    phi=mstopphi), 
                                                            nu=0.1))
              return(results)
        }, error=function(error_message){
          out <- error_message
          return(out)
        }
      )
    }
    modelBoot <- bootingfunc()
    if(any(class(modelBoot) == "FDboostLSS")) {
        thiscondition <- TRUE
    } else {
      temp <- paste(as.character(modelBoot))
      temp <- str_sub(word(temp,-1),end=-2)
      mstopmu <- as.numeric(temp) - 1
    }
  }

  # Keep mstop values
  results[20,] <- mstop(modelBoot)[[1]]
  results[21,] <- mstop(modelBoot)[[2]]
    
  #  #B1 on link scale (logit)
  results[1,] <- coef(modelBoot, raw=TRUE, which=2, parameter="mu")[[1]]
  #  #B2 on link scale (logit)
  results[2,] <- coef(modelBoot, raw=TRUE, which=3, parameter="mu")[[1]]
  #  #B3 on link scale (logit)
  results[3,] <- coef(modelBoot, raw=TRUE, which=4, parameter="mu")[[1]]
  #  #B4 on link scale (logit)
  results[4,] <- coef(modelBoot, raw=TRUE, which=5, parameter="mu")[[1]]
  #  #B5 on link scale (logit)
  x <- coef(modelBoot)$mu$smterms$`bbs(effort)`$x
  y <- coef(modelBoot)$mu$smterms$`bbs(effort)`$value[,1]
  thisdat <- data.frame("x"=x, "y"=y)
  thismod <- lm(y ~ I(x^2), data=thisdat)
  results[5,] <- as.numeric(thismod$coefficients[2])
  #  #B6 on link scale (logit)
  results[6,] <- coef(modelBoot, raw=TRUE, which=7, parameter="mu")[[1]]
  # #Phi intercept
  results[11,] <- 1 / exp(modelBoot$phi$offset +
                        as.vector(coef(modelBoot)$phi$intercept))
  #B0 on link scale (logit)
  results[10,] <- modelBoot$mu$offset + as.vector(coef(modelBoot)$mu$intercept)
  # yr.sd
  results[7,] <- sd(as.vector(coef(modelBoot, raw=TRUE)$mu$`brandom(year)`))
  
  # prediction on test dataset using bootstrapped model
  testset <- testdata
  oldw <- getOption("warn")
  options(warn=-1)
  temp <- predict(modelBoot, newdata=testset, type="response")
  options(warn=oldw)
  testset$predy <- temp$mu
  testset$predphi <- 1 / temp$phi
  temp2 <- as.vector(coef(modelBoot, raw=TRUE)$mu$`brandom(year)`)
  for(m in 1:nrow(testset)){
    testset$raneff[m] <- temp2[testset$year[m]]
    testset$predy2[m] <- inverselogit(logit(testset$predy[m]) - testset$raneff[m])
    testset$yhat[m] <-rbeta(1, testset$predy2[m] / testset$predphi,
                            (1 - testset$predy2[m]) / testset$predphi)
  }
    
  results[100:199,] <- testset$yhat
  results[200:299,] <- testset$predphi
  return(results)
}


####
####----9.5 Function to parallelize Bootstrapping---------------------------####
####                                                                        

doboot <- function(numcores, pop, fun, ...){
  require(doParallel)
  require(foreach)
  require(gamboostLSS) 
  require(FDboost)
  require(stringr)
  cl <- makeCluster(numcores)
  registerDoParallel(cl)
  results <- foreach(pop, .combine=cbind,
                     .packages=c("gamboostLSS", 
                                   "FDboost", 
                                   "stringr")) %dopar% fun(...)
  stopCluster(cl)
  return(results)
}

####
####----10. Inverse-logit function------------------------------------------####

inverselogit <- function(alpha) {
  exp(alpha) / ( exp(alpha) + 1 )
}

####
####----11. logit function------------------------------------------####

logit <- function(p) {
  log(p/(1-p))
}


####
####----12. Parallelize cvrisk() function-----------------------------------####

myApply <- function(dataset, FDboostLSS, ...){
  myFun <- function(...) {
    library("mboost")
    FDboostLSS(...)
  }
  parLapply(cl=cl, dataset, myFun, ...)
}


####
####----13. Customized predictor function (source: iml package at ----------####
####        https://github.com/giuseppec/iml/blob/main/R/Predictor.R

# This adds "regression" as the task for the iml:Predictor function. Otherwise
# the task is "unknown".

inferTaskFromModel.ranger <- function(model) {
  return("regression")
}