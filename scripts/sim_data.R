####
####----1. ABOUT------------------------------------------------------------####
####
####   Evaluating alternative methods for modeling trap         
####   efficiencies of outmigrating juvenile salmonids          
####
####   Authors:  Walden et al.                               
####
####   Content: R-Code for simulations
####
####   This is the file for Supplement X of the manuscript.
####

####
####----2. SET UP WORKSPACE-------------------------------------------------####
####

##------2.1. Load libraries-------------------------------------------------####
install.packages(c("FDboost", "gamboostLSS", "tuneRanger", "stringr"))

# FDboostLSS()
library(FDboost)
# BetaLSS() [families=...]
library(gamboostLSS)
# tuneRanger()
library(tuneRanger)
# word()
library(stringr)

##------2.2. Set working directory------------------------------------------####

# setwd("C:/Users/mwalden/OneDrive - DOI/CMR_recapture_rate_modeling/analyses")

options(warn=1)       # print all warnings
rm(list=ls())         # clear global environment

##------2.3. Source functions, global variables-----------------------------####

source("scripts/sim_functions.R")
source("scripts/sim_variables.R")


## Check vars
mutuneiter
phituneiter
parallel::detectCores()
kfoldB

####
####----3. TUNING ----------------------------------------------------------####
####

##----3.1. Set up matrix

sim_tunemat <- matrix(data=NA, nrow=length(nset)*length(theseefficiencies),
                      ncol=13)
sim_tunemat <- matrix(data=NA, nrow=length(nset)*length(theseefficiencies),
                      ncol=14)
colnames(sim_tunemat) <- c("samplesize", "efficiency", "tunedat.seed1",
                           "tunedat.seed2",  "tunedat.seed3",  "tunedat.seed4",
                           "tunedat.seed5",  "tunedat.seed6",  "tunedat.seed7",
                           "tuneboost.mustop", "tuneboost.phistop",
                           "tuneforest.mtry", "tuneforest.nodesize", "errorCatch")
sim_tunemat[1:8,1] <- as.vector(rep(nset, each=length(theseefficiencies)))
sim_tunemat[,2] <- rep(theseefficiencies, length(nset))

##----3.2. Set cores

# overwrite global var for speed
# xtracores <- 20

##----3.3. Get dataset seeds

set.seed(seed = NULL)
for(i in 1:nrow(sim_tunemat)) {
  thisss <- as.numeric(sim_tunemat[i,1])
  thiseff <- sim_tunemat[i,2]
  tunedat <- simulateDataset(nobs=thisss, type="tune", efficiency=thiseff)
  st <- thisss-6
  sim_tunemat[i,3:9] <- tunedat[st:thisss,10]
  rm(tunedat)
}

# ----3.4. Tuning

# i <- 1
for(i in 1:nrow(sim_tunemat)) {
  # Start every loop by calling global variables for mstop
  # to ensure reset happens
  mu.ntune <- mutuneiter
  phi.ntune <- phituneiter
  print(paste(mu.ntune, phi.ntune))
  # Get this loop's values
  thisss <- as.numeric(sim_tunemat[i,1])
  thiseff <- sim_tunemat[i,2]
  thisseed <- sim_tunemat[i,3:9]
  print(paste("start_",Sys.time()))
  print(thisseed)
  # Get number of folds so that sample size = number of folds
  # Warning: this will take a very long time to run!!!
  kfoldB <- thisss
  print(kfoldB)
  # Get raw dataset
  tunedat <- simulateDataset2(nobs=thisss, type="tune", efficiency=thiseff,
                              seed=thisseed)
  # Format dataset
  tunedat_lab <- data.frame("discharge"=tunedat[,2],
                            "temperature"=tunedat[,3],
                            "effort"=tunedat[,4],
                            "precip"=as.factor(tunedat[,5]),
                            "year"=as.factor(tunedat[,6]), 
                            "y"=tunedat[,1])
  # Remove raw dataset
  rm(tunedat)
  # Get raw dataset
  print(head(tunedat_lab, 1))
  
  # Tune boosted regression
  thiscondition <- TRUE
  while(thiscondition==TRUE){
    tuningfunc <- function() {
      tryCatch(
        expr={
          modeltune <-  FDboostLSS(list(mu=fmu, phi=fphi), timeformula=NULL, 
                                   families=BetaLSS(), 
                                   data=tunedat_lab, 
                                   control=boost_control(mstop=c(mu=mu.ntune, 
                                                                 phi=phi.ntune),
                                                         nu=nu_hp))
          # print(Sys.time())
          cl <- makeCluster(xtracores)
          myApply <- function(X, FUN, ...){
            myFun <- function(...) {
              library("FDboost")
              tryCatch({
                FUN(...)
              }, error=function(e) {
                out <- paste(e)
                return(out)
              })
            }
            parLapply(cl, X, myFun, ...)
          }      
          set.seed(seed.cvrisk)
          result <- cvrisk(modeltune, papply=myApply, 
                           folds=mboost::cv(model.weights(modeltune), type="kfold", 
                                            B=kfoldB),
                           trace=FALSE)
          stopCluster(cl)
          return(result)
        }, error=function(e) {
          result <- paste(e)
          return(result)
        }
      )
    }
    cvn <- tuningfunc()
    if(any(class(cvn) == "cvriskLSS")) {
      thiscondition <- FALSE
    } else {
      print(cvn)
      if(is.na(sim_tunemat[i,14])) {
        sim_tunemat[i,14] <- cvn
      } else {
        sim_tunemat[i,14] <- paste(sim_tunemat[i,14], cvn, sep="_")
      }
      temp <- paste(as.character(cvn))
      temp <- min(as.numeric(unlist(strsplit(temp, split=" |\n"))), na.rm=TRUE)
      mu.ntune <- temp - 1
      print(mu.ntune)
    }
  }
  # Save hyperparameters - mstop values for mu and phi
  sim_tunemat[i,10] <- mstop(cvn)[[1]]
  sim_tunemat[i,11] <- mstop(cvn)[[2]]
  print(paste(sim_tunemat[i,10], sim_tunemat[i,11]))
  # Remove cvrisk object to make sure next loop starts over correctly
  rm(cvn)
  print(paste("Completed boosted regression tuning_", i, Sys.time()))
  ## Tune random forest
  forest.task <- makeRegrTask(id="fortask", 
                              data=cbind(tunedat_lab[, c("discharge", 
                                                         "temperature", 
                                                         "effort", "precip", 
                                                         "y")]),
                              target="y")
  set.seed(seed.forest)
  res <- tuneRanger(forest.task, num.trees=rftrees,
                    num.threads=xtracores, show.info=FALSE,
                    build.final.model=FALSE)
  # Remove task linked to this loop's dataset
  rm(forest.task)
  # Remove this loop's dataset
  rm(tunedat_lab)
  # Save hyperparameters - mtry and min.node.size
  sim_tunemat[i,12] <- res$recommended.pars$mtry
  sim_tunemat[i,13] <- res$recommended.pars$min.node.size
  print(paste(sim_tunemat[i,12], sim_tunemat[i,13]))
  # Remove tuneRanger object
  rm(res)
  print(paste("Completed random forest tuning_", i, Sys.time()))
}

## Check
sim_tunemat

## Save matrix
saveRDS(sim_tunemat, file="data/sim_tunemat.rds")

####
####----4. TRAINING AND TESTING DATASETS------------------------------------####
####

####
####----4. TRAINING AND TESTING DATASETS------------------------------------####
####

rm(list=ls())         # clear global environment

# Source functions and variables

source("sim_functions.R")
source("sim_variables.R")

##----4.1. Set up matrix
sim_datamat <- matrix(data=NA, nrow=length(nset)*length(theseefficiencies)*nsim,
                      ncol=17)
colnames(sim_datamat) <- c("samplesize", "efficiency", "simulation",
                           "traindat.seed1", "traindat.seed2",  
                           "traindat.seed3", "traindat.seed4",
                           "traindat.seed5", "traindat.seed6",  
                           "traindat.seed7",
                           "testdat.seed1", "testdat.seed2",  
                           "testdat.seed3", "testdat.seed4",
                           "testdat.seed5", "testdat.seed6",  
                           "testdat.seed7")

sim_datamat[,1] <- as.vector(rep(nset, each=length(theseefficiencies)*nsim))
sim_datamat[,2] <- rep(rep(theseefficiencies, each=nsim), length(nset))
sim_datamat[,3] <- rep(seq.int(nsim), length(nset)*length(theseefficiencies))
head(sim_datamat,2)

# i <- 1
for(i in 1:nrow(sim_datamat)) {
  thisss <- as.numeric(sim_datamat[i,1])
  thiseff <- sim_datamat[i,2]
  thissim <- as.numeric(sim_datamat[i,3])
  
  ## Get dataset, save seeds
  traindat <- simulateDataset(nobs=thisss, type="train", efficiency=thiseff)
  sim_datamat[i,4:10] <- traindat[1:7,9]
  rm(traindat)
  
  ## Get dataset, save seeds
  testdat <- simulateDataset(nobs=ntest, type="test", efficiency=thiseff)
  sim_datamat[i,11:17] <- testdat[1:7,9]
  rm(testdat)
}

# head(traindat,8)
# head(testdat,8)

head(sim_datamat,5)
tail(sim_datamat,5)

## Save matrix
saveRDS(sim_datamat, file="data/sim_datamat.rds")