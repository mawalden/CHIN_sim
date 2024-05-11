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
### sim_template ###
### 
# sessionInfo()-------------------------------------------------------------####
# R version 4.2.3 (2023-03-15 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# attached base packages:
#   [1] parallel  tools     stats4    grid      compiler  stats     graphics  grDevices utils     datasets 
# [11] methods   base     
# 
# other attached packages:
#   [1] MCMCvis_0.16.3      iml_0.11.2          betareg_3.1-4       pdp_0.8.1           doParallel_1.0.17  
# [6] iterators_1.0.14    foreach_1.5.2       nimble_1.1.0        tuneRanger_0.7      lhs_1.1.6          
# [11] mlrMBO_1.1.5.1      smoof_1.6.0.3       checkmate_2.3.1     mlr_2.19.1          ParamHelpers_1.14.1
# [16] ranger_0.16.0       rpart_4.1.23        stringr_1.5.1       gamboostLSS_2.0-7   FDboost_1.1-2      
# [21] mboost_2.9-7        stabs_0.6-4         zoo_1.8-12          XML_3.99-0.16.1     withr_3.0.0        
# [26] vctrs_0.6.2         tibble_3.2.1        survival_3.5-3      sandwich_3.1-0      rlang_1.1.1        
# [31] pracma_2.4.4        pillar_1.9.0        partykit_1.2-19     mvtnorm_1.1-3       numDeriv_2016.8-1.1
# [36] nnet_7.3-18         munsell_0.5.0       modeltools_0.2-23   mgcv_1.8-42         nlme_3.1-162       
# [41] Metrics_0.1.4       Matrix_1.6-5        MASS_7.3-58.2       magrittr_2.0.3      libcoin_1.0-9      
# [46] lattice_0.22-5      jagsUI_1.6.2        inum_1.0-5          igraph_1.4.2        globals_0.16.3     
# [51] generics_0.1.3      future_1.33.1       fastmatch_1.1-4     dplyr_1.0.10        DiceKriging_1.6.0  
# [56] data.table_1.15.2   colorspace_2.0-3    codetools_0.2-19    coda_0.19-4.1       cli_3.4.1          
# [61] BBmisc_1.13         ggplot2_3.5.0      
# 
# loaded via a namespace (and not attached):
#   [1] splines_4.2.3     Formula_1.2-5     backports_1.4.1   glue_1.6.2        quadprog_1.5-8   
# [6] digest_0.6.31     pkgconfig_2.0.3   listenv_0.9.1     scales_1.3.0      parallelMap_1.5.1
# [11] fansi_1.0.3       parallelly_1.37.1 lifecycle_1.0.4   rstudioapi_0.15.0 gtable_0.3.4     
# [16] flexmix_2.3-19    DBI_1.2.2         R6_2.5.1          nnls_1.4          utf8_1.2.2       
# [21] stringi_1.8.3     Rcpp_1.0.9        tidyselect_1.2.1  lmtest_0.9-40   

####
####----2. SET UP WORKSPACE-------------------------------------------------####
####

####----Primary libraries
# install.packages(c("FDboost","gamboostLSS","stringr","rpart","ranger",
#                    "nimble","doParallel","foreach","pdp",
#                    "betareg","iml","MCMCvis"))

##------2.1. Load libraries-------------------------------------------------####
## Load dependencies for iml Interaction() function
# library(BBmisc);      library(cli);       library(coda);    library(codetools); 
# library(colorspace);  library(compiler);  library(data.table); 
# library(DiceKriging); library(dplyr);     library(fastmatch); 
# library(future);      library(generics);  library(ggplot2); library(globals); 
# library(grid);        library(igraph);    library(inum);    library(jagsUI); 
# library(lattice);     library(libcoin);   library(magrittr);library(MASS); 
# library(Matrix);      library(Metrics);   library(mgcv);    library(modeltools); 
# library(munsell);     library(nlme);      library(nnet);    library(numDeriv);
# library(partykit);    library(pillar);    library(pracma);  library(rlang); 
# library(sandwich);    library(survival);  library(tibble);  library(tools); 
# library(vctrs);       library(withr);     library(XML);     library(zoo);

# Primary libraries
library(FDboost)      # FDboostLSS()
library(gamboostLSS)  # BetaLSS() [families=...]
library(stringr)      # word()
library(rpart) 
library(ranger)       # ranger()
library(nimble)       # readBUGSmodel() etc.
library(doParallel)   # foreach/%dopar% etc.
library(foreach)      # foreach()
library(pdp)          # partial()
library(betareg)      # betareg()
library(iml)          # Predictor$new(); Interaction$new()
library(MCMCvis)      # MCMCsummary()

##------2.2. Set working directory------------------------------------------####

# setwd("C:/Users/mwalden/OneDrive - DOI/CMR_recapture_rate_modeling/analyses")

options(warn=1)       # print all warnings
# rm(list=ls())         # clear global environment
parallel::detectCores()

##------2.3. Source functions, models, global variables, and dataset seeds--####

source("scripts/sim_functions.R")
source("scripts/sim_models.R")
source("scripts/sim_variables.R")

# Tuning matrix
sim_tunemat <- read.csv("data/sim_tunemat.csv", header=TRUE)
summary(sim_tunemat)

# Dataset seeds
sim_datamat <- read.csv("data/sim_datamat.csv", header=TRUE)
summary(sim_datamat)

####
####----5. SIMULATION FOR-LOOP----------------------------------------------####
####

## For each N, run model for each simulated dataset and store 
## parameter estimates

xtracores <- 3
nimblecores <- 3
nimblechains <- 9
nset <- 150
nsim <- 2

arrResults <- array(data=NA, dim=c(length(nset)*2, nsim, 7, 300))

## Storage array
## First dimension: rows == "thiscounter" different scenarios 
                            # (sample size X efficiency)
## Second dimension: columns == each simulation "j"
## Third dimension: 
    # [,,1, 1:300]: model 1
    # [,,2, 1:300]: model 2
    # [,,3, 1:300]: model 3
    # [,,4, 1]: sample size
    # [,,5, 1]: efficiency
    # [,j,6, 1]: date/time that simulation was run
    # [,j,6, 2]: total time for that simulation 
    # [,j,7, ]: error messages from error catching

## Fourth dimension [,,1:3, 1:300]: model results
# Intercept: [,,,1:6] 1=mean, 2=LCI, 3=UCI, 4=sd, 5=var.imp(M2), 
#                     6=f.stat(M1)/overlap0(M2)
#                     * all NA for Model 3
#                     * link (logit) scale
# discharge: [,,,10:16] 1=mean, 2=LCI, 3=UCI, 4=sd(M1,M2)/se(M3), 
#                     5=var.imp(M2,M3), 6=f.stat(M1)/overlap0(M2)/p.val(M3)
# degreesC:  [,,,20:26] 1=mean, 2=LCI, 3=UCI, 4=sd(M1,M2)/se(M3), 
#                       5=var.imp(M2,M3), 6=f.stat(M1)/overlap0(M2)/p.val(M3)
# effort:    [,,,30:36] 1=mean, 2=LCI, 3=UCI, 4=sd(M1,M2)/se(M3), 
#                       5=var.imp(M2,M3), 6=f.stat(M1)/overlap0(M2)/p.val(M3)
# precip:    [,,,40:46] 1=mean, 2=LCI, 3=UCI, 4=sd(M1,M2)/se(M3), 
#                       5=var.imp(M2,M3), 6=f.stat(M1)/overlap0(M2)/p.val(M3)
# effort^2:  [,,,50:56] 1=mean, 2=LCI, 3=UCI, 4=sd(M1,M2)/se(M3), 
#                       5=var.imp(M2), 6=f.stat(M1)/overlap0(M2)
# QXdegreesC:[,,,60:66] 1=mean, 2=LCI, 3=UCI, 4=sd(M1,M2)/se(M3), 
#                       5=var.imp(M2), 6=f.stat(M1)/overlap0(M2)/H.val(M3)
# yr.sd:     [,,,70:73] 1=mean, 2=LCI, 3=UCI, 4=sd, 6=overlap0(M2)
# phi:       [,,,80:73] 1=mean, 2=LCI, 3=UCI, 4=sd
# [,,,200]: failure rate
# [,,,210]: rmse
# [,,,220]: ll null 
# [,,3,221]: ll null random
# [,,,225]: ll full or final 
# [,,3,226]: ll final random
# [,,,230]: mean CI width
# [,,,231]: 1=null model failed to converge (can't calculate % deviance explained)
# [,,1,232]: 1=full model failed to converge (delete model results from simulation run)
# [,,2,240:241]: mstop mu and phi from boosted regression model
# [,,2,245:246]: number of times mstop mu and phi from bootstraps were different
#                from the model's mstops
# [,,3,250]: total # iterations
# [,,,299]: fatal error (delete that model's results from that simulation run)

print(Sys.time())
## For-loop to run simulations

thiscounter <- 0

# i <- 200
# i <- 100
# nset <- 200
for (i in nset) {
  
   theseefficiencies <- "Max efficiency 0.1"
  # thisefficiency <- "Max efficiency 0.1"
  for (thisefficiency in theseefficiencies) {
    
    # for indexing output array
    thiscounter <- thiscounter + 1
    # Save sample size in output array
    arrResults[thiscounter, , 4, 1] <- i
    # Save efficiency in output array
    arrResults[thiscounter, , 5, 1] <- thisefficiency
    # Get hyperparameters from sim_tunemat
    tunemustop <- sim_tunemat[which(sim_tunemat[,1]==i & 
        sim_tunemat[,2]==thisefficiency),10]
    tunephistop <- sim_tunemat[which(sim_tunemat[,1]==i & 
        sim_tunemat[,2]==thisefficiency),11]
    forestmtry <- sim_tunemat[which(sim_tunemat[,1]==i & 
        sim_tunemat[,2]==thisefficiency),12]
    forestnodesize <- sim_tunemat[which(sim_tunemat[,1]==i & 
        sim_tunemat[,2]==thisefficiency),13]
    
    ####
    ####----5.2 Loop for each simulation--------------------------------------####
    ####    
    
    print(Sys.time())
    # # flush.console()
    
    # j <- 1
    # nsim <- 1
    for (j in 1:nsim) {
      
      starttime <- Sys.time()
      print(starttime)
      # # flush.console()
      arrResults[thiscounter,j,6, 1] <- starttime
      
      ####
      ##----5.2.1. Get  datasets for this modeling iteration------------------####
      ####
      
      trainseeds <- sim_datamat[which(sim_datamat[,1]==i & 
          sim_datamat[,2]==thisefficiency &
          sim_datamat[,3]==j),4:10]
      traindat <- simulateDataset(nobs=i, type="train", 
        efficiency=thisefficiency, seeds=trainseeds)
      testseeds <- sim_datamat[which(sim_datamat[,1]==i & 
          sim_datamat[,2]==thisefficiency & 
          sim_datamat[,3]==j),11:17]
      testdat <- simulateDataset(nobs=ntest, type="test", 
        efficiency=thisefficiency, seeds=testseeds)
      traindat_lab <- data.frame("discharge"=traindat[,2], 
        "temperature"=traindat[,3], 
        "effort"=traindat[,4], 
        "precip"=as.factor(traindat[,5]), 
        "year"=as.factor(traindat[,6]), 
        "y"=traindat[,1], 
        "int"=1, 
        "a"=traindat[,7], 
        "b"=traindat[,8])
      testdat_lab <- data.frame("discharge"=testdat[,2], 
        "temperature"=testdat[,3], 
        "effort"=testdat[,4], 
        "precip"=as.factor(testdat[,5]), 
        "year"=as.factor(testdat[,6]), 
        "y"=testdat[,1], 
        "int"=1, 
        "a"=testdat[,7], 
        "b"=testdat[,8])
      print(paste("Beginning Model 1 at ", Sys.time()))
      # flush.console()
      
      ####
      ##----5.2.2. Modeling Approach #1---------------------------------------####
      ####  
      
      ##----Run null Nimble model
      
      my5.2.2.0 <- function() {
        
        tryCatch(
          expr={
            
            win.data <- list(
              y=as.vector(traindat[,1]),
              itsyr=as.numeric(as.factor(traindat[,6])),
              nyr=length(unique(traindat[,6])),
              length_y=nrow(traindat))
            params <- c("p0","b.year","phi","yr.sd","ll")
            inits <- list(p0=runif(1, 0.0001, 0.3), 
              phi=runif(1, 0.000001, 0.01), 
              yr.sd=runif(1, 0, 2), 
              b.year=runif(10, -2, 2))
            modGLM <- do(cores=nimblecores, pop=seq(nimblechains), 
              run_MCMC_allcode,
              model="jagsnull",
              data=win.data,
              monitors=params,
              inits=inits,
              niter=50000,
              nburnin=40000)
            if(any(MCMCsummary(modGLM, 
              params=c("b.year","p0", 
                "phi","yr.sd"))$Rhat > 1.1)) {
              
              modGLM <- do(cores=nimblecores, pop=seq(nimblechains), 
                run_MCMC_allcode,
                model="jagsnull",
                data=win.data,
                monitors=params,
                inits=inits,
                niter=100000,
                nburnin=90000)
            }
            return(modGLM)
          }, error=function(error_message){
            
            out <- error_message
            return(out)
          }
        )}
      modGLM5.2.2.0 <- my5.2.2.0()
      if (any(class(modGLM5.2.2.0) == "error")) {
        
        thisstring <- paste0("Sample ", i, ", ", thisefficiency, 
          "Iteration ", j, "Model 1 null error ", 
          modGLM5.2.2.0, "_", Sys.time(), collapse="")
        print(thisstring)
        # flush.console()
        arrResults[thiscounter, j, 7, 1] <- thisstring
        #fatal error flag for null model
        arrResults[thiscounter, j, 1, 298] <- 1
      } else {
        
        arrResults[thiscounter, j, 1, 231] <- 
          ifelse(any(MCMCsummary(modGLM5.2.2.0, 
            params=c("b.year","p0","phi",
              "yr.sd"))$Rhat > 1.1), 1, 0)
        outmcmc <- do.call(rbind, modGLM5.2.2.0)
        arrResults[thiscounter, j, 1, 220] <- mean(outmcmc[,"ll"], 
          na.rm=TRUE) 
        #OK fatal error flag for null model
        arrResults[thiscounter, j, 1, 298] <- 0
        rm(outmcmc)
        # Ends 5.2.2.0. Run null Nimble model
      }
      rm(modGLM5.2.2.0)
      
      ##----Run full Nimble model
      
      my5.2.2.1 <- function() {
        
        tryCatch(
          expr={
            
            win.data <- list(
              y=as.vector(traindat[,1]),
              itsyr=as.numeric(as.factor(traindat[,6])),
              nyr=length(unique(traindat[,6])),
              covs=cbind(traindat[,2:5], traindat[,4]^2, 
                traindat[,2]*traindat[,3]),
              length_y=nrow(traindat))
            paramsFull <- c("Beta","p0","phi","b.year",
              "yr.sd","ll")
            inits <- list(p0=runif(1, 0.001, 0.999),
              phi=runif(1, 0.000001, 0.01),
              Beta=runif(6, -2, 2),
              yr.sd=runif(1, 0, 5),
              b.year=runif(10, -5, 5))
            modGLM <- do(cores=nimblecores, pop=seq(nimblechains), 
              run_MCMC_allcode,
              model="jagsfull",
              data=win.data,
              monitors=paramsFull,
              inits=inits,
              niter=50000,
              nburnin=40000)
            if(any(MCMCsummary(modGLM, params=c("Beta","b.year",
              "p0","phi",
              "yr.sd"))$Rhat > 1.1)) {
              
              modGLM <- do(cores=nimblecores, pop=seq(nimblechains), 
                run_MCMC_allcode,
                model="jagsfull",
                data=win.data,
                monitors=paramsFull,
                inits=inits,
                niter=100000,
                nburnin=90000)
            }
            return(modGLM)
          }, error=function(error_message){
            
            out <- error_message
            return(out)
          }
        )
        # ends function my5.2.2.1.
      }
      modGLM5.2.2.1 <- my5.2.2.1()
      if(any(class(modGLM5.2.2.1) == "error")) {
        
        thisstring <- paste0("Sample ", i, ", ", thisefficiency, 
          "Iteration ", j, "Model 1 full error ", 
          modGLM5.2.2.1, "_", Sys.time(), collapse="")
        print(thisstring)
        # flush.console()
        arrResults[thiscounter, j, 7, 2] <- thisstring
        #fatal error flag for full model
        arrResults[thiscounter, j, 1, 299] <- 1 
      } else {
        
        #OK fatal error flag for full model
        arrResults[thiscounter, j, 1, 299] <- 0 
        arrResults[thiscounter, j, 1, 232] <- 
          ifelse(any(MCMCsummary(modGLM5.2.2.1, 
            params=c("Beta","b.year",
              "p0","phi",
              "yr.sd"))$Rhat > 1.1), 1, 0)
        outmcmc <- do.call(rbind, modGLM5.2.2.1)
        arrResults[thiscounter, j, 1, 225] <- mean(outmcmc[,"ll"], na.rm=TRUE) 
        
        ##----Marginal prediction failure rate on test dataset
        
        testset <- testdat_lab
        tempdat <- matrix(data=NA, nrow=nrow(testset), ncol=nsample)
        for (m in 1:nrow(testset)) {
          
          for (n in 1:nsample) {
            
            tempmu <- inverselogit(logit(outmcmc[n,"p0"]) + 
                outmcmc[n,"Beta[1]"]*testset$discharge[m] + 
                outmcmc[n,"Beta[2]"]*testset$temperature[m] + 
                outmcmc[n,"Beta[3]"]*testset$effort[m] + 
                outmcmc[n,"Beta[4]"]*
                as.numeric(paste(testset$precip[m])) + 
                outmcmc[n,"Beta[5]"]*testset$effort[m]^2 + 
                outmcmc[n,"Beta[6]"]*testset$discharge[m]*
                testset$temperature[m])
            tempphi <- outmcmc[n,"phi"]
            tempdat[m,n] <- rbeta(1, shape1=tempmu / tempphi, 
              shape2=(1 - tempmu) / tempphi)
          }
        }
        # lowerCI "observed"
        testset$obs025 <- qbeta(0.025, testset$a, testset$b, lower.tail=TRUE)
        # upperCI "observed"
        testset$obs0975 <- qbeta(0.975, testset$a, testset$b, lower.tail=TRUE)
        # lower prediction interval y
        testset$pred025 <- apply(tempdat, 1, quantile, probs=0.025)
        # upper prediction interval y
        testset$pred0975 <- apply(tempdat, 1, quantile, probs=0.975)
        testset$predMean <- rowMeans(tempdat)
        testset$dist <- abs(testset$pred0975 - testset$pred025)
        predfail <- subset(testset,(obs025 > pred0975) | (obs0975 < pred025))
        
        ##----5.2.2.4. Store results -----------------------------------------####
        
        arrResults[thiscounter, j, 1, 10:13] <- c(mean(outmcmc[,"Beta[1]"], 
          na.rm=TRUE), 
          quantile(outmcmc[,"Beta[1]"], 
            0.025, na.rm=TRUE),
          quantile(outmcmc[,"Beta[1]"], 
            0.975, na.rm=TRUE),
          sd(outmcmc[,"Beta[1]"], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 1, 16] <- fstat(outmcmc[,"Beta[1]"])
        arrResults[thiscounter, j, 1, 20:23] <- c(mean(outmcmc[,"Beta[2]"],
          na.rm=TRUE),
          quantile(outmcmc[,"Beta[2]"], 
            0.025, na.rm=TRUE),
          quantile(outmcmc[,"Beta[2]"], 
            0.975, na.rm=TRUE),
          sd(outmcmc[,"Beta[2]"], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 1, 26] <- fstat(outmcmc[,"Beta[2]"])
        arrResults[thiscounter, j, 1, 30:33] <- c(mean(outmcmc[,"Beta[3]"], 
          na.rm=TRUE), 
          quantile(outmcmc[,"Beta[3]"], 
            0.025, na.rm=TRUE),
          quantile(outmcmc[,"Beta[3]"], 
            0.975, na.rm=TRUE),
          sd(outmcmc[,"Beta[3]"], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 1, 36] <- fstat(outmcmc[,"Beta[3]"])
        arrResults[thiscounter, j, 1, 40:43] <- c(mean(outmcmc[,"Beta[4]"], 
          na.rm=TRUE),
          quantile(outmcmc[,"Beta[4]"], 
            0.025, na.rm=TRUE),
          quantile(outmcmc[,"Beta[4]"], 
            0.975, na.rm=TRUE),
          sd(outmcmc[,"Beta[4]"], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 1, 46] <- fstat(outmcmc[,"Beta[4]"])
        arrResults[thiscounter, j, 1, 50:53] <- c(mean(outmcmc[,"Beta[5]"], 
          na.rm=TRUE),
          quantile(outmcmc[,"Beta[5]"], 
            0.025, na.rm=TRUE),
          quantile(outmcmc[,"Beta[5]"], 
            0.975, na.rm=TRUE),
          sd(outmcmc[,"Beta[5]"], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 1, 56] <- fstat(outmcmc[,"Beta[5]"])
        arrResults[thiscounter, j, 1, 60:63] <- c(mean(outmcmc[,"Beta[6]"], 
          na.rm=TRUE),
          quantile(outmcmc[,"Beta[6]"], 
            0.025, na.rm=TRUE),
          quantile(outmcmc[,"Beta[6]"], 
            0.975, na.rm=TRUE),
          sd(outmcmc[,"Beta[6]"], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 1, 66] <- fstat(outmcmc[,"Beta[6]"])
        arrResults[thiscounter, j, 1, 70:73] <- c(mean(outmcmc[,"yr.sd"], 
          na.rm=TRUE),
          quantile(outmcmc[,"yr.sd"], 
            0.025, na.rm=TRUE),
          quantile(outmcmc[,"yr.sd"], 
            0.975, na.rm=TRUE),
          sd(outmcmc[,"yr.sd"], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 1, 76] <- fstat(outmcmc[,"yr.sd"])
        arrResults[thiscounter, j, 1, 80:83] <- c(mean(outmcmc[,"phi"], 
          na.rm=TRUE),
          quantile(outmcmc[,"phi"], 
            0.025, na.rm=TRUE),
          quantile(outmcmc[,"phi"], 
            0.975, na.rm=TRUE),
          sd(outmcmc[,"phi"], na.rm=TRUE))
        arrResults[thiscounter, j, 1, 1:4] <- c(mean(logit(outmcmc[,"p0"]), 
          na.rm=TRUE),
          quantile(logit(outmcmc[,"p0"]), 
            0.025, na.rm=TRUE),
          quantile(logit(outmcmc[,"p0"]), 
            0.975, na.rm=TRUE),
          sd(logit(outmcmc[,"p0"]), 
            na.rm=TRUE))
        arrResults[thiscounter, j, 1, 200] <- nrow(predfail)/nrow(testset)
        arrResults[thiscounter, j, 1, 210] <- rmse(testset$y, testset$predMean, 
          phiTrue)    
        arrResults[thiscounter, j, 1, 230] <- mean(testset$dist, na.rm=TRUE)
        rm(outmcmc)
        rm(testset)
      }
      rm(modGLM5.2.2.1)
      # end 5.2.2. Modeling Approach #1.
      print(paste("Completed Model 1 at ", Sys.time(), 
        ". Beginning Model 2."))
      # flush.console()
      
      ####
      ##-----5.2.3. Modeling Approach #2--------------------------------------####
      ####
      
      ##----5.2.3.1. run Full Boosted model-----------------------------------####
      
      # Same formula for baselearners as above in tuning
      # Run primary model for variable selection
      oldw <- getOption("warn")
      options(warn=-1)
      mutrainiter <- tunemustop
      thiscondition <- FALSE
      while(thiscondition==FALSE){
        
        trainfunc <- function(x) {
          
          tryCatch(
            expr={
              
              results <- FDboostLSS(list(mu=fmu, phi=fphi), 
                timeformula=NULL, 
                families=BetaLSS(), 
                data=traindat_lab, 
                control=boost_control(mstop=c(mu=mutrainiter, 
                  phi=tunephistop), 
                  nu=nu_hp))
              return(results)
            }, error=function(error_message){
              
              out <- error_message
              return(out)
            }
          )}
        modelBoost <- trainfunc()
        if(any(class(modelBoost) == "FDboostLSS")) {
          
          thiscondition <- TRUE
        } else {
          
          temp <- paste(as.character(modelBoost))
          temp <- str_sub(word(temp,-1),end=-2)
          mutrainiter <- as.numeric(temp) - 1
        }
      }
      options(warn=oldw)
      arrResults[thiscounter, j, 2, 240:241] <- c(mstop(modelBoost)[[1]],
        mstop(modelBoost)[[2]])
      
      ##------5.2.3.2. Variable importance------------------------------------####
      
      my5.2.3.2 <- function() {
        
        tryCatch(
          expr={
            
            ## mu
            ## 1. which baselearners were selected in boosting steps
            blearner_names <- names(modelBoost$mu$baselearner)
            blearner_selected <- modelBoost$mu$xselect()
            ## 2. compute risks for each step
            # 2.1 initial risk for intercept model
            risk0 <- modelBoost$mu$family@risk(y=modelBoost$mu$response, 
              f=modelBoost$mu$offset)
            # 2.2 risk after each boosting-steps
            riskstep <- modelBoost$mu$risk()
            # 2.3 risk reduction per step
            riskdiff <- c(risk0, riskstep[-length(riskstep)]) - riskstep
            ## 3. compute empirical risk (according to output in cvrisk)
            riskdiff <- riskdiff / length(modelBoost$mu$response)
            ## 4.  explained Risk attributed to baselearners
            explained <- sapply(seq_along(blearner_names), 
              FUN=function(u) {
                sum(riskdiff[which(blearner_selected == u)])
              })
            # names(explained) <- temp
            names(explained) <- blearner_names
            #percent of total
            explainedmu <- explained / sum(explained)
            return(list(explainedmu, risk0))
          }, error=function(error_message){
            
            out <- error_message
            return(out)
          }
          # ends function my5.2.3.2.
        )} 
      explainedmu <- my5.2.3.2()
      if (any(class(explainedmu) == "error")) {
        
        thisstring <- paste0("Sample ", i, ", ", thisefficiency, 
          "Iteration ", j, "Model 2 explained risk error ", 
          explainedmu, "_", Sys.time(), collapse="")
        print(thisstring)
        # flush.console()
        arrResults[thiscounter, j, 7, 3] <- thisstring
        # fatal error flag
        arrResults[thiscounter, j, 2, 299] <- 1 
      } else {
        
        #intercept
        arrResults[thiscounter, j, 2, 5] <- as.numeric(explainedmu[[1]][1])  
        #discharge
        arrResults[thiscounter, j, 2, 15] <- as.numeric(explainedmu[[1]][2])  
        #temperature linear
        arrResults[thiscounter, j, 2, 25] <- as.numeric(explainedmu[[1]][3])  
        #effort
        arrResults[thiscounter, j, 2, 35] <- as.numeric(explainedmu[[1]][4])  
        #precip
        arrResults[thiscounter, j, 2, 45] <- as.numeric(explainedmu[[1]][5])  
        #effort^2
        arrResults[thiscounter, j, 2, 55] <- as.numeric(explainedmu[[1]][6])  
        #interaction
        arrResults[thiscounter, j, 2, 65] <- as.numeric(explainedmu[[1]][7])  
        #year
        arrResults[thiscounter, j, 2, 75] <- as.numeric(explainedmu[[1]][8])  
        #null LL (risk0)
        arrResults[thiscounter, j, 2, 220] <- as.numeric(explainedmu[[2]][1]) 
        # Ends 5.2.3.2. Compute explained risk.
      } 
      rm(explainedmu)
      
      ##----5.2.3.3. Bootstrapping for CI around parameter estimates and------####
      ##             predictions
      
      my5.2.3.3 <- function() {
        
        tryCatch(
          expr={
            
            results <- doboot(numcores=xtracores, p=1:nboot, run_bootstrap,
              nboot=nboot, 
              mstopmu=mstop(modelBoost)[[1]], 
              mstopphi=mstop(modelBoost)[[2]], 
              maindata=traindat_lab, 
              mainfmu=fmu, 
              mainfphi=fphi, 
              testdata=testdat_lab,
              i=i, 
              nyear=nyear)
            return(results)
          }, error=function(error_message){
            
            out <- error_message
            return(out)
          }
          # ends function my5.2.3.3
        )} 
      resultsBoot <- my5.2.3.3()
      if(any(class(resultsBoot) == "error")) {
        
        thisstring <- paste0("Sample ", i, ", ", thisefficiency, 
          "Iteration ", j, "Model 2 bootstrapping error: ", 
          resultsBoot, "_", Sys.time(), collapse="")
        print(thisstring)
        # flush.console()
        arrResults[thiscounter, j, 7, 4] <- thisstring
        #fatal error flag
        arrResults[thiscounter, j, 2, 299] <- 1
      } else {
        
        ##----5.2.3.4. Failure rate & RMSE -----------------------------------####
        testset <- testdat_lab
        oldw <- getOption("warn")
        options(warn=-1)
        pred <- predict(modelBoost, newdata=testset, type="response")
        options(warn=oldw)
        testset$predy <- pred$mu
        testset$predphi <- 1 / pred$phi
        yhat <- NA
        if(is.null(as.vector(coef(modelBoost, raw=TRUE)$mu$`brandom(year)`))) {
          
          yhat <- rbeta(nrow(testset), testset$predy / testset$predphi,
            (1 - testset$predy) / testset$predphi)
        } else {
          
          temp2 <- as.vector(coef(modelBoost, raw=TRUE)$mu$`brandom(year)`)
          for (v in 1:nrow(testset)) {
            
            testset$raneff[v] <- temp2[testset$year[v]]
            testset$predy2[v] <- inverselogit(logit(testset$predy[v]) - 
                testset$raneff[v])
            yhat[v] <- rbeta(1, testset$predy2[v] / testset$predphi[v],
              (1 - testset$predy2[v]) / testset$predphi[v])
          }
        }
        
        ## failure rate
        
        testset$obs0.025 <- qbeta(0.025, testset$a,
          testset$b,
          lower.tail=TRUE)
        testset$obs0.975 <- qbeta(0.975, testset$a,
          testset$b,
          lower.tail=TRUE)
        yobs <- resultsBoot[100:199,]
        testset$yhat0.025 <- apply(yobs, 1, quantile, probs=0.025, 
          lower.tail=TRUE, na.rm=TRUE)
        testset$yhat0.975 <- apply(yobs, 1, quantile, probs=0.975, 
          lower.tail=TRUE, na.rm=TRUE)
        predfail <- subset(testset,
          (yhat0.025 > obs0.975) | (yhat0.975 < obs0.025))
        
        ##----5.2.3.5. Store results ---------------------------------------####
        
        arrResults[thiscounter, j, 2, 10:13] <- c(mean(resultsBoot[1,], 
          na.rm=TRUE),
          quantile(resultsBoot[1,], 
            0.025, na.rm=TRUE),
          quantile(resultsBoot[1,], 
            0.975, na.rm=TRUE),
          sd(resultsBoot[1,], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 2, 20:23] <- c(mean(resultsBoot[2,], 
          na.rm=TRUE),
          quantile(resultsBoot[2,], 
            0.025, na.rm=TRUE),
          quantile(resultsBoot[2,], 
            0.975, na.rm=TRUE),
          sd(resultsBoot[2,], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 2, 30:33] <- c(mean(resultsBoot[3,], 
          na.rm=TRUE),
          quantile(resultsBoot[3,], 
            0.025, 
            na.rm=TRUE),
          quantile(resultsBoot[3,], 
            0.975, 
            na.rm=TRUE),
          sd(resultsBoot[3,], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 2, 40:43] <- c(mean(resultsBoot[4,], 
          na.rm=TRUE),
          quantile(resultsBoot[4,], 
            0.025, 
            na.rm=TRUE),
          quantile(resultsBoot[4,], 
            0.975, 
            na.rm=TRUE),
          sd(resultsBoot[4,], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 2, 50:53] <- c(mean(resultsBoot[5,], 
          na.rm=TRUE),
          quantile(resultsBoot[5,], 
            0.025, 
            na.rm=TRUE),
          quantile(resultsBoot[5,], 
            0.975, 
            na.rm=TRUE),
          sd(resultsBoot[5,], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 2, 60:63] <- c(mean(resultsBoot[6,], 
          na.rm=TRUE),
          quantile(resultsBoot[6,], 
            0.025, 
            na.rm=TRUE),
          quantile(resultsBoot[6,], 
            0.975, 
            na.rm=TRUE),
          sd(resultsBoot[6,], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 2, 70:73] <- c(mean(resultsBoot[7,], 
          na.rm=TRUE),
          quantile(resultsBoot[7,], 
            0.025, 
            na.rm=TRUE),
          quantile(resultsBoot[7,], 
            0.975, 
            na.rm=TRUE),
          sd(resultsBoot[7,], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 2, 80:83] <- c(mean(resultsBoot[11,]),
          quantile(resultsBoot[11,], 
            0.025, na.rm=TRUE),
          quantile(resultsBoot[11,], 
            0.975, na.rm=TRUE),
          sd(resultsBoot[11,], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 2, 1:4] <- c(mean(resultsBoot[10,], 
          na.rm=TRUE),
          quantile(resultsBoot[10,], 
            0.025, 
            na.rm=TRUE),
          quantile(resultsBoot[10,], 
            0.975, 
            na.rm=TRUE),
          sd(resultsBoot[10,], 
            na.rm=TRUE))
        arrResults[thiscounter, j, 2, 245] <-  
          length(which(resultsBoot[20,]!=mstop(modelBoost)[[1]]))
        arrResults[thiscounter, j, 2, 246] <-  
          length(which(resultsBoot[20,]!=mstop(modelBoost)[[2]]))
        
        arrResults[thiscounter, j, 2, 210] <- rmse(testset$y, yhat, 
          phiTrue)    
        
        arrResults[thiscounter, j, 2, 200] <- nrow(predfail)/nrow(testset)
        
        arrResults[thiscounter, j, 2, 225] <- tail(modelBoost$mu$risk(),1)
        
        arrResults[thiscounter, j, 2, 230] <- mean(abs(testset$yhat0.975 - 
            testset$yhat0.025), 
          na.rm=TRUE)
        rm(testset)
        # end else, end my5.2.3.3
      }
      rm(resultsBoot)
      rm(modelBoost)
      print(paste("Completed Model 2 at ", Sys.time(), 
        ". Beginning Model 3."))
      # flush.console()
      
      ####
      ##----5.2.4. Modeling Approach #3---------------------------------------####
      ####  ~ 5 min (will vary depending on # iterations)
      
      ##----5.2.4.1. While loop-----------------------------------------------####
      
      data <- traindat_lab
      convfail <- NA
      fatalflag <- NA
      my5.2.4.1 <- function() {
        
        tryCatch(
          expr={
            
            MaxIterations <- 20
            #parse formula 
            Target <- data[,6]
            initialRandomEffects=rep(0,length(data[,6])) 
            #set up variables for loop 
            ContinueCondition<-TRUE 
            iterations<-0 
            #initial values 
            AdjustedTarget<-as.numeric(Target)-initialRandomEffects
            oldlik<- Inf
            while(ContinueCondition==TRUE){ 
              
              # Current values of variables 
              iterations <- iterations+1 
              # For first iteration, get null model (intercept-only) 
              if(iterations == 1) {
                
                forest <- ranger(AdjustedTarget ~ int,
                  data=data, splitrule="beta", 
                  num.threads=xtracores, 
                  num.trees=rftrees, replace=FALSE, mtry=1,
                  importance="permutation", keep.inbag=TRUE,
                  min.node.size=forestnodesize, node.stats=TRUE,
                  write.forest=TRUE, seed=seed.forest)
              } else {
                
                forest <- ranger(AdjustedTarget ~ discharge + temperature + 
                    effort + precip, 
                  data=data, splitrule="beta", 
                  num.threads=xtracores, 
                  num.trees=rftrees, replace=FALSE, mtry=forestmtry, 
                  importance="permutation", keep.inbag=TRUE, 
                  min.node.size=forestnodesize, node.stats=TRUE, 
                  quantreg=TRUE, scale.permutation.importance=TRUE, 
                  write.forest=TRUE, seed=seed.forest)
              }
              forestprob <- predict(forest, data=data, 
                num.threads=xtracores)$predictions
              forestprob2 <- as.vector(logit(forestprob))
              if(iterations == 1) {
                
                par_phi <- optim(par=1,fn=betaLoss, mu=forest$predictions,
                  y=data$y, method="Brent", 
                  lower=0, upper=100000)$par
                nullLLfixed <- sum(betaLogLik(y=data$y, mu=forest$predictions,
                  phi=par_phi))
              }
              # Estimate New Random Effects and Errors
              win.data <- list(
                y=data$y,fixedeff=forestprob2,
                itsyr=as.numeric(as.factor(data$year)),
                nyr=length(unique(data$year)),
                length_y=nrow(traindat))
              paramsRF <- c("ll","b.year","yr.sd","phi")
              inits <- list(yr.sd=runif(1, 0, 5),
                b.year=runif(10, -5, 5),
                phi=runif(1, 0.00001, 1))
              modRF <- do(cores=nimblecores, pop=seq(nimblechains), 
                run_MCMC_allcode,
                model="jagsRF",
                data=win.data,
                monitors=paramsRF,
                inits=inits,
                niter=50000,
                nburnin=40000)
              if(any(MCMCsummary(modRF, 
                params=c("b.year",
                  "phi","yr.sd"))$Rhat > 1.1)) {
                
                modRF <- do(cores=nimblecores, pop=seq(nimblechains), 
                  run_MCMC_allcode,
                  model="jagsRF",
                  data=win.data,
                  monitors=paramsRF,
                  inits=inits,
                  niter=100000,
                  nburnin=90000)
              }
              outmcmc <- do.call(rbind, modRF)
              newlik <- mean(outmcmc[,"ll"], na.rm=TRUE)
              if(any(MCMCsummary(modRF, 
                params=c("b.year",
                  "phi","yr.sd"))$Rhat > 1.1)) {
                
                convfail <- 1
                ContinueCondition <- FALSE
              } else if(iterations == 1) {
                
                nullLLrandom <- newlik
                ContinueCondition <- TRUE
              } else if(iterations == 2) {
                
                ErrorTolerance <- sd(outmcmc[,"ll"], na.rm=TRUE)
                ContinueCondition <- TRUE
              } else {
                
                ContinueCondition <- (abs(newlik-oldlik) > ErrorTolerance & 
                    iterations <  MaxIterations)
              }
              oldlik <- newlik
              # Extract random effects to make the new adjusted target 
              newbyear <- colMeans(outmcmc[,grep("b.year", colnames(outmcmc))])
              # n <- 1
              for (n in 1:nrow(data)) {
                
                data$newy[n] <- inverselogit(logit(data$y[n]) - 
                    newbyear[data$year[n]])
              }
              AdjustedTarget <- data$newy
            }
            return(list(forest, iterations, ErrorTolerance, nullLLrandom, newlik,
              convfail, outmcmc, nullLLfixed))
          }, error=function(error_message){
            
            out <- error_message
            return(out)
          }
          # ends function my5.2.4.1.
        )}
      resultsRF <- my5.2.4.1()
      if(any(class(resultsRF) == "error")) {
        
        thisstring <- paste0("Sample ", i, ", ", thisefficiency, 
          "Iteration ", j, "Model 3 while-loop error ", 
          resultsRF, "_", Sys.time(), collapse="")
        print(thisstring)
        # flush.console()
        arrResults[thiscounter, j, 7, 5] <- thisstring
        fatalflag <- 1
        arrResults[thiscounter, j, 3, 299] <- fatalflag
      } else {
        
        forest <- resultsRF[[1]]
        arrResults[thiscounter, j, 3, 252] <- resultsRF[[2]]
        arrResults[thiscounter, j, 3, 253] <- resultsRF[[3]]
        arrResults[thiscounter, j, 3, 221] <- resultsRF[[4]]
        arrResults[thiscounter, j, 3, 226] <- resultsRF[[5]]
        arrResults[thiscounter, j, 3, 232] <- resultsRF[[6]]
        outmcmc <- resultsRF[[7]]
        arrResults[thiscounter, j, 3, 220] <- resultsRF[[8]]
        #end else
      } 
      rm(resultsRF)
      # everything else for model 3 needs to be wrapped in "if" statement for 
      # fatal error check
      
      ##----5.2.4.2. Get parameter estimates and variable importance----------####
      
      if (is.na(fatalflag) & is.na(arrResults[thiscounter, j, 3, 232])) {
        
        my5.2.4.2.1 <- function() {
          
          tryCatch(
            expr={
              
              cl <- makeCluster(xtracores)
              registerDoParallel(cl)
              ice.dat <- pdp::partial(forest, pred.var=c("discharge", 
                "temperature"),
                pred.fun=pred.fun.ranger, train=data, 
                chull=TRUE, parallel=TRUE, 
                paropts=list(.packages="ranger"))
              stopCluster(cl)
              results <- betareg(yhat ~ discharge*temperature, data=ice.dat, 
                link="logit")
              return(results)
            }, error=function(error_message){
              
              out <- error_message
              return(out)
            }
            # ends function my5.2.4.2.1
          )}
        results5.2.4.2.1 <- my5.2.4.2.1()
        if (any(class(results5.2.4.2.1) == "error")) {
          
          thisstring <- paste0("Sample ", i, ", ", thisefficiency, 
            "Iteration ", j, "Model 3 5.2.4.2.1 error ", 
            results5.2.4.2.1, "_", Sys.time(), collapse="")
          print(thisstring)
          # flush.console()
          arrResults[thiscounter, j, 7, 6] <- thisstring
          fatalflag <- 1
          arrResults[thiscounter, j, 3, 299] <- fatalflag
        } else {
          
          arrResults[thiscounter, j, 3, 10] <- as.numeric(results5.2.4.2.1$coefficients$mean[2])
          arrResults[thiscounter, j, 3, 13] <- summary(results5.2.4.2.1, type="response")[1]$coefficients$mean[6]
          arrResults[thiscounter, j, 3, 20] <- as.numeric(results5.2.4.2.1$coefficients$mean[3])
          arrResults[thiscounter, j, 3, 23] <- summary(results5.2.4.2.1, type="response")[1]$coefficients$mean[7]
          arrResults[thiscounter, j, 3, 60] <- as.numeric(results5.2.4.2.1$coefficients$mean[4])
          arrResults[thiscounter, j, 3, 63] <- summary(results5.2.4.2.1, type="response")[1]$coefficients$mean[8]
          # end else 5.2.4.2.1
        } 
        rm(results5.2.4.2.1)
        my5.2.4.2.2 <- function() {
          
          tryCatch(
            expr={
              
              cl <- makeCluster(xtracores)
              registerDoParallel(cl)
              ice.dat <- pdp::partial(forest, pred.var="effort", ice=TRUE,
                train=data, chull=TRUE, parallel=TRUE,
                paropts=list(.packages="ranger"))
              stopCluster(cl)
              results <- betareg(yhat ~ poly(effort, 2, raw=TRUE), data=ice.dat, 
                link="logit")
              return(results)
            }, error=function(error_message){
              
              out <- error_message
              return(out)
            }
            # ends function my5.2.4.2.2
          )} 
        results5.2.4.2.2 <- my5.2.4.2.2()
        if(any(class(results5.2.4.2.2) == "error")) {
          
          thisstring <- paste0("Sample ", i, ", ", thisefficiency, 
            "Iteration ", j, "Model 3 5.2.4.2.2 error ", 
            results5.2.4.2.2, "_", Sys.time(), collapse="")
          print(thisstring)
          # flush.console()
          arrResults[thiscounter, j, 7, 7] <- thisstring
          fatalflag <- 1
          arrResults[thiscounter, j, 3, 299] <- fatalflag 
        } else {
          
          arrResults[thiscounter, j, 3, 30] <- as.numeric(results5.2.4.2.2$coefficients$mean[2])
          arrResults[thiscounter, j, 3, 33] <- as.numeric(summary(results5.2.4.2.2)[1]$coefficients$mean[5])
          arrResults[thiscounter, j, 3, 50] <- as.numeric(results5.2.4.2.2$coefficients$mean[3])
          arrResults[thiscounter, j, 3, 53] <- as.numeric(summary(results5.2.4.2.2)[1]$coefficients$mean[6])
          # end else 5.2.4.2.2
        } 
        rm(results5.2.4.2.2)
        my5.2.4.2.3 <- function() {
          
          tryCatch(
            expr={
              
              cl <- makeCluster(xtracores)
              registerDoParallel(cl)
              ice.dat <- pdp::partial(forest, pred.var="precip", ice=TRUE,
                train=data, chull=TRUE, parallel=TRUE,
                paropts=list(.packages="ranger"))
              stopCluster(cl)
              results <- betareg(yhat ~ precip, data=ice.dat, link="logit")
              return(results)
            }, error=function(error_message){
              
              out <- error_message
              return(out)
            }
            # ends function my5.2.4.2.3
          )} 
        results5.2.4.2.3 <- my5.2.4.2.3()
        if(any(class(results5.2.4.2.3) == "error")) {
          
          thisstring <- paste0("Sample ", i, ", ", thisefficiency, 
            "Iteration ", j, "Model 3 5.2.4.2.3 error ", 
            results5.2.4.2.3, "_", Sys.time(), collapse="")
          print(thisstring)
          # flush.console()
          arrResults[thiscounter, j, 7, 8] <- thisstring
          fatalflag <- 1
          arrResults[thiscounter, j, 3, 299] <- fatalflag
        } else {
          arrResults[thiscounter, j, 3, 40] <- as.numeric(results5.2.4.2.3$coefficients$mean[2])
          arrResults[thiscounter, j, 3, 43] <- summary(results5.2.4.2.3)[1]$coefficients$mean[4]
          # end else 5.2.4.2.3
        } 
        rm(results5.2.4.2.3)
        
        ##----5.2.4.3. Final LL fixed and random effects----------------------####
        
        par_phi <- optim(par=1,fn=betaLoss, mu=forest$predictions, 
          y=data$y, method="Brent", lower=0, upper=100000)$par
        arrResults[thiscounter, j, 3, 225] <- sum(betaLogLik(y=data$y, 
          mu=forest$predictions, 
          phi=par_phi))
        
        ##----5.2.4.4. rmse---------------------------------------------------####
        
        arrResults[thiscounter, j, 3, 80:83] <- c(mean(outmcmc[,"phi"], 
          na.rm=TRUE),
          quantile(outmcmc[,"phi"], 0.025, 
            na.rm=TRUE),
          quantile(outmcmc[,"phi"], 0.975, 
            na.rm=TRUE), 
          sd(outmcmc[,"phi"], na.rm=TRUE))
        testset <- testdat_lab
        predcolsm <- predict(forest, data=testset, type="quantiles", 
          quantiles=0.5, num.threads=xtracores)
        predcolsq <- predict(forest, data=testset, type="quantiles", 
          quantiles=c(0.025, 0.975), num.threads=xtracores)
        
        ##----5.2.4.5. CI width-----------------------------------------------####
        
        dist <- abs(predcolsq$predictions[,2] - predcolsq$predictions[,1])
        arrResults[thiscounter, j, 3, 230] <- mean(dist)
        
        ##----5.2.4.6. Failure rate-------------------------------------------####
        
        testset$obs0.025 <- qbeta(0.025, testset$a,
          testset$b,
          lower.tail=TRUE)
        testset$obs0.975 <- qbeta(0.975, testset$a,
          testset$b,
          lower.tail=TRUE)
        testset$yhat0.025 <- predcolsq$predictions[,1]
        testset$yhat0.975 <- predcolsq$predictions[,2]
        predfail <- subset(testset,
          (yhat0.025 > obs0.975) | (yhat0.975 < obs0.025))
        
        ##----5.2.4.7. Store results------------------------------------------####
        
        # variable importance and p-value
        temp <- importance_pvalues(forest, method="altmann", num.permutations=100, 
          formula=y ~ discharge + temperature + effort + 
            precip, data=data)
        arrResults[thiscounter, j, 3, 15] <- forest$variable.importance[1]
        arrResults[thiscounter, j, 3, 16] <- temp[1,2]
        arrResults[thiscounter, j, 3, 25] <- forest$variable.importance[2]
        arrResults[thiscounter, j, 3, 26] <- temp[2,2]
        arrResults[thiscounter, j, 3, 35] <- forest$variable.importance[3]
        arrResults[thiscounter, j, 3, 36] <- temp[3,2]
        arrResults[thiscounter, j, 3, 45] <- forest$variable.importance[4]
        arrResults[thiscounter, j, 3, 46] <- temp[4,2]
        rm(temp)
        ## Friedman's H statistic
        keep <- c("discharge", "temperature", "effort", "precip")
        X <- data[,keep]
        summary(X)
        thispredictor <- iml::Predictor$new(model=forest, 
          predict.function = pred.fun.ranger,
          data=X,
          y=data$y)
        interact.discharge <- iml::Interaction$new(thispredictor, 
          feature="discharge")
        arrResults[thiscounter, j, 3, 66] <- interact.discharge$results$.interaction[1]
        # parameter estimates
        arrResults[thiscounter, j, 3, 70:73] <- c(mean(outmcmc[,"yr.sd"], 
          na.rm=TRUE),
          quantile(outmcmc[,"yr.sd"], 
            0.025, na.rm=TRUE),
          quantile(outmcmc[,"yr.sd"], 
            0.975, na.rm=TRUE), 
          sd(outmcmc[,"yr.sd"], 
            na.rm=TRUE))
        # prediction
        arrResults[thiscounter, j, 3, 210] <- rmse(testset$y, 
          predcolsm$predictions[,1], 
          phiTrue)
        arrResults[thiscounter, j, 3, 200] <- nrow(predfail)/nrow(testset)
        rm(testset)
        rm(outmcmc)
        rm(forest)
      }
      print(paste("Completed Model 3 at ", Sys.time()))
      # flush.console()
      
      ##----Save simulation results/print to console
      
      if(j %% 5 == 0) {
        
        filename <- paste0('output/R_ss',thiscounter,'_',j,'.RData', collapse='')
        save(arrResults, file=filename)
      }
      thisstring <- paste0("Completed one simulation: Sample size: ", i, 
        ", Efficiency: ", thisefficiency, 
        ", Iteration: ", j, ". ", Sys.time(), collapse="")
      print(thisstring)
      # flush.console()
      endtime <- Sys.time()
      dtime <- difftime(endtime, starttime, units="secs")
      arrResults[thiscounter, j, 6, 2] <- dtime
    }
  }
}


print(Sys.time())

#12 min 20 sec

##----6. Save output----
## Save output array
thisstring <- paste0('output/R_nsim_', nsim, '_nset_', nset, '.RData', collapse='')
save(arrResults, file=thisstring)
# sessionInfo()
