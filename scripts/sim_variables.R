####
####----1. ABOUT------------------------------------------------------------####
####
####   Evaluating alternative methods for modeling trap         
####   efficiencies of outmigrating juvenile salmonids          
####
####   Authors:  Walden et al.                               
####
####   Content: R-Code for sourcing declared variables.
####
####   This is the file for Supplement X of the manuscript.
####

####----2. "True" parameter values for simulated data-----------------------####

## Covariates
Qmin <- 500             # discharge
Qmax <- 7000            # discharge
degCmin <- 6            # temperature
degCmax <- 21           # temperature
effortmin <- 720        # effort
effortmax <- 3840       # effort
precipitation <- 0      # precipitation (yes/no) (nuisance)
nyear <- 10             # groups (random effects)

## Set parameter values
Beta0_maxpt1 <- -3.2    # Trap efficiency max 0.10
Beta0_maxpt25 <- -2.1   # Trap efficiency max 0.25
Beta1 <- -0.3       # discharge
Beta2 <- 0.2        # temperature
Beta3 <- 0.2        # effort
Beta4 <- 0.0        # precipitation yes/no
Beta5 <- -0.4       # effort^2
Beta6 <- 0.4        # discharge x temperature interaction
raneffSD <- 0.1      # year effect standard deviation
phiTrue <- 0.0002

####----3. Simulation variables---------------------------------------------####

# vector of sample sizes to test
nset <- c(50,100,150,200) 
# vector of trap efficiencies to test
theseefficiencies <- c("Max efficiency 0.1","Max efficiency 0.25")
# Number of simulations to complete
nsim <- 1000              
# number of MCMC samples to extract for marginal predictions
nsample <- 1000           
# size of test dataset
ntest <- 100

####----4. Model parameters, hyperparameters, or variables------------------####

## For boosted regression modeling
seed.cvrisk <- 143180916
mutuneiter <- 2000
phituneiter <- 2000
nu_hp=0.1
fmu <- as.formula(y ~ 1 + 
                    bols(discharge, intercept=FALSE) + 
                    bols(temperature, intercept=FALSE) + 
                    bols(effort, intercept=FALSE) + 
                    bols(precip, intercept=FALSE) + 
                    bbs(effort, knots=2, boundary.knots=NULL, degree=2, 
                        df=1, center=TRUE) + 
                    bols(temperature, intercept=FALSE, df=1) %X% 
                    bols(discharge, intercept=FALSE, df=1) + 
                    brandom(year)
)
fphi <- as.formula(y ~ 1)
nboot <- 1000
kfoldB <- 25

## For random forest modeling
seed.forest <- 234672817
rftrees <- 2000