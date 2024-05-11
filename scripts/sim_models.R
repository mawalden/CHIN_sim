####----ABOUT---------------------------------------------------------------####
####   
####   Evaluating alternative methods for modeling trap         
####   efficiencies of outmigrating juvenile salmonids          
####
####   Authors:  Walden et al.                               
####
####   Content: JAGS/Nimble code for text files of models      
####
####   This is the file with all models used in the manuscript.
####

####
####----1. Nimble NULL -----------------------------------------------------####
####

cat(file = "jagsnull.txt","
model{
  ##----PRIORS----
  phi ~ dgamma(0.001, 0.001)

  # Main effects
  p0 ~ dbeta(1,1)                     # Intercept for mean
  b_p0 <- logit(p0)                   # Derived, for random effect estimation

  # Random effect
  yr.sd ~ dunif(0, 10)                # Prior
  for(w in 1:nyr){
    b.year[w] ~ dnorm(b_p0, sd = yr.sd)  # Hierarchical centering around intercept
  }
  
  ##----LIKELIHOOD----
  
  for (q in 1:length_y) {
    logit(p[q]) <- b.year[itsyr[q]]
    a[q] <- p[q] / phi
    b[q] <- (1 - p[q]) / phi
    y[q] ~ dbeta(a[q], b[q])
    
    llobs[q] <- dbeta(y[q], a[q], b[q], log = TRUE)
  }
    ll <- sum(llobs[1:length_y])
}")

####
####----2. Nimble FULL 5 covariates-----------------------------------------####
####

cat(file = "jagsfull.txt","
model{
  ##----PRIORS----
  phi ~ dgamma(0.001, 0.001)

  # Main effects
  p0 ~ dbeta(1,1)
  b_p0 <- logit(p0)

  for(w in 1:6) {
    Beta[w] ~ dnorm(0, 1.75^(-2))
  }
  
  # Random effect
  yr.sd ~ dunif(0, 10)
  for(d in 1:nyr){
    b.year[d] ~ dnorm(b_p0, sd = yr.sd)
  }
  
  ##----LIKELIHOOD----
  
  for (q in 1:length_y) {
    logit(p[q]) <- sum(Beta[1:6] * covs[q,1:6]) + b.year[itsyr[q]]
    a[q] <- p[q] / phi
    b[q] <- (1 - p[q]) / phi
    y[q] ~ dbeta(a[q], b[q])
    
    llobs[q] <- dbeta(y[q], a[q], b[q], log = TRUE)
  }
    ll <- sum(llobs[1:length_y])

}")

####
####----3. Nimble known fixed effects and unknown random effects------------####
####

cat(file = "jagsRF.txt","
model{
  ##----PRIORS----
  phi ~ dgamma(0.001, 0.001)

  # Random effect
  yr.sd ~ dunif(0, 10)
  for(w in 1:nyr){
    b.year[w] ~ dnorm(0, sd = yr.sd)
  }
  
  ##----LIKELIHOOD----
  
  for (q in 1:length_y) {
    logit(p[q]) <- fixedeff[q] + b.year[itsyr[q]]
    a[q] <- p[q] / phi
    b[q] <- (1 - p[q]) / phi
    y[q] ~ dbeta(a[q], b[q])
    
    llobs[q] <- dbeta(y[q], a[q], b[q], log = TRUE)
  }
    ll <- sum(llobs[1:length_y])

}")