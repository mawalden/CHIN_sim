####
####----ABOUT------------------------------------------------------------####
####
####   Evaluating alternative methods for modeling trap         
####   efficiencies of outmigrating juvenile salmonids          
####
####   Authors:  Walden et al.                               
####
####   Content: R-Code for summary metrics
####
####   This is the file for Supplement X of the manuscript.
####
####
####

##------X. Load libraries-------------------------------------------------####
library(ggplot2)

##------X. Set working directory------------------------------------------####

# setwd("C:/Users/mwalden/OneDrive - DOI/CMR_recapture_rate_modeling/analyses")

options(warn=1)       # print all warnings
rm(list=ls())         # clear global environment

##------X. Source variables-------------------------------------------------####

source("scripts/sim_variables.R")

####
####----X.X. Lookup table for true parameter values-------------------------####
####

truevals <- data.frame("efficiency" = c(rep("Max efficiency 0.1", 9), 
                                        rep("Max efficiency 0.25", 9)),
                       "parameter" = rep(c("intercept", "discharge", "degreesC", 
                                           "effort", "precip",  "effort.sqrd", 
                                           "dischargeXdegreesC", "yr.sd", 
                                           "phi"), 2),
                       "value" = c(Beta0_maxpt1, Beta1, Beta2, Beta3, Beta4, 
                         Beta5, Beta6, raneffSD, phiTrue,
                         Beta0_maxpt25, Beta1, Beta2, Beta3, Beta4, Beta5, 
                         Beta6, raneffSD, phiTrue))

head(truevals)

####
####----X.X. Create summary dataframe---------------------------------------####
####

metricsnames <- c("intercept", "discharge", "degreesC", "effort", "precip",
                  "effort.sqrd", "dischargeXdegreesC", "yr.sd", "phi",
                  "LL.null", "ll.final", "LL.null.random", "ll.final.random",
                  "failure.rate", "RMSE", "avg.CI.width", "convergence.fails", 
                  "R2D")
resultsSave <- data.frame("sample size" = c(rep(nset[1], length(metricsnames)*2*3),
  rep(nset[2], length(metricsnames)*2*3),
  rep(nset[3], length(metricsnames)*2*3),
  rep(nset[4], length(metricsnames)*2*3)),
  "efficiency" = c(rep(c(rep("Max efficiency 0.1", length(metricsnames)*3), 
  rep("Max efficiency 0.25", length(metricsnames)*3)), 4)),
  "model" = rep(c(rep("Bayesian GLM", length(metricsnames)),
                  rep("Boosted beta regression", length(metricsnames)),
                  rep("Beta mixed forest", length(metricsnames))), 8),
  "metric"= rep(metricsnames, 24),
  "value"=NA, "std.dev"=NA, "std.err"=NA, "LCI0.025"=NA, "UCI0.975"=NA, 
  "f.stat"=NA, "p.val"=NA,  "incl.0"=NA, "var.imp"=NA, "H.val"=NA,
  "percent.recovered"=NA)
head(resultsSave,20)
resultsSave$modelAlias <- ifelse(resultsSave$model == "Bayesian GLM", 1, 
  ifelse(resultsSave$model == "Boosted beta regression", 2, 3))

####
####----X.X. Load results array(s)------------------------------------------####
####

thearray <- array(data=NA, dim=c(length(nset)*2, nsim, 7, 300))

# thearray <- readRDS("output/first1_200run16042024.rds")
d1 <- get(load('output/jupyter_nsim_20_nset_50.RData'))
dim(d1)
d2 <- get(load('output/jupyter_nsim_20_nset_100.RData'))
dim(d1)
d3 <- get(load('output/jupyter_nsim_50_nset_150.RData'))
dim(d3)
d4 <- get(load('output/jupyter_nsim_50_nset_200.RData'))
dim(d4)
rm(arrResults)

####
####----X.X. Calculate summaries--------------------------------------------####
####

samplesizes <- nset
efficiencies <- theseefficiencies
models <- c(1, 2, 3)
covars <- metricsnames[2:7]

thearray <- d1
thearray <- d2
thearray <- d3
thearray <- d4

# i <- 50
# for(i in samplesizes) {

  # j <- "Max efficiency 0.1"
  # j <- "Max efficiency 0.25"
  for(j in efficiencies) {
    thiscounter <- which(thearray[,1,4,1] == i &
                           thearray[,1,5,1] == j)
    if(length(thiscounter)==0) next
    # k <- 1
    # k <- 2
    # k <- 3
    # thiscounter <- 1
    for(k in models) {
      
      thisset <- thearray[thiscounter, , k,]
      class(thisset) <- "numeric"
      
      ##----X.X.X. LL, r^2, convergence---------------------------------------
      if(k == 1) {
        # Number of convergence failures of full model
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="convergence.fails")] <- 
          sum(thisset[,232], na.rm=TRUE)
        # remove any simulations with convergence failures of full model
        # or fatal errors
        discard <- which(thisset[,232]==1 | thisset[,299]==1)
        if(length(discard) > 0) thisset <- thisset[-discard,]
        # mean null LL where null LL successfully converged
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="LL.null")] <- 
          mean(thisset[which(thisset[,231]!=1),220], na.rm=TRUE)
        # mean full LL
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="ll.final")] <- 
          mean(thisset[,225], na.rm=TRUE)
        # R2D
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="R2D")] <- 
          mean(1 - (-2*thisset[which(thisset[,231]!=1),225]) / 
                 (-2*thisset[which(thisset[,231]!=1),220]), na.rm=TRUE)
      } # end k = 1
      if(k == 2) {
        # remove any simulations with fatal errors
        discard <- which(thisset[,299]==1)
        if(length(discard) > 0) thisset <- thisset[-discard,]
        # mean start LL
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="LL.null")] <- 
          mean(thisset[,220], na.rm=TRUE)
        # mean full LL
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="ll.final")] <- 
          mean(thisset[,225], na.rm=TRUE)
        # R2D
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="R2D")] <- 
          mean(1 - (-2*thisset[,225]) / 
                 (-2*thisset[,220]), na.rm=TRUE)
      } # end k = 2
      if(k == 3) {
        # Number of convergence failures of full model
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="convergence.fails")] <- 
          sum(thisset[,232], na.rm=TRUE)
        # remove any simulations with convergence failures of full model
        # or fatal errors
        discard <- which(thisset[,232]==1 | thisset[,299]==1)
        if(length(discard) > 0) thisset <- thisset[-discard,]
        # mean null fixed LL
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="LL.null")] <- 
          mean(thisset[,220], na.rm=TRUE)
        # mean last fixed LL
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="ll.final")] <- 
          mean(thisset[,225], na.rm=TRUE)
        # mean null random LL
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="LL.null.random")] <- 
          mean(thisset[,221], na.rm=TRUE)
        # mean last fixed LL
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="ll.final.random")] <- 
          mean(thisset[,226], na.rm=TRUE) 
        # R2D
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="R2D")] <- 
          mean(1 - (-2*thisset[,225]) / 
                 (-2*thisset[,220]), na.rm=TRUE)
      } # end k = 3
      
      ##----X.X.X. intercept----------------------------------------------------
      if(k == 1 | k == 2) {
        #get true parameter val
        thistrue <- truevals$value[which(truevals$efficiency==j & 
                                           truevals$parameter=="intercept")]
        resultsSave$value[which(resultsSave$sample.size==i & 
                                           resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="intercept")] <- 
          mean(thisset[,1], na.rm=TRUE)
        resultsSave$std.dev[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                    resultsSave$modelAlias==k & 
                                    resultsSave$metric=="intercept")] <- 
          mean(thisset[,4], na.rm=TRUE)
        resultsSave$LCI0.025[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                    resultsSave$modelAlias==k & 
                                    resultsSave$metric=="intercept")] <- 
          mean(thisset[,2], na.rm=TRUE)
        resultsSave$UCI0.975[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                    resultsSave$modelAlias==k & 
                                    resultsSave$metric=="intercept")] <- 
          mean(thisset[,3], na.rm=TRUE)
        resultsSave$percent.recovered[which(resultsSave$sample.size==i &
                                    resultsSave$efficiency==j &
                                      resultsSave$modelAlias==k & 
                                      resultsSave$metric=="intercept")] <-
          length(which(thistrue >= thisset[,2] & thistrue <= thisset[,3]))/nrow(thisset)
      } #end model == 1 or 2 
      if(k == 2) {
        resultsSave$var.imp[which(resultsSave$sample.size==i & 
                                    resultsSave$efficiency==j & 
                                    resultsSave$modelAlias==k & 
                                    resultsSave$metric=="intercept")] <- 
          mean(thisset[,5], na.rm=TRUE)
        resultsSave$incl.0[which(resultsSave$sample.size==i & 
                                   resultsSave$efficiency==j & 
                                   resultsSave$modelAlias==k & 
                                   resultsSave$metric=="intercept")] <- 
        length(which(thisset[,2] <= 0 & thisset[,3] >= 0))/nrow(thisset)
      } #end model == 2 
      
      
      ##----X.X.X. Betas 1 through 6--------------------------------------------
      
      # m <- 1
      for(m in 1:length(covars)) {
        
        #get true parameter val
        thistrue <- truevals$value[which(truevals$efficiency==j & 
                                           truevals$parameter==covars[m])]
        # all models
        resultsSave$value[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric==covars[m])] <- 
          mean(thisset[,m*10], na.rm=TRUE)

        if(k == 1 | k == 2) {
          resultsSave$std.dev[which(resultsSave$sample.size==i & 
                                      resultsSave$efficiency==j & 
                                      resultsSave$modelAlias==k & 
                                      resultsSave$metric==covars[m])] <- 
            mean(thisset[,m*10+3], na.rm=TRUE)
          resultsSave$LCI0.025[which(resultsSave$sample.size==i & 
                                       resultsSave$efficiency==j & 
                                       resultsSave$modelAlias==k & 
                                       resultsSave$metric==covars[m])] <- 
            mean(thisset[,m*10+1], na.rm=TRUE)
          resultsSave$UCI0.975[which(resultsSave$sample.size==i & 
                                       resultsSave$efficiency==j & 
                                       resultsSave$modelAlias==k & 
                                       resultsSave$metric==covars[m])] <- 
            mean(thisset[,m*10+2], na.rm=TRUE)
          resultsSave$percent.recovered[which(resultsSave$sample.size==i &
                                                resultsSave$efficiency==j &
                                                resultsSave$modelAlias==k & 
                                                resultsSave$metric==covars[m])] <-
            length(which(thistrue >= thisset[,m*10+1] & 
                           thistrue <= thisset[,m*10+2]))/nrow(thisset)
        } # end k = 1 or 2
        
        if(k == 1) {
          resultsSave$f.stat[which(resultsSave$sample.size==i & 
                                     resultsSave$efficiency==j & 
                                     resultsSave$modelAlias==k & 
                                     resultsSave$metric==covars[m])] <- 
            mean(thisset[,m*10+6], na.rm=TRUE)
        } # end k=1
        
        if(k == 2) {
          resultsSave$incl.0[which(resultsSave$sample.size==i & 
                                     resultsSave$efficiency==j & 
                                     resultsSave$modelAlias==k & 
                                     resultsSave$metric==covars[m])] <- 
            length(which(thisset[,m*10+1] <= 0 & thisset[,m*10+2] >= 0))/nrow(thisset)
          resultsSave$var.imp[which(resultsSave$sample.size==i & 
                                      resultsSave$efficiency==j & 
                                      resultsSave$modelAlias==k & 
                                      resultsSave$metric==covars[m])] <- 
            mean(thisset[,m*10+5], na.rm=TRUE)
        } # end k=2
        
        if(k == 3) {
          resultsSave$std.err[which(resultsSave$sample.size==i & 
                                      resultsSave$efficiency==j & 
                                      resultsSave$modelAlias==k & 
                                      resultsSave$metric==covars[m])] <- 
            mean(thisset[,m*10+3], na.rm=TRUE)
          resultsSave$LCI0.025[which(resultsSave$sample.size==i & 
                                       resultsSave$efficiency==j & 
                                       resultsSave$modelAlias==k & 
                                       resultsSave$metric==covars[m])] <- 
            mean(thisset[,m*10] - thisset[,m*10+3]*1.96, na.rm=TRUE)
          resultsSave$UCI0.975[which(resultsSave$sample.size==i & 
                                       resultsSave$efficiency==j & 
                                       resultsSave$modelAlias==k & 
                                       resultsSave$metric==covars[m])] <- 
            mean(thisset[,m*10] + thisset[,m*10+3]*1.96, na.rm=TRUE)
          resultsSave$percent.recovered[which(resultsSave$sample.size==i &
                                                resultsSave$efficiency==j &
                                                resultsSave$modelAlias==k & 
                                                resultsSave$metric==covars[m])] <-
            length(which(thistrue >= (thisset[,m*10] - thisset[,m*10+3]*1.96) & 
                           thistrue <= (thisset[,m*10] + thisset[,m*10+3]*1.96)))/
            nrow(thisset)
          resultsSave$p.val[which(resultsSave$sample.size==i & 
                                    resultsSave$efficiency==j & 
                                    resultsSave$modelAlias==k & 
                                    resultsSave$metric==covars[m])] <- 
            mean(thisset[,m*10+6], na.rm=TRUE)
          resultsSave$var.imp[which(resultsSave$sample.size==i & 
                                      resultsSave$efficiency==j & 
                                      resultsSave$modelAlias==k & 
                                      resultsSave$metric==covars[m])] <- 
            mean(thisset[,m*10+5], na.rm=TRUE)
        } # end k=3
      } # end Betas 1:6

      ##----X.X.X. yr.sd--------------------------------------------------------
      #get true parameter val
      thistrue <- truevals$value[which(truevals$efficiency==j & 
                                         truevals$parameter=="yr.sd")]
      # all models
      resultsSave$value[which(resultsSave$sample.size==i & 
                                resultsSave$efficiency==j & 
                                resultsSave$modelAlias==k & 
                                resultsSave$metric=="yr.sd")] <- 
        mean(thisset[,70], na.rm=TRUE)
      resultsSave$std.dev[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="yr.sd")] <- 
        mean(thisset[,73], na.rm=TRUE)
      resultsSave$LCI0.025[which(resultsSave$sample.size==i & 
                                   resultsSave$efficiency==j & 
                                   resultsSave$modelAlias==k & 
                                   resultsSave$metric=="yr.sd")] <- 
        mean(thisset[,71], na.rm=TRUE)
      resultsSave$UCI0.975[which(resultsSave$sample.size==i & 
                                   resultsSave$efficiency==j & 
                                   resultsSave$modelAlias==k & 
                                   resultsSave$metric=="yr.sd")] <- 
        mean(thisset[,72], na.rm=TRUE)
      resultsSave$percent.recovered[which(resultsSave$sample.size==i &
                                            resultsSave$efficiency==j &
                                            resultsSave$modelAlias==k & 
                                            resultsSave$metric=="yr.sd")] <-
        length(which(thistrue >= thisset[,71] & thistrue <= thisset[,72]))/nrow(thisset)
      if(k == 2) {
        resultsSave$incl.0[which(resultsSave$sample.size==i & 
                                   resultsSave$efficiency==j & 
                                   resultsSave$modelAlias==k & 
                                   resultsSave$metric=="yr.sd")] <- 
          length(which(thisset[,71] <= 0 & thisset[,72] >= 0))/nrow(thisset)
        resultsSave$var.imp[which(resultsSave$sample.size==i & 
                                    resultsSave$efficiency==j & 
                                    resultsSave$modelAlias==k & 
                                    resultsSave$metric=="yr.sd")] <- 
          mean(thisset[,75], na.rm=TRUE)
      } # end k=2
      
      ##----X.X.X. phi----------------------------------------------------------
      #get true parameter val
      thistrue <- truevals$value[which(truevals$efficiency==j & 
                                         truevals$parameter=="phi")]
      # all models
      resultsSave$value[which(resultsSave$sample.size==i & 
                                resultsSave$efficiency==j & 
                                resultsSave$modelAlias==k & 
                                resultsSave$metric=="phi")] <- 
        mean(thisset[,80], na.rm=TRUE)
      resultsSave$std.dev[which(resultsSave$sample.size==i & 
                                  resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="phi")] <- 
        mean(thisset[,83], na.rm=TRUE)
      resultsSave$LCI0.025[which(resultsSave$sample.size==i & 
                                   resultsSave$efficiency==j & 
                                   resultsSave$modelAlias==k & 
                                   resultsSave$metric=="phi")] <- 
        mean(thisset[,81], na.rm=TRUE)
      resultsSave$UCI0.975[which(resultsSave$sample.size==i & 
                                   resultsSave$efficiency==j & 
                                   resultsSave$modelAlias==k & 
                                   resultsSave$metric=="phi")] <- 
        mean(thisset[,82], na.rm=TRUE)
      resultsSave$percent.recovered[which(resultsSave$sample.size==i &
                                            resultsSave$efficiency==j &
                                            resultsSave$modelAlias==k & 
                                            resultsSave$metric=="phi")] <-
        length(which(thistrue >= thisset[,81] & thistrue <= thisset[,82]))/nrow(thisset)
      
      ##----X.X.X. H value, clean up for Model 3 RF-----------------------------
      if(k == 3) {
        resultsSave$H.val[which(resultsSave$sample.size==i & 
                                    resultsSave$efficiency==j & 
                                  resultsSave$modelAlias==k & 
                                  resultsSave$metric=="dischargeXdegreesC")] <- 
          mean(thisset[,66], na.rm=TRUE)
        resultsSave$var.imp[which(resultsSave$sample.size==i &
                                  resultsSave$efficiency==j &
                                    resultsSave$modelAlias==k &
                                    resultsSave$metric=="dischargeXdegreesC")] <- NA
        resultsSave$var.imp[which(resultsSave$sample.size==i &
                                    resultsSave$efficiency==j &
                                    resultsSave$modelAlias==k &
                                    resultsSave$metric=="effort.sqrd")] <- NA
        resultsSave$p.val[which(resultsSave$sample.size==i &
                                    resultsSave$efficiency==j &
                                  resultsSave$modelAlias==k &
                                  resultsSave$metric=="dischargeXdegreesC")] <- NA
        resultsSave$p.val[which(resultsSave$sample.size==i &
                                    resultsSave$efficiency==j &
                                  resultsSave$modelAlias==k &
                                  resultsSave$metric=="effort.sqrd")] <- NA
      } # end k=3
      
      ##----X.X.X. Failure rate, RMSE, mean CI width----------------------------
      resultsSave$value[which(resultsSave$sample.size==i & 
                                resultsSave$efficiency==j & 
                                resultsSave$modelAlias==k & 
                                resultsSave$metric=="failure.rate")] <- 
        mean(thisset[,200], na.rm=TRUE)
      resultsSave$value[which(resultsSave$sample.size==i & 
                                resultsSave$efficiency==j & 
                                resultsSave$modelAlias==k & 
                                resultsSave$metric=="RMSE")] <- 
        mean(thisset[,210], na.rm=TRUE)
      resultsSave$value[which(resultsSave$sample.size==i & 
                                resultsSave$efficiency==j & 
                                resultsSave$modelAlias==k & 
                                resultsSave$metric=="avg.CI.width")] <- 
        mean(thisset[,230], na.rm=TRUE)
      
  } #models
} # efficiencies
# } #sample sizes

resultsSave[which(resultsSave$sample.size==50 & 
                    resultsSave$efficiency=="Max efficiency 0.1" & 
                    resultsSave$modelAlias==3),]

####
####----X.X. Save summaries-------------------------------------------------####
####

write.csv(resultsSave, "output/jupyter_SummaryResults.csv", row.names = F)


####
####----X.X. Figures--------------------------------------------------------####
####

####----Plot X. RMSE------------------------------------------------------------

scale_factor <- 1

rmsedat <- subset(resultsSave, metric=="RMSE")
rmsedat$sample.size <- as.numeric(rmsedat$sample.size)
rmsedat$grp <- paste(rmsedat$modelAlias, "_", rmsedat$efficiency)
# rmsedat$clr <- ifelse(rmsedat$modelAlias==1, "blue", 
#                       ifelse(rmsedat$modelAlias==2, "orange", "purple"))
# rmsedat$clr <- as.factor(rmsedat$clr)

rmseest <- ggplot(rmsedat, aes(x = sample.size, y = value, group = grp,
                               color=grp, shape = grp)) + 
  geom_point(size=3) + 
  labs(x = "Sample size", 
       y = "Root mean square error on test dataset") + 
  scale_color_manual(name = "Model & Efficiency",
                     labels = c("Bayesian GLM, < 0.1",
                                "Bayesian GLM, < 0.25",
                                "Boosted beta regression, < 0.1",
                                "Boosted beta regression, < 0.25",
                                "Beta mixed forest, < 0.1",
                                "Beta mixed forest, < 0.25"),
                     values=c("orange", "orange", 
                                         "purple", "purple", 
                                         "brown","brown")) + 
  scale_shape_manual(name = "Model & Efficiency",
                     labels = c("Bayesian GLM, < 0.1",
                                "Bayesian GLM, < 0.25",
                                "Boosted beta regression, < 0.1",
                                "Boosted beta regression, < 0.25",
                                "Beta mixed forest, < 0.1",
                                "Beta mixed forest, < 0.25"),
                     values=c(15,16,15,16,15,16)) +
  theme(legend.position="top")

  # theme(
  #   text = element_text(size= scale_factor * 8, family="serif" ),
  #   axis.text = element_text(size= scale_factor * 7, family="serif"),
  #   axis.title.x = element_text(hjust=0.5, size= scale_factor * 8, family="serif"),
  #   axis.title.y = element_text(hjust=0.5, size= scale_factor * 8, family="serif"),
  #   legend.title=element_blank(),
  #   legend.text = element_text(size= scale_factor * 7, family="serif"),
  #   legend.position = c(.8, .98),
  #   legend.justification = c("right", "top"),
  #   legend.direction = "horizontal",
  #   legend.key.size = unit(0.2, "cm")
  # )
rmseest

keep <- grid.arrange(rmseest, rmseerr, ncol=2)
ggsave(plot = keep, filename = "output/errorTolerance_rmse.png", 
       scale = scale_factor, width = 6, height = 3, units = "in", dpi = 300)


