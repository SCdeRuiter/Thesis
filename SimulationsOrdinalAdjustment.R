
knitr::opts_chunk$set(echo = TRUE)
library(foreign)
library(dplyr)
library(rms)

set.seed(1) # set seed for reproducibility


## 0. Read in data and fit a model to the data

# Read in data 
MP <- read.spss("MP trial database complete cases n 221.sav", use.value.label = TRUE, to.data.frame = TRUE)

MP$treat <- droplevels(MP$treat, exclude = if(anyNA(levels(MP$treat))) NULL else NA) # Drop empty levels
MP$treat <- recode(MP$treat, IVIg = "IVIg+MP") # Correct labels
MP$treat <- recode(MP$treat, plasmaferese = "IVIg") 

ReverseGBS <- function(x) { 
  # Function to:
  # 1) Combine extreme categories ('death' with 'requiring assisted ventilation'; and 'healthy state' with 'minor symptoms')
  # 2) Reverse the scale such that higher GBSDS = better outcome
  Rec <- x
  Rec[x == 0] <- 5 
  Rec[x == 1] <- 5
  Rec[x == 2] <- 4 
  Rec[x == 4] <- 2
  Rec[x == 5] <- 1 
  Rec[x == 6] <- 1
  return(Rec)
}

# Put data in a data frame 
MPdata <- data.frame(GBSDS = ReverseGBS(MP[,"f9"]),  # GBSDS score at 4 weeks
                     treat = MP$treat, # treatment indicator
                     mrcss = MP$mrcss, # MRC sum score
                     weakonse = MP$weakonse, 
                     age = MP$age, 
                     diarrhea = as.numeric(MP$diarrhea)-1) # days from onset of weakness until randomisation
 
MPdata <- MPdata[-165, ] # the 165th row was a subject with missing weakonse; remove subject


# Fit a model to the data (coefficients will be used for simulation of data later)
## Original covariates
dd <- datadist(MPdata); options(datadist = "dd")
MPfit <- lrm(GBSDS ~ treat + mrcss + weakonse, data = MPdata)
MPcoefficients <- coef(MPfit)

 # OTHER COVS
MPfit.agediar <- lrm(GBSDS ~ treat + age + diarrhea, data = MPdata)
MPcoefficients.agediar <- coef(MPfit.agediar)

MPcoefficients # coefficients on logit scale
MPcoefficients.agediar

dd <- datadist(MPdata); options(datadist = "dd")

## Condition A: Betas for covariate strength as observed, assuming proportional odds (2 covariates: mrcss and days from weakness to randomisation)
### A.1 Function for simulation (Wald)

SimFun <- function(nsim, groupsize, txeffectmort, txeffect) {
  
  # input arguments:
  ## nsim = number of iterations per simulation 
  ## groupsize = the number of subjects per treatment arm 
  ## txeffectmort = the treatment effect for mortality (cut-off GBS DS 01 vs GBSDS 23456)
  ## txeffect = the treatment effect for remaining cut-offs
  ## note: this function can be used to simulate both under satisfaction and violation of the PO assumption
  
  # empty objects for storage
  ORs <- matrix(nrow = nsim, ncol = 4)
  SEs <-  matrix(nrow = nsim, ncol = 4)
  Pvals <- matrix(nrow = nsim, ncol = 4)
  Rej <-  matrix(nrow = nsim, ncol = 4)
  
  for (r in 1:nsim) {
    
    if (r %% 500 == 0) cat(paste0("Iteration: ", r, "\n")) # counter to track progression
    
    # Make X matrix for sample data set
    ind <- sample(1:nrow(MPdata), size = 2 * groupsize, replace = T) # random sampling with replacement (obtain indices of observations that will be included in sample)
    
    SampleXmat <- MPdata[ind, c("mrcss", "weakonse")] # obtain X matrix for sampled observations
    SampleXmat$treat <- c(rep(0, groupsize), rep(1, groupsize)) # assign patients to a treatment group (two arms, balanced)
    
    # To get GBS DS for sampled patients:
    # Calculate linear predictors for each patient, for each cut off  
    Linpreds <- matrix(nrow  = 2*groupsize, ncol = 4) # empty object for storage
    Linpreds[, 1] <- MPcoefficients[1] + txeffectmort * SampleXmat$treat +          
                                          MPcoefficients["mrcss"] * SampleXmat$mrcss  +
                                          MPcoefficients["weakonse"] * SampleXmat$weakonse
    
    Linpreds[, 2:4] <- sapply(2:4, function(x) {MPcoefficients[x] + txeffect * SampleXmat$treat +         
                                          MPcoefficients["mrcss"] * SampleXmat$mrcss  +
                                          MPcoefficients["weakonse"] * SampleXmat$weakonse})
    
    probs <- matrix(nrow = nrow(SampleXmat), ncol = 5) # empty object for storage
    GBSDS <- numeric(nrow(SampleXmat)) # empty object for storage
    
    for (i in 1:nrow(SampleXmat)) { # For each subject, calculate probability for each category
      probs[i, 1] <- 1 - plogis(Linpreds[i, 1])                           # P(GBSDS = 0/1)
      probs[i, 2] <- plogis(Linpreds[i, 1]) - plogis(Linpreds[i, 2])      # P(GBSDS = 2)
      probs[i, 3] <- plogis(Linpreds[i, 2]) - plogis(Linpreds[i, 3])      # P(GBSDS = 3)
      probs[i, 4] <- plogis(Linpreds[i, 3]) - plogis(Linpreds[i, 4])      # P(GBSDS = 4)
      probs[i, 5] <- plogis(Linpreds[i, 4])                               # P(GBSDS = 5/6)
      
      GBSDS[i] <-  which(rmultinom(1, 1, (probs[i, ]))==1)                # sample a GBS DS score from multinomial distribution (using the probabilities)
                                                                          # scale is 1 through 5 (1 standing for 0 and 1; 5 standing for 5 and 6)
    }
    
    # compile data set
    comp <- cbind(SampleXmat, GBSDS)
    dd <<- datadist(comp) # the <<- operator assigns the variable to the global environment
    options(datadist = "dd")
    
    # Fit dichotomous models (binary log. reg.)
    moddich.unadj <- lrm(as.numeric(GBSDS > 3) ~ treat, data = comp)                  # unadjusted
    moddich.adj <- lrm(as.numeric(GBSDS > 3) ~ treat + mrcss + weakonse, data = comp) # adjusted
    
    # Fit ordinal models (proportional odds log. reg.)
    modPO.unadj <- lrm(factor(GBSDS) ~ treat, data = comp)                            # unadjusted 
    modPO.adj <- lrm(factor(GBSDS) ~ treat + mrcss + weakonse, data = comp)           # adjusted

    # For each model: extract OR, SE, p val of tx effect and indication of hypothesis rejection (0/1)    
    ORs[r, ] <- exp(c(summary(moddich.unadj)['treat', 'Effect'],
                  summary(moddich.adj)['treat', 'Effect'],
                  summary(modPO.unadj)['treat', 'Effect'],
                  summary(modPO.adj)['treat', 'Effect']))
    
    SEs[r, ] <- c(summary(moddich.unadj)['treat', 'S.E.'], 
                  summary(moddich.adj)['treat', 'S.E.'],
                  summary(modPO.unadj)['treat', 'S.E.'],
                  summary(modPO.adj)['treat', 'S.E.'])
    
    Pvals[r, ] <- 1 - pchisq(log(ORs[r,])^2 / SEs[r,]^2, 1)
    
    Rej[r, ] <- as.numeric(Pvals[r, ] < 0.05)
  }
  
  # For each model compute average Z value 
  ZU.dich <- mean(log(ORs[,1]) / SEs[,1])
  ZA.dich <- mean(log(ORs[,2]) / SEs[,2])
  ZU.PO <- mean(log(ORs[,3]) / SEs[,3])
  ZA.PO <- mean(log(ORs[,4]) / SEs[,4])
  
  RSS.covadj.dich <- 100 - (100 * (ZU.dich/ZA.dich)^2)      # Reduction in sample size by adjustment (in dich. models)
  RSS.covadj.PO <- 100 - (100 * (ZU.PO/ZA.PO)^2)            # Reduction in sample size by adjustment (in ord. models)
  
  RSS.PO.unadj <- 100 - (100 * (ZU.dich/ZU.PO)^2)           # Reduction in sample size by exploiting ordinality (in unadjusted models)
  RSS.PO.adj <- 100 - (100 * (ZA.dich/ZA.PO)^2)             # Reduction in sample size by exploiting ordinality (in adjusted models)
  
  RSS.adjPO.unadjdich <- 100 - (100 * (ZU.dich / ZA.PO)^2)  # Reduction in sample size by adjusted PO as compared to unadjusted dich
  
  MCE <- sapply(1:4, function(x) {sqrt((colMeans(Rej)[x] * (1-colMeans(Rej)[x])) / nsim)}) # monte carlo error of the power estimate
  
  objects <- c(ORs=ORs, SEs=SEs, Pvals=Pvals, Rej=Rej, RejectionRate=colMeans(Rej),
                    RSS.covadj.PO = RSS.covadj.PO, RSS.covadj.dich = RSS.covadj.dich,
                    RSS.PO.adj = RSS.PO.adj, RSS.PO.unadj = RSS.PO.unadj,
                    RSS.adjPO.unadjdich = RSS.adjPO.unadjdich)
  
  return(list(matrix(data = c(colMeans(Rej)*100, colMeans(ORs), colMeans(SEs), NA, RSS.covadj.dich, RSS.PO.unadj, RSS.adjPO.unadjdich), 
                     ncol = 4, byrow = TRUE, dimnames = list(c("Rejection Rate (%)", "Mean OR", "Mean SE", "Red in SS (%)"),
                                                             c("Unadjusted dich", "Adjusted dich", "Unadjusted PO", "Adjusted PO"))), 
              objects,
              MCE = MCE))
}


### A.2 Perform simulation (Wald)

# 60 n
T1error.60n <- SimFun(nsim = 10000, groupsize = 60, txeffectmort = log(1), txeffect = log(1))

txlog1.5.60n <- SimFun(nsim = 10000, groupsize = 60, txeffectmort = log(1.5), txeffect = log(1.5))

txlog2.1.60n <- SimFun(nsim = 10000, groupsize = 60, txeffectmort = log(2.1), txeffect = log(2.1))

txlog2.7.60n <- SimFun(nsim = 10000, groupsize = 60, txeffectmort = log(2.7), txeffect = log(2.7))

# 100 n
T1error.100n <- SimFun(nsim = 10000, groupsize = 100, txeffectmort = log(1), txeffect = log(1))

txlog1.5.100n <- SimFun(nsim = 10000, groupsize = 100, txeffectmort = log(1.5), txeffect = log(1.5))

txlog2.1.100n <- SimFun(nsim = 10000, groupsize = 100, txeffectmort = log(2.1), txeffect = log(2.1))

txlog2.7.100n <- SimFun(nsim = 10000, groupsize = 100, txeffectmort = log(2.7), txeffect = log(2.7))

### A.3 Print results (Wald)

cat("Sample size:", 120, "(60/arm)\n")
cat("OR = 1 \n")
T1error.60n[c(3, 1)]

cat("OR = 1.5 \n")
txlog1.5.60n[c(3, 1)]

cat("OR = 2.1 \n")
txlog2.1.60n[c(3, 1)]

cat("OR = 2.7 \n")
txlog2.7.60n[c(3, 1)]

cat("Sample size:", 200, "(100/arm)\n")
cat("OR = 1 \n")
T1error.100n[c(3, 1)]

cat("OR = 1.5 \n")
txlog1.5.100n[c(3, 1)]

cat("OR = 2.1 \n")
txlog2.1.100n[c(3, 1)]

cat("OR = 2.7 \n")
txlog2.7.100n[c(3, 1)]


### A.1b Function for simulation (LRT) 

SimFunLRT <- function(nsim, groupsize, txeffectmort, txeffect) {
  
  # input arguments:
  ## nsim = number of iterations per simulation 
  ## groupsize = the number of subjects per treatment arm 
  ## txeffectmort = the treatment effect for mortality (cut-off GBS DS 01 vs GBSDS 23456)
  ## txeffect = the treatment effect for remaining cut-offs
  
  # empty objects for storage
  ORs <- matrix(nrow = nsim, ncol = 4)
  SEs <-  matrix(nrow = nsim, ncol = 4)
  Pvals <- matrix(nrow = nsim, ncol = 4)
  Rej <-  matrix(nrow = nsim, ncol = 4)
  
  for (r in 1:nsim) {
    
    if (r %% 500 == 0) cat(paste0("Iteration: ", r, "\n")) # counter to track progression
    
    # Make X matrix for sample data set
    ind <- sample(1:nrow(MPdata), size = 2 * groupsize, replace = T) # random sampling with raplacement (obtain indices of observations that will be included in sample)
    
    SampleXmat <- MPdata[ind, c("mrcss", "weakonse")] # obtain X matrix for sampled observations
    SampleXmat$treat <- c(rep(0, groupsize), rep(1, groupsize)) # assign patients to a treatment group (two arms, balanced)
    
    # To get GBS DS for sampled patients:
    # Calculate linear predictors for each patient, for each cut off 
    Linpreds <- matrix(nrow  = 2*groupsize, ncol = 4) # empty object for storage
    
    Linpreds[, 1] <- MPcoefficients[1] + txeffectmort * SampleXmat$treat +          
                                          MPcoefficients["mrcss"] * SampleXmat$mrcss  +
                                          MPcoefficients["weakonse"] * SampleXmat$weakonse
    
    Linpreds[, 2:4] <- sapply(2:4, function(x) {MPcoefficients[x] + txeffect * SampleXmat$treat +         
                                          MPcoefficients["mrcss"] * SampleXmat$mrcss  +
                                          MPcoefficients["weakonse"] * SampleXmat$weakonse})
    
    probs <- matrix(nrow = nrow(SampleXmat), ncol = 5) # empty object for storage
    GBSDS <- numeric(nrow(SampleXmat)) # empty object for storage
    
    for (i in 1:nrow(SampleXmat)) { # For each subject, calculate probability for each category
      probs[i, 1] <- 1 - plogis(Linpreds[i, 1])                           # P(GBSDS = 0/1)
      probs[i, 2] <- plogis(Linpreds[i, 1]) - plogis(Linpreds[i, 2])      # P(GBSDS = 2)
      probs[i, 3] <- plogis(Linpreds[i, 2]) - plogis(Linpreds[i, 3])      # P(GBSDS = 3)
      probs[i, 4] <- plogis(Linpreds[i, 3]) - plogis(Linpreds[i, 4])      # P(GBSDS = 4)
      probs[i, 5] <- plogis(Linpreds[i, 4])                               # P(GBSDS = 5/6)
      
      GBSDS[i] <-  which(rmultinom(1, 1, (probs[i, ]))==1)                # sample a GBS DS score from multinomial distribution (using the probabilities)
                                                                          # scale is 1 through 5 (1 standing for 0 and 1; 5 standing for 5 and 6)
    }
    
    # compile data set
    comp <- cbind(SampleXmat, GBSDS)
    dd <<- datadist(comp) # the <<- operator assigns the variable to the global environment
    options(datadist = "dd")
    
    # Fit dichotomous models (binary log. reg.)
    moddich.unadj <- lrm(as.numeric(GBSDS > 3) ~ treat, data = comp)                  # unadjusted
    moddich.adj <- lrm(as.numeric(GBSDS > 3) ~ treat + mrcss + weakonse, data = comp) # adjusted
    
    # Fit ordinal models (proportional odds log. reg.)
    modPO.unadj <- lrm(factor(GBSDS) ~ treat, data = comp)                            # unadjusted 
    modPO.adj <- lrm(factor(GBSDS) ~ treat + mrcss + weakonse, data = comp)           # adjusted
    
    # Fit nested models in order to obtain LRT p value
    ## dich adj
    moddich.adj.nested <- lrm(as.numeric(GBSDS > 3) ~ mrcss + weakonse, data = comp) 
    ll.dich.nested <- logLik(moddich.adj.nested)
    ll.dich.complex <- logLik(moddich.adj)   
    teststat.dich <- -2*(as.numeric(ll.dich.nested)-as.numeric(ll.dich.complex))
    pval.dich <- pchisq(teststat.dich, df = 1, lower.tail = F)
    
    ## po adj
    modPO.adj.nested <- lrm(factor(GBSDS) ~  mrcss + weakonse, data = comp)  
    ll.PO.nested <- logLik(modPO.adj.nested)
    ll.PO.complex <- logLik(modPO.adj)
    teststat.PO <- -2*(as.numeric(ll.PO.nested)-as.numeric(ll.PO.complex))
    pval.PO <- pchisq(teststat.PO, df = 1, lower.tail = F)

    # For each model: extract OR, SE, p val of tx effect and indication of hypothesis rejection (0/1)    
    ORs[r, ] <- exp(c(summary(moddich.unadj)['treat', 'Effect'],
                  summary(moddich.adj)['treat', 'Effect'],
                  summary(modPO.unadj)['treat', 'Effect'],
                  summary(modPO.adj)['treat', 'Effect']))
    
    SEs[r, ] <- c(summary(moddich.unadj)['treat', 'S.E.'], 
                  summary(moddich.adj)['treat', 'S.E.'],
                  summary(modPO.unadj)['treat', 'S.E.'],
                  summary(modPO.adj)['treat', 'S.E.'])
    
    Pvals[r, ] <- c(moddich.unadj$stats[5],
                    pval.dich,
                    modPO.unadj$stats[5], 
                    pval.PO)
    
    Rej[r, ] <- as.numeric(Pvals[r, ] < 0.05)
  }
  
  # For each model compute average Z value
  ZU.dich <- mean(log(ORs[,1]) / SEs[,1])
  ZA.dich <- mean(log(ORs[,2]) / SEs[,2])
  ZU.PO <- mean(log(ORs[,3]) / SEs[,3])
  ZA.PO <- mean(log(ORs[,4]) / SEs[,4])
  
  RSS.covadj.dich <- 100 - (100 * (ZU.dich/ZA.dich)^2)      # Reduction in sample size by adjustment (in dich. models)
  RSS.covadj.PO <- 100 - (100 * (ZU.PO/ZA.PO)^2)            # Reduction in sample size by adjustment (in ord. models)
  
  RSS.PO.unadj <- 100 - (100 * (ZU.dich/ZU.PO)^2)           # Reduction in sample size by exploiting ordinality (in unadjusted models)
  RSS.PO.adj <- 100 - (100 * (ZA.dich/ZA.PO)^2)             # Reduction in sample size by exploiting ordinality (in adjusted models)
  
  RSS.adjPO.unadjdich <- 100 - (100 * (ZU.dich / ZA.PO)^2)  # Reduction in sample size by adjusted PO as compared to unadjusted dich
  
  MCE <- sapply(1:4, function(x) {sqrt((colMeans(Rej)[x] * (1-colMeans(Rej)[x])) / nsim)}) # monte carlo error of the power estimate
  
  objects <- c(ORs=ORs, SEs=SEs, Pvals=Pvals, Rej=Rej, RejectionRate=colMeans(Rej),
                    RSS.covadj.PO = RSS.covadj.PO, RSS.covadj.dich = RSS.covadj.dich,
                    RSS.PO.adj = RSS.PO.adj, RSS.PO.unadj = RSS.PO.unadj,
                    RSS.adjPO.unadjdich = RSS.adjPO.unadjdich)
  
  return(list(matrix(data = c(colMeans(Rej)*100, colMeans(ORs), colMeans(SEs), NA, RSS.covadj.dich, RSS.PO.unadj, RSS.adjPO.unadjdich), 
                     ncol = 4, byrow = TRUE, dimnames = list(c("Rejection Rate (%)", "Mean OR", "Mean SE", "Red in SS (%)"),
                                                             c("Unadjusted dich", "Adjusted dich", "Unadjusted PO", "Adjusted PO"))), 
              objects,
              MCE = MCE))
}



### A.2b Perform simulation (LRT)

# 60 n
T1error.60n.LRT <- SimFunLRT(nsim = 10000, groupsize = 60, txeffectmort = log(1), txeffect = log(1))

txlog1.5.60n.LRT <- SimFunLRT(nsim = 10000, groupsize = 60, txeffectmort = log(1.5), txeffect = log(1.5))

txlog2.1.60n.LRT <- SimFunLRT(nsim = 10000, groupsize = 60, txeffectmort = log(2.1), txeffect = log(2.1))

txlog2.7.60n.LRT <- SimFunLRT(nsim = 10000, groupsize = 60, txeffectmort = log(2.7), txeffect = log(2.7))

# 100 n
T1error.100n.LRT <- SimFunLRT(nsim = 10000, groupsize = 100, txeffectmort = log(1), txeffect = log(1))

txlog1.5.100n.LRT <- SimFunLRT(nsim = 10000, groupsize = 100, txeffectmort = log(1.5), txeffect = log(1.5))

txlog2.1.100n.LRT <- SimFunLRT(nsim = 10000, groupsize = 100, txeffectmort = log(2.1), txeffect = log(2.1))

txlog2.7.100n.LRT <- SimFunLRT(nsim = 10000, groupsize = 100, txeffectmort = log(2.7), txeffect = log(2.7))


### A.3b Print results (LRT)

cat("Sample size:", 120, "(60/arm) (LRT)\n")
cat("OR = 1 \n")
T1error.60n.LRT[c(3, 1)]

cat("OR = 1.5 \n")
txlog1.5.60n.LRT[c(3, 1)]

cat("OR = 2.1 \n")
txlog2.1.60n.LRT[c(3, 1)]

cat("OR = 2.7 \n")
txlog2.7.60n.LRT[c(3, 1)]

cat("Sample size:", 200, "(100/arm) (LRT)\n")
cat("OR = 1 \n")
T1error.100n.LRT[c(3, 1)]

cat("OR = 1.5 \n")
txlog1.5.100n.LRT[c(3, 1)]

cat("OR = 2.1 \n")
txlog2.1.100n.LRT[c(3, 1)]

cat("OR = 2.7 \n")
txlog2.7.100n.LRT[c(3, 1)]

## Condition B: Weaker covariates, assuming proportional odds
### B.1 Function for simulation (Wald)

SimFunWeakerCov <- function(nsim, groupsize, txeffectmort, txeffect) {
  
  # input arguments:
  ## nsim = number of iterations per simulation 
  ## groupsize = the number of subjects per treatment arm 
  ## txeffectmort = the treatment effect for mortality (cut-off GBS DS 01 vs GBSDS 23456)
  ## txeffect = the treatment effect for remaining cut-offs
  
  # empty objects for storage
  ORs <- matrix(nrow = nsim, ncol = 4)
  SEs <-  matrix(nrow = nsim, ncol = 4)
  Pvals <- matrix(nrow = nsim, ncol = 4)
  Rej <-  matrix(nrow = nsim, ncol = 4)
  
  for (r in 1:nsim) {
    
    if (r %% 500 == 0) cat(paste0("Iteration: ", r, "\n")) # counter to track progression
    
    # Make X matrix for sample data set
    ind <- sample(1:nrow(MPdata), size = 2 * groupsize, replace = T) # random sampling with raplacement (obtain indices of observations that will be included in sample)
    
    SampleXmat <- MPdata[ind, c("age", "diarrhea")] # obtain X matrix for sampled observations
    SampleXmat$treat <- c(rep(0, groupsize), rep(1, groupsize)) # assign patients to a treatment group (two arms, balanced)

    # To get GBS DS for sampled patients:
    # Calculate linear predictors for each patient, for each cut off 
    Linpreds <- matrix(nrow  = 2*groupsize, ncol = 4) # empty object for storage
    
    Linpreds[, 1] <- MPcoefficients.agediar[1] + txeffectmort * SampleXmat$treat +          
                                          MPcoefficients.agediar["age"] * SampleXmat$age  +
                                          MPcoefficients.agediar["diarrhea"] * SampleXmat$diarrhea
    
    Linpreds[, 2:4] <- sapply(2:4, function(x) {MPcoefficients.agediar[x] + txeffect * SampleXmat$treat +         
                                          MPcoefficients.agediar["age"] * SampleXmat$age  +
                                          MPcoefficients.agediar["diarrhea"] * SampleXmat$diarrhea})
    
    probs <- matrix(nrow = nrow(SampleXmat), ncol = 5) # empty object for storage
    GBSDS <- numeric(nrow(SampleXmat)) # empty object for storage
    
    for (i in 1:nrow(SampleXmat)) { # For each subject, calculate probability for each category
      probs[i, 1] <- 1 - plogis(Linpreds[i, 1])                           # P(GBSDS = 0/1)
      probs[i, 2] <- plogis(Linpreds[i, 1]) - plogis(Linpreds[i, 2])      # P(GBSDS = 2)
      probs[i, 3] <- plogis(Linpreds[i, 2]) - plogis(Linpreds[i, 3])      # P(GBSDS = 3)
      probs[i, 4] <- plogis(Linpreds[i, 3]) - plogis(Linpreds[i, 4])      # P(GBSDS = 4)
      probs[i, 5] <- plogis(Linpreds[i, 4])                               # P(GBSDS = 5/6)
      
      GBSDS[i] <-  which(rmultinom(1, 1, (probs[i, ]))==1)                # sample a GBS DS score from multinomial distribution (using the probabilities)
                                                                          # scale is 1 through 5 (1 standing for 0 and 1; 5 standing for 5 and 6)
    }
    
    # compile data set
    comp <- cbind(SampleXmat, GBSDS)
    dd <<- datadist(comp) # the <<- operator  assigns the variable to the global environment
    options(datadist = "dd")
    
    # Fit dichotomous models (binary log. reg.)
    moddich.unadj <- lrm(as.numeric(GBSDS > 3) ~ treat, data = comp)                  # unadjusted
    moddich.adj <- lrm(as.numeric(GBSDS > 3) ~ treat + age + diarrhea, data = comp) # adjusted
    
    # Fit ordinal models (proportional odds log. reg.)
    modPO.unadj <- lrm(factor(GBSDS) ~ treat, data = comp)                            # unadjusted 
    modPO.adj <- lrm(factor(GBSDS) ~ treat + age + diarrhea, data = comp)           # adjusted

    # For each model: extract OR, SE, p val of tx effect and indication of hypothesis rejection (0/1)    
    ORs[r, ] <- exp(c(summary(moddich.unadj)['treat', 'Effect'],
                  summary(moddich.adj)['treat', 'Effect'],
                  summary(modPO.unadj)['treat', 'Effect'],
                  summary(modPO.adj)['treat', 'Effect']))
    
    SEs[r, ] <- c(summary(moddich.unadj)['treat', 'S.E.'], 
                  summary(moddich.adj)['treat', 'S.E.'],
                  summary(modPO.unadj)['treat', 'S.E.'],
                  summary(modPO.adj)['treat', 'S.E.'])
    
    Pvals[r, ] <- 1 - pchisq(log(ORs[r,])^2 / SEs[r,]^2, 1)
    
    Rej[r, ] <- as.numeric(Pvals[r, ] < 0.05)
  }
  
  # For each model compute average Z value 
  ZU.dich <- mean(log(ORs[,1]) / SEs[,1])
  ZA.dich <- mean(log(ORs[,2]) / SEs[,2])
  ZU.PO <- mean(log(ORs[,3]) / SEs[,3])
  ZA.PO <- mean(log(ORs[,4]) / SEs[,4])
  
  RSS.covadj.dich <- 100 - (100 * (ZU.dich/ZA.dich)^2)      # Reduction in sample size by adjustment (in dich. models)
  RSS.covadj.PO <- 100 - (100 * (ZU.PO/ZA.PO)^2)            # Reduction in sample size by adjustment (in ord. models)
  
  RSS.PO.unadj <- 100 - (100 * (ZU.dich/ZU.PO)^2)           # Reduction in sample size by exploiting ordinality (in unadjusted models)
  RSS.PO.adj <- 100 - (100 * (ZA.dich/ZA.PO)^2)             # Reduction in sample size by exploiting ordinality (in adjusted models)
  
  RSS.adjPO.unadjdich <- 100 - (100 * (ZU.dich / ZA.PO)^2)  # Reduction in sample size by adjusted PO as compared to unadjusted dich
  
  MCE <- sapply(1:4, function(x) {sqrt((colMeans(Rej)[x] * (1-colMeans(Rej)[x])) / nsim)}) # monte carlo error of the power estimate
  
  objects <- c(ORs=ORs, SEs=SEs, Pvals=Pvals, Rej=Rej, RejectionRate=colMeans(Rej),
                    RSS.covadj.PO = RSS.covadj.PO, RSS.covadj.dich = RSS.covadj.dich,
                    RSS.PO.adj = RSS.PO.adj, RSS.PO.unadj = RSS.PO.unadj,
                    RSS.adjPO.unadjdich = RSS.adjPO.unadjdich)
  
  return(list(matrix(data = c(colMeans(Rej)*100, colMeans(ORs), colMeans(SEs), NA, RSS.covadj.dich, RSS.PO.unadj, RSS.adjPO.unadjdich), 
                     ncol = 4, byrow = TRUE, dimnames = list(c("Rejection Rate (%)", "Mean OR", "Mean SE", "Red in SS (%)"),
                                                             c("Unadjusted dich", "Adjusted dich", "Unadjusted PO", "Adjusted PO"))), 
              objects,
              MCE = MCE))
}


### B.2 Perform simulation (Wald)

# 60 n
T1error.60n.weakcov <- SimFunWeakerCov(nsim = 10000, groupsize = 60, txeffectmort = log(1), txeffect = log(1))

txlog1.5.60n.weakcov <- SimFunWeakerCov(nsim = 10000, groupsize = 60, txeffectmort = log(1.5), txeffect = log(1.5))

txlog2.1.60n.weakcov <- SimFunWeakerCov(nsim = 10000, groupsize = 60, txeffectmort = log(2.1), txeffect = log(2.1))

txlog2.7.60n.weakcov <- SimFunWeakerCov(nsim = 10000, groupsize = 60, txeffectmort = log(2.7), txeffect = log(2.7))

# 100 n
T1error.100n.weakcov <- SimFunWeakerCov(nsim = 10000, groupsize = 100, txeffectmort = log(1), txeffect = log(1))

txlog1.5.100n.weakcov <- SimFunWeakerCov(nsim = 10000, groupsize = 100, txeffectmort = log(1.5), txeffect = log(1.5))

txlog2.1.100n.weakcov <- SimFunWeakerCov(nsim = 10000, groupsize = 100, txeffectmort = log(2.1), txeffect = log(2.1))

txlog2.7.100n.weakcov <- SimFunWeakerCov(nsim = 10000, groupsize = 100, txeffectmort = log(2.7), txeffect = log(2.7))


### B.3 Print results (Wald)

cat("Sample size:", 120, "(60/arm)\n")
cat("OR = 1 \n")
T1error.60n.weakcov[c(3, 1)]

cat("OR = 1.5 \n")
txlog1.5.60n.weakcov[c(3, 1)]

cat("OR = 2.1 \n")
txlog2.1.60n.weakcov[c(3, 1)]

cat("OR = 2.7 \n")
txlog2.7.60n.weakcov[c(3, 1)]

cat("Sample size:", 200, "(100/arm)\n")
cat("OR = 1 \n")
T1error.100n.weakcov[c(3, 1)]

cat("OR = 1.5 \n")
txlog1.5.100n.weakcov[c(3, 1)]

cat("OR = 2.1 \n")
txlog2.1.100n.weakcov[c(3, 1)]

cat("OR = 2.7 \n")
txlog2.7.100n.weakcov[c(3, 1)]


### B.1b Function for simulation (LRT)

SimFunWeakerCovLRT <- function(nsim, groupsize, txeffectmort, txeffect) {
  
  # input arguments:
  ## nsim = number of iterations per simulation 
  ## groupsize = the number of subjects per treatment arm 
  ## txeffectmort = the treatment effect for mortality (cut-off GBS DS 01 vs GBSDS 23456)
  ## txeffect = the treatment effect for remaining cut-offs
  
  # empty objects for storage
  ORs <- matrix(nrow = nsim, ncol = 4)
  SEs <-  matrix(nrow = nsim, ncol = 4)
  Pvals <- matrix(nrow = nsim, ncol = 4)
  Rej <-  matrix(nrow = nsim, ncol = 4)
  
  for (r in 1:nsim) {
    
    if (r %% 500 == 0) cat(paste0("Iteration: ", r, "\n")) # counter to track progression
    
    # Make X matrix for sample data set
    ind <- sample(1:nrow(MPdata), size = 2 * groupsize, replace = T) # random sampling with raplacement (obtain indices of observations that will be included in sample)
    
    SampleXmat <- MPdata[ind, c("age", "diarrhea")] # obtain X matrix for sampled observations
    SampleXmat$treat <- c(rep(0, groupsize), rep(1, groupsize)) # assign patients to a treatment group (two arms, balanced)
    
    # To get GBS DS for sampled patients:
    # Calculate linear predictors for each patient, for each cut off  
    Linpreds <- matrix(nrow  = 2*groupsize, ncol = 4) # empty object for storage
     
    Linpreds[, 1] <- MPcoefficients.agediar[1] + txeffectmort * SampleXmat$treat +          
                                          MPcoefficients.agediar["age"] * SampleXmat$age  +
                                          MPcoefficients.agediar["diarrhea"] * SampleXmat$diarrhea 
    
    Linpreds[, 2:4] <- sapply(2:4, function(x) {MPcoefficients.agediar[x] + txeffect * SampleXmat$treat +         
                                          MPcoefficients.agediar["age"] * SampleXmat$age  +
                                          MPcoefficients.agediar["diarrhea"] * SampleXmat$diarrhea})
    
    probs <- matrix(nrow = nrow(SampleXmat), ncol = 5) # empty object for storage
    GBSDS <- numeric(nrow(SampleXmat)) # empty object for storage
    
    for (i in 1:nrow(SampleXmat)) { # For each subject, calculate probability for each category
      probs[i, 1] <- 1 - plogis(Linpreds[i, 1])                           # P(GBSDS = 0/1)
      probs[i, 2] <- plogis(Linpreds[i, 1]) - plogis(Linpreds[i, 2])      # P(GBSDS = 2)
      probs[i, 3] <- plogis(Linpreds[i, 2]) - plogis(Linpreds[i, 3])      # P(GBSDS = 3)
      probs[i, 4] <- plogis(Linpreds[i, 3]) - plogis(Linpreds[i, 4])      # P(GBSDS = 4)
      probs[i, 5] <- plogis(Linpreds[i, 4])                               # P(GBSDS = 5/6)
      
      GBSDS[i] <-  which(rmultinom(1, 1, (probs[i, ]))==1)                # sample a GBS DS score from multinomial distribution (using the probabilities)
                                                                          # scale is 1 through 5 (1 standing for 0 and 1; 5 standing for 5 and 6)
    }
    
    # compile data set
    comp <- cbind(SampleXmat, GBSDS)
    dd <<- datadist(comp) # the <<- operator  assigns the variable to the global environment
    options(datadist = "dd")
    
    # Fit dichotomous models (binary log. reg.)
    moddich.unadj <- lrm(as.numeric(GBSDS > 3) ~ treat, data = comp)                  # unadjusted
    moddich.adj <- lrm(as.numeric(GBSDS > 3) ~ treat + age + diarrhea , data = comp) # adjusted
    
    # Fit ordinal models (proportional odds log. reg.)
    modPO.unadj <- lrm(factor(GBSDS) ~ treat, data = comp)                            # unadjusted 
    modPO.adj <- lrm(factor(GBSDS) ~ treat + age + diarrhea , data = comp)           # adjusted
    
    # Fit nested models in order to obtain LRT p value
    ## dich adj
    moddich.adj.nested <- lrm(as.numeric(GBSDS > 3) ~ age + diarrhea , data = comp) 
    ll.dich.nested <- logLik(moddich.adj.nested)
    ll.dich.complex <- logLik(moddich.adj)   
    teststat.dich <- -2*(as.numeric(ll.dich.nested)-as.numeric(ll.dich.complex))
    pval.dich <- pchisq(teststat.dich, df = 1, lower.tail = F)
    
    ## po adj
    modPO.adj.nested <- lrm(factor(GBSDS) ~  age + diarrhea , data = comp)  
    ll.PO.nested <- logLik(modPO.adj.nested)
    ll.PO.complex <- logLik(modPO.adj)
    teststat.PO <- -2*(as.numeric(ll.PO.nested)-as.numeric(ll.PO.complex))
    pval.PO <- pchisq(teststat.PO, df = 1, lower.tail = F)

    # For each model: extract OR, SE, p val of tx effect and indication of hypothesis rejection (0/1)    
    ORs[r, ] <- exp(c(summary(moddich.unadj)['treat', 'Effect'],
                  summary(moddich.adj)['treat', 'Effect'],
                  summary(modPO.unadj)['treat', 'Effect'],
                  summary(modPO.adj)['treat', 'Effect']))
    
    SEs[r, ] <- c(summary(moddich.unadj)['treat', 'S.E.'], 
                  summary(moddich.adj)['treat', 'S.E.'],
                  summary(modPO.unadj)['treat', 'S.E.'],
                  summary(modPO.adj)['treat', 'S.E.'])
    
    Pvals[r, ] <- c(moddich.unadj$stats[5],
                    pval.dich,
                    modPO.unadj$stats[5], 
                    pval.PO)
    
    Rej[r, ] <- as.numeric(Pvals[r, ] < 0.05)
  }
  
  # For each model compute average Z value 
  ZU.dich <- mean(log(ORs[,1]) / SEs[,1])
  ZA.dich <- mean(log(ORs[,2]) / SEs[,2])
  ZU.PO <- mean(log(ORs[,3]) / SEs[,3])
  ZA.PO <- mean(log(ORs[,4]) / SEs[,4])
  
  RSS.covadj.dich <- 100 - (100 * (ZU.dich/ZA.dich)^2)      # Reduction in sample size by adjustment (in dich. models)
  RSS.covadj.PO <- 100 - (100 * (ZU.PO/ZA.PO)^2)            # Reduction in sample size by adjustment (in ord. models)
  
  RSS.PO.unadj <- 100 - (100 * (ZU.dich/ZU.PO)^2)           # Reduction in sample size by exploiting ordinality (in unadjusted models)
  RSS.PO.adj <- 100 - (100 * (ZA.dich/ZA.PO)^2)             # Reduction in sample size by exploiting ordinality (in adjusted models)
  
  RSS.adjPO.unadjdich <- 100 - (100 * (ZU.dich / ZA.PO)^2)  # Reduction in sample size by adjusted PO as compared to unadjusted dich
  
  MCE <- sapply(1:4, function(x) {sqrt((colMeans(Rej)[x] * (1-colMeans(Rej)[x])) / nsim)}) # monte carlo error of the power estimate
  
  objects <- c(ORs=ORs, SEs=SEs, Pvals=Pvals, Rej=Rej, RejectionRate=colMeans(Rej),
                    RSS.covadj.PO = RSS.covadj.PO, RSS.covadj.dich = RSS.covadj.dich,
                    RSS.PO.adj = RSS.PO.adj, RSS.PO.unadj = RSS.PO.unadj,
                    RSS.adjPO.unadjdich = RSS.adjPO.unadjdich)
  
  return(list(matrix(data = c(colMeans(Rej)*100, colMeans(ORs), colMeans(SEs), NA, RSS.covadj.dich, RSS.PO.unadj, RSS.adjPO.unadjdich), 
                     ncol = 4, byrow = TRUE, dimnames = list(c("Rejection Rate (%)", "Mean OR", "Mean SE", "Red in SS (%)"),
                                                             c("Unadjusted dich", "Adjusted dich", "Unadjusted PO", "Adjusted PO"))), 
              objects,
              MCE = MCE))
}


### B.2b Perform simulation (LRT)

# 60 n
T1error.60n.weakcov.LRT <- SimFunWeakerCovLRT(nsim = 10000, groupsize = 60, txeffectmort = log(1), txeffect = log(1))

txlog1.5.60n.weakcov.LRT <- SimFunWeakerCovLRT(nsim = 10000, groupsize = 60, txeffectmort = log(1.5), txeffect = log(1.5))

txlog2.1.60n.weakcov.LRT <- SimFunWeakerCovLRT(nsim = 10000, groupsize = 60, txeffectmort = log(2.1), txeffect = log(2.1))

txlog2.7.60n.weakcov.LRT <- SimFunWeakerCovLRT(nsim = 10000, groupsize = 60, txeffectmort = log(2.7), txeffect = log(2.7))

# 100 n
T1error.100n.weakcov.LRT <- SimFunWeakerCovLRT(nsim = 10000, groupsize = 100, txeffectmort = log(1), txeffect = log(1))

txlog1.5.100n.weakcov.LRT <- SimFunWeakerCovLRT(nsim = 10000, groupsize = 100, txeffectmort = log(1.5), txeffect = log(1.5))

txlog2.1.100n.weakcov.LRT <- SimFunWeakerCovLRT(nsim = 10000, groupsize = 100, txeffectmort = log(2.1), txeffect = log(2.1))

txlog2.7.100n.weakcov.LRT <- SimFunWeakerCovLRT(nsim = 10000, groupsize = 100, txeffectmort = log(2.7), txeffect = log(2.7))


### B.3b Print results (LRT)

cat("Sample size:", 120, "(60/arm) (LRT)\n")
cat("OR = 1 \n")
T1error.60n.weakcov.LRT[c(3, 1)]

cat("OR = 1.5 \n")
txlog1.5.60n.weakcov.LRT[c(3, 1)]

cat("OR = 2.1 \n")
txlog2.1.60n.weakcov.LRT[c(3, 1)]

cat("OR = 2.7 \n")
txlog2.7.60n.weakcov.LRT[c(3, 1)]

cat("Sample size:", 200, "(100/arm) (LRT)\n")
cat("OR = 1 \n")
T1error.100n.weakcov.LRT[c(3, 1)]

cat("OR = 1.5 \n")
txlog1.5.100n.weakcov.LRT[c(3, 1)]

cat("OR = 2.1 \n")
txlog2.1.100n.weakcov.LRT[c(3, 1)]

cat("OR = 2.7 \n")
txlog2.7.100n.weakcov.LRT[c(3, 1)]


## Condition C: Betas for covariate strength as observed, but quantitative non-proportionality introduced for treatment OR 

# Grid search for combination of odds ratios (QUALITATIVE / SEVERE NPO)

a <- seq(1.05, 3.4, 0.05) # values for treatment effect at all but mortality cut offs
b <- seq(0.80, 0.95, 0.05) # values for treatment effect mortality 
mygrid <- expand.grid(ORs = a, ORmort = b)
cOR <- numeric(length=nrow(mygrid))

for (i in 1:nrow(mygrid)) {
  a <- mygrid[i, "ORs"]
  b <- mygrid[i, "ORmort"]
  sim <- SimFun(txeffect = log(a), txeffectmort = log(b), nsim = 10, groupsize = 10000)
  cOR[i] <- data.frame(sim[1])[2,4]
  print(c(i, cOR[i], mygrid[i, ]))
}

(txeffects1 <- mygrid[which.min(abs(cOR - 1)), ])      
cOR[which.min(abs(cOR - 1))]          

(txeffects1.5 <- mygrid[which.min(abs(cOR - 1.5)), ])   
cOR[which.min(abs(cOR - 1.5))]        

(txeffects2.1 <- mygrid[which.min(abs(cOR - 2.1)), ])  
cOR[which.min(abs(cOR - 2.1))]       

(txeffects2.7 <- mygrid[which.min(abs(cOR - 2.7)), ])   
cOR[which.min(abs(cOR - 2.7))]       


### C.1 Function for simulation 
Same as in proportionality, but now we specify different betas for mortality and the other cut-offs

### C.2 Perform simulation

# 60 n
NPOT1error.60n <- SimFun(nsim = 10000, groupsize = 60, txeffectmort = as.numeric(log(txeffects1$ORmort)), txeffect = as.numeric(log(txeffects1$ORs)))

NPOtxlog1.5.60n <- SimFun(nsim = 10000, groupsize = 60, txeffectmort = as.numeric(log(txeffects1.5$ORmort)), txeffect = as.numeric(log(txeffects1.5$ORs)))

NPOtxlog2.1.60n <- SimFun(nsim = 10000, groupsize = 60, txeffectmort = as.numeric(log(txeffects2.1$ORmort)), txeffect = as.numeric(log(txeffects2.1$ORs)))

NPOtxlog2.7.60n <- SimFun(nsim = 10000, groupsize = 60, txeffectmort = as.numeric(log(txeffects2.7$ORmort)), txeffect = as.numeric(log(txeffects2.7$ORs)))

# 100 n
NPOT1error.100n <- SimFun(nsim = 10000, groupsize = 100, txeffectmort = as.numeric(log(txeffects1$ORmort)), txeffect = as.numeric(log(txeffects1$ORs)))

NPOtxlog1.5.100n <- SimFun(nsim = 10000, groupsize = 100, txeffectmort = as.numeric(log(txeffects1.5$ORmort)), txeffect = as.numeric(log(txeffects1.5$ORs)))

NPOtxlog2.1.100n <- SimFun(nsim = 10000, groupsize = 100, txeffectmort = as.numeric(log(txeffects2.1$ORmort)), txeffect = as.numeric(log(txeffects2.1$ORs)))

NPOtxlog2.7.100n <- SimFun(nsim = 10000, groupsize = 100, txeffectmort = as.numeric(log(txeffects2.7$ORmort)), txeffect = as.numeric(log(txeffects2.7$ORs)))


### C.3 Print results

cat("Sample size:", 100, "(60/arm)\n")
cat("OR = 1 \n")
NPOT1error.60n[c(3, 1)]

cat("OR = 1.5 \n")
NPOtxlog1.5.60n[c(3, 1)]

cat("OR = 2.1 \n")
NPOtxlog2.1.60n[c(3, 1)]

cat("OR = 2.7 \n")
NPOtxlog2.7.60n[c(3, 1)]

cat("Sample size:", 200, "(100/arm)\n")
cat("OR = 1 \n")
NPOT1error.100n[c(3, 1)]

cat("OR = 1.5 \n")
NPOtxlog1.5.100n[c(3, 1)]

cat("OR = 2.1 \n")
NPOtxlog2.1.100n[c(3, 1)]

cat("OR = 2.7 \n")
NPOtxlog2.7.100n[c(3, 1)]


### C.2b Perform simulation (LRT)

# 60 n
NPOT1error.60n.LRT <- SimFunLRT(nsim = 10000, groupsize = 60, txeffectmort = as.numeric(log(txeffects1$ORmort)), txeffect = as.numeric(log(txeffects1$ORs)))

NPOtxlog1.5.60n.LRT <- SimFunLRT(nsim = 10000, groupsize = 60, txeffectmort = as.numeric(log(txeffects1.5$ORmort)), txeffect = as.numeric(log(txeffects1.5$ORs)))

NPOtxlog2.1.60n.LRT <- SimFunLRT(nsim = 10000, groupsize = 60, txeffectmort = as.numeric(log(txeffects2.1$ORmort)), txeffect = as.numeric(log(txeffects2.1$ORs)))

NPOtxlog2.7.60n.LRT <- SimFunLRT(nsim = 10000, groupsize = 60, txeffectmort = as.numeric(log(txeffects2.7$ORmort)), txeffect = as.numeric(log(txeffects2.7$ORs)))

# 100 n
NPOT1error.100n.LRT <- SimFunLRT(nsim = 10000, groupsize = 100, txeffectmort = as.numeric(log(txeffects1$ORmort)), txeffect = as.numeric(log(txeffects1$ORs)))

NPOtxlog1.5.100n.LRT <- SimFunLRT(nsim = 10000, groupsize = 100, txeffectmort = as.numeric(log(txeffects1.5$ORmort)), txeffect = as.numeric(log(txeffects1.5$ORs)))

NPOtxlog2.1.100n.LRT <- SimFunLRT(nsim = 10000, groupsize = 100, txeffectmort = as.numeric(log(txeffects2.1$ORmort)), txeffect = as.numeric(log(txeffects2.1$ORs)))

NPOtxlog2.7.100n.LRT <- SimFunLRT(nsim = 10000, groupsize = 100, txeffectmort = as.numeric(log(txeffects2.7$ORmort)), txeffect = as.numeric(log(txeffects2.7$ORs)))


### C.3b Print results (LRT)

cat("Sample size:", 100, "(60/arm) (LRT) \n")
cat("OR = 1 \n")
NPOT1error.60n.LRT[c(3, 1)]

cat("OR = 1.5 \n")
NPOtxlog1.5.60n.LRT[c(3, 1)]

cat("OR = 2.1 \n")
NPOtxlog2.1.60n.LRT[c(3, 1)]

cat("OR = 2.7 \n")
NPOtxlog2.7.60n.LRT[c(3, 1)]

cat("Sample size:", 200, "(100/arm) (LRT) \n")
cat("OR = 1 \n")
NPOT1error.100n.LRT[c(3, 1)]

cat("OR = 1.5 \n")
NPOtxlog1.5.100n.LRT[c(3, 1)]

cat("OR = 2.1 \n")
NPOtxlog2.1.100n.LRT[c(3, 1)]

cat("OR = 2.7 \n")
NPOtxlog2.7.100n.LRT[c(3, 1)]


## Condition D: Betas for covariate strength as observed, but quantitative non-proportionality introduced for treatment OR 

# Grid search for combination of odds ratios (QUANTITATIVE / MILD NPO)

a <- seq(1.05, 3.5, 0.05) # values for treatment effect for other cut offs
b.2 <- rep(1, 50) # values for treatment effect mortality 
mygrid.2 <- data.frame(ORs = a, ORmort = b.2)
cOR.2 <- numeric(length=nrow(mygrid.2))

for (i in 1:nrow(mygrid.2)) {
  a <- mygrid.2[i, 1]
  b <- mygrid.2[i, 2]
  sim <- SimFun(txeffect = log(a), txeffectmort = log(b), nsim = 10, groupsize = 10000)
  cOR.2[i] <- data.frame(sim[1])[2,4]
  print(c(i, cOR.2[i], mygrid.2[i, ]))
}

(txeffects1.mnpo <- mygrid.2[which.min(abs(cOR.2 - 1)), ]  )   
cOR.2[which.min(abs(cOR.2 - 1))]          

(txeffects1.5.mnpo <- mygrid.2[which.min(abs(cOR.2 - 1.5)), ]   )
cOR.2[which.min(abs(cOR.2 - 1.5))]        

(txeffects2.1.mnpo <- mygrid.2[which.min(abs(cOR.2 - 2.1)), ]   )
cOR.2[which.min(abs(cOR.2 - 2.1))]        

(txeffects2.7.mnpo <- mygrid.2[which.min(abs(cOR.2 - 2.7)), ]   )
cOR.2[which.min(abs(cOR.2 - 2.7))]        


### D.1 Function for simulation 
Same as in proportionality, but now we specify different betas for mortality and the other cut-offs

### D.2 Perform simulation


# 60 n
mNPOT1error.60n <- SimFun(nsim = 10000, groupsize = 60, txeffectmort = log(txeffects1.mnpo$ORmort), txeffect = log(txeffects1.mnpo$ORs))

mNPOtxlog1.5.60n <- SimFun(nsim = 10000, groupsize = 60, txeffectmort = log(txeffects1.5.mnpo$ORmort), txeffect = log(txeffects1.5.mnpo$ORs))

mNPOtxlog2.1.60n <- SimFun(nsim = 10000, groupsize = 60, txeffectmort = log(txeffects2.1.mnpo$ORmort), txeffect = log(txeffects2.1.mnpo$ORs))

mNPOtxlog2.7.60n <- SimFun(nsim = 10000, groupsize = 60, txeffectmort = log(txeffects2.7.mnpo$ORmort), txeffect = log(txeffects2.7.mnpo$ORs))

# 100 n
mNPOT1error.100n <- SimFun(nsim = 10000, groupsize = 100, txeffectmort = log(txeffects1.mnpo$ORmort), txeffect = log(txeffects1.mnpo$ORs))

mNPOtxlog1.5.100n <- SimFun(nsim = 10000, groupsize = 100, txeffectmort = log(txeffects1.5.mnpo$ORmort), txeffect = log(txeffects1.5.mnpo$ORs))

mNPOtxlog2.1.100n <- SimFun(nsim = 10000, groupsize = 100, txeffectmort = log(txeffects2.1.mnpo$ORmort), txeffect = log(txeffects2.1.mnpo$ORs))

mNPOtxlog2.7.100n <- SimFun(nsim = 10000, groupsize = 100, txeffectmort = log(txeffects2.7.mnpo$ORmort), txeffect = log(txeffects2.7.mnpo$ORs))


### D.3 Print results

cat("Sample size:", 100, "(60/arm)\n")
cat("OR = 1 \n")
mNPOT1error.60n[c(3, 1)]

cat("OR = 1.5 \n")
mNPOtxlog1.5.60n[c(3, 1)]

cat("OR = 2.1 \n")
mNPOtxlog2.1.60n[c(3, 1)]

cat("OR = 2.7 \n")
mNPOtxlog2.7.60n[c(3, 1)]

cat("Sample size:", 200, "(100/arm)\n")
cat("OR = 1 \n")
mNPOT1error.100n[c(3, 1)]

cat("OR = 1.5 \n")
mNPOtxlog1.5.100n[c(3, 1)]

cat("OR = 2.1 \n")
mNPOtxlog2.1.100n[c(3, 1)]

cat("OR = 2.7 \n")
mNPOtxlog2.7.100n[c(3, 1)]


### D.1b Function for simulation  (LRT)
Same as in proportionality, but now we specify different betas for mortality and the other cut-offs

### D.2b Perform simulation (LRT)

# 60 n
mNPOT1error.60n.LRT <- SimFunLRT(nsim = 10000, groupsize = 60, txeffectmort = log(txeffects1.mnpo$ORmort), txeffect = log(txeffects1.mnpo$ORs))

mNPOtxlog1.5.60n.LRT <- SimFunLRT(nsim = 10000, groupsize = 60, txeffectmort = log(txeffects1.5.mnpo$ORmort), txeffect = log(txeffects1.5.mnpo$ORs))

mNPOtxlog2.1.60n.LRT <- SimFunLRT(nsim = 10000, groupsize = 60, txeffectmort = log(txeffects2.1.mnpo$ORmort), txeffect = log(txeffects2.1.mnpo$ORs))

mNPOtxlog2.7.60n.LRT <- SimFunLRT(nsim = 10000, groupsize = 60, txeffectmort = log(txeffects2.7.mnpo$ORmort), txeffect = log(txeffects2.7.mnpo$ORs))

# 100 n
mNPOT1error.100n.LRT <- SimFunLRT(nsim = 10000, groupsize = 100, txeffectmort = log(txeffects1.mnpo$ORmort), txeffect = log(txeffects1.mnpo$ORs))

mNPOtxlog1.5.100n.LRT <- SimFunLRT(nsim = 10000, groupsize = 100, txeffectmort = log(txeffects1.5.mnpo$ORmort), txeffect = log(txeffects1.5.mnpo$ORs))

mNPOtxlog2.1.100n.LRT <- SimFunLRT(nsim = 10000, groupsize = 100, txeffectmort = log(txeffects2.1.mnpo$ORmort), txeffect = log(txeffects2.1.mnpo$ORs))

mNPOtxlog2.7.100n.LRT <- SimFunLRT(nsim = 10000, groupsize = 100, txeffectmort = log(txeffects2.7.mnpo$ORmort), txeffect = log(txeffects2.7.mnpo$ORs))


### D.3b Print results (LRT)

cat("Sample size:", 100, "(60/arm) (LRT) \n")
cat("OR = 1 \n")
mNPOT1error.60n.LRT[c(3, 1)]

cat("OR = 1.5 \n")
mNPOtxlog1.5.60n.LRT[c(3, 1)]

cat("OR = 2.1 \n")
mNPOtxlog2.1.60n.LRT[c(3, 1)]

cat("OR = 2.7 \n")
mNPOtxlog2.7.60n.LRT[c(3, 1)]

cat("Sample size:", 200, "(100/arm) (LRT) \n")
cat("OR = 1 \n")
mNPOT1error.100n.LRT[c(3, 1)]

cat("OR = 1.5 \n")
mNPOtxlog1.5.100n.LRT[c(3, 1)]

cat("OR = 2.1 \n")
mNPOtxlog2.1.100n.LRT[c(3, 1)]

cat("OR = 2.7 \n")
mNPOtxlog2.7.100n.LRT[c(3, 1)]

