---
title: "R Notebook"
output: html_notebook
---

```{r}
library(foreign)
library(dplyr)
library(rms)
library(lme4)
library(reshape2)
library(tidyr)
library(cowplot)
library(ggplot2)
library(MASS)
library(Hmisc)
library(forestplot)
library(brms)
library(rmsb)
library(rstanarm)
library(rstan)
```

## Read in data
```{r}
MP.original <- read.spss("MP trial database complete cases n 221.sav", use.value.label=TRUE, to.data.frame=TRUE)

MP.original$treat <- droplevels(MP.original$treat, exclude = if(anyNA(levels(MP.original$treat))) NULL else NA) # empty levels dropped
MP.original$treat <- recode(MP.original$treat, IVIg = "IVIg+MP") # correction for labels
MP.original$treat <- recode(MP.original$treat, plasmaferese = "IVIg")

ReverseGBS <- function(x) { # reverse the scale, such that higher GBSDS = better outcome
  Rec <- x
  Rec[x == 0] <- 5 
  Rec[x == 1] <- 5
  Rec[x == 2] <- 4 
  Rec[x == 4] <- 2
  Rec[x == 5] <- 1 
  Rec[x == 6] <- 1
  return(Rec)
}

MPdata <- data.frame(id = MP.original$trial,
                     treat = MP.original$treat,
                     ReverseGBS(MP.original[,c(25, 28:38)]), # GBSDS scores at each time point without baseline 
                     mrcss = MP.original$mrcss,
                     weakonse = MP.original$weakonse,
                     age = MP.original$age, 
                     diarrhea = MP.original$diarrhea)

visits.labels <- c("day 6", "day 13", "week 3", "week 4", "week 5", "week 6", "week 8", "week 10", "week 14", "week 18", "week 22", "week 26")
```

## Data in long format
```{r MP data to long format MP}
# MP data set to long format
# a row for each observation within a patient
mplong <- gather(MPdata, key = time, value = GBSDS, f4:f17) %>% arrange(id) 
dd <- datadist(mplong); options(datadist='dd')

tmp <-factor(mplong$time, levels = mplong$time[1:12], labels = c(6/7, 13/7, 3, 4, 5, 6, 8, 10, 14, 18, 22, 26))
mplong$weeks <- as.numeric(levels(tmp)[tmp])

mplong$timefac <- as.factor(mplong$time)

mplong$sqrtweeks <- sqrt(mplong$weeks)
```

## 1. Fit a PO model for each visit in the follow-up
```{r}
# UNADJUSTED PO

dd <- datadist(MPdata); options(datadist = "dd") 


my.po.unadj <- matrix(nrow = 12, ncol = 4, dimnames = list(visits.labels, c("Effect", "S.E.", "Lower 0.95", "Upper 0.95")))

for (i in 1:12) {
  mod <- lrm(MPdata[, 2+i] ~ treat, data = MPdata)
  my.po.unadj[i, ] <- exp(summary(mod)[1, 4:7])
  my.po.unadj[i, 2] <- summary(mod)[1, 5]
}

round(my.po.unadj, 2)

# ADJUSTED PO
my.po.adj <- matrix(nrow = 12, ncol = 4, dimnames = list(visits.labels, c("Effect", "S.E.", "Lower 0.95", "Upper 0.95")))

for (i in 1:12) {
  mod <- lrm(MPdata[, 2+i] ~ treat + age + diarrhea, data = MPdata)
  my.po.adj[i, ] <- exp(summary(mod)["treat - IVIg+MP:IVIg", 4:7])
  my.po.adj[i, 2] <- summary(mod)["treat - IVIg+MP:IVIg", 5]
}

round(my.po.adj, 2)
```

## 2. Fit a binary logistic regression for each possible dichotomy at each time point in the follow-up
```{r}
BinLogRegUnadj <- function(sc, data, cutoffs = 1:4) {
  # For a given time point, the function loops over each cumulative cut-off, performs
  # binary logistic regression and returns an OR for the unadjusted treatment effect, 
  # with its corresponding confidence interval and standard error of the beta coefficient.
  mat <- matrix(nrow = 4, ncol = 4)
  colnames(mat) <- c("OR", "sd", "Lower 0.95", "Upper 0.95")
  for (i in cutoffs) {
    dd <- datadist(data); options(datadist = "dd")
    fit <- lrm(as.numeric(sc>i) ~ treat, data = data)
    mat[i, 1] <- exp(summary(fit)["treat - IVIg+MP:IVIg", "Effect"]) #OR
    mat[i, 2] <- summary(fit)["treat - IVIg+MP:IVIg", "S.E."] #sd
    mat[i, c(3,4)] <- exp(summary(fit)["treat - IVIg+MP:IVIg", c("Lower 0.95", "Upper 0.95")]) #lower and upper bound
  }
  return(round(mat, 2))
}

BinLogRegAdj <- function(sc, data, cutoffs = 1:4) {
  # For a given time point, the function loops over each cumulative cut-off, performs
  # binary logistic regression and returns an OR for the adjusted treatment effect, 
  # with its corresponding confidence interval and standard error of the beta coefficient.
  mat1 <- matrix(nrow = 4, ncol = 4)
  colnames(mat1) <- c("OR", "sd", "Lower 0.95", "Upper 0.95")
  for (i in cutoffs) {
    dd <- datadist(data); options(datadist = "dd")
    fit <- lrm(as.numeric(sc>i) ~ treat + mrcss + weakonse, data = data)
    mat1[i, 1] <- exp(summary(fit)["treat - IVIg+MP:IVIg", "Effect"]) #OR
    mat1[i, 2] <- summary(fit)["treat - IVIg+MP:IVIg", "S.E."] #sd
    mat1[i, c(3,4)] <- exp(summary(fit)["treat - IVIg+MP:IVIg", c("Lower 0.95", "Upper 0.95")]) #lower and upper bound
  }
  return(round(mat1, 2))
}

# Apply functions to each time point
ORs.unadj <- ORs.adj <- vector(mode = "list", length = 12)
names(ORs.adj) <- names(ORs.unadj) <- visits.labels
for (i in 1:12) {
  ORs.adj[[i]] <- BinLogRegAdj(sc = MPdata[,2+i], data = MPdata)
  ORs.unadj[[i]] <- BinLogRegUnadj(sc = MPdata[,2+i], data = MPdata)
}

ORs.unadj
ORs.adj
```

## Plot of treatment effect over time
```{r}
# Plots of treatment effect over time

logORslabels <- c("log(0.1)", "log(0.2)", "log(0.5)", "log(1)", "log(2)", "log(5)", "log(10)")

visits <- c("day 6", "day 13", "week 3", "week 4", "week 5", "week 6", "week 8", "week 10", "week 14", "week 18", "week 22", "week 26")
visits <- factor(visits, levels = visits)

cORunadj <- data.frame(logOR = log(c(my.po.unadj[, 1], 
                            sapply(ORs.unadj, function(x) {x[1, 1]}), sapply(ORs.unadj, function(x) {x[2, 1]}), 
                            sapply(ORs.unadj, function(x) {x[3, 1]}), sapply(ORs.unadj, function(x) {x[4, 1]}))),
                      Dichotomy = c(rep("common OR (PO model)", 12), rep("01 vs 23456", 12), rep("012 vs 3456", 12), rep("0123 vs 456", 12), rep("01234 vs 56", 12)),
                      visits = rep(visits, 5),
                      time = rep(c(round(6/7, 1), round(13/7, 1), 3, 4, 5, 6, 8, 10, 14, 18, 22, 26), 5),
                      low = log(c(my.po.unadj[, 3], 
                            sapply(ORs.unadj, function(x) {x[1, 3]}), sapply(ORs.unadj, function(x) {x[2, 3]}), 
                            sapply(ORs.unadj, function(x) {x[3, 3]}), sapply(ORs.unadj, function(x) {x[4, 3]}))),
                      high = log(c(my.po.unadj[, 4], 
                            sapply(ORs.unadj, function(x) {x[1, 4]}), sapply(ORs.unadj, function(x) {x[2, 4]}), 
                            sapply(ORs.unadj, function(x) {x[3, 4]}), sapply(ORs.unadj, function(x) {x[4, 4]}))))

p1 <- ggplot(cORunadj, aes(x=time, y=logOR, group = Dichotomy, color = Dichotomy)) + 
  theme_cowplot(font_family = "serif",
                font_size = 14) +
  geom_hline(yintercept=0, color = "grey") +
  scale_x_discrete(limits = (c(round(6/7, 1), round(13/7, 1), 3, 4, 5, 6, 8, 10, 14, 18, 22, 26)),
                   labels = visits) +
  scale_y_discrete(limits = log(c(0.1, 0.2, 0.5, 1, 2, 5, 10)),
                   labels = logORslabels)  +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 0.85),
        legend.position = "none") +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.05, linetype = "dotted") +  
  labs(title="Treatment effect over time (unadjusted)", x="Time", y = "log(OR)") +
  scale_fill_manual(values=c('#999999','#E69F00'))  

cORadj <- data.frame(logOR = log(c(my.po.adj[, 1], 
                            sapply(ORs.adj, function(x) {x[1, 1]}), sapply(ORs.adj, function(x) {x[2, 1]}), 
                            sapply(ORs.adj, function(x) {x[3, 1]}), sapply(ORs.adj, function(x) {x[4, 1]}))),
                     Dichotomy = c(rep("commmon OR (PO model)", 12), rep("01 vs 23456", 12), rep("012 vs 3456", 12), rep("0123 vs 456", 12), rep("01234 vs 56", 12)),
                     visits = rep(visits, 5),
                     time = rep(c(round(6/7, 1), round(13/7, 1), 3, 4, 5, 6, 8, 10, 14, 18, 22, 26), 5),
                     low = log(c(my.po.adj[, 3], 
                            sapply(ORs.adj, function(x) {x[1, 3]}), sapply(ORs.adj, function(x) {x[2, 3]}), 
                            sapply(ORs.adj, function(x) {x[3, 3]}), sapply(ORs.adj, function(x) {x[4, 3]}))),
                     high = log(c(my.po.adj[, 4], 
                            sapply(ORs.adj, function(x) {x[1, 4]}), sapply(ORs.adj, function(x) {x[2, 4]}), 
                            sapply(ORs.adj, function(x) {x[3, 4]}), sapply(ORs.adj, function(x) {x[4, 4]}))))

p2 <- ggplot(cORadj, aes(x=time, y=logOR, group = Dichotomy, color = Dichotomy)) + 
  theme_cowplot(font_family = "serif",
                font_size = 14) +
  geom_hline(yintercept = 0, color = "grey") +
  scale_x_discrete(limits = (c(round(6/7, 1), round(13/7, 1), 3, 4, 5, 6, 8, 10, 14, 18, 22, 26)),
                   labels = visits) +
  scale_y_discrete(limits = log(c(0.1, 0.2, 0.5, 1, 2, 5, 10)),
                     labels = logORslabels) +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 0.85),
        legend.position = "none") +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.05, linetype = "dotted") + 
  labs(title="Treatment effect over time (adjusted)", x="Time", y = "log(OR)")+
  scale_fill_manual(values=c('#999999','#E69F00'))  
  
prow <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1))
legend <- get_legend(p1 + 
                       guides(color = guide_legend(nrow = 1, title = "")) +
                       theme(legend.position = c(0.25, 0.8))
)

ptot <- plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))
ptot
```
```{r}
idplot <- xyplot(GBSDS ~ sqrtweeks | id, data = mplong, type = "p", pch = 20, strip = F, xlab = "Time points (sqrt(weeks))", ylab = "GBS DS (reversed)", scales = list(list(cex = 0.6)), par.settings = my.theme, layout = c(13, 17))
```


## Longitudinal binary logistic models (mixed effects model)
### Optimal number of gaussian quadrature points (unadjusted analysis)
```{r}
# UNADJUSTED

quadraturepoints <- c(0, 1, 2, 5, 10, 25, 50)
modAGQres1 <- modAGQres2 <- modAGQres3 <- modAGQres4 <- matrix(nrow = length(quadraturepoints), ncol = 10, dimnames = list(rep(NA, length(quadraturepoints)), c("nAGQ", "(Intercept)", "treatIVIg+MP", "weeks", "treatIVIg+MP:weeks", "seInt", "seTreat", "seWeeks", "seTreat*Weeks", "varre")))

# GBSDS 01 vs 23456 (40 is sufficient)
# modAGQres1 <- matrix(nrow = length(quadraturepoints), ncol = 10)
# colnames(modAGQres1) <- c("nAGQ", "(Intercept)", "treatIVIg+MP", "weeks", "treatIVIg+MP:weeks", "seInt", "seTreat", "seWeeks", "seTreat*Weeks", "varre")

for (i in 1:length(quadraturepoints)) {
  modAGQ <- glmer((GBSDS>1) ~ treat * weeks + (1 | id),  family = "binomial", nAGQ = quadraturepoints[i], data = mplong, control = glmerControl(optimizer="bobyqa"))
  modAGQres1[i, ] <- c(quadraturepoints[i], fixef(modAGQ), sqrt(diag(vcov(modAGQ))), sqrt(VarCorr(modAGQ)[[1]][1]))
}

melt1 <- cbind(melt(modAGQres1[,-1], id.vars = "nAGQ"), p = rep(quadraturepoints, 9))
ggplot(melt1, aes(p, value, color = Var2)) + geom_line(aes(colour = Var2)) 

# GBSDS 012 vs 3456 (20 is sufficient)
modAGQres2 <- matrix(nrow = length(quadraturepoints), ncol = 10)
colnames(modAGQres2) <- c("nAGQ", "(Intercept)", "treatIVIg+MP", "weeks", "treatIVIg+MP:weeks", "seInt", "seTreat", "seWeeks", "seTreat*Weeks", "varre")

for (i in 1:length(quadraturepoints)) {
  modAGQ <- glmer((GBSDS>2) ~ treat * weeks + (1 | id),  family = "binomial", nAGQ = quadraturepoints[i], data = mplong, control = glmerControl(optimizer="bobyqa"))
  modAGQres2[i, ] <- c(quadraturepoints[i], fixef(modAGQ), sqrt(diag(vcov(modAGQ))), sqrt(VarCorr(modAGQ)[[1]][1]))
}
melt2 <- cbind(melt(modAGQres2[,-1], id.vars = "nAGQ"), p = rep(quadraturepoints, 9))
ggplot(melt2, aes(p, value, color = Var2)) + geom_line(aes(colour = Var2)) 

# GBSDS 0123 vs 456 (20 is sufficient)
modAGQres3 <- matrix(nrow = length(quadraturepoints), ncol = 10)
colnames(modAGQres3) <- c("nAGQ", "(Intercept)", "treatIVIg+MP", "weeks", "treatIVIg+MP:weeks", "seInt", "seTreat", "seWeeks", "seTreat*Weeks", "varre")

for (i in 1:length(quadraturepoints)) {
  modAGQ <- glmer((GBSDS>3) ~ treat * weeks + (1 | id),  family = "binomial", nAGQ = quadraturepoints[i], data = mplong, control = glmerControl(optimizer="bobyqa"))
  modAGQres3[i, ] <- c(quadraturepoints[i], fixef(modAGQ), sqrt(diag(vcov(modAGQ))), sqrt(VarCorr(modAGQ)[[1]][1]))
}
melt3 <- cbind(melt(modAGQres3[,-1], id.vars = "nAGQ"), p = rep(quadraturepoints, 9))
ggplot(melt3, aes(p, value, color = Var2)) + geom_line(aes(colour = Var2)) 

# GBSDS 01234 vs 56 (20 is sufficient)
modAGQres4 <- matrix(nrow = length(quadraturepoints), ncol = 10)
colnames(modAGQres4) <- c("nAGQ", "(Intercept)", "treatIVIg+MP", "weeks", "treatIVIg+MP:weeks", "seInt", "seTreat", "seWeeks", "seTreat*Weeks", "varre")

for (i in 1:length(quadraturepoints)) {
  modAGQ <- glmer((GBSDS>4) ~ treat * weeks + (1 | id),  family = "binomial", nAGQ = quadraturepoints[i], data = mplong, control = glmerControl(optimizer="bobyqa"))
  modAGQres4[i, ] <- c(quadraturepoints[i], fixef(modAGQ), sqrt(diag(vcov(modAGQ))), sqrt(VarCorr(modAGQ)[[1]][1]))
}
melt4 <- cbind(melt(modAGQres4[,-1], id.vars = "nAGQ"), p = rep(quadraturepoints, 9))
ggplot(melt4, aes(p, value, color = Var2)) + geom_line(aes(colour = Var2)) 

#We compare the parameter estimates (fixed effects and random effects variance parameters) and standard errors with increasing number of quadrature points. 
#The estimates stabilize after 40 points. We would proceed with 40 points since it is no computationally intensive. 
#However, nAGQ > 1 is only available for models with a single, scalar random-effects term 
# So if random slopes are needed, we should use nAGQ = 0 or 1. From the plots, it is clear that nAGQ = 0 is closer to the estimates after stabilization than nAGQ = 1.
```
### Model building, backward selection (unadjusted analysis)
```{r}
# UNADJUSTED

# GBSDS 01 vs 23456
model1.1 <- glmer((GBSDS>1) ~ treat * weeks + (1 | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
model1.2 <- glmer((GBSDS>1) ~ treat * weeks + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(model1.1, model1.2) # based on LRT, the random slopes term is statistically significant
summary(model1.2) # try without interaction term

model1.3 <- glmer((GBSDS>1) ~ treat + weeks + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(model1.2, model1.3) # choose the simplest model: without interaction

summary(model1.3) # the final model (random intercepts and random slopes, no interaction terms)

model1.4 <- glmer((GBSDS>1) ~ treat + weeks + (1 | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
anova(model1.4, model1.3)


c(exp(summary(model1.3)$coef[2, 1]), exp(summary(model1.3)$coef[2, 1]+1.96 * summary(model1.3)$coef[2, 2]), exp(summary(model1.3)$coef[2, 1]-1.96 * summary(model1.3)$coef[2, 2]), summary(model1.3)$coef[2, 2])

# model1.4 <- glmer((GBSDS>1) ~ treat + weeks + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
# 
# anova(model1.4, model1.3) # choose the simplest model: without interaction
# 
# summary(model1.3)
# plot(model1.3, resid(., type = "pearson") ~ fitted(.), type = c("p", "smooth"))
# plot(model1.4, resid(., type = "pearson") ~ fitted(.), type = c("p", "smooth"))
# ggplot(data.frame(eta = predict(model1.4, type = "link"), pearson = residuals(model1.4, type = "pearson"), treat = mplong$treat[!is.na(mplong$GBSDS)]), aes(x=eta, y=pearson, group = treat)) + geom_point(aes(colour = treat))
# ggplot(data.frame(eta = predict(model1.3, type = "link"), pearson = residuals(model1.4, type = "pearson")), aes(x=eta, y=pearson)) + geom_point()


# GBSDS 012 vs 3456
model2.1 <- glmer((GBSDS>2) ~ treat * weeks + (1 | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
model2.2 <- glmer((GBSDS>2) ~ treat * weeks + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(model2.1, model2.2) # based on LRT, the random slopes term is statistically significant
summary(model2.2) # try without interaction term

model2.3 <- glmer((GBSDS>2) ~ treat + weeks + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
anova(model2.2, model2.3) # choose the simplest model: without interaction

summary(model2.3) # the final model (random intercepts and random slopes, no interaction terms) 
c(exp(summary(model2.3)$coef[2, 1]), exp(summary(model2.3)$coef[2, 1]+1.96 * summary(model2.3)$coef[2, 2]), exp(summary(model2.3)$coef[2, 1]-1.96 * summary(model2.3)$coef[2, 2]), summary(model2.3)$coef[2, 2])

# GBSDS 0123 vs 456
model3.1 <- glmer((GBSDS>3) ~ treat * weeks + (1 | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
model3.2 <- glmer((GBSDS>3) ~ treat * weeks + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(model3.1, model3.2) # based on LRT, the random slopes term is statistically significant
summary(model3.2) # try without interaction term

model3.3 <- glmer((GBSDS>3) ~ treat + weeks + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
anova(model3.2, model3.3) # choose the simplest model: without interaction

summary(model3.3) # the final model (random intercepts and random slopes, no interaction terms)
c(exp(summary(model3.3)$coef[2, 1]), exp(summary(model3.3)$coef[2, 1]+1.96 * summary(model3.3)$coef[2, 2]), exp(summary(model3.3)$coef[2, 1]-1.96 * summary(model3.3)$coef[2, 2]), summary(model3.3)$coef[2, 2])

# GBSDS 01234 vs 56
model4.1 <- glmer((GBSDS>4) ~ treat * weeks + (1 | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
model4.2 <- glmer((GBSDS>4) ~ treat * weeks + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(model4.1, model4.2) # based on LRT, the random slopes term is statistically significant
summary(model4.2) # try without interaction term

model4.3 <- glmer((GBSDS>4) ~ treat + weeks + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
anova(model4.2, model4.3) # choose the simplest model: without interaction

summary(model4.3) # the final model (random intercepts and random slopes, no interaction terms)
c(exp(summary(model4.3)$coef[2, 1]), exp(summary(model4.3)$coef[2, 1]+1.96 * summary(model4.3)$coef[2, 2]), exp(summary(model4.3)$coef[2, 1]-1.96 * summary(model4.3)$coef[2, 2]), summary(model4.3)$coef[2, 2])
```

### Optimal number of gaussian quadrature points (adjusted analysis)
```{r}
# --------------------- ADJUSTED

adjmodAGQres1 <- adjmodAGQres2 <- adjmodAGQres3 <- adjmodAGQres4 <- matrix(nrow = length(quadraturepoints), ncol = 18, dimnames = list(rep(NA, length(quadraturepoints)), c("nAGQ", "(Intercept)", "treatIVIg+MP", "weeks", "mrcss", "weakonse", "treatIVIg+MP:weeks", "treatIVIg+MP:mrcss", "treatIVIg+MP:weakonse ", "seInt", "seTreat", "seWeeks", "seMrcss", "seWeakonse", "seTreatIVIg+MP:weeks", "treatIVIg+MP:mrcss", "setreatIVIg+MP:weakonse", "varre")))

# GBSDS 01 vs 23456 (40 is sufficient)
for (i in 1:length(quadraturepoints)) {
  adjmodAGQ <- glmer((GBSDS>1) ~ treat * weeks + treat * mrcss + treat * weakonse +  (1 | id),  family = "binomial", nAGQ = quadraturepoints[i], data = mplong, control = glmerControl(optimizer="bobyqa"))
  adjmodAGQres1[i, ] <- c(quadraturepoints[i], fixef(adjmodAGQ), sqrt(diag(vcov(adjmodAGQ))), sqrt(VarCorr(adjmodAGQ)[[1]][1]))
}

adjmelt1 <- cbind(melt(adjmodAGQres1[,-1], id.vars = "nAGQ"), p = rep(quadraturepoints, 17))
ggplot(adjmelt1, aes(p, value, color = Var2)) + geom_line(aes(colour = Var2)) 

# GBSDS 012 vs 3456 (20 is sufficient)
# modAGQres2 <- matrix(nrow = length(quadraturepoints), ncol = 10)
# colnames(modAGQres2) <- c("nAGQ", "(Intercept)", "treatIVIg+MP", "weeks", "treatIVIg+MP:weeks", "seInt", "seTreat", "seWeeks", "seTreat*Weeks", "varre")

for (i in 1:length(quadraturepoints)) {
  adjmodAGQ <- glmer((GBSDS>2) ~ treat * weeks + treat * mrcss + treat * weakonse +   (1 | id),  family = "binomial", nAGQ = quadraturepoints[i], data = mplong, control = glmerControl(optimizer="bobyqa"))
  adjmodAGQres2[i, ] <- c(quadraturepoints[i], fixef(adjmodAGQ), sqrt(diag(vcov(adjmodAGQ))), sqrt(VarCorr(adjmodAGQ)[[1]][1]))
}
adjmelt2 <- cbind(melt(adjmodAGQres2[,-1], id.vars = "nAGQ"), p = rep(quadraturepoints, 17))
ggplot(adjmelt2, aes(p, value, color = Var2)) + geom_line(aes(colour = Var2)) 

# GBSDS 0123 vs 456 (20 is sufficient)
# modAGQres3 <- matrix(nrow = length(quadraturepoints), ncol = 10)
# colnames(modAGQres3) <- c("nAGQ", "(Intercept)", "treatIVIg+MP", "weeks", "treatIVIg+MP:weeks", "seInt", "seTreat", "seWeeks", "seTreat*Weeks", "varre")

for (i in 1:length(quadraturepoints)) {
  adjmodAGQ <- glmer((GBSDS>3) ~ treat * weeks + treat * mrcss + treat * weakonse +  (1 | id),  family = "binomial", nAGQ = quadraturepoints[i], data = mplong, control = glmerControl(optimizer="bobyqa"))
  adjmodAGQres3[i, ] <- c(quadraturepoints[i], fixef(adjmodAGQ), sqrt(diag(vcov(adjmodAGQ))), sqrt(VarCorr(adjmodAGQ)[[1]][1]))
}
adjmelt3 <- cbind(melt(adjmodAGQres3[,-1], id.vars = "nAGQ"), p = rep(quadraturepoints, 17))
ggplot(adjmelt3, aes(p, value, color = Var2)) + geom_line(aes(colour = Var2)) 

# GBSDS 01234 vs 56 (20 is sufficient)
# modAGQres4 <- matrix(nrow = length(quadraturepoints), ncol = 10)
# colnames(modAGQres4) <- c("nAGQ", "(Intercept)", "treatIVIg+MP", "weeks", "treatIVIg+MP:weeks", "seInt", "seTreat", "seWeeks", "seTreat*Weeks", "varre")

for (i in 1:length(quadraturepoints)) {
  adjmodAGQ <- glmer((GBSDS>4) ~ treat * weeks + treat * mrcss + treat * weakonse +   (1 | id),  family = "binomial", nAGQ = quadraturepoints[i], data = mplong, control = glmerControl(optimizer="bobyqa"))
  adjmodAGQres4[i, ] <- c(quadraturepoints[i], fixef(adjmodAGQ), sqrt(diag(vcov(adjmodAGQ))), sqrt(VarCorr(adjmodAGQ)[[1]][1]))
}
adjmelt4 <- cbind(melt(adjmodAGQres4[,-1], id.vars = "nAGQ"), p = rep(quadraturepoints, 17))
ggplot(adjmelt4, aes(p, value, color = Var2)) + geom_line(aes(colour = Var2)) 

#We compare the parameter estimates (fixed effects and random effects variance parameters) and standard errors with increasing number of quadrature points. 
#The estimates stabilize after 40 points. We would proceed with 40 points since it is no computationally intensive. 
#However, nAGQ > 1 is only available for models with a single, scalar random-effects term 
# So if random slopes are needed, we should use nAGQ = 0 or 1. From the plots, it is clear that nAGQ = 0 is closer to the estimates after stabilization than nAGQ = 1.
```
### Model building, backward selection (adjusted analysis)
```{r}
# ----------- ADJUSTED
# the models are fitted using the laplace approach (nAGQ = 1)

# GBSDS 01 vs 23456
adj.model1.1 <- glmer((GBSDS>1) ~ treat * weeks + treat * mrcss + treat * weakonse +  (1 | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
adj.model1.2 <- glmer((GBSDS>1) ~ treat * weeks + treat * mrcss + treat * weakonse +  (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model1.1, adj.model1.2) # based on LRT, the random slopes term is statistically significant
summary(adj.model1.2) # try without interaction treatment by weeks term

adj.model1.3 <- glmer((GBSDS>1) ~ treat + weeks + treat * mrcss + treat * weakonse +  (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model1.2, adj.model1.3) # choose the simplest model: without interaction treatment by weeks
summary(adj.model1.3) # try without interaction treatment by mrcss term

adj.model1.4 <- glmer((GBSDS>1) ~ treat + weeks + mrcss + treat * weakonse +  (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model1.3, adj.model1.4) # choose the simplest model: without interaction treatment by mrcss
summary(adj.model1.4) # try without interaction treatment by weakonset term

adj.model1.5 <- glmer((GBSDS>1) ~ treat + weeks + mrcss + weakonse +(weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model1.4, adj.model1.5) # choose the simplest model: without interaction treatment by weakonset

summary(adj.model1.5)  # the final model (random intercepts and random slopes, no interaction terms)
c(exp(summary(adj.model1.5)$coef[2, 1]), exp(summary(adj.model1.5)$coef[2, 1]+1.96 * summary(adj.model1.5)$coef[2, 2]), exp(summary(adj.model1.5)$coef[2, 1]-1.96 * summary(adj.model1.5)$coef[2, 2]), summary(adj.model1.5)$coef[2, 2])



# GBSDS 012 vs 3456
adj.model2.1 <- glmer((GBSDS>2) ~ treat * weeks + treat * mrcss + treat * weakonse + (1 | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
adj.model2.2 <- glmer((GBSDS>2) ~ treat * weeks + treat * mrcss + treat * weakonse + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model2.1, adj.model2.2) # based on LRT, the random slopes term is statistically significant
summary(adj.model2.2) # try without interaction treatment by weeks term

adj.model2.3 <- glmer((GBSDS>2) ~ treat + weeks + treat * mrcss + treat * weakonse + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model2.2, adj.model2.3) # choose the simplest model: without interaction treatment by weeks
summary(adj.model2.3) # try without interaction treatment by weakonset term

adj.model2.4 <- glmer((GBSDS>2) ~ treat + weeks + treat * mrcss + weakonse + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model2.3, adj.model2.4) # choose the simplest model: without interaction treatment by weakonset
summary(adj.model2.4) # try without interaction treatment by mrcss term

adj.model2.5 <- glmer((GBSDS>2) ~ treat + weeks + mrcss + weakonse +(weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model2.4, adj.model2.5) # choose the simplest model: without interaction treatment by mrcss

summary(adj.model2.5)  # the final model (random intercepts and random slopes, no interaction terms)
c(exp(summary(adj.model2.5)$coef[2, 1]), exp(summary(adj.model2.5)$coef[2, 1]+1.96 * summary(adj.model2.5)$coef[2, 2]), exp(summary(adj.model2.5)$coef[2, 1]-1.96 * summary(adj.model2.5)$coef[2, 2]), summary(adj.model2.5)$coef[2, 2])



# GBSDS 0123 vs 456
adj.model3.1 <- glmer((GBSDS>3) ~ treat * weeks + treat * mrcss + treat * weakonse +  (1 | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
adj.model3.2 <- glmer((GBSDS>3) ~ treat * weeks + treat * mrcss + treat * weakonse +  (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model3.1, adj.model3.2) # based on LRT, the random slopes term is statistically significant
summary(adj.model3.2) # try without interaction treatment by weeks term

adj.model3.3 <- glmer((GBSDS>3) ~ treat + weeks + treat * mrcss + treat * weakonse + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model3.2, adj.model3.3) # choose the simplest model: without interaction treatment by weeks
summary(adj.model3.3) # try without interaction treatment by weakonset term

adj.model3.4 <- glmer((GBSDS>3) ~ treat + weeks + treat * mrcss + weakonse + (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model3.3, adj.model3.4) # choose the simplest model: without interaction treatment by weakonset
summary(adj.model3.4) # try without interaction treatment by mrcss term

adj.model3.5 <- glmer((GBSDS>3) ~ treat + weeks + mrcss + weakonse +(weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model3.4, adj.model3.5) # choose the simplest model: without interaction treatment by mrcss

summary(adj.model3.5)  # the final model (random intercepts and random slopes, no interaction terms)
c(exp(summary(adj.model3.5)$coef[2, 1]), exp(summary(adj.model3.5)$coef[2, 1]+1.96 * summary(adj.model3.5)$coef[2, 2]), exp(summary(adj.model3.5)$coef[2, 1]-1.96 * summary(adj.model3.5)$coef[2, 2]), summary(adj.model3.5)$coef[2, 2])

# GBSDS 01234 vs 56
adj.model4.1 <- glmer((GBSDS>4) ~ treat * weeks + treat * mrcss + treat * weakonse +  (1 | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))
adj.model4.2 <- glmer((GBSDS>4) ~ treat * weeks + treat * mrcss + treat * weakonse +  (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model4.1, adj.model4.2) # based on LRT, the random slopes term is statistically significant
summary(adj.model4.2) #  try without interaction treatment by weakonset term

adj.model4.3 <- glmer((GBSDS>4) ~  treat * weeks + treat * mrcss + weakonse +  (weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model4.2, adj.model4.3) # choose the simplest model: without interaction treatment by weakonset
summary(adj.model4.3) # try without interaction treatment by weeks

adj.model4.4 <- glmer((GBSDS>4) ~ treat + weeks + treat * mrcss + weakonse +(weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model4.3, adj.model4.4) # choose the simplest model: without interaction treatment by weeks
summary(adj.model4.4)  # try without interaction treatment by mrcss

adj.model4.5 <- glmer((GBSDS>4) ~ treat + weeks + mrcss + weakonse +(weeks | id),  family = "binomial", nAGQ = 0, data = mplong, control = glmerControl(optimizer="bobyqa"))

anova(adj.model4.4, adj.model4.5) # choose the simplest model: without interaction treatment by mrcss

summary(adj.model4.5)  # the final model (random intercepts and random slopes, no interaction terms)
c(exp(summary(adj.model4.5)$coef[2, 1]), exp(summary(adj.model4.5)$coef[2, 1]+1.96 * summary(adj.model4.5)$coef[2, 2]), exp(summary(adj.model4.5)$coef[2, 1]-1.96 * summary(adj.model4.5)$coef[2, 2]), summary(adj.model4.5)$coef[2, 2])
```

## Longitudinal proportional odds models 
### Model building unadjusted LPO using compareBmods() 
```{r}
# Start with an elaborate fixed effects part
# LPO with two random components vs LPO with one random component 
 
# sqrt of weeks seems to be appropriate based on the plots

# start with an elaborate fixed part, main effects and interaction effects
# test for one versus two random effects
unadj.LPO1 <- blrm(formula = GBSDS ~ sqrtweeks * treat + cluster(id), 
                   data = mplong,
                   iter = 5000, # by default half of iterations as burn in
                   file = 'unadj.LPO1.Rdata',
                   loo = T,
                   seed = 1) 

unadj.LPO2 <- blrm(formula = GBSDS ~ time(sqrtweeks) * treat + cluster(id), 
                   data = mplong,
                   iter = 5000, # by default half of iterations as burn in
                   file = 'unadj.LPO2.Rdata',
                   loo = T,
                   seed = 1) 

compareBmods(unadj.LPO1, unadj.LPO2)

# Model with both random slopes and random intercepts has a larger weight, indicating a better fit (McElreath's book)
# Continue with random intercept only model 

# test for various correlation structures not possible

# test for different variance covariance structures for our random intercepts and random slopes not possible

# test whether the fixed part can be simplified
unadj.LPO3 <- blrm(formula = GBSDS ~ sqrtweeks + treat + cluster(id), loo = TRUE, 
                                   data = mplong, 
                                   iter = 5000, # by default half of iterations as burn in
                                   file = 'unadj.LPO3.Rdata',
                                   seed = 1) 

compareBmods(unadj.LPO1, unadj.LPO2, unadj.LPO3) 
# model 1 most likely to fit, interaction term is left in 

# Final model:
# Interaction term weeks * treat
# random intercepts only
# variance covariance matrix and correlation structures not possible

# Whereas final model based on lme was:
# no interaction terms
# random intercepts and random slopes
# autocorrelation structure 
# unstructured variance covariance matrix
```
 
### Model building unadjusted LPO (as if continuous)
```{r}
# sqrt of weeks seems to be appropriate based on the plots

# start with an elaborate fixed part, main effects and interaction effects
# test for one versus two random effects
unadj.lme1 <- nlme::lme(GBSDS ~ sqrtweeks * treat, 
                  data = mplong, 
                  random = ~ 1 | id,
                  na.action = na.omit)

unadj.lme2 <- nlme::lme(GBSDS ~ sqrtweeks * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  na.action = na.omit)

anova(unadj.lme1, unadj.lme2)
# model with both random slopes and random intercepts has a larger loglikelihood, indicating a better fit
# The LRT shows that we can reject the null hypothesis of no difference in fit between the two models. 
# Note that anova method uses the chi-squared asymptotic distribution. Here you need to consider the proper mixture of Chi-squared distributions.
# The model with random slopes and random intercepts fitted the data significantly better than the model with only random intercepts.
# Continue with random slope and random intercept model 

# test for various correlation structures
unadj.lme3 <- nlme::lme(GBSDS ~ sqrtweeks * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corAR1(),
                  na.action = na.omit)

unadj.lme4 <- nlme::lme(GBSDS ~ sqrtweeks * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corCAR1(),
                  na.action = na.omit)

unadj.lme5 <- nlme::lme(GBSDS ~ sqrtweeks * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corCompSymm(),
                  na.action = na.omit)

unadj.lme6 <- nlme::lme(GBSDS ~ sqrtweeks * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corExp(),
                  na.action = na.omit)

unadj.lme7 <- nlme::lme(GBSDS ~ sqrtweeks * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corRatio(),
                  na.action = na.omit)

# Correlation structures corLin, CorSpher, CorSymm do not converge 
# and correlation structure corARMA is not available with our formula

anova(unadj.lme2, unadj.lme3, unadj.lme4, unadj.lme5, unadj.lme6, unadj.lme7)

# corCAR has the best fit and is better than the default correlation structure (corresponding to no within subject correlation)
# Continue with corCAR1, as it allows for unequally spaced continuous time intervals (Singer and Willet 2003) and is simpler to implement than corExp.

# Next, try different variance covariance structures for our random intercepts and random slopes 
# Compare unstructured variance covariance matrix with a diagonal structure, identical structure and compound symmetry structure 

unadj.lme8 <- nlme::lme(GBSDS ~ sqrtweeks * treat, 
                  data = mplong, 
                  random = list(id = pdDiag(form = ~ sqrtweeks)),
                  correlation = corCAR1(),
                  method = "REML",
                  na.action = na.omit)

unadj.lme9 <- nlme::lme(GBSDS ~ sqrtweeks * treat, 
                  data = mplong, 
                  random = list(id = pdIdent(form = ~ sqrtweeks)),
                  correlation = corCAR1(),
                  method = "REML",
                  na.action = na.omit)

unadj.lme10 <- nlme::lme(GBSDS ~ sqrtweeks * treat, 
                  data = mplong, 
                  random = list(id = pdCompSymm(form = ~ sqrtweeks)),
                  correlation = corCAR1(),
                  method = "REML",
                  na.action = na.omit)

anova(unadj.lme4, unadj.lme8, unadj.lme9, unadj.lme10) 

# Different variance covariance structtures do not fit better than the unstructured variance covariance matrix
# Continue with unstructured variance covariance matrix (unadj.lme4 is the winner)

plot(unadj.lme4, resid(., type= "pearson") ~ fitted(.) | treat,
     type = c("p", "smooth"),
     lwd = 3)

# standardized residuals against the fitted values
# standardized resisuals increase as fitted values increase
# log and square root transformation

unadj.lme4log <- nlme::lme(log(GBSDS, base = 10) ~ sqrtweeks * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corCAR1(),
                  na.action = na.omit)
plot(unadj.lme4log, resid(., type= "pearson") ~ fitted(.) | treat,
     type = c("p", "smooth"),
     lwd = 3)


unadj.lme4sqrt <- nlme::lme(sqrt(GBSDS) ~ sqrtweeks * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corCAR1(),
                  na.action = na.omit)
plot(unadj.lme4sqrt, resid(., type= "pearson") ~ fitted(.) | treat,
     type = c("p", "smooth"),
     lwd = 3)

# The transformations do not improve the fit
# Continue with untransformed dependent variable (lme4 still the winner)

# Now, simplify the fixed part
anova(unadj.lme4)
# remove interaction term sqrtweeks:treat

unadj.lme4_2 <- nlme::lme(GBSDS ~ sqrtweeks + treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corCAR1(),
                  na.action = na.omit)
anova(unadj.lme4_2)

plot(unadj.lme4_2, resid(., type= "pearson") ~ fitted(.) | treat,
     type = c("p", "smooth"),
     lwd = 3)b

# Final model:
# no interaction terms
# random intercepts and random slopes
# autocorrelation structure 
# unstructured variance covariance matrix
```
 
 
### Model building adjusted LPO using compareBmods()
```{r}
# Start with an elaborate fixed effects part
# LPO with two random components vs LPO with one random component 
 
# sqrt of weeks seems to be appropriate based on the plots

# start with an elaborate fixed part, main effects and interaction effects
# test for one versus two random effects
adj.LPO.1 <- blrm(formula = GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse * treat + cluster(id), 
                  data = mplong,
                  iter = 5000, # by default half of iterations as burn in
                  file = 'adj.LPO.1.Rdata',
                  loo = TRUE,
                  seed = 1) 

memory.limit(size=56000)
adj.LPO.2 <- blrm(formula = GBSDS ~ time(sqrtweeks) * treat + mrcss * treat + weakonse * treat + cluster(id), 
                  data = mplong,
                  iter = 5000, # by default half of iterations as burn in
                  file = 'adj.LPO.2.Rdata',
                  loo = TRUE,
                  seed = 1) 

compareBmods(adj.LPO.1, adj.LPO.2)

# model with random intercepts only has a larger weight, indicating a better fit (McElreath's book)
# Continue with random intercept model 

# test for various correlation structures not possible

# test for different variance covariance structures for our random intercepts and random slopes not possible

# test whether the fixed part can be simplified
adj.LPO.1 # leave interaction term treat=IVIg+MP * mrcss out

adj.LPO.3 <- blrm(formula = GBSDS ~ sqrtweeks + treat + mrcss * treat + weakonse * treat +  cluster(id), loo = TRUE, 
                                   data = mplong, 
                                   iter = 5000, # by default half of iterations as burn in
                                   file = 'adj.LPO3.Rdata',
                                   seed = 1) 

adj.LPO.4 <- blrm(formula = GBSDS ~ sqrtweeks * treat + mrcss + treat + weakonse * treat +  cluster(id), loo = TRUE, 
                                   data = mplong, 
                                   iter = 5000, # by default half of iterations as burn in
                                   file = 'adj.LPO4.Rdata',
                                   seed = 1) 

adj.LPO.5 <- blrm(formula = GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse + treat +  cluster(id), loo = TRUE, 
                                   data = mplong, 
                                   iter = 5000, # by default half of iterations as burn in
                                   file = 'adj.LPO5.Rdata',
                                   seed = 1) 

compareBmods(adj.LPO.1, adj.LPO.4) # which (if any) interaction term can be left out?



# Final model:
# all interaction terms
# random intercepts only
```
 
 
### Model building adjusted LPO (as if continuous)
```{r}
# sqrt of weeks seem to be appropriate based on the plots

library(nlme)

# start with an elaborate fixed part, main effects and interaction effects
# test for one versus two random effects
lme1 <- nlme::lme(GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse * treat, 
                  data = mplong, 
                  random = ~ 1 | id,
                  na.action = na.omit)

lme2 <- nlme::lme(GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  na.action = na.omit)

anova(lme1, lme2)
# model with both random slopes and random intercepts has a larger loglikelihood, indicating a better fit
# The LRT shows that we can reject the null hypothesis of no difference in fit between the two models. 
# Note that anova method uses the chi-squared asymptotic distribution. Here you need to consider the proper mixture of Chi-squared distributions.
# The model with random slopes and random intercepts fitted the data significantly better than the model with only random intercepts.
# Continue with random slope and random intercept model 

# test for various correlation structures
lme3 <- nlme::lme(GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corAR1(),
                  na.action = na.omit)

lme4 <- nlme::lme(GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corCAR1(),
                  na.action = na.omit)

lme5 <- nlme::lme(GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corCompSymm(),
                  na.action = na.omit)

lme6 <- nlme::lme(GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corExp(),
                  na.action = na.omit)

lme7 <- nlme::lme(GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corRatio(),
                  na.action = na.omit)

# Correlation structures corLin, CorSpher, CorSymm do not converge 
# and correlation structure corARMA is not available with our formula

anova(lme2, lme3, lme4, lme5, lme6, lme7)

# corAR, corCAR and corExp have an equal fit and are better than the default correlation structure (corresponding to no within subject correlation)
# Continue with corCAR1, as it allows for unequally spaced continuous time intervals (Singer and Willet 2003) and is simpler to implement than corExp.

# Next, try different variance covariance structures for our random intercepts and random slopes 
# Compare unstructured variance covariance matrix with a diagonal structure, identical structure and compound symmetry structure 

lme8 <- nlme::lme(GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse * treat, 
                  data = mplong, 
                  random = list(id = pdDiag(form = ~ sqrtweeks)),
                  correlation = corCAR1(),
                  method = "REML",
                  na.action = na.omit)

lme9 <- nlme::lme(GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse * treat, 
                  data = mplong, 
                  random = list(id = pdIdent(form = ~ sqrtweeks)),
                  correlation = corCAR1(),
                  method = "REML",
                  na.action = na.omit)

lme10 <- nlme::lme(GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse * treat, 
                  data = mplong, 
                  random = list(id = pdCompSymm(form = ~ sqrtweeks)),
                  correlation = corCAR1(),
                  method = "REML",
                  na.action = na.omit)

anova(lme4, lme8, lme9, lme10) 

# Different variance covariance structures do not fit better than the unstructured variance covariance matrix
# Continue with unstructured variance covariance matrix (lme4 is the winner)

plot(lme4, resid(., type= "pearson") ~ fitted(.) | treat,
     type = c("p", "smooth"),
     lwd = 3)

# standardized residuals against the fitted values
# standardized resisuals increase as fitted values increase
# log and square root transformation

lme4log <- nlme::lme(log(GBSDS, base = 10) ~ sqrtweeks * treat + mrcss * treat + weakonse * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corCAR1(),
                  na.action = na.omit)
plot(lme4log, resid(., type= "pearson") ~ fitted(.) | treat,
     type = c("p", "smooth"),
     lwd = 3)


lme4sqrt <- nlme::lme(sqrt(GBSDS) ~ sqrtweeks * treat + mrcss * treat + weakonse * treat, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corCAR1(),
                  na.action = na.omit)
plot(lme4sqrt, resid(., type= "pearson") ~ fitted(.) | treat,
     type = c("p", "smooth"),
     lwd = 3)

# The transformations do not improve the fit
# Continue with untransformed dependent variable (lme4 still the winner)

# Now, simplify the fixed part
anova(lme4)
# remove treatIVIg+MP:weakonse

lme4_2 <- nlme::lme(GBSDS ~ sqrtweeks * treat + mrcss * treat + weakonse, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corCAR1(),
                  na.action = na.omit)

anova(lme4_2)
# remove treatIVIg+MP:mrcss

lme4_3 <- nlme::lme(GBSDS ~ sqrtweeks * treat + mrcss  + weakonse, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corCAR1(),
                  na.action = na.omit)

anova(lme4_3)
# remove sqrtweeks:treatIVIg+MP ( so none of the interactions were sign.)

lme4_4 <- nlme::lme(GBSDS ~ treat + sqrtweeks + mrcss + weakonse, 
                  data = mplong, 
                  random = ~ sqrtweeks | id,
                  correlation = corCAR1(),
                  na.action = na.omit)

anova(lme4_4)

plot(lme4_4, resid(., type= "pearson") ~ fitted(.) | treat,
     type = c("p", "smooth"),
     lwd = 3)

# Final model:
# no interactions
# random intercepts and random slopes
# autocorrelation structure 
# unstructured variance covariance matrix
```

## The final LPOs unadjusted and adjusted
```{r}
LPO.unadj <- blrm(formula = GBSDS ~ time(sqrtweeks) + treat + cluster(id), 
                  data = mplong,
                  iter = 5000, # by default half of iterations as burn in
                  file = 'LPO.unadj.Rdata',
                  loo = T,
                  seed = 1)

LPO.adj <- blrm(formula = GBSDS ~ time(sqrtweeks) + treat + mrcss + weakonse + cluster(id), 
                data = mplong,
                iter = 5000, # by default half of iterations as burn in
                file = 'LPO.adj.Rdata',
                loo = T,
                seed = 1)
```
