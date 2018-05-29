######################################################################
########## ----- PREMIERS ESSAIS SUR LA METHODE
########## ----- (22/05/2018)
##########
######################################################################

rm(list = ls())
setwd("C:/Users/marie/Desktop/Samuel/aCoxStory/Code R/Premières simulations")
print(getwd())
source("cmprskFrail.R")
library(survival)
library(coxme)
library(MASS)
library(cmprsk)

load("datCrPack.rda")
summary(datCrPack)

######################################################################
########## Modèle avec risque compétitif

### Méthode 'crr' (package 'cmprsk')
crr(datCrPack$time, datCrPack$cens, datCrPack$group, failcode = 2)

# convergence:  TRUE
# coefficients:
#   datCrPack$group1
# -0.703
# standard errors:
#   [1] 0.3368
# two-sided p-values:
#   datCrPack$group1
# 0.037

######################################################################
########## Modèle avec fragilités

### Méthode 'coxph' (package 'survival')
coxph(formula = Surv(time, status) ~ ph.ecog + frailty.gaussian(inst),
      data = lung)

# Call:
#   coxph(formula = Surv(time, status) ~ ph.ecog + frailty.gaussian(inst),
#         data = lung)
#
# coef se(coef)    se2  Chisq   DF       p
# ph.ecog                 0.507    0.117  0.115 18.835 1.00 1.4e-05
# frailty.gaussian(inst)                         3.253 2.59    0.29
#
# Iterations: 8 outer, 31 Newton-Raphson
# Variance of random effect= 0.0221
# Degrees of freedom for terms= 1.0 2.6
# Likelihood ratio test=23.3  on 3.56 df, p=6.82e-05
# n=226 (2 observations deleted due to missingness)

### Méthode 'coxme' (package 'coxme')
coxme(formula = Surv(time, status) ~ ph.ecog + (1 | inst),
      data = lung)

# Cox mixed-effects model fit by maximum likelihood
# Data: lung
# events, n = 163, 226 (2 observations deleted due to missingness)
# Iterations= 11 47
# NULL Integrated    Fitted
# Log-likelihood -739.375    -730.47 -727.6774
#
# Chisq  df          p   AIC  BIC
# Integrated loglik 17.81 2.0 1.3571e-04 13.81 7.62
# Penalized loglik 23.40 3.6 6.7858e-05 16.20 5.08
#
# Model:  Surv(time, status) ~ ph.ecog + (1 | inst)
# Fixed coefficients
# coef exp(coef)  se(coef)    z       p
# ph.ecog 0.5073993  1.660966 0.1167809 4.34 1.4e-05
#
# Random effects
# Group Variable  Std Dev    Variance
# inst  Intercept 0.15058058 0.02267451

######################################################################
########## Modèle avec pénalisation L2 (Ridge)

coxph(formula = Surv(time, status) ~ ridge(ph.ecog, theta = 1),
      data = lung)

# Call:
#   coxph(formula = Surv(time, status) ~ ridge(ph.ecog, theta = 1),
#         data = lung)
#
# coef se(coef)    se2  Chisq DF       p
# ridge(ph.ecog)  0.473    0.113  0.113 17.509  1 2.9e-05
#
# Iterations: 1 outer, 3 Newton-Raphson
# Degrees of freedom for terms= 1
# Likelihood ratio test=17.6  on 0.99 df, p=2.72e-05
# n=227 (1 observation deleted due to missingness)

######################################################################
########## Modèle combiné frailty / rique compétitif

modelgauss = cmprskFrail(
  formula = Surv(time, cens) ~ group + frailty(CENTRE),
  data = datCrPack,
  distrib = "gaussian"
)
modelgam = cmprskFrail(
  formula = Surv(time, cens) ~ group + frailty(CENTRE),
  data = datCrPack,
  distrib = "gamma"
)

datCrPack$CENTRE2 = ifelse(as.numeric(as.character(datCrPack$CENTRE)) >
                             4, 1, 0)
modelgam2 = cmprskFrail(
  formula = Surv(time, cens) ~ group + frailty(CENTRE2),
  data = datCrPack,
  distrib = "gamma"
)

#Print et plot
print(modelgauss)
print(modelgam)
plot(modelgam, addCR = F)
plot(modelgam, addCR = T)
