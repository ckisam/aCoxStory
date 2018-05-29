######################################################################
########## ----- COMPARAISON DES EXISTANTS PENALISES
########## ----- (29/05/2018)
##########
######################################################################

rm(list = ls())
setwd("C:/Users/marie/Desktop/Samuel/aCoxStory/Code R/Premières simulations")
print(getwd())

# source("cmprskFrail.R")
# library(survival)
# library(coxme)
# library(MASS)
# library(cmprsk)

load("datCrPack.rda")
lung <- lung
summary(datCrPack)

library(MASS)

### coxph
library(survival)
### coxme
library(coxme)

### Park et Hastie (2007)
library(glmpath)
### Goeman (2010)
library(penalized)
### Simon, Friedman, Hastie et Tibshirani (2011)
library(coxnet)
### Groll, Hastie et Tutz (2017)
library(PenCoxFrail)


######################################################################
########## Simulation d'un modèle naïf (sans corrélations)
########## -- on considère une fonction de base exponentielle

# Paramètres structurants
n <- 100
p <- 100
censureRate <- .4
nbNonSparse <- 10
possibleBeta <- c(-10:-1, 1:10)
lambda <- 10

# Structure des données
covariateMeanVector <- rep(0, p)
covariateVarianceMatrix <- diag(1, p, p)
beta <- rep(0, p)
nonSparseIndex <- sample(1:length(beta), nbNonSparse)
beta[nonSparseIndex] <-
  sample(possibleBeta, nbNonSparse, replace = TRUE)

# Simulation des données
covariateValues <- lapply(1:n, function(i) {
  return(mvrnorm(1, covariateMeanVector, covariateVarianceMatrix))
})
cdfValues <- runif(n)
timeValues <- lapply(1:n, function(i) {
  return(-log(cdfValues[i]) / (lambda * exp(sum(
    beta * covariateValues[[i]]
  ))))
})
censureValues <- rbinom(n, 1, censureRate)
