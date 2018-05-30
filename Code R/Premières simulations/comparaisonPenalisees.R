######################################################################
########## ----- COMPARAISON DES EXISTANTS PENALISES
########## ----- (29/05/2018)
##########
######################################################################

rm(list = ls())
gc()
setwd("C:/Users/marie/Desktop/Samuel/aCoxStory/Code R/Premières simulations")
print(getwd())

# source("cmprskFrail.R")
# library(survival)
# library(coxme)
# library(MASS)
# library(cmprsk)

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
if (FALSE) {
  #url <- "https://cran.r-project.org/src/contrib/Archive/Coxnet/Coxnet_0.2.tar.gz"
  #pkgFile <- "Coxnet_0.2.tar.gz"
  url <-
    "https://cran.r-project.org/src/contrib/Archive/Coxnet/Coxnet_0.1-1.tar.gz"
  pkgFile <- "Coxnet_0.1-1.tar.gz"
  download.file(url = url, destfile = pkgFile)
  
  # Install dependencies
  
  install.packages(c("ada", "ipred", "evd"))
  
  # Install package
  install.packages(pkgs = pkgFile,
                   type = "source",
                   repos = NULL)
  
  # Delete package tarball
  unlink(pkgFile)
}
### Groll, Hastie et Tutz (2017)
library(PenCoxFrail)

### Données
load("datCrPack.rda")
lung <- lung
summary(datCrPack)

######################################################################
########## Simulation d'un modèle naïf (sans corrélations)
########## -- on considère une fonction de base exponentielle

# Paramètres structurants
n <- 1000
p <- 100
censureRate <- .4
nbNonSparse <- 10
possibleBeta <- c(-10:-1, 1:10)
lambda <- 10

# Structure des données
covariateNames <- paste("X", 1:p, sep = "_")
covariateMeanVector <- rep(0, p)
covariateVarianceMatrix <- diag(1, p, p)
beta <- rep(0, p)
nonSparseIndex <- sample(1:length(beta), nbNonSparse)
beta[nonSparseIndex] <-
  sample(possibleBeta, nbNonSparse, replace = TRUE)

# Simulation des données brutes
covariateValues <- lapply(1:n, function(i) {
  return(mvrnorm(1, covariateMeanVector, covariateVarianceMatrix))
})
cdfValues <- runif(n)
realTimeValues <- lapply(1:n, function(i) {
  return(-log(cdfValues[i]) / (lambda * exp(sum(
    beta * covariateValues[[i]]
  ))))
})
censureValues <- rbinom(n, 1, 1 - censureRate)
observedTimeValues <- lapply(1:n, function(i) {
  if (censureValues[i] == 0) {
    return(runif(1, 0, realTimeValues[[i]]))
  }
  return(realTimeValues[[i]])
})

# Mise en forme du jeu de données
donnees <-
  data.frame(
    time_real = do.call(c, realTimeValues),
    time_obs = do.call(c, observedTimeValues),
    status = censureValues
  )
covariateValues <- as.data.frame(do.call(rbind, covariateValues))
colnames(covariateValues) <- covariateNames
donnees <- cbind(donnees, covariateValues)

# Nettoyage
rm(
  list = c(
    "cdfValues",
    "censureValues",
    "covariateValues",
    "observedTimeValues",
    "realTimeValues"
  )
)

######################################################################
########## On teste sur la base 'donnees'

formule <-
  as.formula(paste(
    "Surv(time_obs, status) ~ ",
    paste(covariateNames, collapse = " + "),
    sep = ""
  ))
res.coxph <- coxph(formula = formule, data = donnees)
### 'coxme' pose problème parce qu'il lui faut des effets mixtes (fixes et aléatoires)
# res.coxme <- coxme(formula = formule, data = donnees)


coeffFound <-
  data.frame(true_beta = beta, coxph = res.coxph$coefficients)

