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
library(glmnet)
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
### Mittal, Madigan et al. (2013)
library(Cyclops)
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
p <- 1000
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
# covariateValues <- lapply(1:n, function(i) {
#   return(mvrnorm(1, covariateMeanVector, covariateVarianceMatrix))
# })
# cdfValues <- runif(n)
# realTimeValues <- lapply(1:n, function(i) {
#   return(-log(cdfValues[i]) / (lambda * exp(sum(
#     beta * covariateValues[[i]]
#   ))))
# })
# censureValues <- rbinom(n, 1, 1 - censureRate)
# observedTimeValues <- lapply(1:n, function(i) {
#   if (censureValues[i] == 0) {
#     return(runif(1, 0, realTimeValues[[i]]))
#   }
#   return(realTimeValues[[i]])
# })
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

### 'coxpath met très longtemps à tourner ...
# res.coxpath <-
#   coxpath(data = list(
#     x = as.matrix(donnees[, covariateNames]),
#     time = donnees$time_obs,
#     status = donnees$status
#   ))

res.penalized.nopen <- penalized(formule, data = donnees)
### Très long la suite :
res.penalized.lasso <-
  penalized(formule, data = donnees, lambda1 = 1)
# res.penalized.ridge <-
#   penalized(formule, data = donnees, lambda2 = 1)
# res.penalized.elasticNet <-
#   penalized(formule,
#             data = donnees,
#             lambda1 = 1,
#             lambda2 = 1)

# res.pencoxfrail <- pencoxfrail(fix = formule, data = donnees, xi=10)
res.pencoxfrail <-
  pencoxfrail(
    fix = formule,
    vary.coef = "~",
    rnd = list(),
    data = donnees,
    xi = 10
  )

res.cv.glmnet <- cv.glmnet(
  x = as.matrix(donnees[, covariateNames]),
  y = Surv(time = donnees$time_obs, donnees$status),
  family = "cox"
)
test <- cv.glmnet(
  x = as.matrix(donnees[, covariateNames]),
  y = Surv(time = donnees$time_obs, donnees$status),
  family = "cox",
  alpha = 0
)
res.glmnet.min <-
  glmnet(
    x = as.matrix(donnees[, covariateNames]),
    y = Surv(time = donnees$time_obs, donnees$status),
    family = "cox",
    lamdda = c(1, .5)
  )
res.glmnet.1se <-
  glmnet(
    x = as.matrix(donnees[, covariateNames]),
    y = Surv(time = donnees$time_obs, donnees$status),
    family = "cox",
    lamdda = res.cv.glmnet$lambda.1se
  )

donnees.cyclops <- createCyclopsData(formula = formule,
                                     data = donnees,
                                     modelType = "cox")
res.cyclops.nopen <-
  fitCyclopsModel(donnees.cyclops)
res.cyclops.lasso <-
  fitCyclopsModel(donnees.cyclops,
                  prior = createPrior("laplace", useCrossValidation = TRUE))
res.cyclops.ridge <-
  fitCyclopsModel(donnees.cyclops,
                  prior = createPrior("normal", useCrossValidation = TRUE))

### Comparaison des coefficients obtenus :
coeffFound <-
  data.frame(
    true_beta = beta,
    coxph = res.coxph$coefficients,
    penalized_nopen = res.penalized.nopen@penalized,
    penalized_lasso = res.penalized.lasso@penalized,
    # penalized_ridge = res.penalized.ridge@penalized,
    # penalized_elastic_net = res.penalized.elasticNet@penalized,
    # cyclops_nopen = res.cyclops.nopen$estimation$estimate,
    cyclops_lasso = res.cyclops.lasso$estimation$estimate
    #cyclops_ridge = res.cyclops.ridge$estimation$estimate
  )
