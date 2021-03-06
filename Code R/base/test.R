##################################################################
############################# TESTS ##############################
##################################################################


setwd("C:/Users/Marie DA/Dropbox/a_flux_marie/th�se/Package cmpFrail/figure crcox")
library(survival); library(coxme); library(MASS); library(cmprsk)

load("datCrPack.rda")
summary(datCrPack)

#Mod�le avec risque comp�titif
crr(datCrPack$time,datCrPack$cens,datCrPack$group,failcode=2)

#Mod�le avec Frailty
coxph(formula = Surv(time, status) ~ ph.ecog + frailty.gaussian(inst), data = lung)
coxme(formula = Surv(time, status) ~ ph.ecog+ (1|inst), data=lung)

#Mod�le avec Ridge
coxph(formula = Surv(time, status) ~ ridge(ph.ecog,theta=1), data = lung)

#Mod�le combin� frailty/rique comp�titif
modelgauss=cmprskFrail(formula=Surv(time,cens)~ group + frailty(CENTRE),data=datCrPack,distrib="gaussian")
modelgam=cmprskFrail(formula=Surv(time,cens)~ group + frailty(CENTRE),data=datCrPack,distrib="gamma")

datCrPack$CENTRE2=ifelse(as.numeric(as.character(datCrPack$CENTRE)) >4,1,0)
modelgam2=cmprskFrail(formula=Surv(time,cens)~ group + frailty(CENTRE2),data=datCrPack,distrib="gamma")

#Print et plot
print(modelgauss)
print(modelgam)
plot(modelgam,addCR=F)
plot(modelgam,addCR=T)
 
 