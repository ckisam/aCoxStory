# cmprskFrail: R function for Katsahian et al. (2011) and for gamma frailties in a competing risk model
# =====================

# Authors / Version:
# =====================
# M. De Antonio, C. Boudreau, H. Varet
# Last updated: Feb 20, 2018

# Packages:           
# =====================
# make sure that the following packages are loaded
# require(survival); require(coxme); require(MASS); require(data.table)

# @todo V2
# =====================
# différente manière d'entrer la survie (surv(u~d) ou start=stop=..) ou  style crr plutôt que coxph ?
# vérifier les ...
# regarder les différences de calculs likelihood
# covariates can NOT be named weights and/or Gj and/or tmp (à revoir)
 
# The cmprskFrail function:           
# =====================


#@todo V2
#différente manière d'entrer la survie (surv(u~d) ou start=stop=..)
#data manquante omit
#regarder les ...
#regarder les différences de calculs likelihood
cmprskFrail <- function(formula, data, distrib="gaussian"){
  # DESCRIPTION:
  # Fits a competing risk Cox mixed-effects model with normal and gamma frailties
  #
  # Author: M De Antonio, C. Boudreau, K. Liu & H. Varet
  # Last updated/version: Jul 2017
  #
  # Developed on R version 3.0.1, with coxme version 2.2-3 (15/05/2012)
  #
  # ARGUMENTS:
  # formula = a formula object, with the response on the left of a '~'
  #            operator and the terms on the right.  The response must be a
  #            survival object as returned by the 'Surv' function. The terms
  #            on the right must include a frailty at the end
  # data = data frame of observations in which to interpret the variables
  #        named in the 'formula'
  # distrib = distribution "gaussian" or "gamma"
  # REFERENCES:
  #   Katsahian, Resche-Rigon, Chevret & Porcher (2006), "Analysing multicentre
  #         competing risks data with a mixed proportional hazards model for the
  #         subdistributionz", Statistics in Medicine, Vol. 25, 4267-4278.
  #   Ripatti & Palmgren (2000), "Estimation of multivariate frailty models using
  #         penalized partial likelihood", Biometrics, Vol. 56, 1016-1022
  #   Therneau, Grambsh & Pankratz (2003), "Penalized survival models with frailty",
  #         J Computational & Graphical Statistics, Vol. 12, 156-175
  #   Kinship package: mixed-effects Cox models, sparse matrices, and model
  #         ing data from large pedigrees, http://CRAN.R-project.org/package=kinship

  # NOTE: covariates can NOT be named weights and/or Gj and/or tmp (à revoir)
  # modifier l'entrée de la fonction ? style crr plutôt que coxph ?

  # make sure that the survival, kinship & MASS packages are loaded
  # require(survival); require(coxme); require(MASS); 050515

  # call matching     
  print(match.call())
  input=NULL
  ridge=NULL
  formula <- deparse(substitute(formula))

  # input verifications
  input$data=substitute(data)
  name_data=substitute(data)
  if (!is.data.frame(data)) stop("'", name_data, "' must be a data.frame or a list")
  input$distrib=distrib
  if (!distrib %in% c("gaussian","gamma")) stop("distrib must be 'gaussian' or 'gamma'")

  # extracting the frailty name from formula
  formula <- paste(formula, collapse=" ")
  formula <- gsub("[[:space:]]+", " ", formula)
  frailname <- strsplit(formula, "frailty")
  frailname <- gsub("[[:space:]]+|\\(|\\)", "", frailname[[1]][2])
  input$frailname= frailname
  
  #in case of ridge penalisation, extracting the formula
  if(length(strsplit(formula, "ridge")[[1]]>1)){
     ridge=strsplit(strsplit(gsub("[[:space:]]+", "", formula), "ridge")[[1]][2],")")[[1]][1]
     ridge=substr(ridge,2,nchar(ridge))
     thetapenal=strsplit(ridge,"theta=")[[1]][2]
     # extracting the variable names from formula
     varnames <- gsub(paste("\\<Surv\\>|\\(|\\)|\\<ridge\\>|\\(|\\<theta","=",thetapenal,"\\>|\\)|\\<frailty\\> ?\\(.*\\)",sep=""), "", gsub("[[:space:]]+", "", formula))
  }else{
     varnames <- gsub("\\<Surv\\>|\\(|\\)|\\<frailty\\> ?\\(.*\\)", "", formula)
  }
  varnames <- unlist(strsplit(varnames,  "[[:space:]]*(\\+|,|~)[[:space:]]*"))
  varnames <- gsub("[[:space:]]+", "", varnames[which(varnames!="")])
  varnames <- c(varnames, frailname)
  
  # extracting the covariate names from varnames
  varnamestc <- varnames[c(1:2)]
  covnames <- varnames[-c(1:2, length(varnames))]
  input$covnames= covnames
  input$varnamestc= varnamestc
  # deleting "1" from varnames if no fixed covariate in the model
  if (any(covnames=="1")){
    varnames <- varnames[-which(varnames=="1")]
  }
  # formula in input
  if(length(strsplit(formula, "ridge")[[1]])>1){
   input$formula=paste("Surv(",varnames[1],", ",varnames[2],") ~ ridge(",paste(covnames,collapse=", "),", theta=",thetapenal,") + frailty(",frailname,")",sep="")
  }else{
  input$formula=paste("Surv(",varnames[1],", ",varnames[2],") ~ ",paste(covnames,collapse=" + ")," + frailty(",frailname,")",sep="") }
 
  # keeping only the required columns of data and sorting
  data <- data[ , varnames]
  data <- data[order(data[ , varnames[1]]), ] # sorting (not required)

  # checking for missing values (this function doesn't support missing values)
  if (any(is.na(data)))
     stop("one or more missing values in '", name_data,
          "'; \n\t  this fct. doesn't support missing values")

  # creating the uniqid, fail = I(epsilon=1), cens = I(epsilon=0), etc variables
  othernames <- make.unique(c(varnames,
                "fail", "cens", "uniqid", "start", "stop", "kmcens", "exp.eta"))
  othernames <- tail(othernames, 7)
  othernames <- list(fail=othernames[1], cens=othernames[2], uniqid=othernames[3],
                     start=othernames[4], stop=othernames[5], kmcens=othernames[6], exp.eta=othernames[7])
  data[ , othernames$fail] <- ifelse(data[ , varnames[2]]==1, 1, 0)
  data[ , othernames$cens] <- ifelse(data[ , varnames[2]]==0, 1, 0)

  # computing the Kaplan-Meier of the censoring distribution
  formula.survfit=paste("Surv(",varnames[1],",",othernames$cens,")~1")
  km.cens <- survfit(as.formula(formula.survfit), data, conf.type="none")
  intervals <- summary(km.cens); tmp <- length(intervals$time)
  intervals <- cbind(c(0, intervals$time[-tmp]), intervals$time, c(1, intervals$surv[-tmp]))
  colnames(intervals) <- c(othernames$start, othernames$stop, othernames$kmcens)
  maxcause1 <- max(data[data[ , varnames[2]]==1, varnames[1]])
  if (maxcause1 > max(intervals[ , othernames$stop])){
    intervals <- rbind(intervals, c(max(intervals[ , othernames$stop]), maxcause1, min(intervals[ , othernames$kmcens])))
  } else{
    intervals <- intervals[intervals[ , othernames$start] < maxcause1, ]
    if (!is.matrix(intervals)){intervals <- as.matrix(t(intervals))}
    intervals[dim(intervals)[1], othernames$stop] <- maxcause1
  }
  intervals <- rbind(intervals, rep(intervals[dim(intervals)[1], 2:3], c(2, 1)))

  # assigning an uniqid # to each subject
  data[ , othernames$uniqid] <- 1:nrow(data)

  # creating dataset with failure from other causes
  data2 <- data[data[ , varnames[2]]>1, ]
  tmp <- summary(km.cens, data2[ , varnames[1]]*(1-.Machine$double.eps))$surv
  if(length(tmp) < dim(data2)[1]){
    tmp <- c(tmp, rep(tmp[length(tmp)], dim(data2)[1] - length(tmp)))
  }
  data2$Gj <- tmp

  # putting data in the required format for coxme
  data[ , othernames$start] <- 0
  data[ , othernames$stop] <- data[ , varnames[1]]
  data$weights <- 1
  data[ , othernames$fail] <- ifelse(data[ , varnames[2]]==1, 1, 0)

  # putting data2 in the required format for coxme
  tmp <- matrix(intervals[ , othernames$start], dim(data2)[1], dim(intervals)[1], byrow=T)
  tmp2 <- matrix(data2[ , varnames[1]], dim(data2)[1], dim(intervals)[1])
  tmp3 <- matrix(data2[ , othernames$uniqid], dim(data2)[1], dim(intervals)[1])
  data3 <- cbind(as.vector(tmp2), as.vector(tmp), as.vector(tmp3))
  colnames(data3) <- c(othernames$start, othernames$stop, othernames$uniqid)
  data3 <- data3[data3[ , othernames$start] < data3[ , othernames$stop], ]
  data3 <- merge(data3[ , -1], intervals, by=othernames$stop, sort=F)
  data2 <- merge(data2, data3, by=othernames$uniqid)
  tmp <- NULL; tmp2 <- NULL; tmp3 <- NULL; data3 <- NULL # clearing some memory
  tmp <- data2[ , varnames[1]] > data2[ , othernames$start]
  data2[ , othernames$start][tmp] <- data2[ , varnames[1]][tmp]
  data2[ , othernames$fail] <- 0

  # removing lines with start = stop
  tmp <- data2[ , othernames$start]==data2[ , othernames$stop]
  data2 <- data2[!tmp, ]

  # computing the weights
  data2$weights <- data2[ , othernames$kmcens] / data2$Gj

  # putting all the data together for coxme
  tmp <- unlist(c(varnames, othernames[1:5], "weights"), use.names=F)
  data.final <- rbind(data[ , tmp], data2[ , tmp])
  input$export_weights=data.final
  
  # output
  if (distrib=="gaussian" & length(strsplit(formula, "ridge")[[1]])==1){
    # calling coxph and coxme for normal frailties
      formula.random <- paste("+(1 |", frailname,")")
      formula.cox <- paste("Surv(",othernames$start,",", othernames$stop,",",othernames$fail,")~",paste(covnames,collapse="+"))
      fit.coxme <- coxme(as.formula(paste(formula.cox,formula.random)), data=data.final, weights=data.final[,"weights"], x=FALSE, y=FALSE)
      res=list(input=input,fit.model=fit.coxme) 
      #fit.coxph <- coxph(as.formula(formula.cox), data=data.final, weights=data.final[,"weights"], x=FALSE, y=FALSE)
    }    
  if (distrib=="gamma" | length(strsplit(formula, "ridge")[[1]])>1){
    # calling coxph for gamma frailties
    formula.cox <- paste(othernames$start, ",", othernames$stop)
    formula.cox <- sub(varnames[1], formula.cox, formula)
    formula.cox <- sub(varnames[2], othernames$fail, formula.cox)
    formula.cox <- sub("frailty\\(.*\\)", paste("frailty.",distrib,"(", frailname, ")",sep=""), formula.cox)
  
    # model for gamma frailties
    fit.gamma <- coxph(formula=as.formula(formula.cox), data=data.final, weights=data.final[,"weights"], x=FALSE, y=FALSE)
    res=list(input=input,fit.model=fit.gamma) 
  }
  attr(res, "class") <- "cmprskFrail"
  return(res)
}

##################################################################                                                                              
######################### Fonction PRINT #########################
##################################################################

print.cmprskFrail <- function(obj){
    cat("  Formula:", obj$input$formula,"\n")   
    cat("  Data:", deparse(obj$input$data),"\n")
        omit <- obj$fit.model$na.action
    if (length(omit)){
    cat("  Subject number: ",nrow(eval(obj$input$data))," (",naprint(omit)," are deleted due to missing data)","\n",sep="")
    }else{cat("  Subject number: ",length(unique(obj$input$export_weights$uniqid)),"\n",sep="")}
    ncause=length(names(table(eval(obj$input$data)[,obj$input$varnamestc[2]]!=0)))
    cat("  Censored: ",table(eval(obj$input$data)[,obj$input$varnamestc[2]])[1], paste(" Cause ",names(table(eval(obj$input$data)[,obj$input$varnamestc[2]])[-1]),": ",table(eval(obj$input$data)[,obj$input$varnamestc[2]])[-1], sep=""),"\n",sep="")
    nvar <- length(obj$input$covnames)
    cat("  Covariates number: ",nvar,"\n",sep="")
    nfrail <- length(table(obj$input$export_weights[,obj$input$frailname]))
    cat("  Frailties -number: ",nfrail,"   -distribution: ", obj$input$distrib, "\n",sep="")
    cat("  Iterations:", obj$fit.model$iter, "\n \n",sep=" ")
    
    cat("Competing risk Cox mixed-effects model with ", obj$input$distrib, " frailties \n \n",sep="")
       
    if(obj$input$distrib=="gaussian" & length(strsplit(obj$input$formula, "ridge")[[1]])==1){    
       ### output for normal frailties 
       tmp <- nfrail+1
       if (length(obj$input$covnames)==1 && obj$input$covnames=="1"){
          var.naive.coef=NULL
       } else{
          var.naive.coef <- as.matrix(obj$fit.model$var)
          var.naive.coef <- as.matrix(var.naive.coef[tmp:dim(var.naive.coef)[1], tmp:dim(var.naive.coef)[1]])
       
       dimnames(var.naive.coef) <- rep(list(names(obj$fit.model$coefficients)), 2)}
       var=diag(var.naive.coef)
       var.random <- obj$fit.model$vcoef[[1]]
    
       temp <- matrix(obj$fit.model$loglik[-3], nrow = 1) #delete fitted
       dimnames(temp) <- list("  Log-likelihood", c("NULL", "Integrated"))#, "Fitted"))
       print(temp)
       cat("\n")
       chi1 <- 2 * diff(obj$fit.model$loglik[c(1, 2)])
       chi2 <- 2 * diff(obj$fit.model$loglik[c(1, 3)])      
       temp <- rbind(c(round(chi1, 2), round(obj$fit.model$df[1], 2), signif(1 - 
       pchisq(chi1, obj$fit.model$df[1]), 5), round(chi1 - 2 * obj$fit.model$df[1], 
       2), round(chi1 - log(obj$fit.model$n[1]) * obj$fit.model$df[1], 2)), c(round(chi2, 
       2), round(obj$fit.model$df[2], 2), signif(1 - pchisq(chi2, obj$fit.model$df[2]), 
       5), round(chi2 - 2 * obj$fit.model$df[2], 2), round(chi2 - log(obj$fit.model$n[1]) * 
       obj$fit.model$df[2], 2)))
       dimnames(temp) <- list(c("  Integrated loglik", "  Penalized loglik"), c("Chisq", "df", "p", "AIC", "BIC"))
       print(temp, quote = F) #, digits = digits) 
       
       if (nvar > 0) {
           beta <- obj$fit.model$coefficients
           se <- sqrt(diag(obj$fit.model$var)[nfrail + 1:nvar])
           tmp <- cbind(beta, exp(beta), se, round(beta/se, 2), 
           signif(1 - pchisq((beta/se)^2, 1), 2))
           dimnames(tmp) <- list(paste(" ",names(beta)," "),c("coef", "exp(coef)", "se(coef)", "z", "p"))
        cat("\n  Fixed coefficients\n")
        print(tmp)} #, digits = digits) 
        
        cat("\n  Random effects\n")
        random <- VarCorr(obj$fit.model)
        nrow <- sapply(random, function(x) if (is.matrix(x))
             nrow(x)
        else length(x))
    maxcol <- max(sapply(random, function(x) if (is.matrix(x)) 1 + ncol(x) else 2))
    temp1 <- matrix(NA, nrow = sum(nrow), ncol = maxcol)
    indx <- 0        
    for (term in random) {
        if (is.matrix(term)) {
            k <- nrow(term)
            nc <- ncol(term)
            for (j in 1:k) {
                temp1[j + indx, 1] <- sqrt(term[j, j])
                temp1[j + indx, 2] <- term[j, j]
                if (nc > j) {
                  indx2 <- (j + 1):nc
                  temp1[j + indx, 1 + indx2] <- term[j, indx2]
                }
            }
        }
        else {
            k <- length(term)
            temp1[1:k + indx, 1] <- sqrt(term)
            temp1[1:k + indx, 2] <- term
        }
        indx <- indx + k
    }   
    indx <- cumsum(c(1, nrow))
    temp3 <- rep("", nrow(temp1))
    temp3[indx[-length(indx)]] <- names(random)
    xname <- unlist(lapply(random, function(x) if (is.matrix(x)) 
        dimnames(x)[[1]]
    else names(x)))
    temp <- cbind("",temp3, xname, ifelse(is.na(temp1), "", format(temp1, digits = 6)))    #format(temp1, digits = digits)  
    if (maxcol == 2) { 
        temp4 <- c("","Group", "Variable", "Std Dev", "Variance")  }
    else {temp4 <- c("","Group", "Variable", "Std Dev", "Variance", "Corr", rep("", maxcol - 3))  }
    dimnames(temp) <- list(rep("", nrow(temp)), temp4)   
    print(temp, quote = F)
    invisible(obj$fit.model)
    cat("\n")}
### output for gamma frailties
    if(obj$input$distrib=="gamma" | length(strsplit(obj$input$formula, "ridge")[[1]])>1){
        temp <- matrix(obj$fit.model$loglik, nrow = 1)
        dimnames(temp) <- list("  Log-likelihood", c("NULL", "Integrated"))
        print(temp)
        cat("\n")
        chi1 <- 2 * diff(obj$fit.model$loglik[c(1, 2)])
        temp <- rbind(c(round(chi1, 2), round(obj$fit.model$df[1], 2), round(signif(1 -pchisq(chi1, obj$fit.model$df[1]), 5),6), 
                      round(chi1 - 2 * obj$fit.model$df[1], 2), round(chi1 - log(obj$fit.model$n[1]) * obj$fit.model$df[1], 2)))
        dimnames(temp) <- list(c("  Integrated loglik"), c("Chisq", "df", "p", "AIC", "BIC"))
        print(temp, quote = F) #, digits = digits)
        if (nvar > 0) {
            beta <- obj$fit.model$coefficients
            var <- diag(obj$fit.model$var)
            se <- round(sqrt(var),6)
            if (is.null(obj$fit.model$naive.var)) {
               tmp <- cbind(beta, exp(beta), se, round(beta/se,2), round(signif(1 - pchisq((beta/se)^2, 1), 2),6))
               dimnames(tmp) <- list(paste(" ",names(beta)," "), c("coef", "exp(coef)", "se(coef)", "z", "p"))
            } else {
               nse <- sqrt(diag(obj$fit.model$naive.var))
               tmp <- cbind(beta, exp(beta), nse, se, round(beta/se,2), round(signif(1 - pchisq((beta/se)^2, 1), 2),6))
               dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", "se(coef)", "robust se", "z", "p"))
            }

            cat("\n  Fixed coefficients\n")
            print(tmp[1:nvar,])} #, digits = digits)

            cat("\n  Random effects\n")
            theta <- obj$fit.model$history[[1]]$theta
            temp <- cbind("",obj$input$frailname,"", round(sqrt(theta),6),round(theta,6))
            dimnames(temp) <- list("",c("","Group", "Variable", "Std Dev", "Variance"))
            print(temp, quote = F)
            cat("\n")
}}




##################################################################
######################### Fonction PLOT ##########################
##################################################################

plot.cmprskFrail <- function(obj,addCR=F,export=F,plot_name=NULL,xlab = obj$input$varnamestc[1],ylab=NULL,ylim=c(0,1.3),main=NULL,conf.type="none",lwd = c(2,1),lty=c(1,2), col=NULL,cex.lab=1,cex.axis=1,leg.pos=NULL,leg.txt=NULL,bty="l",mar=c(5, 4.5, 4, 2) + 0.1,...){
    database=eval(obj$input$data)
    ncause=length(names(table(database[,obj$input$varnamestc[2]]!=0)))
    nfrail=length(table(obj$input$export_weights[,obj$input$frailname]))
    #Default parameters varying among the addCR option
    if(is.null(col) && addCR==F) col=c(2,1)
    if(is.null(col) && addCR==T) col=c(1:ncause)
    if(is.null(plot_name) && addCR==F) plot_name="centreeffect_glob"
    if(is.null(plot_name) && addCR==T) plot_name="centreeffect_CR"
    if(is.null(main) && addCR==F) main="Global survival with frailties"
    if(is.null(main) && addCR==T) main="Cumulative rate"
    if(is.null(ylab) && addCR==F) ylab=expression(paste("1 -",widehat(S), "(t)"))
    if(is.null(ylab) && addCR==T) ylab="Cumulative incidence"
    if(is.null(leg.pos) && addCR==F) leg.pos="bottomright"
    if(is.null(leg.pos) && addCR==T) leg.pos="topright"
    if(is.null(leg.txt) && addCR==F) leg.txt=c("All causes","Mean curve","Centers curves")
    if(is.null(leg.txt) && addCR==T) leg.txt=c("Cause ","Mean curve C","Centers curves C")
    #Dot Parameters
    dotnames <- names(list(...))
    if (any(dotnames == "type")) stop("The graphical argument 'type' is not allowed")
    if (export==T) png(paste(plot_name,file=".png",sep="")) #,width = 800, height = 800)
    #Graphical Parameters
    par(bty=bty,mar=mar)
    if(addCR==F){  
        kmg_centre <- survfit(Surv(database[,obj$input$varnamestc[1]], database[,obj$input$varnamestc[2]]!=0) ~ database[,obj$input$frailname],                       conf.type=conf.type)
        plot(kmg_centre, xlab=xlab, ylab=ylab,lwd=lwd[2],lty=lty[2],col=col[2],fun="event",main=main)
        kmg_mean <- survfit(Surv(database[,obj$input$varnamestc[1]], database[,obj$input$varnamestc[2]]!=0)~1, conf.type=conf.type)
        lines(kmg_mean, lwd=lwd[1],lty=lty[1], col=col[1], fun="event")
        legend(leg.pos,legend=leg.txt,lty=c(NA,lty),col=c(NA,col),box.lty=0,lwd=c(NA,lwd))
    }
    if(addCR==T){
        plot(cuminc(database[,obj$input$varnamestc[1]],database[,obj$input$varnamestc[2]],database[,obj$input$frailname]),
        main=main,lty=lty[2], col=rep(1:ncause,each=nfrail), ylab=ylab,lwd=lwd[2],cex.lab=cex.lab,
        cex.axis=cex.axis,ylim=ylim,wh=c(nfrail,nfrail)) 
        cuminctemp <- cuminc(database[,obj$input$varnamestc[1]],database[,obj$input$varnamestc[2]])
        for (i in (1:ncause)){
             lines(cuminctemp[[i]]$time,cuminctemp[[i]]$est,col=col[i],lwd=lwd[1])
             }
        legend(leg.pos,legend=paste(leg.txt,rep(1:ncause,each=3),sep=""),lty=rep(c(NA,lty),ncause),col=rep(col,each=ncause+1),
        lwd=rep(c(NA,lwd),ncause),box.lty=0)
    }
    if (export==T) dev.off()
 }
 