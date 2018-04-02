###################
###Main Function###
###################
isoVCT <- function(method = 'Emp', nper = 5000, 
                   y, p, label){
  
  ##method: character, statistics distribution, default setting is 'Emp'; 
  ##nper: numeric, number of permutation, default setting is '1000';
  ##y: vector, expression level of isoforms;
  ##p: numeric, numeber of isoforms;
  ##label: vector, label fo samples.
  
  ##Fit NB model
  ranGroup <- rep(c(1: p), length(label))
  fitdata <- data.frame(y=y, ranGroup=ranGroup)
  fitNB <- try(glmer.nb(y ~ (1|ranGroup), fitdata), silent=T)

  if(method == 'Emp'){
    
    if(inherits(fitNB, 'try-error')){
      
      pval <- UTestEmpPois(nper=nper, y=y, p=p, label=label)
      
    }else{
      
      if(getME(fitNB, 'theta') !=0 & getME(fitNB, 'glmer.nb.theta') >= 0.01){
        
         pval <- UTestEmpNB(nper=nper, y=y, p=p, label=label)
        
      }else{
        
        pval <- UTestEmpPois(nper=nper, y=y, p=p, label=label)
        
      }
      
    }
    
  }else{
      
    if(inherits(fitNB, 'try-error')){
      
      pval <- UTestThePois(y=y, p=p, label=label)
      
    }else{
      
      if(getME(fitNB, 'theta') != 0 & getME(fitNB, 'glmer.nb.theta') >= 0.01){
        
        pval <- UTestTheNB(y=y, p=p, label=label)
        
      }else{
        
        pval <- UTestThePois(y=y, p=p, label=label)
        
      }
      
    }
  
  }   
  
  return(pval)
  
}


###Trace
tr <- function(mat){
  
  trace <- sum(diag(mat))
  return(trace)
  
}



##############
###NB Assmp###
##############
###Estimation of Statistics U
UEstNB <- function(label, p, fit, dp, gamma){
  
  require(plyr)
  X <- as.matrix(adply(label, 1, function(a){a*diag(1, nrow=p)})[,-1])
  residVec <- resid(fit, type='response')
  phi <- 1/(1 + dp*gamma)
  xMat <- (residVec * phi) %*% X
  W <- gamma/(1 + gamma*dp)
  e <- (gamma*dp)/(1+gamma*dp)^2
  W0 <- W + e*residVec 
  U <- 0.5*(xMat%*%t(xMat) - tr(W0*(X%*%t(X))))
  return(list(U, X))
  
}


###Empirical distribution
##Permutation
UPerNB <- function(label, p, fit, dp, gamma){
  
  index <- sample(label)
  X <- as.matrix(adply(index, 1, function(a){a*diag(1, nrow=p)})[,-1])
  residVec <- resid(fit, type='response')
  phi <- 1/(1 + dp*gamma)
  xMat <-  (residVec * phi) %*% X
  W <- gamma/(1 + gamma*dp)
  e <- (gamma*d)/(1+gamma*dp)^2
  W0 <- W + e*residVec 
  U <- 0.5*(xMat%*%t(xMat)-tr(W0*X%*%t(X)))
  return(U)
  
}

##Test
UTestEmpNB <- function(nper, y, p, label){
  
  require(lme4)
  require(plyr)
  #Fit model
  ranGroup <- rep(c(1: p), length(label))
  fitdata <- data.frame(y = y, ranGroup = ranGroup)
  fitNB <- glmer.nb(y ~ (1|ranGroup), data = fitdata)
  
  if(getME(fitNB, 'theta')==0){
    
    stop('This situation is perhaps unsuitable to mixed model!')
    
  }else{  
    
    gamma <- getME(fitNB, 'mu')
    dp <- 1/getME(fitNB, 'glmer.nb.theta')
    U <- UEstNB(label=label, p=p, fit=fitNB, dp=dp, gamma=gamma)[[1]]
    
    #Permutation and test
    UDis <- replicate(nper, UPerNB(label=label, p=p, fit=fitNB, dp=dp, gamma=gamma))
    pval <- sum(as.numeric(U) <= UDis)/nper
    
    return(as.numeric(pval))
    
  }
  
}



###Theoritical distribution
UTestTheNB <-  function(y, p, label){
  
  require(plyr)
  
  ##Fit model
  ranGroup <- rep(c(1: p), length(label))
  fitdata <- data.frame(y=y, ranGroup=ranGroup)
  fitNB <- glmer.nb(y ~ (1|ranGroup), fitdata)
  
  if(getME(fitNB, 'theta')==0){
    
    stop('This situation is perhaps unsuitable to mixed model!')
    
  }else{ 
    
    gamma <- getME(fitNB, 'mu')[1:p]
    dp <- 1/getME(fitNB, 'glmer.nb.theta')
    
    ##Statistics U
    UEst <- UEstNB(label=label, p=p, fit=fitNB, dp=dp, gamma=gamma)
    U <- UEst[[1]]
    X <- UEst[[2]]
    
    ##Infomation Matrix
    rij <- 2*(gamma/(1+dp*gamma))^2
    rii <- (3*gamma^3*dp^2 + 2*gamma^3*dp + 2*gamma^2*dp + 2*gamma^2 + gamma)/(1+gamma*dp)^3
    
    I.tautau <- 0.25*((sampleSize*p/2)/2*sum(rii)+
                        (sampleSize*p/2)*((sampleSize*p/2)-2)/4*sum(rij))
    I.alphatau <- 0.5*((sampleSize*p/2)*sum(c((gamma+gamma^2*dp)/((1+gamma*dp)^2))))
    I.alphaalpha <- (sampleSize*p/2)*sum(gamma/(1+dp*gamma))
    I.eff <- I.tautau - I.alphatau^2/I.alphaalpha
    
    pval <- 1-as.numeric(pchisq(U^2/I.eff,df=1))
    
    return(pval)
    
  }
  
}




###################
###Poisson Assmp###
###################
###Estimation of Statistics U
UEstPois <- function(label, p, fit, W){
  
  require(plyr)
  X <- as.matrix(adply(label, 1, function(a){a*diag(1, nrow=p)})[,-1])
  xMat <- resid(fit, type='response') %*% X
  W0 <- W
  U <- 0.5*(xMat%*%t(xMat) - tr(W0*X%*%t(X)))
  return(list(as.numeric(U), X))

}


###Empirical distribution
##Permutation
UPerPois <- function(label, p, fit, W){

  index <- sample(label)
  X <- as.matrix(adply(index, 1, function(a){a*diag(1, nrow=p)})[,-1])
  xMat <- resid(fit, type='response') %*% X
  W0 <- W
  U <- 0.5*(xMat%*%t(xMat)-tr(W0*X%*%t(X)))
  return(U)

}

##Test
UTestEmpPois <- function(nper, y, p, label){

  require(lme4)
  #Fit model
  ranGroup <- rep(c(1: p), length(label))
  fitPois <- try(glmer(y ~ (1|ranGroup), family = 'poisson'), silent=T)
  
  if(getME(fitPois, 'theta')==0 | inherits(fitPois, 'try-error')){
 
    stop('This situation is perhaps unsuitable to mixed model!')
  
  }else{
    
    W <- getME(fitPois, 'mu')
    U <- UEstPois(label=label, p=p, fit=fitPois, W=W)[[1]]
  
    #Permutation and test
    UDis <- replicate(nper, UPerPois(label=label, p=p, fit=fitPois, W=diag(W)))
    pval <- sum(as.numeric(U) <= UDis)/nper
    return(as.numeric(pval))
  
  }

}



###Theoritical distribution
UTestThePois <- function(y, label, p){
  
  require(plyr)
  require(lme4)
  require(MASS)
  
  #Fit model
  ranGroup <- rep(c(1: p), length(label))
  sampleSize <- length(label)
  fitPois <- try(glmer(y ~ (1|ranGroup), family = 'poisson'), silent=T)
  
  if(getME(fitPois, 'theta')==0 | inherits(fitPois, 'try-error')){
 
    stop('This situation is perhaps unsuitable to mixed model!')
  
  }else{
    
    ##Statistic U
    W <- getME(fitPois, 'mu')
    UEst <- UEstPois(label=label, p=p, fit=fitPois, W=W)
    U <- UEst[[1]]
    X <- UEst[[2]]

    ##Information Matrix I
    gamma <- W[1:p]
    rij <- 2*gamma*gamma
    rii <- gamma + 2*gamma^2
    I.tautau <- 0.25*((sampleSize*p/2)/2*sum(rii)+
                      (sampleSize*p/2)*((sampleSize*p/2)-2)/4*sum(rij))
   
    C <- diag(rep(gamma, length(label)))
    K <- do.call('rbind', rlply(length(label), cbind(1, diag(p))))
    a <- diag(X%*%t(X))
    I.alphatau <- 0.5*(t(K)%*%C%*%a)
    I.alphaalpha <- t(K)%*%C%*%K
    
    I.eff <- I.tautau - t(I.alphatau)%*%ginv(I.alphaalpha)%*%I.alphatau
    pval <- 1 - pchisq(U*(1/as.numeric(I.eff))*U, df=1)
   
    return(as.numeric(pval))
    
  }

}

