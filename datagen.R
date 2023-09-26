datagen<-function(n,CL,sdu){
  y <-0
  yx<-0
  delta <- 0                        #Censorship indicator (0 if data is censored and 1 otherwise)
  c     <- matrix(999,n,1)          #Censoring variable independent-identical distributed with Y
  
  x<--2.5*runif(n,0,3)              #without error
  x<-sort(x)                        #sorted values
  u<-stats::rnorm(n,mean=0,sd=sdu)  #U errors (EIV)
  e<-stats::rnorm(n,mean=0,sd=0.5)  #\varepsilon error (model)
  W<-x+u                            #W (X+u)
  m<-W*sin(W)+20                       #function with W
  mx<-x*sin(x)+20                      #function without error
  z<-m+e                            #Response with EIV
  zx<-mx+e                          #Response without EIV
  #Censoring Part for EIV---------------------------------------------------------
  delta <-1-rbern(n,CL)             #Censorship indicator (0):censored, (1):observed
  for (i in 1:n){
    if (delta[i]==0){
      while (z[i]<=c[i]){
        c[i]<-rnorm(1,mean(z),sd=sd(z))
      }
    }
    else{
      c[i]<-z[i]
    }
  }                           #generating censoring variable
  for (j in 1:n){
    if (z[j]<=c[j]){
      y[j]<-z[j]
    }
    else{
      y[j]<-c[j]
    }
  }          
  #Censoring Part for non-EIV-----------------------------------------------------
  library(Rlab)
  delta <-1-rbern(n,CL)      #Censorship indicator (0):censored, (1):observed
  for (i in 1:n){
    if (delta[i]==0){
      while (zx[i]<=c[i]){
        c[i]<-rnorm(1,mean(zx),sd=sd(zx))
      }
    }
    else{
      c[i]<-zx[i]
    }
  }        #generating censoring variable
  for (j in 1:n){
    if (zx[j]<=c[j]){
      yx[j]<-zx[j]
    }
    else{
      yx[j]<-c[j]
    }
  }    
  #-------------------------------------------------------------------------------
  bjNP <- function(x,y,delta){
    obj   <- smooth.spline(x,y)
    fit   <- obj$y  
    resid <- y-fit
    cond  <- mean(resid^2)#
    #-------------------------------------------------------------------------------
    
    ry <- y
    KME<-function(y,delta){
      library(pracma)
      # where y: right-censored observations, delta: censorship indicator
      n<-length(y)
      M<-0
      yg<-0
      M<-0
      delta2<-0
      #Synthetic data transformation
      y1<-cbind(y,delta)
      y2<-sortrows(y1,k=2)
      delta2<-y2[,2]
      delta2[n]<-1
      sy<-sort(y)
      for (i1 in 1:n){
        Ma<-1
        for (i2 in 1:n){
          mGt<-((n-i2)/(n-i2+1))^((delta2[i2]==0)*(sy[i2]<=sy[i1]))
          Ma<-Ma*mGt
        }
        M[i1]=Ma
        yg[i1]<-(delta[i1]*y[i1])/M[i1]
      }
      return(M)
    }
    
    M <- KME(y,delta)
    a <- 1
    Mr <- 0
    for (i in 1:n){
      if (delta[i]==0){
        Mr[i] <-M[a]
        a <- a+1
      }
      else{
        Mr[i] <- 1
      }
    }
    ybj <- y/Mr
    tol <- 0
    count <- 1
    while (tol<cond){
      obj2   <- smooth.spline(x,ybj)
      fit2   <- obj2$y  
      ec <- (y-fit2)*(1-delta)
      for (j in 1:n){
        ybj[j] <- y[j]*delta[j]+(y[j]+abs(ec[j]))*(1-delta[j]) 
      }
      resid2 <- ybj-fit2
      tol <-mean(resid2^2) 
      #plot(ybj,type="l",ylim=c(min(y),max(y)))
      #par(new=TRUE)
      #plot(y,type="l",col="red",lwd=1)
      #par(new=TRUE)
      count <- count+1
      if (count>10){
        break
      }
    }
    lm3 <- smooth.spline(x,ybj)
    fit3 <- lm3$y
    for (l in 1:n){
      if (delta[l]==0){
        ybj[l] <- fit3[l]
      }
      else{
        ybj[l] <- y[l]
      }
    }
    #plot(ybj,type="l",ylim=c(min(y),max(y)))
    #par(new=TRUE)
    #plot(y,type="l",col="red",lwd=1)
    return(ybj)
  }
  ystar<-bjNP(x,y,delta)
  #-------------------------------------------------------------------------------
  res<-new.env()
  res$z  <-z     #Completely observed EIV response variable
  res$zx <-zx    #Completely observed non-EIV response variable
  res$x  <-x     #non-EIV covariate 
  res$y  <-y     #right-censored EIV response variable
  res$W  <-W     #EIV covariate 
  res$m  <-m     #EIV function
  res$mx <-mx    #non-EIV function
  res$yx <-yx    #non-EIV right-censored response variable
  res$u  <-u     #error of x (W=x+u)
  res$e  <-e     #random error terms of the model (y=m(x)+e)
  res$yg<-ystar  #Synthetic response data
  res$delta <- delta
  return(res)
} 