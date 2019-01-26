library(mvtnorm)
library(geoR)
library(matrixStats)
library(doParallel)
library(foreach)
library(fields)
registerDoParallel(20)


setwd("/mnt/home/aroy2")

# Transform that to matern parameters
theta2mat <- function(theta,nugget=TRUE){
  c(exp(theta[1]),
    ifelse(nugget,1/(1+exp(-theta[2])),1),
    exp(theta[3]),
    exp(theta[4]))
}

# Transform matern parameters to theta
mat2theta <- function(mat,nugget=TRUE){
  c(log(mat[1]),
    ifelse(nugget,log(mat[2]/(1-mat[2])),Inf),
    log(mat[3]),
    log(mat[4]))
}

simulateMatern <- function(theta, d){
  thetan <- rep(0, 4)
  temp <- theta2mat(theta)
  thetan[1] <- temp[1]*(1-temp[2])
  thetan[2] <- temp[1]*temp[2]
  thetan[3:4] <- temp[3:4]
  Var <- thetan[2] * matern(d, thetan[3], thetan[4])
  diag(Var) <- diag(Var) + thetan[1]
  ret <- rmvnorm(1, sigma = Var)
  return(ret)
}

Maternvar <- function(theta, d){
  thetan <- rep(0, 4)
  temp <- theta2mat(theta)
  thetan[1] <- temp[1] * (1 -temp[2])
  thetan[2] <- temp[2] * temp[1]
  thetan[3:4] <- temp[3:4]
  Var <- thetan[2] * matern(d, thetan[3], thetan[4])
  diag(Var) <- diag(Var) + thetan[1]
  return(Var)
}

updateBeta <- function(l, y, IVarc, Beta, IVare, X, B =  NULL){
  if(is.null(B)){B <- rep(1, n)}
  mean <- apply(X, 3, function(x){rowSums(Beta[, c(2:11)[-l]]*x[, -l])})
  Beta.mean <- rowSums(sapply(1:m, function(k){diag((X[, l,k])*(B))%*%IVare%*%(y[, k]-mean[, k])})) #rowSums((y - mean)*X[, l,])*(B)
  Beta.ivar <- lapply(1:m, function(k){diag((X[,l,k])*B)%*%IVare%*%diag((X[,l,k])*B)})
  Beta.ivar <- Reduce('+', Beta.ivar)+ IVarc
  Beta.var <- solve(Beta.ivar)
  Beta.var <- (Beta.var + t(Beta.var))/2
  Beta.mean <- Beta.var %*% Beta.mean
  gen <- rmvnorm(1, Beta.mean, Beta.var)
  return(gen)
}

ADJ1D <- function(m){
  A <- 0*diag(m)
  for(j in 2:m){
    A[j,j-1]<-A[j-1,j]<-1
  }
  A}

conprb <- function(a, b, r){
  J <- length(a)
  t1 <- min((r*a) / b)
  t2 <- max(a / (r * b))
  p_conditional <- ((r / (r^2 - 1))^J) * (prod(1/a)) * (t1 ^ J - t2^J) / J 
  return(p_conditional)
}

q=1
foreach(repli = 1:10) %dopar% {
  foreach(vind = 1:3)  %dopar% {
    
    vx <- c(3, 6)
    var <- c(0.1, 1, 1.5)
    m        <- 30
    n1 <- n2 <- m
    p        <- round(m/2)
    m2       <- m^2
    s        <- expand.grid(1:m,1:m)
    n <- 100
    set.seed(8)
    A1 <- rep(1:n2, each = n1)
    A2 <- rep(1:n1, n2)
    tempA <- cbind(A2, A1)
    #Ap <- matrix(rep(array(t(tempA)), n3), ncol=2, byrow = T)
    Ap <- tempA #cbind(Ap, rep(1:n3, each=n1*n2))
    
    
    RE      <- TRUE  # Generate data with random effects?
    
    
    pri.mn=c(0,0,0,0)
    pri.sd=c(10,2,10,1)
    
    L=1
    MHY=.01
    
    X <- matrix(0, n, n1*n2)
    
    loc <- Ap
    dis <- as.matrix(dist(loc))
    nux <- 3 #vx[nuind]
    xvar <- Exponential(dis, range = nux)
    
    for(i in 1:n){
      X[i, ] <- rmvnorm(1,sigma = xvar) 
    }
    
    h <- 2 #round(runif(1,0,2)) + 1
    u <- matrix(runif(h*2), nrow=2)
    
    # d <- .4*exp(-5*rowSums((Ap-matrix(c(n1*u[1,1], n2*u[1,2]), nrow = nrow(Ap), 2, byrow=T) )^2)/50)
    if(h == 2){
      d <- 1*exp(-5*rowSums((Ap-matrix(c(n1*u[1,1], n2*u[1,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50) + 1*exp(-5*rowSums((Ap-matrix(c(n1*u[2,1], n2*u[2,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50)
    }
    # if(h == 3){
    #   d <- .4*exp(-5*rowSums((Ap-matrix(c(n1*u[1,1], n2*u[1,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50) + 
    #     .4*exp(-5*rowSums((Ap-matrix(c(n1*u[2,1], n2*u[2,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50) + .4*exp(-5*rowSums((Ap-matrix(c(n1*u[3,1], n2*u[3,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50)
    # }
    
    h <- 5 #round(runif(1,0,2)) + 1
    u <- matrix(c(.2,.8,.2,.8,.5,.8,.2,.2,.8,.5), nrow=5)
    
    # d <- .4*exp(-5*rowSums((Ap-matrix(c(n1*u[1,1], n2*u[1,2]), nrow = nrow(Ap), 2, byrow=T) )^2)/50)
    # if(h == 2){
    #   d <- .4*exp(-5*rowSums((Ap-matrix(c(n1*u[1,1], n2*u[1,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50) + .4*exp(-5*rowSums((Ap-matrix(c(n1*u[2,1], n2*u[2,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50)
    # }
    # if(h == 3){
    #   d <- .4*exp(-5*rowSums((Ap-matrix(c(n1*u[1,1], n2*u[1,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50) +
    #     .4*exp(-5*rowSums((Ap-matrix(c(n1*u[2,1], n2*u[2,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50) + .4*exp(-5*rowSums((Ap-matrix(c(n1*u[3,1], n2*u[3,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50)
    # }
    d <- 0
    for(i in 1:5){
      d <- d + 2*exp(-20*rowSums((Ap-matrix(c(n1*u[i,1], n2*u[i,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50)
    }
    
    B0  <- d
    B0[which(B0<1e-1)]<- 0 
    sigma0 <- var[vind]
    set.seed(repli)
    
    Y <- rnorm(n, mean = X%*%B0, sd=sigma0)
    
    knots <- expand.grid(seq(0,m+1,length=p),
                         seq(0,m+1,length=p))
    bw    <- min(dist(knots))
    
    A     <- ADJ1D(p)
    ADJ   <- kronecker(A,diag(p))+kronecker(diag(p),A)
    
    D     <- rdist(s,knots)
    B     <- exp(-0.5*(D/bw)^2)
    B     <- ifelse(D>3*bw,0,B)
    
    S     <- sqrt(diag(B%*%solve(diag(colSums(ADJ))-0.99*ADJ)%*%t(B)))    
    B     <- diag(1/S)%*%B
    S     <- sqrt(diag(B%*%solve(diag(colSums(ADJ))-0.99*ADJ)%*%t(B)))   
    
    init.theta= c(1,.9)#c(0.3216614,  0.1695090, -1.2245358,  1.0982590)#c(0,2,2,0)
    theta <- init.theta
    
    sigma = 1
    
    M <- diag(rowSums(ADJ))
    
    isdmat <- (1/theta[1])*(M - theta[2]*ADJ)
    
    Xt <- X %*% B
    Beta.mean <- t(Xt) %*% Y / sigma^2
    Beta.ivar <- t(Xt) %*% Xt / sigma^2 + isdmat
    Beta.var <- solve(Beta.ivar)
    Beta.var <- (Beta.var + t(Beta.var))/2
    Beta.mean <- Beta.var %*% Beta.mean
    temp <- array(rmvnorm(1, Beta.mean, Beta.var))
    
    # psd <- isdmat    #Maternvar(thetaA, dis)
    # BB       <- (theta[1])*(t(temp)%*%(psd)%*%(temp))/2+.1
    # 
    # theta[1]    <- 1/(rgamma(1,(n1*n2)/2+.1,BB))
    # 
    # isdmat <- (1/theta[1])*(M - theta[2]*ADJ)
    # 
    # Beta.ivar <- t(Xt) %*% Xt / sigma^2 + isdmat
    # Beta.var <- solve(Beta.ivar)
    # Beta.var <- (Beta.var + t(Beta.var))/2
    # Beta.mean <- Beta.var %*% Beta.mean
    # temp <- array(rmvnorm(1, Beta.mean, Beta.var))
    # 
    # psd <- isdmat    #Maternvar(thetaA, dis)
    # BB       <- (theta[1])*(t(temp)%*%(psd)%*%(temp))/2+.1
    # 
    # theta[1]    <- 1/(rgamma(1,(n1*n2)/2+.1,BB))
    # 
    # isdmat <- (1/theta[1])*(M - theta[2]*ADJ)
    
    temp <- B %*% array(temp)
    temp1 <- sign(temp)*(abs(temp))^(1/q)
    Betac <- matrix(rep(temp1, q), ncol = q)
    Beta <- array(temp)
    coef <- matrix(0, ncol=q, nrow = ncol(B))
    
    
    Total_itr <- 2000
    itr = 0
    theta_p <- list()
    Beta_p <- list()
    sigma_p <- rep(0, Total_itr)
    acceptedthetano_p <- rep(0, Total_itr)
    
    tol=0.000001
    
    sdl <- 1e1
    alpha0 <- 0.1
    beta0 <- 0.1
    
    while(itr < Total_itr){
      itr <- itr + 1
      al <- alpha0 + n / 2
      be <- beta0 + sum((Y-X%*%Beta) ^ 2) / 2
      sigma <- sqrt(1 / rgamma(1, al, be))
      sigma_p[itr] <- sigma
      
      if(q>1){
        for(k in 1:q){
          Ba <- Betac[, -k]
          if(q > 2){
            Ba <- rowProds(Betac[, -k])
          }
          ivar <- isdmat*(k==1) + isdmat*(theta[1])*(k>1)
          X1 <- X*matrix(Ba, n, n1*n2, byrow = T)
          X1t <- X1 %*% B
          Beta.mean <- t(X1t) %*% Y / sigma^2
          Beta.ivar <- t(X1t) %*% X1t / sigma^2 + ivar
          Beta.var <- solve(Beta.ivar)
          Beta.var <- (Beta.var + t(Beta.var))/2
          Beta.mean <- Beta.var %*% Beta.mean
          coef[, k] <- array(rmvnorm(1, Beta.mean, Beta.var))
          Betac[, k] <- B %*% coef[, k] 
        }
        temp <- array(rowProds(Betac))
        Beta <- temp
        
        u <- rbeta(1, sdl, sdl)
        thetaA <- theta
        cant <- rep(1, 2)
        cant[2]     <- thetaA[2] * u / (thetaA[2] * u+ (1 - thetaA[2]) * (1-u))
        cansd    <- M - cant[2] * ADJ
        psd <- isdmat    #Maternvar(thetaA, dis)
        y   <-     coef
        bb       <- (t(y[, 1])%*%(cansd)%*%(y[, 1]))/2+.1
        cant1    <- 1/(rgamma(1,(n1*n2)/2+.1,bb))
        cant[1]  <- cant1
        cansd    <- cansd/(cant1)
        BB       <- (thetaA[1])*(t(y[, 1])%*%(psd)%*%(y[, 1]))/2+.1
        
        term1    <- t(y[, 2])%*%psd%*%y[, 2]*(thetaA[1]) / 2
        
        if(q > 2){
          term1    <- sum(apply(y[, 2:q], 2, function(x){t(x)%*%psd%*%x*(thetaA[1])}))/2
        }
        
        term2    <- t(y[, 2])%*%cansd%*%y[, 2]*(cant[1]) / 2
        
        if(q > 2){
          term2    <- sum(apply(y[, 2:q], 2, function(x){t(x)%*%cansd%*%x*(cant[1])}))/2
        }
        a <- theta[2]
        b <- cant[2]
        curll    <- 0.5*as.numeric(determinant(psd)$modulus) + 0.5*(q-1)*as.numeric(determinant(psd*(thetaA[1]))$modulus)-
          (t(y[, 1])%*%(psd)%*%(y[, 1]))/2 - term1 +
          dgamma(1/(thetaA[1]),.1,.1,log=TRUE) + dbeta(thetaA[2], 9, 1, log = T)
        canll    <- 0.5*as.numeric(determinant(cansd)$modulus) + 0.5*(q-1)*as.numeric(determinant(cansd*(cant1))$modulus)-
          (t(y[, 1])%*%(cansd)%*%(y[, 1]))/2 - term2 +
          dgamma(1/(cant[1]),.1,.1,log=TRUE) + + dbeta(cant[2], 9, 1, log = T)
        Q1       <- dgamma(1/(thetaA[1]),(n1*n2)/2+.1,BB,log=TRUE)
        Q2       <- dgamma(1/(cant[1]),(n1*n2)/2+.1,bb,log=TRUE)
        R        <- exp(canll-curll+Q1-Q2)
        
        p_con <- dbeta((a*(1-b))/(a+b-2*a*b), sdl, sdl) / dbeta((b*(1-a))/(a+b-2*a*b),  sdl, sdl)
        p_coditional <- (b*(1 - b)) / (a*(1 - a))
        R <- R*p_coditional * p_con
        if(!is.na(R)){if(log(runif(1))< R){
          acceptedthetano_p[itr] <- 1
          theta  <- cant
          isdmat  <- cansd
        }}
        
      }
      
      if(q==1){
        Beta.mean <- t(Xt) %*% Y / sigma^2
        Beta.ivar <- t(Xt) %*% Xt / sigma^2 + isdmat
        Beta.var <- solve(Beta.ivar)
        Beta.var <- (Beta.var + t(Beta.var))/2
        Beta.mean <- Beta.var %*% Beta.mean
        coef <- array(rmvnorm(1, Beta.mean, Beta.var))
        Betac <- B %*% coef 
        Beta <- Betac
        
        u <- rbeta(1, sdl, sdl)
        thetaA <- theta
        cant <- rep(1, 2)
        cant[2]     <- theta[2] * u / (theta[2] * u+ (1 - theta[2]) * (1-u))
        cansd    <- M - cant[2] * ADJ
        psd <- isdmat    #Maternvar(thetaA, dis)
        y   <-    coef
        bb       <- (t(y)%*%(cansd)%*%(y))/2+.1
        cant1    <- 1/(rgamma(1,(n1*n2)/2+.1,bb))
        cant[1]  <- cant1
        cansd    <- cansd/(cant1)
        BB       <- (theta[1])*(t(y)%*%(psd)%*%(y))/2+.1
        
        psd <- isdmat*(theta[1]/cant1)
        theta[1] <- cant1
        isdmat <- psd
        
        a <- theta[2]
        b <- cant[2]
        curll    <- 0.5*as.numeric(determinant(psd)$modulus)-
          (t(y)%*%(psd)%*%(y))/2 + dbeta(thetaA[2], 9, 1, log = T)
        #dgamma(1/(thetaA[1]),.1,.1,log=TRUE) + 
        
        canll    <- 0.5*as.numeric(determinant(cansd)$modulus) -
          (t(y)%*%(cansd)%*%(y))/2 + dbeta(cant[2], 9, 1, log = T)
        #dgamma(1/(cant[1]),.1,.1,log=TRUE) + 
        
        Q1       <- dgamma(1/(thetaA[1]),(n1*n2)/2+.1,BB,log=TRUE)
        Q2       <- dgamma(1/(cant[1]),(n1*n2)/2+.1,bb,log=TRUE)
        R        <- exp(canll-curll+Q1-Q2)
        
        p_con <- dbeta((a*(1-b))/(a+b-2*a*b), sdl, sdl) / dbeta((b*(1-a))/(a+b-2*a*b),  sdl, sdl)
        p_coditional <- (b*(1 - b)) / (a*(1 - a))
        R <- R*p_coditional * p_con
        
        if(!is.na(R)){if((runif(1))< R){
          acceptedthetano_p[itr] <- 1
          theta  <- cant
          isdmat  <- cansd
        }}
        
      }
      
      Beta_p[[itr]] <- Beta
      
      if(itr %% 100 == 0){
        if(mean(acceptedthetano_p[1:itr]) > 0.45){sdl <- sdl*0.8}
        if(mean(acceptedthetano_p[1:itr]) < 0.3){sdl <- sdl*1.8}
        if(sdl<1){sdl=1}
        print(itr)
        print(mean((Beta-B0)^2))
      }
      
    }
    Beta_post <- Beta_p[501:2000]
    save(Beta_post, file = paste("3rdwork3rdsim", simtype,"rep", repli,".rda", sep =""))
  }
}