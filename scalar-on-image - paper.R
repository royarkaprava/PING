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


foreach(repli = 1:10) %dopar% {
  foreach(nuind = 1:2)  %dopar% {
    foreach(vind = 1:2)  %dopar% {
      q=3
      vx <- c(3, 6)
      var <- c(.1, 2)
      n1 <- 20
      n2 <- 20
      n <- 100
      set.seed(8)
      A1 <- rep(1:n2, each = n1)
      A2 <- rep(1:n1, n2)
      tempA <- cbind(A2, A1)
      #Ap <- matrix(rep(array(t(tempA)), n3), ncol=2, byrow = T)
      Ap <- tempA #cbind(Ap, rep(1:n3, each=n1*n2))
      
      
      m       <- 20    # Number of subjects
      RE      <- TRUE  # Generate data with random effects?
      
      
      pri.mn=c(0,0,0,0)
      pri.sd=c(10,2,10,1)
      
      L=1
      MHY=.01
      
      X <- matrix(0, n, n1*n2)
      
      loc <- Ap
      dis <- as.matrix(dist(loc))
      nux <- vx[nuind]
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
      
      init.theta= c(0,2,2,0)#c(0.3216614,  0.1695090, -1.2245358,  1.0982590)#c(0,2,2,0)
      theta <- init.theta
      
      sigma = 1
      
      isdmat <- solve(Maternvar(theta, dis))
      
      Beta.mean <- t(X) %*% Y / sigma^2
      Beta.ivar <- t(X) %*% X / sigma^2 + isdmat
      Beta.var <- solve(Beta.ivar)
      Beta.var <- (Beta.var + t(Beta.var))/2
      Beta.mean <- Beta.var %*% Beta.mean
      temp <- rmvnorm(1, Beta.mean, Beta.var)
      temp1 <- sign(temp)*(abs(temp))^(1/q)
      Betac <- matrix(rep(temp1, q), ncol = q)
      Beta <- array(temp)
      
      Total_itr <- 2000
      itr = 0
      theta_p <- list()
      Beta_p <- list()
      sigma_p <- rep(0, Total_itr)
      acceptedthetano_p <- rep(0, Total_itr)
      
      tol=0.000001
      
      sdl <- 1e-1
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
            B <- Betac[, -k]
            if(q > 2){
              B <- rowProds(Betac[, -k])
            }
            ivar <- isdmat*(k==1) + isdmat*exp(theta[1])*(k>1)
            X1 <- X*matrix(B, n, n1*n2, byrow = T)
            Beta.mean <- t(X1) %*% Y / sigma^2
            Beta.ivar <- t(X1) %*% X1 / sigma^2 + ivar
            Beta.var <- solve(Beta.ivar)
            Beta.var <- (Beta.var + t(Beta.var))/2
            Beta.mean <- Beta.var %*% Beta.mean
            Betac[, k] <- rmvnorm(1, Beta.mean, Beta.var) 
          }
          temp <- array(rowProds(Betac))
          Beta <- temp
          
          thetaA <- theta
          cant <- rep(0, 4)
          cant[-1]     <- thetaA[-1] + rnorm(3,sd = sdl) #MH[2]*tCthetaA%*%
          cansd    <- solve(Maternvar(cant, dis))
          psd <- isdmat    #Maternvar(thetaA, dis)
          y   <-     Betac
          bb       <- (t(y[, 1])%*%(cansd)%*%(y[, 1]))/2+.1
          cant1    <- -log(rgamma(1,(n1*n2)/2+.1,bb))
          cant[1]  <- cant1
          cansd    <- cansd/exp(cant1)
          BB       <- exp(thetaA[1])*(t(y[, 1])%*%(psd)%*%(y[, 1]))/2+.1
          
          term1    <- t(y[, 2])%*%psd%*%y[, 2]*exp(thetaA[1]) / 2
          
          if(q > 2){
            term1    <- sum(apply(y[, 2:q], 2, function(x){t(x)%*%psd%*%x*exp(thetaA[1])}))/2
          }
          
          term2    <- t(y[, 2])%*%cansd%*%y[, 2]*exp(cant[1]) / 2
          
          if(q > 2){
            term2    <- sum(apply(y[, 2:q], 2, function(x){t(x)%*%cansd%*%x*exp(cant[1])}))/2
          }
          
          curll    <- 0.5*as.numeric(determinant(psd)$modulus) + 0.5*(q-1)*as.numeric(determinant(psd*exp(thetaA[1]))$modulus)-
            (t(y[, 1])%*%(psd)%*%(y[, 1]))/2 - term1 +  
            sum(dnorm(thetaA[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
            dgamma(exp(-thetaA[1]),.1,.1,log=TRUE)
          canll    <- 0.5*as.numeric(determinant(cansd)$modulus) + 0.5*(q-1)*as.numeric(determinant(cansd*exp(cant1))$modulus)-
            (t(y[, 1])%*%(cansd)%*%(y[, 1]))/2 - term2 +
            sum(dnorm(cant[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
            dgamma(exp(-cant[1]),.1,.1,log=TRUE)
          Q1       <- dgamma(exp(-thetaA[1]),(n1*n2)/2+.1,BB,log=TRUE)
          Q2       <- dgamma(exp(-cant[1]),(n1*n2)/2+.1,bb,log=TRUE)
          R        <- canll-curll+Q1-Q2
          if(!is.na(R)){if(log(runif(1))< R){
            acceptedthetano_p[itr] <- 1
            theta  <- cant
            isdmat  <- cansd
          }}
        }
        
        if(q==1){
          Beta.mean <- t(X) %*% Y / sigma^2
          Beta.ivar <- t(X) %*% X / sigma^2 + isdmat
          Beta.var <- solve(Beta.ivar)
          Beta.var <- (Beta.var + t(Beta.var))/2
          Beta.mean <- Beta.var %*% Beta.mean
          temp <- array(rmvnorm(1, Beta.mean, Beta.var))
          Beta <- temp
          
          thetaA <- theta
          cant <- rep(0, 4)
          cant[-1]     <- thetaA[-1] + rnorm(3,sd = sdl) #MH[2]*tCthetaA%*%
          cansd    <- solve(Maternvar(cant, dis))
          psd <- isdmat    #Maternvar(thetaA, dis)
          y <-     Beta
          bb       <- (t(y)%*%(cansd)%*%(y))/2+.1
          cant1    <- -log(rgamma(1,(n1*n2)/2+.1,bb))
          cant[1]  <- cant1
          cansd    <- cansd/exp(cant1)
          BB       <- exp(thetaA[1])*(t(y)%*%(psd)%*%(y))/2+.1
          
          curll    <- 0.5*as.numeric(determinant(psd)$modulus) - (t(y)%*%(psd)%*%(y))/2 +
            sum(dnorm(thetaA[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
            dgamma(exp(-thetaA[1]),.1,.1,log=TRUE)
          canll    <- 0.5*as.numeric(determinant(cansd)$modulus) - (t(y)%*%(cansd)%*%(y))/2 +
            sum(dnorm(cant[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
            dgamma(exp(-cant[1]),.1,.1,log=TRUE)
          Q1       <- dgamma(exp(-thetaA[1]),(n1*n2)/2+.1,BB,log=TRUE)
          Q2       <- dgamma(exp(-cant[1]),(n1*n2)/2+.1,bb,log=TRUE)
          R        <- canll-curll+Q1-Q2
          if(!is.na(R)){if(log(runif(1))< R){
            acceptedthetano_p[itr] <- 1
            theta  <- cant
            isdmat  <- cansd
          }}
        }
        
        Beta_p[[itr]] <- Beta
        
        if(itr %% 100 == 0){
          if(mean(acceptedthetano_p[1:itr]) > 0.45){sdl <- sdl*1.2}
          if(mean(acceptedthetano_p[1:itr]) < 0.3){sdl <- sdl*0.8}
          print(itr)
          print(mean((Beta-B0)^2))
        }
        
      }
      
      Beta_post <- Beta_p[501:2000]
      save(Beta_post, file = paste("3rdwork3rdsim", q,"g", sigma0,"_",nux,"_", repli,".rda", sep =""))
    }
  }
}
