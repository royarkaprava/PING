library(mvtnorm)
library(geoR)
library(matrixStats)
library(fields)

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
  thetan[1] <- temp[1]*(1-temp[2])
  thetan[2] <- temp[1]*temp[2]
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
# 0.9, 0.18, 0.09
set.seed(8)

n1 <- 50
n2 <- 50
n <- 100
c1       <- round(runif(n, 1, n1))
c2       <- round(runif(n, 1, n1))
#dimvec <- c(n1, n2, n3)
buff     <- n1/8

thetaA0 <- c(0,log(0.95/0.05),log(10),0)

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

X <- matrix(rnorm(n1*n2), n1)
M<-m

#MU_b  <- matrix(0,J,M)
tauQ  <- 1

loc <- cbind(c1,c2)
dis <- as.matrix(dist(loc/50))

B0  <- simulateMatern(thetaA0,dis)

Beta0 <- matrix(0, n, 10)

for(i in 1:5){
  Beta0[, i] <- rep(0, n)
}
layout(matrix(1:6, nrow = 2, byrow = T), respect = T)
uh <- round(runif(5,0,2)) + 1
for(i in 6:10){
  h <- uh[i-5]
  u <- matrix(runif(h*2), ncol = 2, nrow=h)
  d <- 2*exp(-3*rowSums((Ap-matrix(c(n1*u[1,1], n2*u[1,2]), nrow = nrow(Ap), 2, byrow=T) )^2)/50)
  if(h == 2){
    d <- 2*exp(-3*rowSums((Ap-matrix(c(n1*u[1,1], n2*u[1,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50) + 2*exp(-3*rowSums((Ap-matrix(c(n1*u[2,1], n2*u[2,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50)
  }
  if(h == 3){
    d <- 2*exp(-3*rowSums((Ap-matrix(c(n1*u[1,1], n2*u[1,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50) + 
      2*exp(-3*rowSums((Ap-matrix(c(n1*u[2,1], n2*u[2,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50) + 2*exp(-3*rowSums((Ap-matrix(c(n1*u[3,1], n2*u[3,2]), nrow = nrow(Ap), 2, byrow=T))^2)/50)
  }
  
  B1  <- d
  B1[which(B1<1e-1)]<- 0  
  
  B1 <- matrix(B1, n1, n2)
  #image(B1, col = tim.colors(), main = paste("Beta_",i,sep = ""))
  Beta0[, i] <- B1[loc]
}

X <- array(0, c(n, 10, m))

for(i in 1:10){
  thetat <- rnorm(4)
  for(k in 1:m){
    X[,i,k] <- simulateMatern(thetat, dis) 
  }
}

for(q in 3){
  for(repli in 9:10){
    for(SNR in c(1, 5, 10)){
      
      set.seed(repli)
      
      #grd <- grid(1:n1, 1:n2) 
      
      
      mean <- apply(X, 3, function(x){rowSums(Beta0[, 1:10]*x)})
      temp <- matrix(rep(B0, m), ncol = m) + mean
      temp <- mean((apply(temp, 1, sd))^2)
      
      theta0  <- c(log(temp/SNR),log(0.90/0.10),log(10),0)
      
      Ylist <- matrix(0, n, m) 
      
      
      for(t in 1:m){
        Ylist[, t] <- B0 + rowSums(X[,,t]*Beta0) + simulateMatern(theta0, dis)
      }
      
      Beta <- matrix(0, n, 11)
      
      Beta[, 1] <- B0#cbind(colMeans(Ylist), Beta) #solve(crossprod(design))%*%(design*rowMeans(Ylist))
      Beta[, 2:11] <- Beta0
      #Beta <- lm(Ylist~design-1)$coefficients
      
      init.theta=c(0,2,2,0)
      theta <- matrix(init.theta, ncol = 4,nrow = (12),byrow=T)
      
      isdmat <- list()
      Betacom <- list()
      IVare <- solve(Maternvar(theta[1, ], dis))
      for(i in 1:12){
        isdmat[[i]] <- solve(Maternvar(theta[i, ], dis))
        if(i > 2){
          Ivar <- isdmat[[i]]
          mean <- apply(X, 3, function(x){rowSums(Beta[, c(2:11)[-i+2]]*x[, -i+2])})
          y <- Ylist - mean - matrix(rep(Beta[, 1], m), ncol = m)
          Beta.mean <- rowSums(sapply(1:m, function(k){diag((X[, i-2,k]))%*%IVare%*%(y[, k])}))
          Beta.ivar <- lapply(1:m, function(k){diag((X[,i-2,k]))%*%IVare%*%diag((X[,i-2,k]))})
          Beta.ivar <- Reduce('+', Beta.ivar) + Ivar
          Beta.var <- solve(Beta.ivar)
          Beta.var <- (Beta.var + t(Beta.var))/2
          Beta.mean <- Beta.var %*% Beta.mean
          Beta[, i-1] <- rmvnorm(1, Beta.mean, Beta.var)
          temp <- Beta[, i-1] #sqrt(abs(Betasp[, 2]))
          temp1 <- sign(temp)*(abs(temp))^(1/q)
          Betac <- matrix(rep(temp1, q), ncol = q)
          Betacom[[i-2]] <- Betac
        }
      }
      
      Total_itr <- 2000
      itr = 0
      Be0_p <- list()
      Be_p <- list()
      theta_p <- list()
      Beta_p <- list()
      acceptedthetano_p <- matrix(0, Total_itr, 12)
      accepteddeltamatno_p <- rep(0, Total_itr)
      accepteddeltano_p <- rep(0, Total_itr)
      
      tol=0.000001
      
      sdl <- rep(0.001, 12)
      sdl[1] <- 0.1
      
      Y <- Ylist
      while(itr < Total_itr){
        itr <- itr + 1
        flag <- rep(0, 12) 
        IVare <- isdmat[[1]]
        Ivar <- isdmat[[2]]
        
        mean <- apply(X, 3, function(x){rowSums(Beta[, 2:11]*x)})
        Beta.mean <- rowSums(IVare%*%(Ylist - mean))
        Beta.ivar <- m*IVare + Ivar
        Beta.var <- solve(Beta.ivar)
        Beta.var <- (Beta.var + t(Beta.var))/2
        Beta.mean <- (Beta.var %*% Beta.mean)
        Beta[, 1] <- rmvnorm(1, Beta.mean, Beta.var)
        
        if(q>1){
          for(i in 1:10){
            Betac <- Betacom[[i]]
            Ivar <- (isdmat[[i+2]])
            for(k in 1:q){
              B <- Betac[, - k]
              if(q>2){
                B <- rowProds(Betac[, - k]) 
              }
              Ivar <- Ivar*(k==1) + (Ivar/exp(theta[i+2, 1]))*(k>1)
              Betac[, k] <- updateBeta(i, Ylist-matrix(rep(Beta[, 1], m), ncol = m), Ivar, Beta, IVare, X, B)
            }
            Betacom[[i]] <- Betac
            Beta[, i+1] <- rowProds(Betac)
          }
          
          for(i in 3:12){
            thetaA <- theta[i, ]
            cant <- rep(0, 4)
            cant[-1]     <- thetaA[-1] + rnorm(3,sd = sdl[i]) #MH[2]*tCthetaA%*%
            cansd    <- solve(Maternvar(cant, dis))
            psd <- isdmat[[i]]    #Maternvar(thetaA, dis)
            y   <-     Betacom[[i-2]]
            bb       <- (t(y[, 1])%*%(cansd)%*%(y[, 1]))/2+.1
            cant1    <- -log(rgamma(1,n/2+.1,bb))
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
            Q1       <- dgamma(exp(-thetaA[1]),n/2+.1,BB,log=TRUE)
            Q2       <- dgamma(exp(-cant[1]),n/2+.1,bb,log=TRUE)
            R        <- canll-curll+Q1-Q2
            if(!is.na(R)){if(log(runif(1))< R){
              flag[i] <- 1
              theta[i,]  <- cant
              isdmat[[i]]  <- cansd
            }}
          }
        }
        
        if(q == 1){
          for(i in 1:10){
            Ivar <- isdmat[[i+2]]
            mean <- apply(X, 3, function(x){rowSums(Beta[, c(2:11)[-i]]*x[, -i])})
            y <- Ylist - matrix(rep(Beta[, 1], m), ncol = m) - mean
            Beta.mean <- rowSums(sapply(1:m, function(k){diag((X[,i,k]))%*%IVare%*%(y[, k])}))
            Beta.ivar <- lapply(1:m, function(k){diag((X[,i,k]))%*%IVare%*%diag((X[,i,k]))})
            Beta.ivar <- Reduce('+', Beta.ivar) +  Ivar
            Beta.var <- solve(Beta.ivar)
            Beta.var <- (Beta.var + t(Beta.var))/2
            Beta.mean <- Beta.var %*% Beta.mean
            Beta[, i+1] <- rmvnorm(1, Beta.mean, Beta.var)
          }
          
          for(i in 3:12){
            thetaA <- theta[i, ]
            cant <- rep(0, 4)
            cant[-1]     <- thetaA[-1] + rnorm(3,sd = sdl[i]) #MH[2]*tCthetaA%*%
            cansd    <- solve(Maternvar(cant, dis))
            psd <- isdmat[[i]]    #Maternvar(thetaA, dis)
            y <-     Beta[, i-1]
            bb       <- (t(y)%*%solve(cansd)%*%(y))/2+.1
            cant1    <- -log(rgamma(1,n/2+.1,bb))
            cant[1]  <- cant1
            cansd    <- cansd/exp(cant1)
            BB       <- exp(thetaA[1])*(t(y)%*%(psd)%*%(y))/2+.1
            
            curll    <- 0.5*as.numeric(determinant(psd)$modulus) - (t(y)%*%(psd)%*%(y))/2 +
              sum(dnorm(thetaA[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
              dgamma(exp(-thetaA[1]),.1,.1,log=TRUE)
            canll    <- 0.5*as.numeric(determinant(cansd)$modulus) - (t(y)%*%(cansd)%*%(y))/2 +
              sum(dnorm(cant[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
              dgamma(exp(-cant[1]),.1,.1,log=TRUE)
            Q1       <- dgamma(exp(-thetaA[1]),n/2+.1,BB,log=TRUE)
            Q2       <- dgamma(exp(-cant[1]),n/2+.1,bb,log=TRUE)
            R        <- canll-curll+Q1-Q2
            if(!is.na(R)){if(log(runif(1))< R){
              flag[i] <- 1
              theta[i,]  <- cant
              isdmat[[i]]  <- cansd
            }}
          }
          
        }
        
        thetaA <- theta[1, ]
        cant <- rep(0, 4)
        cant[-1]     <- thetaA[-1] + rnorm(3,sd = sdl[1]) #MH[2]*tCthetaA%*%
        cansd    <- solve(Maternvar(cant, dis))
        psd <- isdmat[[1]]    #Maternvar(thetaA, dis)
        mean <- apply(X, 3, function(x){rowSums(Beta[, 2:11]*x)})
        y <-     Ylist - matrix(rep(Beta[, 1], m), ncol = m) - mean
        bb       <- sum(apply(y, 2, function(x){t(x)%*%cansd%*%x}))/2+.1
        cant1    <- -log(rgamma(1,m*n/2+.1,bb))
        cant[1]  <- cant1
        cansd    <- cansd/exp(cant1)
        BB       <- exp(thetaA[1])*sum(apply(y, 2, function(x){t(x)%*%psd%*%x}))/2+.1
        
        curll    <- 0.5*m*as.numeric(determinant(psd)$modulus) - sum(apply(y, 2, function(x){t(x)%*%psd%*%x}))/2+
          sum(dnorm(thetaA[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
          dgamma(exp(-thetaA[1]),.1,.1,log=TRUE)
        canll    <- 0.5*m*as.numeric(determinant(cansd)$modulus) - sum(apply(y, 2, function(x){t(x)%*%cansd%*%x}))/2+
          sum(dnorm(cant[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
          dgamma(exp(-cant[1]),.1,.1,log=TRUE)
        Q1       <- dgamma(exp(-thetaA[1]),m*n/2+.1,BB,log=TRUE)
        Q2       <- dgamma(exp(-cant[1]),m*n/2+.1,bb,log=TRUE)
        R        <- canll-curll+Q1-Q2
        if(!is.na(R)){if(log(runif(1))< R){
          flag[1] <- 1
          theta[1,]  <- cant
          isdmat[[1]]  <- cansd
        }}
        
        thetaA <- theta[2, ]
        cant <- rep(0, 4)
        cant[-1]     <- thetaA[-1] + rnorm(3,sd = sdl[2]) #MH[2]*tCthetaA%*%
        cansd    <- solve(Maternvar(cant, dis))
        psd <- isdmat[[2]]    #Maternvar(thetaA, dis)
        y <-     Beta[, 1]
        bb       <- (t(y)%*%solve(cansd)%*%(y))/2+.1
        cant1    <- -log(rgamma(1,n/2+.1,bb))
        cant[1]  <- cant1
        cansd    <- cansd/exp(cant1)
        BB       <- exp(thetaA[1])*(t(y)%*%(psd)%*%(y))/2+.1
        
        curll    <- 0.5*as.numeric(determinant(psd)$modulus) - (t(y)%*%(psd)%*%(y))/2 +
          sum(dnorm(thetaA[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
          dgamma(exp(-thetaA[1]),.1,.1,log=TRUE)
        canll    <- 0.5*as.numeric(determinant(cansd)$modulus) - (t(y)%*%(cansd)%*%(y))/2 +
          sum(dnorm(cant[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
          dgamma(exp(-cant[1]),.1,.1,log=TRUE)
        Q1       <- dgamma(exp(-thetaA[1]),n/2+.1,BB,log=TRUE)
        Q2       <- dgamma(exp(-cant[1]),n/2+.1,bb,log=TRUE)
        R        <- canll-curll+Q1-Q2
        if(!is.na(R)){if(log(runif(1))< R){
          flag[2] <- 1
          theta[2,]  <- cant
          isdmat[[2]]  <- cansd
        }}
        acceptedthetano_p[itr, ] <- flag
        Beta_p[[itr]] <- Beta
        
        if(itr %% 100 == 0){
          for(l in 1:12){
            if(mean(acceptedthetano_p[1:itr,l]) > 0.45){sdl[l] <- sdl[l]*1.2}
            if(mean(acceptedthetano_p[1:itr,l]) < 0.3){sdl[l] <- sdl[l]*0.8}
          }
        }
        if(itr %% 100 == 0) {
          print(itr)
          print(mean((Beta[, 2:11]-Beta0)^2))}
        #theta_p[[itr]] <- theta
      }
      Beta_post <- Beta_p[501:2000]
      save(Beta_post, file = paste("3rdwork2ndsim", q,"g", SNR,"_", repli,".rda", sep =""))
    }}}

