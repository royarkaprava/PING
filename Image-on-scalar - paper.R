library(mvtnorm)
library(splines)
library(parallel)
library(MCMCpack)
library(foreach)
#library(doMC)
library(matrixStats)
library(fields)

greal <- function(dim){ #get_real
  B = 0
  ind=0
  yf <- c(1,2)
  tol <- 1e-9
  m <- 1 - dim %% 2
  
  while(sum(yf)!=(prod(dim) - prod(2^m))/2 | length(ind)!=sum(yf)){
    yf <- rep(0, prod(dim))
    #while(length(unique(B))<prod(dim)){
    B <- array(rnorm(prod(dim),sd=length(dim)*exp(mean(log(dim)))), dim)
    #}
    
    Bfft <- fft(B)/sqrt(prod(dim))
    #bfa <- bfa/4
    bfa <- array(Re(Bfft))
    tem <- cbind(bfa, 1:prod(dim))
    tem <- tem[order(tem[,2]),]
    tem <- tem[order(tem[,1]),]
    tem1 <- tem[2:prod(dim),1] - tem[1:(prod(dim)-1),1]
    c <- which(abs(tem1)<tol)
    c1 <- (tem[c+1,2] < tem[c,2])
    yf[tem[c[c1], 2]] <- 1
    yf[tem[c[!c1]+1,2]]<-1
    
    yf <- as.logical(yf)
    
    clust <- rep(0, prod(dim))
    clust[!yf] <- 1:sum(!yf)
    
    c2 <- bfa[yf]
    ind <- sort(union(tem[c[c1]+1, 2], tem[c[!c1],2]))
    tol <- tol/10
  }
  
  temp1 <- cbind(bfa[ind], clust[ind])
  temp2 <- cbind(c2, which(yf==T))
  temp1 <- temp1[order(temp1[,1]),]
  temp2 <- temp2[order(temp2[,1]),]
  clust[temp2[, 2]] <- temp1[,2]
  
  return(list(realindex = 1 - yf, clustervec = clust))
}




tilda <- function(grealout, xf){
  id <- grealout$realindex
  xtilda <- array(0, prod(dim(xf)))
  xtilda[as.logical(id)] <- array(Re(xf))[as.logical(id)]
  xtilda[!id] <- array(Im(xf))[!id]
  xtilda <- array(xtilda, dim(xf))
}

recon1 <- function(grealout, xtilda){
  id <- grealout$realindex
  clust <- grealout$clustervec
  id <- as.logical(id)
  y <- xtilda[id]
  im <- xtilda[!id]
  imtotal <- array(0, prod(dim(xtilda)))
  realtotal <- y[clust]
  imtotal[!id] <- im
  c1 <- id*clust
  c2 <- clust[!id]
  ind <- which(c1 %in% c2)
  c2y <- cbind(c2, im)
  c2y <- c2y[order(c2y[,1]),]
  imtotal[ind] <- - c2y[,2]
  return(array(complex(real=realtotal,imaginary=imtotal), dim(xtilda)))
}


convert1 <- function(x,gout=NULL,inv=FALSE){
  if(is.null(gout)){
    gout <- greal(dim(x))
  }
  if(inv==FALSE){
    y <- fft(x)/sqrt(prod(dim(x)))
    y <- tilda(gout, y)
  }
  if(inv==TRUE){
    y <- recon1(gout,x)
    y <- Re(fft(y,inv=T))/sqrt(prod(dim(x)))
  }
  return(y)}

FFT <- function(Y, dimvec = dim(Y)){
  Y <- array(Y, dimvec)
  fft(Y)/sqrt(prod(dim(Y)))
}

IFFT <- function(Y, dimvec = dim(Y)){
  Y <- array(Y, dimvec)
  fft(Y,inverse=TRUE)/sqrt(prod(dim(Y)))
}

homega3D <- function(dim){ #make.delta2D
  n1 <- dim[1]
  n2 <- dim[2]
  n3 <- dim[3]
  omega1 <- 2*pi*(seq(0,n1-1,1))/n1        
  omega2 <- 2*pi*(seq(0,n2-1,1))/n2
  omega3 <- 2*pi*(seq(0,n3-1,1))/n3
  
  #create matrix each column represent co-ordinate of col-major 3D matrix
  A1 <- rep(1:n2, each = n1)
  A2 <- rep(1:n1, n2)
  tempA <- cbind(A2, A1)
  A <- matrix(rep(array(t(tempA)), n3), ncol=2, byrow = T)
  A <- cbind(A, rep(1:n3, each=n1*n2))
  
  omegacombine <- rbind(omega1[A[,1]],omega2[A[,2]], omega3[A[,3]])
  
  h <- colSums((sin(omegacombine/2))^2)
  return(h)
} 

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

# Matern spectral density
spec_d_matern <- function(delta,theta,dim,d=2,nugget=TRUE, whichreal){
  pars     <- theta2mat(theta,nugget)
  num      <- -(pars[4]+d/2)
  den      <- 1/pars[3]^2+delta/4 #why delta/4?
  f        <- den^num
  f        <- prod(dim)*f/sum(f) #why this normalization?
  f        <- (1-pars[2]) + pars[2]*f
  f            <- pars[1]*f
  f[whichreal] <- f[whichreal]/2
  return(f)}

simulateMatern <- function(specd,dim){
  Z     <- array(rnorm(prod(dim)),dim)
  y     <- sqrt(specd)*FFT(Z)
  Y     <- Re(IFFT(y))
  return(Y)}

impute <- function(Y.o, mu,miss,specd, dimvec,tol,gout,maxiter){
  Y.o <- array(Y.o, dimvec)
  mu <- array(mu, dimvec)
  miss <- array(miss, dimvec)
  Y.o  <- mu+update.Y(Y.o-mu,miss,specd, dimvec,tol=tol,maxiter=maxiter)
  dat  <- convert1(array(Y.o,dimvec),gout)
  return(dat)}

condExp <- function(Arr,A,miss, dimvec,specd,maxiter=50,tol=0.1){
  xArr         <- 0*specd
  xArr[!miss]  <- pcg(A,Arr,miss, dimvec,specd,maxiter=maxiter,tol=tol)
  ArrImp       <- A(xArr,miss,specd, dimvec, all=TRUE)
  return(ArrImp)}     



A <- function(y,miss,specd, dimvec, all=FALSE){
  y[miss] <- 0
  x       <- Re(IFFT(sqrt(specd)*FFT(y,dimvec),dimvec))
  if(!all){x<-x[!miss]}
  return(x)}



pcg <-function(A,yArr,miss,dimvec,specd,maxiter=1e4,tol=1e-1){
  
  # A is forward multiplication, takes in an array
  
  # we define the preconditioner to be A with 1/specden   
  
  # if * denotes a vector, *Arr denotes corresp. array with zeros or NAs
  # in place of missing value locations
  
  y     = yArr[!miss]
  
  x     = rep(0,sum(!miss))
  xArr  = array(0,dimvec)
  xArr[!miss] <- x;
  r     = y-A(xArr,miss,specd, dimvec)
  z     = A(yArr,miss,1/specd, dimvec)
  p     = z
  iter  = 0
  pArr <- array(0,dimvec)
  r1Arr <- array(0,dimvec)
  while (mean(r^2)>tol & iter<maxiter){
    iter  = iter+1
    pArr[!miss] <- p  
    Ap    = A(pArr,miss,specd, dimvec)
    a     = sum(r*z)/sum(p*Ap)
    x     = x+a*p
    r1    = r-a*Ap
    r1Arr[!miss] <- r1   
    z1    = A(r1Arr,miss,1/specd, dimvec)
    bet   = sum(z1*r1)/sum(z*r)
    p     = z1+bet*p
    z     = z1
    r     = r1
  }
  return(x)}

update.Y <- function(Y,miss,specd, dimvec,maxiter=50,tol=0.1){
  YcondExp     <- condExp(Y,A,miss, dimvec,specd,maxiter=maxiter,tol=tol)
  eArr         <- simulateMatern(specd,dimvec)
  econdExp     <- condExp(eArr,A,miss, dimvec,specd,maxiter=maxiter,tol=tol)
  Y[miss]      <- YcondExp[miss] + eArr[miss]-econdExp[miss] 
  return(Y)}

##############################################################################


updatetheta <- function(l, y, theta, sdl1){
  n <- prod(dimvec)
  y <- y[,l]
  speccur <- SPECD
  if(l>1)
  {speccur <- SPECDC[,l-1]}
  flag = 0
  thetaA <- theta[l, ]
  cant     <- thetaA[-1] + rnorm(3,sd = sdl1) #MH[2]*tCthetaA%*%
  cansd    <- spec_d_matern(h.omega,c(0,cant),dimvec, whichreal = realind)
  #bb       <- sum(y/cansd)/2+.1
  #cant1    <- -log(rgamma(1,n/2+.1,bb))
  aaa       <- 1/sdl1
  bbb       <- 1/sdl1
  cant1    <- -log(rgamma(1,aaa,bbb/exp(-thetaA[1])))
  cant     <- c(cant1,cant)
  #cansd    <- exp(cant1)*cansd
  #BB       <- exp(thetaA[1])*sum(y/speccur)/2+.1
  
  curll    <- -0.5*sum(log(speccur)) - 0.5*sum(y/speccur)+
    sum(dnorm(thetaA[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
    dcauchy(exp(-thetaA[1]),scale = delta[l-2],log=TRUE)#dgamma(exp(-thetaA[1]),.1,.1,log=TRUE)
  canll    <- -0.5*sum(log(cansd)) - 0.5*sum(y/cansd)+
    sum(dnorm(cant[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
    dcauchy(exp(-cant[1]),scale = delta[l-2],log=TRUE)#dgamma(exp(-cant[1]),.1,.1,log=TRUE)
  Q1       <- dgamma(exp(-cant1),aaa,bbb/exp(-thetaA[1]),log=TRUE)#dgamma(exp(-thetaA[1]),n/2+.1,BB,log=TRUE)
  Q2       <- dgamma(exp(-thetaA[1]),   aaa,bbb/exp(-cant1),log=TRUE)#dgamma(exp(-cant[1]),n/2+.1,bb,log=TRUE)
  R        <- canll-curll+Q1-Q2
  if(!is.na(R)){if(log(runif(1))< R){
    flag <- 1
    theta[l, ]  <- cant
    SPECDA  <- cansd
    curll <- canll
  }}
  if(l>1){
    SPECDC[,l-1] <- cansd
  }
  else{
    SPECD <- cansd
  }
  
  return(list(flag=flag, theta=thetaA))
}

updatetheta1 <- function(l, y, theta, sdl1){
  n <- prod(dimvec)
  y <- y[,l]
  speccur <- SPECD
  if(l>1)
  {speccur <- SPECDC[,l-1]}
  flag = 0
  thetaA <- theta[l, ]
  cant     <- thetaA[-1] + rnorm(3,sd = sdl1) #MH[2]*tCthetaA%*%
  cansd    <- spec_d_matern(h.omega,c(0,cant),dimvec, whichreal = realind)
  bb       <- sum(y/cansd)/2+.1
  cant1    <- -log(rgamma(1,n/2+.1,bb))
  cant     <- c(cant1,cant)
  cansd    <- exp(cant1)*cansd
  BB       <- exp(thetaA[1])*sum(y/speccur)/2+.1
  
  curll    <- -0.5*sum(log(speccur)) - 0.5*sum(y/speccur)+
    sum(dnorm(thetaA[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
    dgamma(exp(-thetaA[1]),.1,.1,log=TRUE)
  canll    <- -0.5*sum(log(cansd)) - 0.5*sum(y/cansd)+
    sum(dnorm(cant[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
    dgamma(exp(-cant[1]),.1,.1,log=TRUE)
  Q1       <- dgamma(exp(-thetaA[1]),n/2+.1,BB,log=TRUE)
  Q2       <- dgamma(exp(-cant[1]),n/2+.1,bb,log=TRUE)
  R        <- canll-curll+Q1-Q2
  if(!is.na(R)){if(log(runif(1))< R){
    flag <- 1
    thetaA  <- cant
    SPECDA  <- cansd
    curll <- canll
  }}
  if(l>1){
    SPECDC[,l-1] <- cansd
  }
  if(l==1){
    SPECD <- cansd
  }
  
  return(list(flag=flag, theta=thetaA))
}


updatebeta <- function(l, Beta){
  mean1     <- Beta[ch.omega, ]
  mean1 <- mean1[, l]
  SPECDl <- SPECDC[ch.omega, ]
  SPECDl <- SPECDl[, l]
  mean <- mean1/SPECDl
  mean <- (t(B[ch.omega,])%*%mean)
  betavar <- solve(diag(1/deltamat[, l]) + crossprod(B[ch.omega, ])*sum(1/SPECDl))
  mean <- betavar%*%array(mean)
  betal <- rmvnorm(1, mean, betavar)
  return(betal)
}

updateBetaomega <- function(omegaid, y){
  Betaomega.mean1 <- colSums(matrix(y[omegaid, ],nrow = nrow(X), ncol = ncol(X), byrow = T)*X) / SPECD[omegaid]
  Betaomega.mean2 <- B[omegaid, ] %*% beta /SPECDC[omegaid, ]
  Betaomega.var <- solve(CX/SPECD[omegaid] + sum(1/SPECDC[omegaid, ]))
  Betaomega.mean <- Betaomega.var%*%array(Betaomega.mean1+Betaomega.mean2)
  Betaomega <- rmvnorm(1, Betaomega.mean, Betaomega.var)
  return(Betaomega)
}

updateBeta <- function(l, y, SPECDCa, Beta, SPECD, X, B, n){
  Beta.mean1 <- rowSums((y - Beta[,-l]%*%t(X[,-l]))*matrix(X[,l], nrow=nrow(y),ncol=ncol(y),byrow=T)) / SPECD
  #Beta.mean2 <- B %*% beta[, l] /SPECDC[, l]
  Beta.var <- array(sum(X[,l]^2) / SPECD + 1 / SPECDCa)
  Beta.mean <- array(Beta.mean1) #+ array(Beta.mean2)
  gen <- rnorm(n, Beta.mean/Beta.var, sqrt(1/Beta.var))
  return(gen)
}


updatebeta1 <- function(l, y, Beta){
  mean1     <- rowSums((y - Beta[,-l]%*%t(X[,-l]))*matrix(X[,l], nrow=nrow(y),ncol=ncol(y),byrow=T)) / SPECD
  mean <- (t(B)%*%mean1)
  betavar <- solve(diag(1/deltamat[, l]) + sum(X[,l]^2)*crossprod(B/sqrt(SPECD)))
  mean <- betavar%*%array(mean)
  betal <- rmvnorm(1, mean, betavar)
  return(betal)
}

updatedeltamat <- function(l, deltamat, delta, sdl2){
  flag      <- 0
  sig       <- deltamat[, l]
  #s0        <- 1/sqrt(taub[q+1])   
  aaa       <- 1/sdl2
  bbb       <- 1/sdl2
  cansig    <- rgamma(J,aaa,bbb/sig)
  lp1       <- dnorm(beta[, l],0,cansig,log=TRUE) + dcauchy(cansig,scale = delta[l],log=TRUE)
  lp2       <- dnorm(beta[, l],0,sig,   log=TRUE) + dcauchy(   sig,scale = delta[l],log=TRUE)
  lq1       <- dgamma(cansig,aaa,bbb/sig,log=TRUE)
  lq2       <- dgamma(sig,   aaa,bbb/cansig,log=TRUE)
  R         <- lp1-lp2+lq2-lq1
  R         <- (log(runif(J))<R)
  sig       <- ifelse(R,cansig,sig)
  flag      <- (flag+sum(R))/length(R)
  deltamat[, l] <- array(sig)
  return(flag)
}

updatedelta <- function(l, deltamat, delta, sdl3){
  flag      <- 0
  sig       <- delta[l]
  aaa       <- 1/sdl3
  bbb       <- 1/sdl3
  cansig    <- rgamma(1,aaa,bbb/sig)
  lp1       <- sum(dcauchy(deltamat[l],0,cansig,log=TRUE) + dgamma(1/cansig,0.1,0.1,log=TRUE))
  lp2       <- sum(dcauchy(deltamat[l],0,   sig,log=TRUE) + dgamma(1/sig,0.1,0.1,log=TRUE))
  lq1       <- sum(dgamma(cansig,aaa,bbb/sig,log=TRUE))
  lq2       <- sum(dgamma(sig,   aaa,bbb/cansig,log=TRUE))
  R         <- lp1-lp2+lq2-lq1
  if(log(runif(1))<R)
  {
    delta[l-1] <- cansig
    flag <- 1
  }
  return(flag)
}


HMC = function (U, grad_U, epsilon, L, current_q,x,y,z=Beta[,1],a=SPECD,b=fixedy)
{
  q = current_q
  p = rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p
  flag = 0
  # Make a half step for momentum at the beginning
  
  p = p - epsilon * grad_U(q) / 2
  
  # Alternate full steps for position and momentum
  
  for (k in 1:L)
  {
    # Make a full step for the position
    
    q = q + epsilon * p
    
    # Make a full step for the momentum, except at end of trajectory
    
    if (k!=L) p = p - epsilon * grad_U(q)
  }
  
  # Make a half step for momentum at the end.
  
  p = p - epsilon * grad_U(q) / 2
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  
  current_U = U(current_q,x,y,z=Beta[,1],a=SPECD,b=fixedy)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q,x,y,z=Beta[,1],a=SPECD,b=fixedy)
  proposed_K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  
  
  R <- exp(current_U-proposed_U+current_K-proposed_K)
  if(is.nan(R)){R=1}
  ret <- list(up = current_q, flag = flag)
  {if (runif(1) < R)
  {
    # accept
    flag <- flag + 1
    ret <- list(up = q, flag = flag) 
  }}
  return(ret)
}

U <- function(can,x,y,z=Beta[,1],a=SPECD,b=fixedy){
  Betacan <- array(convert1(array(can*(x),dimvec), gout))
  fcan <- array(convert1(array(can,dimvec), gout))
  term1 <- apply((b - tcrossprod(cbind(z,Betacan), X)), 2, function(x){sum(dnorm(x, sd= sqrt(a), log=T))})
  return(- sum(term1) - sum(dnorm(fcan, sd = sqrt(y), log = T)))
}

grad_U <- function(can){
  canf <- array(convert1(array(can,dimvec), gout, inv = T))
  Betacan <- array(convert1(array(can*rowProds(Betac[, -i]),dimvec), gout))
  #fcan <- array(convert1(array(can,dimvec), gout))
  marg <- rowSums((fixedy - tcrossprod(cbind(Beta[,1],Betacan), X))*matrix(X[,2], nrow = nrow(Y), M, byrow = T))/SPECD
  temp <- array(rowProds(Betac[, -i]), dimvec)
  temp <- convert1(temp, gout)
  temp[21 - Ap] <- temp
  temp <- array(convert1(temp, gout, inv = T))
  temp1 <- (array(convert1(array(marg, dimvec), gout, inv=T)))*temp
  temp1 <- array(convert1(array(temp1, dimvec), gout))
  grad <- -temp1 + canf / SPECDC[, i+1]
  gradf <- array(convert1(array(grad,dimvec), gout, inv = T))
  return(gradf )
}

n1       <- 20
n2       <- 20
n3       <- 20  
n        <- n1*n2*n3
dimvec <- c(n1, n2, n3)
buff     <- 0

theta  <- c(log(.02),log(0.90/0.10),log(10),0)
thetaA <- c(0,log(0.95/0.05),log(10),0)
thetaB <- c(0,log(0.95/0.05),log(10),0)
Sigma  <- matrix(c(0.5,-0.5,-0.5,1),2,2)


delta <- homega3D(dimvec)
#id    <- greal(dim)$realindex

#id <- greal3(n1, n2, n3)

mat    <- theta2mat(theta)
matA   <- theta2mat(thetaA)
matB   <- theta2mat(thetaB)

specd0  <- spec_d_matern(delta,theta,dimvec)
specdA0 <- spec_d_matern(delta,thetaA,dimvec)
specdB <- spec_d_matern(delta,thetaB,dimvec)

A1 <- rep(1:n2, each = n1)
A2 <- rep(1:n1, n2)
tempA <- cbind(A2, A1)
Ap <- matrix(rep(array(t(tempA)), n3), ncol=2, byrow = T)
Ap <- cbind(Ap, rep(1:n3, each=n1*n2))

Ap2 <- Ap-matrix(c(n1*.3, n2*.5,n3*.7), nrow = nrow(Ap), 3, byrow=T) 
Ap1 <- Ap-matrix(c(n1*.3, n2*.7,n3*.3), nrow = nrow(Ap), 3, byrow=T)
Ap3 <- Ap-matrix(c(n1*.7, n2*.3,n3*.7), nrow = nrow(Ap), 3, byrow=T)


d1 <- 2*exp(-4*rowSums(Ap1^2)/20)
d2 <- 2*exp(-1.5*rowSums(Ap2^2)/20)
d3 <- 2*exp(-4*rowSums(Ap3^2)/20)
d4 <- 2*exp(-4*rowSums((Ap-n1*.7)^2)/20)
d5 <- 2*exp(-4*rowSums((Ap-n1*.3)^2)/20)
Z1     <- ifelse(Ap[, 3]<n3/2,1,0)
Z2     <- 1-Z1

# Missing data locations

# miss  <- array(FALSE, dimvec)
# miss[n1-1:buff+1,,]<- TRUE
# miss[,n2-1:buff+1,]<- TRUE
# miss[,,n3-1:buff+1]<- TRUE
# miss[,1:buff,]<- TRUE
# miss[1:buff,,]<- TRUE 
# miss[,,1:buff]<- TRUE 

m       <- 20    # Number of subjects
RE      <- F  # Generate data with random effects?


dimvec <- c(n1,n2,n3)
n_spline <- 6
J <- n_spline^3

gout <- greal(c(n1, n2, n3))

rid <- as.logical(gout$realindex)

h.omega <- homega3D(c(n1, n2, n3))
ch.omega <- which(h.omega <= quantile(h.omega, J/(n1*n2*n3))) 

realind <- array(0, c(n1, n2, n3))

realind[1, 1, 1]  <- 1
if(n1%%2==0){
  realind[(n1/2+1), 1, 1] <- 1
  if(n2%%2==0){
    realind[(n1/2+1), (n2/2+1), 1] <- 1
    realind[1, (n2/2+1), 1] <- 1
  }
  if(n3%%2==0){
    realind[(n1/2+1), (n2/2+1), (n3/2+1)] <- 1
    realind[1, (n2/2+1), (n3/2+1)] <- 1
    realind[(n1/2+1), 1, (n3/2+1)] <- 1
    realind[1, 1, (n3/2+1)]  <- 1
  }}
realind <- as.logical(array(realind))

#index of a 3D array

A1 <- rep(1:n2, each = n1)
A2 <- rep(1:n1, n2)
tempA <- cbind(A2, A1)
Ap <- matrix(rep(array(t(tempA)), n3), ncol=2, byrow = T)
Ap <- cbind(Ap, rep(1:n3, each=n1*n2))

if(n_spline>2){
  B1 <- bs(1:n1,df=n_spline,intercept=TRUE)
  B2 <- bs(1:n2,df=n_spline,intercept=TRUE)
  B3 <- bs(1:n3,df=n_spline,intercept=TRUE)
  B  <- NULL
  #whichrB <- allw[activeB]
  count<-1
  for(j3 in 1:n_spline){for(j2 in 1:n_spline){
    for(j1 in 1:n_spline){
      b <- outer(B1[,j1],B2[,j2])
      b <- array(outer(b, B3[,j3]))
      B <- cbind(B,convert1(array(b,dimvec),gout))
      count <- count + 1 
    }
  }}
  Q  <- diag(J)
}

B[-ch.omega, ] <- 0

Br <- apply(B, 2 ,function(x){convert1(array(x,dimvec), gout, inv = T)})

pri.mn=c(0,0,0,0)
pri.sd=c(10,2,10,1)
L=1
MHY=.01

X   <- 1:m
X <- scale(X)
X <- cbind(rep(1,m),X)
Beta <- matrix(0, n1*n2*n3, ncol(X))
M<-m

MU_b  <- matrix(0,J,M)
tauQ  <- 1
init.theta=c(0,2,2,0)
theta <- matrix(init.theta, ncol = 4,nrow = (7),byrow=T)


B0  <- simulateMatern(specdA0,dimvec)
B1  <- d1 + d2 + d3 + d4 + d5
B1[which(B1<1e-1)]<- 0 
B10 <- array(B1, dimvec)
image.plot(B10[,12,])

#set.seed(Repli[repli])

Ylist <- matrix(0, n, m)
for(t in 1:m){
  re         <- t(chol(Sigma))%*%rnorm(2)
  Ylist[, t] <- B0 + X[t, 2]*B1 + simulateMatern(specd0,dimvec)
  #Ylist[, t] <- ifelse(miss,NA,Ylist[, t])
}

Y.o <- Ylist
miss1    <- is.na(Y.o)
miss2  <- is.na(Y.o[, 1])
X <- X
# number of individuals/time points M


Y.o <- apply(Y.o, 2, function(x){c<-which(is.na(x));x[c]<-mean(x, na.rm = T);return(x)})
#theta[1,1] <- log(.004)
SPECD <- spec_d_matern(h.omega,theta[1, ],dimvec, whichreal = realind)

SPECDC <- mcmapply(2:(7), FUN=function(i){spec_d_matern(h.omega,theta[i, ],dimvec, whichreal = realind)})
#SPECDC[, 2] <- spec_d_matern(h.omega,theta[3, ],dimvec,nugget = F, whichreal = realind)
#SPECDC[, 3] <- spec_d_matern(h.omega,theta[4, ],dimvec,nugget = F, whichreal = realind)
specd <- SPECD

specd[realind] <- 2*specd[realind]

Y <- apply(Y.o, 2 ,function(x){convert1(array(x,dimvec), gout)})

#no of components
p= 1

M<-m

#################################################################################################################################

####starting values

deltamat <- matrix(1/rgamma(p*J, 10, 10), J, p)

delta <- 1/rgamma(p, 10, 10)

deltamat <- matrix(1,ncol=p+1,nrow = J)

#beta <- matrix(rnorm((p+1)*J, sd= deltamat), nrow = J)

Betasp <- matrix(0, n1*n2*n3, ncol(X))  

CX <- crossprod(X)
Betasp <- t(solve(CX)%*%crossprod(X,t(Y.o)))

#resid <- array(convert1(array(rowMeans(Y[,5:10])-Beta1p, dimvec), gout, inv = T))/mean(X[5:10, 2])
temp <- Betasp[, 2] #sqrt(abs(Betasp[, 2]))
temp1 <- sign(temp)*(abs(temp))^(1/5)
Betac <- cbind(temp1, temp1, temp1, temp1, temp1)
Beta <- t(solve(CX)%*%crossprod(X,t(Y)))
#Betac <- matrix(rnorm(2*nrow(Y)),ncol=2)

#Beta[, 1] <- rnorm(nrow(Y), sd = SPECDC[,1])
#Beta[, 2] <- array(convert1(array(Betac[,1]*Betac[,2],dimvec), gout))

#can1 <- Beta[,2]/2
#can2 <- Beta[,2]*2

#beta <- solve(crossprod(B))%*% crossprod(B,Beta)

Q <- diag(J)
itr = 0
Betacf <- apply(Betac, 2 ,function(x){convert1(array(x,dimvec), gout)})
#MCMC
Total_itr <- 7000
itr = 2000
Beta_p <- list()
theta_p <- list()
Betac_p <- list()
acceptedthetano_p <- matrix(0, Total_itr, 3)
accepteddeltamatno_p <- rep(0, Total_itr)
accepteddeltano_p <- rep(0, Total_itr)

tol=0.000001

sdl1 <- sdl2 <- sdl3 <- 0.05
tauQ <- 1
MU <- rep(0, length(ch.omega))
flag <- rep(0, 3)
sd = .05
flagbeta2 <- 0
c=0
while(itr < Total_itr)
{
  flag <- rep(0, 3) 
  itr <- itr + 1
  
  fixedy <- matrix(0, nrow(Y), ncol(Y))
  
  fixedy[ - ch.omega, ] <- Y[ - ch.omega, ]
  
  fixedy[ch.omega, ] <- Y[ch.omega, ] - MU
  
  #Beta <- mcmapply(1:nrow(Y), FUN=updateBetaomega, MoreArgs = list(fixedy))
  
  #Beta <- t(Beta)
  
  #Beta[, 1] <- array(updateBeta(1, fixedy, SPECDC[,1], Beta, SPECD, X, B, n))
  #Beta1_p[[itr]] <- Beta[, 1]
  
  if(itr<=2000){
    Beta[, 2] <- array(updateBeta(2, fixedy, SPECDC[,2], Beta, SPECD, X, B, n))
    
    #Beta_p[[itr]] <- Beta
    #Betac_p[[itr]] <- Betac
    
    #beta <- mcmapply(1:ncol(X), FUN=updatebeta, MoreArgs = list(Beta))
    
    residy <- rowSums((fixedy - Beta%*%t(X))^2)
    
    ty <- cbind(residy, Beta^2)
    
    l=1
    n <- prod(dimvec)
    y <- ty[,l]
    speccur <- SPECD
    SPECDA <- speccur
    thetaA <- theta[l, ]
    cant     <- thetaA[-1] + rnorm(3,sd = sdl1) #MH[2]*tCthetaA%*%
    cansd    <- spec_d_matern(h.omega,c(0,cant),dimvec, whichreal = realind)
    bb       <- sum(y/cansd)/2+.1
    cant1    <- -log(rgamma(1,M*n/2+.1,bb))
    cant     <- c(cant1,cant)
    cansd    <- exp(cant1)*cansd
    BB       <- exp(thetaA[1])*sum(y/speccur)/2+.1
    
    curll    <- -0.5*M*sum(log(speccur)) - 0.5*sum(y/speccur)+
      sum(dnorm(thetaA[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
      dgamma(exp(-thetaA[1]),.1,.1,log=TRUE)
    canll    <- -0.5*M*sum(log(cansd)) - 0.5*sum(y/cansd)+
      sum(dnorm(cant[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
      dgamma(exp(-cant[1]),.1,.1,log=TRUE)
    Q1       <- dgamma(exp(-thetaA[1]),M*n/2+.1,BB,log=TRUE)
    Q2       <- dgamma(exp(-cant[1]),M*n/2+.1,bb,log=TRUE)
    R        <- canll-curll+Q1-Q2
    if(!is.na(R)){if(log(runif(1))< R){
      flag[l] <- 1
      theta[l,]  <- cant
      SPECDA  <- cansd
      curll <- canll
    }}
    if(l>1){
      SPECDC[,l-1] <- SPECDA
    }
    if(l==1){
      SPECD <- SPECDA
    }
    
    for(l in 2){
      n <- prod(dimvec)
      y <- ty[,l]
      speccur <- SPECD
      if(l>1)
      {speccur <- SPECDC[,l-1]}
      SPECDA <- speccur
      thetaA <- theta[l, ]
      cant     <- thetaA[-1] + rnorm(3,sd = sdl2) #MH[2]*tCthetaA%*%
      cansd    <- spec_d_matern(h.omega,c(0,cant),dimvec, whichreal = realind)
      bb       <- sum(y/cansd)/2+.1
      cant1    <- -log(rgamma(1,n/2+.1,bb))
      cant     <- c(cant1,cant)
      cansd    <- exp(cant1)*cansd
      BB       <- exp(thetaA[1])*sum(y/speccur)/2+.1
      
      curll    <- -0.5*sum(log(speccur)) - 0.5*sum(y/speccur)+
        sum(dnorm(thetaA[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
        dgamma(exp(-thetaA[1]),.1,.1,log=TRUE)
      canll    <- -0.5*sum(log(cansd)) - 0.5*sum(y/cansd)+
        sum(dnorm(cant[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
        dgamma(exp(-cant[1]),.1,.1,log=TRUE)
      Q1       <- dgamma(exp(-thetaA[1]),n/2+.1,BB,log=TRUE)
      Q2       <- dgamma(exp(-cant[1]),n/2+.1,bb,log=TRUE)
      R        <- canll-curll+Q1-Q2
      if(!is.na(R)){if(log(runif(1))< R){
        flag[l] <- 1
        theta[l,]  <- cant
        SPECDA  <- cansd
        curll <- canll
      }}
      if(l>1){
        SPECDC[,l-1] <- SPECDA
      }
    }
    for(l in 3){
      n <- prod(dimvec)
      y <- ty[,l]
      speccur <- SPECD
      if(l>1)
      {speccur <- SPECDC[,l-1]}
      SPECDA <- speccur
      thetaA <- theta[l, ]
      cant     <- thetaA[-1] + rnorm(3,sd = sdl3) #MH[2]*tCthetaA%*%
      cansd    <- spec_d_matern(h.omega,c(0,cant),dimvec,  whichreal = realind)
      bb       <- sum(y/cansd)/2+.1
      cant1    <- -log(rgamma(1,n/2+.1,bb))
      cant     <- c(cant1,cant)
      cansd    <- exp(cant1)*cansd
      BB       <- exp(thetaA[1])*sum(y/speccur)/2+.1
      
      curll    <- -0.5*sum(log(speccur)) - 0.5*sum(y/speccur)+
        sum(dnorm(thetaA[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
        dgamma(exp(-thetaA[1]),.1,.1,log=TRUE)
      canll    <- -0.5*sum(log(cansd)) - 0.5*sum(y/cansd)+
        sum(dnorm(cant[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
        dgamma(exp(-cant[1]),.1,.1,log=TRUE)
      Q1       <- dgamma(exp(-thetaA[1]),n/2+.1,BB,log=TRUE)
      Q2       <- dgamma(exp(-cant[1]),n/2+.1,bb,log=TRUE)
      R        <- canll-curll+Q1-Q2
      if(!is.na(R)){if(log(runif(1))< R){
        flag[l] <- 1
        theta[l,]  <- cant
        SPECDA  <- cansd
        curll <- canll
      }}
      if(l>1){
        SPECDC[,l-1] <- SPECDA
      }
    }
    
    Betac <- array(convert1(array(Beta[,2], dimvec), gout, inv = T))
    b11 <- array((Betac), dimvec)
    image.plot(b11[,,1])
    
    #Betac <- matrix(rnorm(2*nrow(Y)),ncol=2)
    
    #Beta[, 1] <- rnorm(nrow(Y), sd = SPECDC[,1])
    #Beta[, 2] <- array(convert1(array(Betac[,1]*Betac[,2],dimvec), gout))
    
    #can1 <- Beta[,2]/2
    #can2 <- Beta[,2]*2
    
    theta_p[[itr]] <- theta
    
    #acceptedthetano <- mean(flag)
    
    acceptedthetano_p[itr, ] <- flag
    
    if(itr>1500)
      Beta2_p[[itr-1500]] <- Beta[, 2]
  }
  
  if(itr==2000){
    Beta2 <- Reduce('+', Beta2_p)/500
    Betac <- array(convert1(array(Beta2, dimvec), gout, inv = T))
    #Beta1p <- Reduce('+', Beta1_p[501:1000])/500
    #resid <- array(convert1(array(Y[,10]-Beta1p, dimvec), gout, inv = T))
    temp1 <- sign(Betac)*(abs(Betac))^(1/5)
    Betac <- cbind(temp1, temp1, temp1, temp1, temp1)
    theta[3, 1] <- (theta[3, 1])
    theta[4, ] <- theta[3, ]
    theta[5, ] <- theta[3, ]
    theta[6, ] <- theta[3, ]
    theta[7, ] <- theta[3, ]
    sdl3 <- 0.01
    SPECDC <- mcmapply(2:(7), FUN=function(i){spec_d_matern(h.omega,theta[i, ],dimvec, whichreal = realind)})
  }
  
  
  if(itr>2000){
    
    for(i in 1:5){
      Beta[, 1] <- array(updateBeta(1, fixedy, SPECDC[,1], Beta, SPECD, X, B, n))
      Beta1sp <- array(convert1(array(Beta[, 1],dimvec), gout, inv = T))
      Ysp <- apply(Y, 2 ,function(x){convert1(array(x,dimvec), gout, inv = T)})
      Yf <- apply(Ysp, 2 ,function(x){convert1(array(x/rowProds(Betac[, -i]),dimvec), gout)})
      Beta1f <- array(convert1(array(Beta1sp/rowProds(Betac[, -i]),dimvec), gout))
      
      esp <- apply(Y-Beta%*%t(X), 2 ,function(x){convert1(array(x,dimvec), gout, inv = T)})  
      #- matrix(rowProds(Betac[,-i]),n,m)*matrix(X[,2], n,m,byrow=T)*matrix(tempsp, n, m)
      espf <- apply(esp, 2 ,function(x){convert1(array((x/rowProds(Betac[, -i])-x),dimvec), gout)})
      inter <- Yf - matrix(Beta1f, nrow = n, ncol= m) - espf
      
      # Beta.mean1 <- rowSums((inter)*matrix(X[,2], n,m,byrow=T)) / SPECD
      # #Beta.mean2 <- B %*% beta[, l] /SPECDC[, l]
      # Beta.var <- array(sum(X[,2]^2) / SPECD + 1 / SPECDC[, i+1])
      # temp <- array(Beta.mean1)/Beta.var #+ array(Beta.mean2)
      # tempsp <- array(convert1(array(temp, dimvec), gout, inv=T))
      #
      # esp <- Ysp - cbind(Beta1sp, rowProds(Betac)) %*% t(X) - matrix(rowProds(Betac[,-i]),n,m)*matrix(X[,2], n,m,byrow=T)*matrix(tempsp-Betac[, i], n, m)
      # espf <- apply(esp, 2 ,function(x){convert1(array((x/rowProds(Betac[, -i])-x),dimvec), gout)})
      # inter <- Yf - matrix(Beta1f, nrow = n, ncol= 10) - espf
      
      #temp <- temp + rnorm(1, 0, sd=sd)
      gen <- array(updateBeta(2, inter, SPECDC[,i+1], matrix(0, n, 2), SPECD, X, B, n))
      Q1 <- U(Betac[, i], rowProds(Betac[, -i]),SPECDC[,i+1], Beta[,1],SPECD, fixedy)
      up <- array(convert1(array(gen, dimvec), gout, inv = T))
      prop <- sd
      nrm <- as.numeric(sqrt(crossprod(up - Betac[, i])))
      up <- Betac[, i] + min(nrm, sd)*(up - Betac[, i])/nrm
      Q2 <- U(up, rowProds(Betac[, -i]),SPECDC[,i+1], Beta[,1],SPECD, fixedy)
      #Betac[, i] <- up
      if(rbinom(1,1,min(exp(Q1-Q2),1))){
        c=c+1
        Betac[, i] <- up
      }
      
      
      Beta[, 2] <- array(convert1(array(rowProds(Betac), dimvec), gout))
    }
    #c <- c+(c1>0)
    if(itr%%100==0){
      it=itr-2000
      print(c/(5*it))
      if(c/(5*it) > 0.70){sd <- sd * (1+5e-1)}
      if(c/(5*it) < 0.55){sd <- sd * (1-5e-1)}
      # if(sd>1){sd=1}
    }      
    Beta_p[[itr-2000]] <- Beta
    Betac_p[[itr-2000]] <- Betac
    #b11 <- array(rowProds(Betac), dimvec)
    #image.plot(b11[,14,])
    
    #beta <- mcmapply(1:ncol(X), FUN=updatebeta, MoreArgs = list(Beta))
    
    residy <- rowSums((fixedy - Beta%*%t(X))^2)
    
    ty <- cbind(residy, Beta[,1]^2, array(convert1(array(Betac[,1],dimvec), gout))^2, array(convert1(array(Betac[,2],dimvec), gout))^2
                ,array(convert1(array(Betac[,3],dimvec), gout))^2,array(convert1(array(Betac[,4],dimvec), gout))^2,array(convert1(array(Betac[,5],dimvec), gout))^2)
    
    l=1
    n <- prod(dimvec)
    y <- ty[,l]
    speccur <- SPECD
    SPECDA <- speccur
    thetaA <- theta[l, ]
    cant     <- thetaA[-1] + rnorm(3,sd = sdl1) #MH[2]*tCthetaA%*%
    cansd    <- spec_d_matern(h.omega,c(0,cant),dimvec, whichreal = realind)
    bb       <- sum(y/cansd)/2+.1
    cant1    <- -log(rgamma(1,M*n/2+.1,bb))
    cant     <- c(cant1,cant)
    cansd    <- exp(cant1)*cansd
    BB       <- exp(thetaA[1])*sum(y/speccur)/2+.1
    
    curll    <- -0.5*M*sum(log(speccur)) - 0.5*sum(y/speccur)+
      sum(dnorm(thetaA[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
      dgamma(exp(-thetaA[1]),.1,.1,log=TRUE)
    canll    <- -0.5*M*sum(log(cansd)) - 0.5*sum(y/cansd)+
      sum(dnorm(cant[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+ 
      dgamma(exp(-cant[1]),.1,.1,log=TRUE)
    Q1       <- dgamma(exp(-thetaA[1]),M*n/2+.1,BB,log=TRUE)
    Q2       <- dgamma(exp(-cant[1]),M*n/2+.1,bb,log=TRUE)
    R        <- canll-curll+Q1-Q2
    if(!is.na(R)){if(log(runif(1))< R){
      flag[l] <- 1
      theta[l,]  <- cant
      SPECDA  <- cansd
      curll <- canll
    }}
    if(l>1){
      SPECDC[,l-1] <- SPECDA
    }
    if(l==1){
      SPECD <- SPECDA
    }
    
    l=2
    n <- prod(dimvec)
    y <- ty[,l]
    speccur <- SPECD
    if(l>1)
    {speccur <- SPECDC[,l-1]}
    SPECDA <- speccur
    thetaA <- theta[l, ]
    cant     <- thetaA[-1] + rnorm(3,sd = sdl2) #MH[2]*tCthetaA%*%
    cansd    <- spec_d_matern(h.omega,c(0,cant),dimvec, whichreal = realind)
    bb       <- sum(y/cansd)/2+.1
    cant1    <- -log(rgamma(1,n/2+.1,bb))
    cant     <- c(cant1,cant)
    cansd    <- exp(cant1)*cansd
    BB       <- exp(thetaA[1])*sum(y/speccur)/2+.1
    
    curll    <- -0.5*sum(log(speccur)) - 0.5*sum(y/speccur)+
      sum(dnorm(thetaA[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
      dgamma(exp(-thetaA[1]),.1,.1,log=TRUE)
    canll    <- -0.5*sum(log(cansd)) - 0.5*sum(y/cansd)+
      sum(dnorm(cant[-1],pri.mn[-1],pri.sd[-1],log=TRUE))+
      dgamma(exp(-cant[1]),.1,.1,log=TRUE)
    Q1       <- dgamma(exp(-thetaA[1]),n/2+.1,BB,log=TRUE)
    Q2       <- dgamma(exp(-cant[1]),n/2+.1,bb,log=TRUE)
    R        <- canll-curll+Q1-Q2
    if(!is.na(R)){if(log(runif(1))< R){
      flag[l] <- 1
      theta[l,]  <- cant
      SPECDA  <- cansd
      curll <- canll
    }}
    if(l>1){
      SPECDC[,l-1] <- SPECDA
    }
    
    l=3
    n <- prod(dimvec)
    y <- ty[,3:7]
    speccur <- SPECDC[,2:6]
    thetaA <- theta[3:7, ]
    cantc     <- thetaA[1,-1] + rnorm(3,sd = sdl3) #MH[2]*tCthetaA%*%
    cansd1 <- matrix(0, n, 5)
    cansd    <- spec_d_matern(h.omega,c(0,cantc),dimvec, whichreal = realind)
    #cansd2    <- spec_d_matern(h.omega,c(0,cant),dimvec, nugget = F, whichreal = realind)
    bb1       <- sum(y[,1]/cansd)/2+.1
    bb2       <- sum(y[,2]/cansd)/2+.1
    bb3       <- sum(y[,3]/cansd)/2+.1
    bb4       <- sum(y[,4]/cansd)/2+.1
    bb5       <- sum(y[,5]/cansd)/2+.1
    
    cant1    <- -log(rgamma(1,n/2+.1,bb1))
    cant11     <- c(cant1,cantc)
    cansd1[, 1]    <- exp(cant1)*cansd
    
    cant2    <- 0#-log(rgamma(1,n/2+.1,bb2))
    cant22     <- c(cant2,cantc)
    cansd1[, 2]    <- exp(cant2)*cansd
    
    cant3    <- 0#-log(rgamma(1,n/2+.1,bb3))
    cant33     <- c(cant3,cantc)
    cansd1[, 3]    <- exp(cant3)*cansd
    
    cant4    <- 0#-log(rgamma(1,n/2+.1,bb4))
    cant44     <- c(cant4,cantc)
    cansd1[, 4]    <- exp(cant4)*cansd
    
    cant5    <- 0#-log(rgamma(1,n/2+.1,bb5))
    cant55     <- c(cant5,cantc)
    cansd1[, 5]    <- exp(cant5)*cansd
    
    BB1      <- exp(thetaA[1, 1])*sum(y[,1]/speccur[,1])/2+.1
    BB2      <- exp(thetaA[2, 1])*sum(y[,2]/speccur[,2])/2+.1
    BB3      <- exp(thetaA[3, 1])*sum(y[,3]/speccur[,3])/2+.1
    BB4      <- exp(thetaA[4, 1])*sum(y[,3]/speccur[,3])/2+.1
    BB5      <- exp(thetaA[5, 1])*sum(y[,3]/speccur[,3])/2+.1
    
    curll1    <- -0.5*sum(log(speccur[,1])) - 0.5*sum(y[,1]/speccur[,1])+
      sum(dnorm(thetaA[1,-c(1)],pri.mn[-c(1)],pri.sd[-c(1)],log=TRUE))+
      dgamma(exp(-thetaA[1,1]),.1,.1,log=TRUE)
    
    curll2    <- -0.5*sum(log(speccur[,2])) - 0.5*sum(y[,2]/speccur[,2])+
      dgamma(exp(-thetaA[2,1]),.1,.1,log=TRUE)
    
    curll3    <- -0.5*sum(log(speccur[,3])) - 0.5*sum(y[,3]/speccur[,3])+
      dgamma(exp(-thetaA[3,1]),.1,.1,log=TRUE)
    
    curll4    <- -0.5*sum(log(speccur[,4])) - 0.5*sum(y[,4]/speccur[,4])+
      dgamma(exp(-thetaA[4,1]),.1,.1,log=TRUE)
    
    curll5    <- -0.5*sum(log(speccur[,5])) - 0.5*sum(y[,5]/speccur[,5])+
      dgamma(exp(-thetaA[5,1]),.1,.1,log=TRUE)
    
    canll1    <- -0.5*sum(log(cansd1[,1])) - 0.5*sum(y[,1]/cansd1[,1])+
      sum(dnorm(cant11[-c(1)],pri.mn[-c(1)],pri.sd[-c(1)],log=TRUE))+
      dgamma(exp(-cant1),.1,.1,log=TRUE)
    
    canll2    <- -0.5*sum(log(cansd1[,2])) - 0.5*sum(y[,2]/cansd1[,2])+
      dgamma(exp(-cant2),.1,.1,log=TRUE)
    
    canll3    <- -0.5*sum(log(cansd1[,3])) - 0.5*sum(y[,3]/cansd1[,3])+
      dgamma(exp(-cant3),.1,.1,log=TRUE)
    
    canll4    <- -0.5*sum(log(cansd1[,4])) - 0.5*sum(y[,4]/cansd1[,4])+
      dgamma(exp(-cant4),.1,.1,log=TRUE)
    
    canll5    <- -0.5*sum(log(cansd1[,5])) - 0.5*sum(y[,5]/cansd1[,5])+
      dgamma(exp(-cant5),.1,.1,log=TRUE)
    
    Q1       <- (dgamma(exp(-thetaA[1,1]),n/2+.1,BB1,log=TRUE))#+(dgamma(exp(-thetaA[2,1]),n/2+.1,BB2,log=TRUE))+(dgamma(exp(-thetaA[3,1]),n/2+.1,BB3,log=TRUE))+
    #(dgamma(exp(-thetaA[4,1]),n/2+.1,BB4,log=TRUE))+(dgamma(exp(-thetaA[5,1]),n/2+.1,BB5,log=TRUE))
    
    Q2       <- dgamma(exp(-cant1),n/2+.1,bb1,log=TRUE)# +dgamma(exp(-cant2),n/2+.1,bb2,log=TRUE)+dgamma(exp(-cant3),n/2+.1,bb3,log=TRUE)+
    # dgamma(exp(-cant4),n/2+.1,bb4,log=TRUE) + dgamma(exp(-cant5),n/2+.1,bb5,log=TRUE)
    
    R        <- canll1-curll1+canll2-curll2+canll3-curll3+canll4-curll4+canll5-curll5+Q1-Q2
    if(!is.na(R)){if(log(runif(1))< R){
      flag[l] <- 1
      theta[3, ]  <- cant11
      theta[4, ] <- cant22
      theta[5, ] <- cant33
      theta[6, ] <- cant44
      theta[7, ] <- cant55
      SPECDC[,2]  <- cansd1[,1]
      SPECDC[,3]  <- cansd1[,2]
      SPECDC[,4]  <-  cansd1[,3]
      SPECDC[,5]  <-  cansd1[,4]
      SPECDC[,6]  <-  cansd1[,5]
    }}    
    
  }
  theta_p[[itr]] <- theta
  
  #acceptedthetano <- mean(flag)
  
  acceptedthetano_p[itr, ] <- flag
  
  specd <- SPECD
  
  specd[realind] <- 2*specd[realind]
  
  # Betaac <- apply(Beta, 2, function(x){convert1(array(x,dimvec), gout, inv = T)})
  # 
  # 
  # mu <- tcrossprod(Betaac,(X)) #+ Br %*% MU_b
  # 
  # if(itr%%10==0){
  #   for(t in 1:M){
  #     Y[,t] <- array(impute(Y.o[,t], mu[,t],miss=miss1[,t],specd=specd, dimvec, tol=tol,gout=gout,maxiter=50))
  #   }
  # }
  # 
  b11 <- array(rowProds(Betac), dimvec)
  image.plot(b11[,12,])
  
  if(itr%%100 == 0){
    if(mean(acceptedthetano_p[1:itr,1]) > 0.45){sdl1 <- sdl1*1.2}
    if(mean(acceptedthetano_p[1:itr,1]) < 0.15){sdl1 <- sdl1*0.8}
    if(mean(acceptedthetano_p[1:itr,2]) > 0.45){sdl2 <- sdl2*1.2}
    if(mean(acceptedthetano_p[1:itr,2]) < 0.15){sdl2 <- sdl2*0.8}
    if(itr <= 2000){
      if(mean(acceptedthetano_p[1:itr,3]) > 0.45){sdl3 <- sdl3*1.2}
      if(mean(acceptedthetano_p[1:itr,3]) < 0.3){sdl3 <- sdl3*0.8}
    }
    
    if(itr >2000){
      if(mean(acceptedthetano_p[2001:itr,3]) > 0.4){sdl3 <- sdl3*1.2}
      if(mean(acceptedthetano_p[2001:itr,3]) < 0.2){sdl3 <- sdl3*0.8}
    }
    # if(itr>1000){
    #   it <- itr-1000
    #   if(flagbeta2/(2*it) > 0.70){sd <- sd * 2}
    #   if(flagbeta2/(2*it) < 0.55){sd <- sd / 2}
    #   print(flagbeta2/(2*it))
    # }
  }
  #print(itr)
}
Betac_post<-Betac_p[(1:5000)*2]