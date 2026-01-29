###Test case to verify the likelihoods and log-posteriors
###Let the four rivers be point masses at 5,10,15, and 20 MYA, respectively.
###Let the mixing coefficients be 0.25 for all loess horizons
###Let the horizon ages sample equally from the four modes of F
###Let the obvious class assignments and class parameters be assumed
set.seed(123)
n_x <- 20
X <- list(rep(5,n_x),rep(10,n_x),rep(15,n_x),rep(20,n_x))
x <- seq(0,25,,1e3+1)
F <- 0.25/n_x*dnorm(x,5,0.01)+0.25/n_x*dnorm(x,10,0.01)+0.25/n_x*dnorm(x,15,0.01)+0.25/n_x*dnorm(x,20,0.01)
Y <- list(c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),
          c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)))
delta <- list(rep(0.25,4),rep(0.25,4),rep(0.25,4),rep(0.25,4))
C <- list(rep(1,n_x),rep(1,n_x),rep(1,n_x),rep(1,n_x))
phi <- list(matrix(c(5,1),2,1),matrix(c(10,1),2,1),matrix(c(15,1),2,1),matrix(c(20,1),2,1))

nu <- 1
mu <- 0
tau <- 0.5
a <- 3
b <- 1

.pi_computer <- function(l,j,pro,num){
  .phi <- if(num) phi[[l]][,j] else pro
  phi[[l]][,j] <- if(num) pro else phi[[l]][,j]
  .pi_1 <- prod((2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(X[[1]]-phi[[1]][1,1])^2/2/phi[[1]][2,1]))*
    prod((2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(X[[2]]-phi[[2]][1,1])^2/2/phi[[2]][2,1]))*
    prod((2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(X[[3]]-phi[[3]][1,1])^2/2/phi[[3]][2,1]))*
    prod((2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(X[[4]]-phi[[4]][1,1])^2/2/phi[[4]][2,1]))
  .pi_2 <- prod(delta[[1]][1]/n_x*(2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(Y[[1]]-phi[[1]][1,])^2/2/phi[[1]][2,1])+
                delta[[1]][2]/n_x*(2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(Y[[1]]-phi[[2]][1,])^2/2/phi[[2]][2,1])+
                delta[[1]][3]/n_x*(2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(Y[[1]]-phi[[3]][1,])^2/2/phi[[3]][2,1])+
                delta[[1]][3]/n_x*(2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(Y[[1]]-phi[[4]][1,])^2/2/phi[[4]][2,1]))*
    prod(delta[[2]][1]/n_x*(2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(Y[[2]]-phi[[1]][1,])^2/2/phi[[1]][2,1])+
        delta[[2]][2]/n_x*(2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(Y[[2]]-phi[[2]][1,])^2/2/phi[[2]][2,1])+
        delta[[2]][3]/n_x*(2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(Y[[2]]-phi[[3]][1,])^2/2/phi[[3]][2,1])+
        delta[[2]][4]/n_x*(2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(Y[[2]]-phi[[4]][1,])^2/2/phi[[4]][2,1]))*
    prod(delta[[3]][1]/n_x*(2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(Y[[3]]-phi[[1]][1,])^2/2/phi[[1]][2,1])+
          delta[[3]][2]/n_x*(2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(Y[[3]]-phi[[2]][1,])^2/2/phi[[2]][2,1])+
          delta[[3]][3]/n_x*(2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(Y[[3]]-phi[[3]][1,])^2/2/phi[[3]][2,1])+
          delta[[3]][4]/n_x*(2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(Y[[3]]-phi[[4]][1,])^2/2/phi[[4]][2,1]))*
    prod(delta[[4]][1]/n_x*(2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(Y[[4]]-phi[[1]][1,])^2/2/phi[[1]][2,1])+
          delta[[4]][2]/n_x*(2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(Y[[4]]-phi[[2]][1,])^2/2/phi[[2]][2,1])+
          delta[[4]][3]/n_x*(2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(Y[[4]]-phi[[3]][1,])^2/2/phi[[3]][2,1])+
          delta[[4]][4]/n_x*(2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(Y[[4]]-phi[[4]][1,])^2/2/phi[[4]][2,1]))
  .pi_3 <- (2*pi*phi[[1]][2,1]/nu)^(-0.5)*exp(-(phi[[1]][1,1]-mu)^2/(2*phi[[1]][2,1]/nu))
  .pi_4 <- b^a/gamma(a)*(phi[[1]][2,1])^(-a-1)*exp(-b/phi[[1]][2,1])
  .pi_1*.pi_2*.pi_3*.pi_4/((2*pi*phi[[l]][2,j]/nu)^(-0.5)*exp(-0.5*(phi[[l]][1,j]-.phi[1])^2/(phi[[l]][2,j]/nu))*
                              (.phi[2]*(.phi[2]^2/tau+1))^(.phi[2]^2/tau+2)/gamma(.phi[2]^2/tau+2)*
                              phi[[l]][2,j]^(-(.phi[2]^2/tau+2)-1)*exp(-.phi[2]*(.phi[2]^2/tau+1)/phi[[1]][2,j]))
}
l <- 1
j <- 1
.pro <- vector("numeric",2) #Draw...
.pro[2] <- rinvgamma(1,phi[[l]][2,j]^2/tau+2,rate=phi[[l]][2,j]*(phi[[l]][2,j]^2/tau+1)) #proposed...
.pro[1] <- rnorm(1,phi[[l]][1,j],sqrt(pro[2])) #proper class parameters.
.pi_computer(1,1,.pro,num=T)/.pi_computer(1,1,.pro,num=FALSE)


set.seed(123)
n_x <- 20
X <- list(rep(5,n_x),rep(10,n_x),rep(15,n_x),rep(20,n_x))
x <- seq(0,25,,1e3+1)
F <- 0.25/n_x*dnorm(x,5,0.01)+0.25/n_x*dnorm(x,10,0.01)+0.25/n_x*dnorm(x,15,0.01)+0.25/n_x*dnorm(x,20,0.01)
Y <- list(c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),
          c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)))
delta <- list(rep(0.25,4),rep(0.25,4),rep(0.25,4),rep(0.25,4))
C <- list(rep(1,n_x),rep(1,n_x),rep(1,n_x),rep(1,n_x))
phi <- list(matrix(c(5,1),2,1),matrix(c(10,1),2,1),matrix(c(15,1),2,1),matrix(c(20,1),2,1))

nu_phi <- 1
mu_phi <- 0
tau <- 0.5
a_phi <- 3
b_phi <- 1

l <- 1
#for(l in 1:4){
  k <- length(na.omit(unique(C[[l]])))
#  for(j in 1:k){
  j <- 1
    pro <- vector("numeric",2) #Draw...
    pro[2] <- rinvgamma(1,phi[[l]][2,j]^2/tau+2,rate=phi[[l]][2,j]*(phi[[l]][2,j]^2/tau+1)) #proposed...
    pro[1] <- rnorm(1,phi[[l]][1,j],sqrt(pro[2])) #proper class parameters.
    .pi_n <- loglikelihood_X(deltas=delta,phis=list_replace(phi,pro,l,j),C=C) +
      loglikelihood_Y(deltas=delta,phis=list_replace(phi,pro,l,j),C=C) +
      dnorm(pro[1],mu_phi,sqrt(pro[2]/nu_phi),log=T) +
      dinvgamma(pro[2],a_phi,rate=b_phi,log=T) #Calculate...
    .pi_d <- loglikelihood_X(deltas=delta,phis=phi,C=C) +
      loglikelihood_Y(deltas=delta,phis=phi,C=C) +
      dnorm(phi[[l]][1,j],mu_phi,sqrt(phi[[l]][2,j]/nu_phi),log=T) +
      dinvgamma(phi[[l]][2,j],a_phi,rate=b_phi,log=T) #acceptance ratio...
    .p <- exp(.pi_n - .pi_d + dnorm(phi[[l]][1,j],pro[1],sqrt(pro[2]/nu_phi),log=T) +
                dinvgamma(phi[[l]][2,j],pro[2]^2/tau+2,rate=pro[2]*(pro[2]^2/tau+1),log=T) -
                dnorm(pro[1],phi[[l]][1,j],sqrt(phi[[l]][2,j]/nu_phi),log=T) -
                dinvgamma(pro[2],phi[[l]][2,j]^2/tau+2,rate=phi[[l]][2,j]*(phi[[l]][2,j]^2/tau+1),log=T)) #and...
#    if(runif(1)<.p){ #accept proposed parameters...
#      phi[[l]][,j] <- pro #with that probability...
#      if(j %in% 1:3) acpt[(l-1)*3 + j] <- 1 #and...
#    } else { #track...
#      if(j %in% 1:3) acpt[(l-1)*3 + j] <- 0 #acceptance rate.
#    }
#  }
#  phi_prop[[l]] <- phi[[l]][,1:k,drop=FALSE] #Record changes to phi in phi_prop.
#}
.p

##########
#loglikelihood_Y is where the disagreement occurs

set.seed(123)
n_x <- 20
X <- list(rep(5,n_x),rep(10,n_x),rep(15,n_x),rep(20,n_x))
x <- seq(0,25,,1e3+1)
F <- 0.25/n_x*dnorm(x,5,0.01)+0.25/n_x*dnorm(x,10,0.01)+0.25/n_x*dnorm(x,15,0.01)+0.25/n_x*dnorm(x,20,0.01)
Y <- list(c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),
          c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)))
delta <- list(rep(0.25,4),rep(0.25,4),rep(0.25,4),rep(0.25,4))
C <- list(rep(1,n_x),rep(1,n_x),rep(1,n_x),rep(1,n_x))
phi <- list(matrix(c(5,1),2,1),matrix(c(10,1),2,1),matrix(c(15,1),2,1),matrix(c(20,1),2,1))

nu <- 1
mu2 <- 0
tau <- 0.5
a <- 3
b <- 1

l <- 1
j <- 1
num <- TRUE
.pro <- vector("numeric",2) #Draw...
.pro[2] <- rinvgamma(1,phi[[l]][2,j]^2/tau+2,rate=phi[[l]][2,j]*(phi[[l]][2,j]^2/tau+1)) #proposed...
.pro[1] <- rnorm(1,phi[[l]][1,j],sqrt(pro[2])) #proper class parameters.

.phi <- if(num) phi[[l]][,j] else pro
phi[[l]][,j] <- if(num) pro else phi[[l]][,j]

prod(delta[[1]][1]/n_x*(2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(Y[[1]]-phi[[1]][1,])^2/2/phi[[1]][2,1])+
       delta[[1]][2]/n_x*(2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(Y[[1]]-phi[[2]][1,])^2/2/phi[[2]][2,1])+
       delta[[1]][3]/n_x*(2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(Y[[1]]-phi[[3]][1,])^2/2/phi[[3]][2,1])+
       delta[[1]][3]/n_x*(2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(Y[[1]]-phi[[4]][1,])^2/2/phi[[4]][2,1]))*
  prod(delta[[2]][1]/n_x*(2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(Y[[2]]-phi[[1]][1,])^2/2/phi[[1]][2,1])+
         delta[[2]][2]/n_x*(2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(Y[[2]]-phi[[2]][1,])^2/2/phi[[2]][2,1])+
         delta[[2]][3]/n_x*(2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(Y[[2]]-phi[[3]][1,])^2/2/phi[[3]][2,1])+
         delta[[2]][4]/n_x*(2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(Y[[2]]-phi[[4]][1,])^2/2/phi[[4]][2,1]))*
  prod(delta[[3]][1]/n_x*(2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(Y[[3]]-phi[[1]][1,])^2/2/phi[[1]][2,1])+
         delta[[3]][2]/n_x*(2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(Y[[3]]-phi[[2]][1,])^2/2/phi[[2]][2,1])+
         delta[[3]][3]/n_x*(2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(Y[[3]]-phi[[3]][1,])^2/2/phi[[3]][2,1])+
         delta[[3]][4]/n_x*(2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(Y[[3]]-phi[[4]][1,])^2/2/phi[[4]][2,1]))*
  prod(delta[[4]][1]/n_x*(2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(Y[[4]]-phi[[1]][1,])^2/2/phi[[1]][2,1])+
         delta[[4]][2]/n_x*(2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(Y[[4]]-phi[[2]][1,])^2/2/phi[[2]][2,1])+
         delta[[4]][3]/n_x*(2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(Y[[4]]-phi[[3]][1,])^2/2/phi[[3]][2,1])+
         delta[[4]][4]/n_x*(2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(Y[[4]]-phi[[4]][1,])^2/2/phi[[4]][2,1]))
deltas <- delta
  cont <- vector("numeric",sum(ns <- unlist(lapply(Y,FUN=length))))
  for(i in 1:4){
    n_i <- length(Y[[i]])
    for(J in 1:n_i){
      Y_like <- 0
      for(l in 1:4){
        n_l <- length(X[[l]])
        sigma2 <- phis[[l]][2,unique(C[[l]])]
        mu <- phis[[l]][1,unique(C[[l]])]
        Y_like <- Y_like + deltas[[i]][l]/n_l*sum(dnorm(Y[[i]][J],mu,sqrt(sigma2)))
      }
      cont[ifelse(i>1,cumsum(ns)[i-1],0)+J] <- Y_like
    }
  }
  exp(sum(log(cont)))
exp(loglikelihood_Y(deltas=delta,phis=phi,C=C)) #mu and sigma2 defined as above, with unique()
##########
# These 3 matched with num=TRUE, l=1:4; J=1:20; i=1:4; j=1
# The first and third matched with num=FALSE

set.seed(123)
n_x <- 20
X <- list(rep(5,n_x),rep(10,n_x),rep(15,n_x),rep(20,n_x))
x <- seq(0,25,,1e3+1)
F <- 0.25/n_x*dnorm(x,5,0.01)+0.25/n_x*dnorm(x,10,0.01)+0.25/n_x*dnorm(x,15,0.01)+0.25/n_x*dnorm(x,20,0.01)
Y <- list(c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),
          c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)))
delta <- list(rep(0.25,4),rep(0.25,4),rep(0.25,4),rep(0.25,4))
C <- list(rep(1,n_x),rep(1,n_x),rep(1,n_x),rep(1,n_x))
phi <- list(matrix(c(5,1),2,1),matrix(c(10,1),2,1),matrix(c(15,1),2,1),matrix(c(20,1),2,1))

nu <- 1
mu <- 0
tau <- 0.5
a <- 3
b <- 1

.pi_computer <- function(l,j,pro,num){
  .phi <- if(num) phi[[l]][,j] else pro
  phi[[l]][,j] <- if(num) pro else phi[[l]][,j]
  .pi_1 <- prod((2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(X[[1]]-phi[[1]][1,1])^2/2/phi[[1]][2,1]))*
    prod((2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(X[[2]]-phi[[2]][1,1])^2/2/phi[[2]][2,1]))*
    prod((2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(X[[3]]-phi[[3]][1,1])^2/2/phi[[3]][2,1]))*
    prod((2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(X[[4]]-phi[[4]][1,1])^2/2/phi[[4]][2,1]))
  .pi_2 <- prod(delta[[1]][1]/n_x*(2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(Y[[1]]-phi[[1]][1,])^2/2/phi[[1]][2,1])+
                  delta[[1]][2]/n_x*(2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(Y[[1]]-phi[[2]][1,])^2/2/phi[[2]][2,1])+
                  delta[[1]][3]/n_x*(2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(Y[[1]]-phi[[3]][1,])^2/2/phi[[3]][2,1])+
                  delta[[1]][3]/n_x*(2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(Y[[1]]-phi[[4]][1,])^2/2/phi[[4]][2,1]))*
    prod(delta[[2]][1]/n_x*(2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(Y[[2]]-phi[[1]][1,])^2/2/phi[[1]][2,1])+
           delta[[2]][2]/n_x*(2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(Y[[2]]-phi[[2]][1,])^2/2/phi[[2]][2,1])+
           delta[[2]][3]/n_x*(2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(Y[[2]]-phi[[3]][1,])^2/2/phi[[3]][2,1])+
           delta[[2]][4]/n_x*(2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(Y[[2]]-phi[[4]][1,])^2/2/phi[[4]][2,1]))*
    prod(delta[[3]][1]/n_x*(2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(Y[[3]]-phi[[1]][1,])^2/2/phi[[1]][2,1])+
           delta[[3]][2]/n_x*(2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(Y[[3]]-phi[[2]][1,])^2/2/phi[[2]][2,1])+
           delta[[3]][3]/n_x*(2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(Y[[3]]-phi[[3]][1,])^2/2/phi[[3]][2,1])+
           delta[[3]][4]/n_x*(2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(Y[[3]]-phi[[4]][1,])^2/2/phi[[4]][2,1]))*
    prod(delta[[4]][1]/n_x*(2*pi*rep(phi[[1]][2,1],n_x))^(-0.5)*exp(-(Y[[4]]-phi[[1]][1,])^2/2/phi[[1]][2,1])+
           delta[[4]][2]/n_x*(2*pi*rep(phi[[2]][2,1],n_x))^(-0.5)*exp(-(Y[[4]]-phi[[2]][1,])^2/2/phi[[2]][2,1])+
           delta[[4]][3]/n_x*(2*pi*rep(phi[[3]][2,1],n_x))^(-0.5)*exp(-(Y[[4]]-phi[[3]][1,])^2/2/phi[[3]][2,1])+
           delta[[4]][4]/n_x*(2*pi*rep(phi[[4]][2,1],n_x))^(-0.5)*exp(-(Y[[4]]-phi[[4]][1,])^2/2/phi[[4]][2,1]))
  .pi_3 <- (2*pi*phi[[1]][2,1]/nu)^(-0.5)*exp(-(phi[[1]][1,1]-mu)^2/(2*phi[[1]][2,1]/nu))
  .pi_4 <- b^a/gamma(a)*(phi[[1]][2,1])^(-a-1)*exp(-b/phi[[1]][2,1])
  .pi_1*.pi_2*.pi_3*.pi_4/((2*pi*phi[[l]][2,j]/nu)^(-0.5)*exp(-0.5*(phi[[l]][1,j]-.phi[1])^2/(phi[[l]][2,j]/nu))*
                             (.phi[2]*(.phi[2]^2/tau+1))^(.phi[2]^2/tau+2)/gamma(.phi[2]^2/tau+2)*
                             phi[[l]][2,j]^(-(.phi[2]^2/tau+2)-1)*exp(-.phi[2]*(.phi[2]^2/tau+1)/phi[[1]][2,j]))
}
.pro <- vector("numeric",2) #Draw...
.pro[2] <- rinvgamma(1,phi[[l]][2,j]^2/tau+2,rate=phi[[l]][2,j]*(phi[[l]][2,j]^2/tau+1)) #proposed...
.pro[1] <- rnorm(1,phi[[l]][1,j],sqrt(pro[2])) #proper class parameters.
.pi_computer(1,1,.pro,num=T)/.pi_computer(1,1,.pro,num=FALSE)
#3.585396e-35

set.seed(123)
n_x <- 20
X <- list(rep(5,n_x),rep(10,n_x),rep(15,n_x),rep(20,n_x))
x <- seq(0,25,,1e3+1)
F <- 0.25/n_x*dnorm(x,5,0.01)+0.25/n_x*dnorm(x,10,0.01)+0.25/n_x*dnorm(x,15,0.01)+0.25/n_x*dnorm(x,20,0.01)
Y <- list(c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),
          c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)),c(rep(5,n_x/4),rep(10,n_x/4),rep(15,n_x/4),rep(20,n_x/4)))
delta <- list(rep(0.25,4),rep(0.25,4),rep(0.25,4),rep(0.25,4))
C <- list(rep(1,n_x),rep(1,n_x),rep(1,n_x),rep(1,n_x))
phi <- list(matrix(c(5,1),2,1),matrix(c(10,1),2,1),matrix(c(15,1),2,1),matrix(c(20,1),2,1))

nu_phi <- 1
mu_phi <- 0
tau <- 0.5
a_phi <- 3
b_phi <- 1

l <- 1
#for(l in 1:4){
  k <- length(na.omit(unique(C[[l]])))
j <- 1
#    for(j in 1:k){
    pro <- vector("numeric",2) #Draw...
    pro[2] <- rinvgamma(1,phi[[l]][2,j]^2/tau+2,rate=phi[[l]][2,j]*(phi[[l]][2,j]^2/tau+1)) #proposed...
    pro[1] <- rnorm(1,phi[[l]][1,j],sqrt(pro[2])) #proper class parameters.
    .pi_n <- loglikelihood_X(deltas=delta,phis=list_replace(phi,pro,l,j),C=C) +
      loglikelihood_Y(deltas=delta,phis=list_replace(phi,pro,l,j),C=C) +
      dnorm(pro[1],mu_phi,sqrt(pro[2]/nu_phi),log=T) +
      dinvgamma(pro[2],a_phi,rate=b_phi,log=T) #Calculate...
    .pi_d <- loglikelihood_X(deltas=delta,phis=phi,C=C) +
      loglikelihood_Y(deltas=delta,phis=phi,C=C) +
      dnorm(phi[[l]][1,j],mu_phi,sqrt(phi[[l]][2,j]/nu_phi),log=T) +
      dinvgamma(phi[[l]][2,j],a_phi,rate=b_phi,log=T) #acceptance ratio...
    .p <- exp(.pi_n - .pi_d + dnorm(phi[[l]][1,j],pro[1],sqrt(phi[[l]][2,j]/nu_phi),log=T) +
                dinvgamma(phi[[l]][2,j],pro[2]^2/tau+2,rate=pro[2]*(pro[2]^2/tau+1),log=T) -
                dnorm(pro[1],phi[[l]][1,j],sqrt(pro[2]/nu_phi),log=T) -
                dinvgamma(pro[2],phi[[l]][2,j]^2/tau+2,rate=phi[[l]][2,j]*(phi[[l]][2,j]^2/tau+1),log=T)) #and...
#    if(runif(1)<.p){ #accept proposed parameters...
#      phi[[l]][,j] <- pro #with that probability...
#      if(j %in% 1:3) acpt[(l-1)*3 + j] <- 1 #and...
#    } else { #track...
#      if(j %in% 1:3) acpt[(l-1)*3 + j] <- 0 #acceptance rate.
#    }
#  }
#  phi_prop[[l]] <- phi[[l]][,1:k,drop=FALSE] #Record changes to phi in phi_prop.
#}
.p
