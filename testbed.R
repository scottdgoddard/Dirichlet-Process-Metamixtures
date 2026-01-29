#This testbed is designed to allow for global testing of the whole algorithm using data
#designed to mimic the real data, but in a simpler way.
rm(list=ls(all=TRUE))
set.seed(123)
source("mcmc.R")
delta <- list("B"=rep(1/4,4),"C"=c(1/3,1/3,1/6,1/6),
              "E"=c(1/2,1/4,1/8,1/8),"I"=c(1/4,1/3,1/4,1/6))
del_mat <- matrix(c(delta[[1]],delta[[2]],delta[[3]],delta[[4]]),4,4,byrow=T)
X <- list(c(rnorm(33,-4,1),rnorm(33,0,1),rnorm(34,4,1)),c(rnorm(33,-6,1),rnorm(33,0,1),rnorm(34,6,1)),
          c(rnorm(33,-2,1),rnorm(33,0,1),rnorm(34,2,1)),c(rnorm(10,-4,1),rnorm(80,0,1),rnorm(10,4,1)))
X_mat <- matrix(c(X[[1]],X[[2]],X[[3]],X[[4]]),4,100,byrow=TRUE)
F <- del_mat %*% X_mat
Y <- list(quantile(F[1,],runif(100)),quantile(F[2,],runif(100)),
          quantile(F[3,],runif(100)),quantile(F[4,],runif(100)))

freqs <- function(x) {
  out <- NULL
  for(j in 1:length(na.omit(unique(x)))) {
    out <- c(out,sum(x==j))
  }
  out
}

list_replace <- function(x,new,l,j=NULL) {#either replaces the jth column of the lth matrix in list x (if j specified)
  if(is.null(j)) x[[l]] <- new #or the entire lth vector (j=NULL)
  else x[[l]][,j] <- new 
  x
}

tightener <- function(Cee,phee,jay=NULL) {#New lines #tightens up C as needed and phi congruently, optionally removing the jth element beforehand
  if(!is.null(jay)) {
    phee <- phee[,sort(unique(Cee)),drop=FALSE]
    if((all(diff(sort(unique(Cee))) == 1) & all(Cee[jay] != Cee[-jay]) & any(Cee[jay] < Cee[-jay])) | 
       (any(diff(sort(unique(Cee))) != 1) & all(Cee[jay] != Cee[-jay]) & any(Cee[jay] < Cee[-jay]))) {
      elim <- setdiff(as.numeric(as.factor(Cee[-jay])),Cee[-jay])
      phee <- cbind(phee[,-(Cee[jay]-max(as.numeric(Cee[jay]>elim)))],phee[,Cee[jay]-max(as.numeric(Cee[jay]>elim))])
      Cee[-jay] <- as.numeric(as.factor(Cee[-jay]))
      Cee[jay] <- max(Cee[-jay]) + 1
    } else {
      Cee <- as.numeric(as.factor(Cee))
    }
  } else {
    phee <- phee[,sort(unique(Cee))]
    Cee <- as.numeric(as.factor(Cee))
  }
  list(Cee,phee)
}

loglikelihood_Y <- function(deltas,phis,C){
  cont <- vector("numeric",sum(ns <- unlist(lapply(Y,FUN=length))))
  for(i in 1:4){
    n_i <- length(Y[[i]])
    for(J in 1:n_i){
      Y_like <- 0
      for(l in 1:4){ #sums over rivers
        n_l <- length(X[[l]])
        sigma2 <- phis[[l]][2,unique(C[[l]])]
        mu <- phis[[l]][1,unique(C[[l]])]
        Y_like <- Y_like + deltas[[i]][l]/n_l*sum(dnorm(Y[[i]][J],mu,sqrt(sigma2))) #sums over river grains
      }
      cont[ifelse(i>1,cumsum(ns)[i-1],0)+J] <- Y_like
    }
  }
  sum(log(cont)) #sums over horizon grains and horizons
}

loglikelihood_X <- function(deltas,phis,C){
  cont <- vector("numeric",4)
  for(l in 1:4){
    sigma2 <- phis[[l]][2,C[[l]]] #should there also be a unique() here as in loglikelihood_Y?
    mu <- phis[[l]][1,C[[l]]] #should there also be a unique() here as in loglikelihood_Y?
    cont[l] <- sum(dnorm(X[[l]],mu,sqrt(sigma2),log=TRUE))
  }
  sum(cont)
}

#Initialization
B <- 1e3
m=2
g0Priors=c(0,1,3,1)
MPrior=c(2,4)
deltaPrior=c(0.5,0.5,0.5,0.5)
tuning=c(2,50)
mu_phi <- g0Priors[1]
nu_phi <- g0Priors[2]
a_phi <- g0Priors[3]
b_phi <- g0Priors[4]
a_M <- MPrior[1]
b_M <- MPrior[2]
alpha_1 <- deltaPrior[1]
alpha_2 <- deltaPrior[2]
alpha_3 <- deltaPrior[3]
alpha_4 <- deltaPrior[4]
tau <- tuning[1]
lambda <- tuning[2]
n_x <- unlist(lapply(X,FUN=\(x) length(x)))
n_y <- unlist(lapply(Y,FUN=\(x) length(x)))
.C1 <- c(rep(1,n_x[1]))
.C2 <- c(rep(1,n_x[2]))
.C3 <- c(rep(1,n_x[3]))
.C4 <- c(rep(1,n_x[4]))
C <- list(.C1,.C2,.C3,.C4) #Each element of C is a map of the components in that river's grains
.k <- lapply(C,FUN=\(x) length(na.omit(unique(x))))
phi_prop <- list(rbind(rep(mu_phi,.k[1]),rep((b_phi/(a_phi-1))/nu_phi,.k[1])),rbind(rep(mu_phi,.k[2]),rep((b_phi/(a_phi-1))/nu_phi,.k[2])),
                 rbind(rep(mu_phi,.k[3]),rep((b_phi/(a_phi-1))/nu_phi,.k[3])),rbind(rep(mu_phi,.k[4]),rep((b_phi/(a_phi-1))/nu_phi,.k[4])))
#Each element of phi_prop is a map of the parameters in that river's components
phi <- list(matrix(c(0,0.5),2,1),matrix(c(0,0.5),2,1),matrix(c(0,0.5),2,1),matrix(c(0,0.5),2,1))
delta <- list("B"=c(1/4,1/4,1/4,1/4),"C"=c(1/4,1/4,1/4,1/4),#Nenana,Delta,
              "E"=c(1/4,1/4,1/4,1/4),"I"=c(1/4,1/4,1/4,1/4))#Tanana,Yukon
phi_chain <- list()
C_chain <- list()
M_chain <- list()
delta_chain <- list()
acpt_chain <- list()
.ps <- list() #This can be deleted
acpt <- vector("numeric",16) #paste0("l=1 phi[,",1:3,"]"),paste0("l=2 phi[,",1:3,"]"),
#paste0("l=3 phi[,",1:3,"]"),paste0("l=4 phi[,",1:3,"]"),"l=1 delta","l=2 delta","l=3 delta","l=4 delta")
M <- 1
progbar <- TRUE
pb <- txtProgressBar(min = 0,max=B,style = 3)

#Algorithm
for(b in 1:B) {
  
  #Step 1--draw new auxiliary mus and sigma^2s and new class labels
  for(l in 1:4) {
    for(j in 1:n_x[l]) {
      k_ <- length(na.omit(unique(C[[l]][-j])))
      h <- k_ + m
      .m <- ifelse(C[[l]][j] %in% C[[l]][-j],m,m-1) #You can't tighten up class labels in .C before this conditional or it will break.
      ttnd <- tightener(C[[l]],phi[[l]],j)
      C[[l]] <- ttnd[[1]]
      phi[[l]] <- ttnd[[2]]
      .k <- length(na.omit(unique(C[[l]]))) #Record changes made to...
      phi_prop[[l]] <- phi[[l]][,1:.k,drop=FALSE] #phi in phi_prop.
      .phi_sigs <- rinvgamma(.m,a_phi,rate=b_phi) #Draw new...
      phi_aux <- matrix(c(rnorm(.m,mu_phi,sqrt(.phi_sigs/nu_phi)),.phi_sigs),2,.m,byrow=T) #auxiliary...
      phi[[l]] <- cbind(phi_prop[[l]],phi_aux) #class parameters.
      p <- c(freqs(C[[l]][-j])/(n_x[l]-1+M)*dnorm(X[[l]][j],phi[[l]][1,1:k_],sqrt(phi[[l]][2,1:k_])), #Draw...
             (M/m)/(n_x[l]-1+M)*dnorm(X[[l]][j],phi[[l]][1,(k_+1):h],sqrt(phi[[l]][2,(k_+1):h]))) #a new...
      C[[l]][j] <- sample(1:h,1,prob=(if(sum(p)==0) rep(1/h,h) else p/sum(p))) #assignment.
      ttnd <- tightener(C[[l]],phi[[l]],j)
      C[[l]] <- ttnd[[1]]
      phi[[l]] <- ttnd[[2]]
      .k <- length(na.omit(unique(C[[l]]))) #Record changes made to...
      phi_prop[[l]] <- phi[[l]][,1:.k,drop=FALSE] #phi in phi_prop.
    }
  }
  
  #Step 2--draw new proper mus and sigma^2s
  .ps[[b]] <- matrix(NA,4,100) #This can be deleted
  for(l in 1:4){
    k <- length(na.omit(unique(C[[l]])))
    for(j in 1:k){
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
      .ps[[b]][l,j] <- .p #This can be deleted
      if(runif(1)<.p){ #accept proposed parameters...
        phi[[l]][,j] <- pro #with that probability...
        if(j %in% 1:3) acpt[(l-1)*3 + j] <- 1 #and...
      } else { #track...
        if(j %in% 1:3) acpt[(l-1)*3 + j] <- 0 #acceptance rate.
      }
    }
    phi_prop[[l]] <- phi[[l]][,1:k,drop=FALSE] #Record changes to phi in phi_prop.
  }
  
  #Step 3--draw new deltas
  for(i in 1:4){
    pro <- MCMCpack::rdirichlet(1,lambda*delta[[i]]) #Draw proposed deltas.
    .pi_n <- loglikelihood_X(deltas=list_replace(delta,pro,i),phis=phi,C=C) + 
      loglikelihood_Y(deltas=list_replace(delta,pro,i),phis=phi,C=C) +
      log(MCMCpack::ddirichlet(pro,deltaPrior)) #Calculate...
    .pi_d <- loglikelihood_X(deltas=delta,phis=phi,C=C) + 
      loglikelihood_Y(deltas=delta,phis=phi,C=C) +
      log(MCMCpack::ddirichlet(delta[[i]],deltaPrior)) #acceptance ratio...
    .p <- exp(.pi_n - .pi_d + MCMCpack::ddirichlet(delta[[i]],lambda*pro) - 
                MCMCpack::ddirichlet(pro,lambda*delta[[i]])) #and...
    if(runif(1)<.p){ #accept proposed parameters...
      delta[[i]] <- pro #with that probability...
      acpt[12+i] <- 1 #and...
    } else { #track...
      acpt[12+i] <- 0 #acceptance rate.
    }
  }
  
  #Step 4--draw new M
  .k <- sum(unlist(lapply(C,FUN=\(x) length(na.omit(unique(x)))))) #Draw...
  z <- rbeta(1,M+1,sum(n_x)) #new...
  pi1 <- a_M+.k+1 #M...
  pi2 <- sum(n_x)*(b_M-log(z)) #from...
  pi <- pi1/(pi1+pi2) #its...
  M <- sample(c(rgamma(1,a_M+.k,scale=b_M-log(z)),rgamma(1,a_M+.k-1,scale=b_M-log(z))),
              1,prob=c(pi,1-pi)) #full conditional.
  
  #Write out iteration
  phi_chain[[b]] <- phi_prop
  C_chain[[b]] <- C
  delta_chain[[b]] <- delta
  M_chain[[b]] <- M
  acpt_chain[[b]] <- acpt
  
  if(progbar) setTxtProgressBar(pb, b)
}

#Return draws
invisible(list("phi"=phi_chain,"c"=C_chain,"delta"=delta_chain,"M"=M_chain,
               "acpt"=acpt_chain,".ps"=.ps))
