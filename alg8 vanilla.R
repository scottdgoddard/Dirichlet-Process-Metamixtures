library(invgamma)

freqs <- function(x) {
  out <- NULL
  for(j in 1:length(unique(x))) {
    out <- c(out,sum(x==j))
  }
  out
}

alg8 <- function(x,m,g0Priors=c(0,1,1,1),alphaPrior=c(2,4),B) {
  #Initialization
  #set.seed(1)
  mu_phi <- g0Priors[1]
  nu_phi <- g0Priors[2]
  a_phi <- g0Priors[3]
  b_phi <- g0Priors[4]
  a_alpha <- alphaPrior[1]
  b_alpha <- alphaPrior[2]
  n <- length(x)
  c <- rep(1,n)
  .k <- length(unique(c))
  phi_prop <- rbind(rep(mu_phi,.k),rep((b_phi/(a_phi-1))/nu_phi,.k)) #top row is mean; bottom row is variance
  phi_chain <- list()
  c_chain <- list()
  alpha_chain <- list()
  alpha <- 1
  
  #Algorithm
  for(b in 1:B){
    for(i in 1:n){
      
      #Bullet point 1
      k_ <- length(unique(c[-i]))
      h <- k_ + m
      c <- as.numeric(as.factor(c))
      if(c[i] %in% c[-i]){
        phi <- cbind(phi_prop,matrix(c(rep(NA,m),#Do we draw directly from G_0 because it is the full conditional
                                       rinvgamma(m,a_phi,scale=b_phi)),2,m,byrow=TRUE))#of phi for these aux. pars?
        phi[1,(k_+1):h] <- rnorm(m,mu_phi,sqrt(phi[2,(k_+1):h]/nu_phi))#Is this because the data give us no info
      } else {#about them? I think so.
          c[i] <- k_ + 1
          phi <- cbind(phi_prop,matrix(c(rep(NA,m-1),
                                         rinvgamma(m-1,a_phi,scale=b_phi)),2,m-1,byrow=TRUE))
          phi[1,(k_+2):h] <- rnorm(m-1,mu_phi,sqrt(phi[2,(k_+2):h]/nu_phi))
      }
      p <- c(freqs(c[-i])/(n-1+alpha)*dnorm(x[i],phi[1,1:k_],sqrt(phi[2,1:k_])),
             (alpha/m)/(n-1+alpha)*dnorm(x[i],phi[1,(k_+1):h],sqrt(phi[2,(k_+1):h])))
      .ci <- c[i]
      c[i] <- sample(1:h,1,prob=p/sum(p))
      c <- as.numeric(as.factor(c))
      if(!(.ci %in% c)) phi <- phi[,-.ci]
      .k_ <- length(unique(c))
      phi_prop <- phi[,unique(c(1:.k_,c[i])),drop=FALSE]
    }
    
    #Bullet point 2
    k <- length(unique(c))
    for(j in 1:k){
      .x <- x[c==j]
      .n <- length(.x)
      phi_prop[2,j] <- rinvgamma(1,a_phi+n/2,scale=b_phi+0.5*sum((.x-mean(.x))^2))
      phi_prop[1,j] <- rnorm(1,(sum(.x)+mu_phi*nu_phi)/(.n+nu_phi),sqrt(phi_prop[2,j]/(.n+nu_phi)))
    }
    
    #Inference on alpha
    z <- rbeta(1,alpha+1,n)
    pi1 <- a_alpha+k+1
    pi2 <- n*(b_alpha-log(z))
    pi <- pi1/(pi1+pi2)
    alpha <- sample(c(rgamma(1,a_alpha+k,scale=b_alpha-log(z)),rgamma(1,a_alpha+k-1,scale=b-log(z))),
                    1,prob=c(pi,1-pi))
    
    #Write out iteration
    phi_chain[[b]] <- phi_prop
    c_chain[[b]] <- c
    alpha_chain[[b]] <- alpha
  }
  
  #Return draws
  invisible(list("phi"=phi_chain,"c"=c_chain,"alpha"=alpha_chain))
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

y <- c(-1.48,-1.4,-1.16,-1.08,-1.02,0.14,0.51,0.53,0.78)
y_ <- scale(y)
plot(density(y_),main="Data set y from Neal (2000)",ylim=c(0,0.4))

ans <- alg8(y_,3,c(0,1,1,1),c(2,4),1e4)
color <- unlist(Mode(ans$c))
segments(y_,0,y_,0.045,col=color,lwd=2)
text(2,0.35,labels=bquote("E("*alpha*"|y)"==.(round(mean(unlist(ans$alpha)),3))))

library(dirichletprocess)
model <- DirichletProcessGaussian(y_, g0Priors = c(0, 1, 1, 1), alphaPriors = c(2, 4))
fit <- Fit(model,its=1e4)
color2 <- unlist(Mode(fit$labelsChain))+5
segments(y_,0.05,y_,0.10,col=color2,lwd=2)
text(2,0.32,labels=bquote("E("*alpha*"|y)"==.(round(mean(unlist(fit$alphaChain)),3))),col="#CD0BBC")
curve(dnorm(x,fit$clusterParameters[[1]][1],sqrt(fit$clusterParameters[[2]][1])),add=TRUE,col="#CD0BBC")
curve(dnorm(x,fit$clusterParameters[[1]][2],sqrt(fit$clusterParameters[[2]][2])),add=TRUE,col="#CD0BBC")
curve(dnorm(x,fit$clusterParameters[[1]][3],sqrt(fit$clusterParameters[[2]][3])),add=TRUE,col="#CD0BBC")
curve(dnorm(x,fit$clusterParameters[[1]][4],sqrt(fit$clusterParameters[[2]][4])),add=TRUE,col="#CD0BBC")
curve(dnorm(x,fit$clusterParameters[[1]][1],sqrt(fit$clusterParameters[[2]][1]))*fit$weights[1]+
        dnorm(x,fit$clusterParameters[[1]][2],sqrt(fit$clusterParameters[[2]][2]))*fit$weights[2]+
        dnorm(x,fit$clusterParameters[[1]][3],sqrt(fit$clusterParameters[[2]][3]))*fit$weights[3]+
        dnorm(x,fit$clusterParameters[[1]][4],sqrt(fit$clusterParameters[[2]][4]))*fit$weights[4],
      add=TRUE,col="blue")

library(BNPmix)
