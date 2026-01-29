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

results <- mcmc(B=1e3,tuning=c(2,25))
 
#Are the answers reasonable?
matrix(apply(matrix(unlist(results$delta),length(results$delta),16,byrow=TRUE),MARGIN=2,FUN=mean),4,4,byrow=TRUE)
#0.25,0.25,0.25,0.25
#0.33,0.33,0.17,0.17
#0.5,0.25,0.13,0.13
#0.25,0.33,0.25,0.17

#Acceptance rates
apply(matrix(unlist(results$acpt),length(results$acpt),16,byrow=T),MARGIN=2,FUN=mean)

#Trace plots
png("/Users/sdgoddard/Desktop/testplot.png")
matplot(matrix(unlist(results$delta),length(results$delta),16,byrow=TRUE)[,1:4],type='l',ylab=expression(delta))
dev.off()
matplot(matrix(unlist(results$delta),length(results$delta),16,byrow=TRUE)[,5:8],type='l',ylab=expression(delta))
matplot(matrix(unlist(results$delta),length(results$delta),16,byrow=TRUE)[,9:12],type='l',ylab=expression(delta))
matplot(matrix(unlist(results$delta),length(results$delta),16,byrow=TRUE)[,13:16],type='l',ylab=expression(delta))

plot(unlist(results$M),type='l',ylab=expression(italic(M)),xlab="")

phi_abr <- lapply(paste0("phi_abr",1:2),FUN=assign,value=matrix(NA,length(results$phi),12))
for(i in seq_along(results$phi)){
  for(j in 1:4){
    phi_abr[[1]][i,((j-1)*3+1):((j-1)*3+3)] <- results$phi[[i]][[j]][1,1:3]
    phi_abr[[2]][i,((j-1)*3+1):((j-1)*3+3)] <- results$phi[[i]][[j]][2,1:3]
  }
}
matplot(phi_abr[[1]][,1:3],ylab=expression(mu),type='l')
matplot(phi_abr[[1]][,4:6],ylab=expression(mu),type='l')
matplot(phi_abr[[1]][,7:9],ylab=expression(mu),type='l')
matplot(phi_abr[[1]][,10:12],ylab=expression(mu),type='l')
matplot(phi_abr[[2]][,1:3],ylab=expression(sigma^2),type='l')
matplot(phi_abr[[2]][,4:6],ylab=expression(sigma^2),type='l')
matplot(phi_abr[[2]][,7:9],ylab=expression(sigma^2),type='l')
matplot(phi_abr[[2]][,10:12],ylab=expression(sigma^2),type='l')

c_abr <- matrix(NA,length(results$c),4)
for(i in seq_along(results$c)){
  for(j in 1:4){
    c_abr[i,j] <- results$c[[i]][[j]][99]
  }
}
matplot(c_abr,type='l',ylab=expression(italic(c)))


#Remove burn-in
bi <- 100
Result <- list()
Result$phi <- results$phi[(bi+1):length(results$phi)]
Result$delta <- results$delta[(bi+1):length(results$delta)]
Result$c <- results$c[(bi+1):length(results$c)]
Result$M <- results$M[(bi+1):length(results$M)]
Result$acpt <- results$acpt[(bi+1):length(results$acpt)]

#Where to go from here?
#-Run more case studies to double-check sampler ☑️
#-Run on cluster for 1e4-1e5 iterations ☑️
#-Search out better convergence/mixing diagnostics
#-Validate my algorithm 8 code with BNP to see what I learn there
#-Work out the mathematics of the sampler for this problem ☑️
#-Design a corrupted sampler with pressure to move correctly to see if it works better
