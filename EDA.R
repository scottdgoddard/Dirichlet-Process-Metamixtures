sampleB <- read.csv("DZMix B.csv")
sampleC <- read.csv("DZMix C.csv")
sampleE <- read.csv("DZMix E.csv")
sampleI <- read.csv("DZMix I.csv")

postscript("/Users/sdgoddard/Library/CloudStorage/GoogleDrive-sdgoddard@alaska.edu/My Drive/Project support/Loess/Stat article/Loess Fall 2025/B.eps",width=6,height=2)
par(mar=c(4,4,0.2,0),cex=0.8)
plot(density(sampleB$Mean.age,na.rm=TRUE,bw="SJ"),main="",ylim=c(0,.0035))
lines(density(sampleB$Nenana.age,na.rm=TRUE,bw="SJ"),col="blue",lty="dashed")
lines(density(sampleB$Delta.age,na.rm=TRUE,bw="SJ"),col="darkgreen",lty="dashed")
lines(density(sampleB$Tanana.age,na.rm=TRUE,bw="SJ"),col="orange",lty="dashed")
lines(density(sampleB$Yukon.age,na.rm=TRUE,bw="SJ"),col="pink3",lty="dashed")
legend("topright",legend=c("Loess","Nenana","Delta","Tanana","Yukon"),
       col=c("black","blue","darkgreen","orange","pink3"),
       lty=c("solid",rep("dashed",4)))
dev.off()

postscript("/Users/sdgoddard/Library/CloudStorage/GoogleDrive-sdgoddard@alaska.edu/My Drive/Project support/Loess/Stat article/Loess Fall 2025/C.eps",width=6,height=2)
par(mar=c(4,4,0.2,0),cex=0.8)
plot(density(sampleC$Mean.age,na.rm=TRUE,bw="SJ"),main="",ylim=c(0,.0035))
lines(density(sampleC$Nenana.age,na.rm=TRUE,bw="SJ"),col="blue",lty="dashed")
lines(density(sampleC$Delta.age,na.rm=TRUE,bw="SJ"),col="darkgreen",lty="dashed")
lines(density(sampleC$Tanana.age,na.rm=TRUE,bw="SJ"),col="orange",lty="dashed")
lines(density(sampleC$Yukon.age,na.rm=TRUE,bw="SJ"),col="pink3",lty="dashed")
legend("topright",legend=c("Loess","Nenana","Delta","Tanana","Yukon"),
       col=c("black","blue","darkgreen","orange","pink3"),
       lty=c("solid",rep("dashed",4)))
dev.off()

postscript("/Users/sdgoddard/Library/CloudStorage/GoogleDrive-sdgoddard@alaska.edu/My Drive/Project support/Loess/Stat article/Loess Fall 2025/E.eps",width=6,height=2)
par(mar=c(4,4,0.2,0),cex=0.8)
plot(density(sampleE$Mean.age,na.rm=TRUE,bw="SJ"),main="",ylim=c(0,.0035))
lines(density(sampleE$Nenana.age,na.rm=TRUE,bw="SJ"),col="blue",lty="dashed")
lines(density(sampleE$Delta.age,na.rm=TRUE,bw="SJ"),col="darkgreen",lty="dashed")
lines(density(sampleE$Tanana.age,na.rm=TRUE,bw="SJ"),col="orange",lty="dashed")
lines(density(sampleE$Yukon.age,na.rm=TRUE,bw="SJ"),col="pink3",lty="dashed")
legend("topright",legend=c("Loess","Nenana","Delta","Tanana","Yukon"),
       col=c("black","blue","darkgreen","orange","pink3"),
       lty=c("solid",rep("dashed",4)))
dev.off()

postscript("/Users/sdgoddard/Library/CloudStorage/GoogleDrive-sdgoddard@alaska.edu/My Drive/Project support/Loess/Stat article/Loess Fall 2025/I.eps",width=6,height=2)
par(mar=c(4,4,0.2,0),cex=0.8)
plot(density(sampleI$Mean.age,na.rm=TRUE,bw="SJ"),main="",ylim=c(0,.0035))
lines(density(sampleI$Nenana.age,na.rm=TRUE,bw="SJ"),col="blue",lty="dashed")
lines(density(sampleI$Delta.age,na.rm=TRUE,bw="SJ"),col="darkgreen",lty="dashed")
lines(density(sampleI$Tanana.age,na.rm=TRUE,bw="SJ"),col="orange",lty="dashed")
lines(density(sampleI$Yukon.age,na.rm=TRUE,bw="SJ"),col="pink3",lty="dashed")
legend("topright",legend=c("Loess","Nenana","Delta","Tanana","Yukon"),
       col=c("black","blue","darkgreen","orange","pink3"),
       lty=c("solid",rep("dashed",4)))
dev.off()