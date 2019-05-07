#5.1 ##############################################################
rm(list = ls())
alpha=.05
mu = matrix(c(7,11), byrow=TRUE, ncol = 1)
X=matrix(c(2,12,
           8,9,
           6,9,
           8,10), byrow=TRUE, ncol=2)
n <- nrow(X)
p <- ncol(X)
(xbar=matrix(c(mean(X[,1]),
              mean(X[,2])), byrow=TRUE, ncol = 1))
S <- matrix(c(0,0,
               0,0),ncol=2)
for (i in 1:2){
  for (j in 1:2){
    S[i,j]=sum((X[,i]-xbar[i,])*(X[,j]-xbar[i,]))/(n-1)
  }
}
(S)
(Sinv=solve(S))
(T2=n*(t(xbar)-t(mu))%*%Sinv%*%(xbar-mu))

(Tdist <- qf(1-alpha,p,n-p)*(n-1)*p/(n-p))

#5.5 ##############################################################
rm(list = ls())
alpha=.05
mu = matrix(c(0.55,0.60), byrow=TRUE, ncol = 1)
xbar = matrix(c(0.564,0.603), byrow=TRUE, ncol = 1)
Sinv<- matrix(c(203.018,-163.391,
                -163.391,200.228),ncol=2)
n <- 42
p <- 2
(T2=n*(t(xbar)-t(mu))%*%Sinv%*%(xbar-mu))

(Tdist <- qf(1-alpha,p,n-p)*(n-1)*p/(n-p))

#5.7 ##############################################################
rm(list = ls())
alpha=.05
dat <- read.table("T5-1.DAT", header=FALSE)
names(dat) <- c("Sweat Rate", "Sodium","Potassium")
xbar = matrix(c(mean(dat$`Sweat Rate`),mean(dat$Sodium),mean(dat$Potassium)), byrow=TRUE, ncol = 1)
n <- nrow(dat)
p <- ncol(dat)
S <- matrix(c(0,0,0,
              0,0,0,
              0,0,0),ncol=3)
for (i in 1:p){
  for (j in 1:p){
    S[i,j]=sum((dat[,i]-xbar[i,])*(dat[,j]-xbar[i,]))/(n-1)
  }
}
(S)
(Fstar<-qf(1-alpha,p,n-p))

(ISweat=matrix(c(mean(dat$`Sweat Rate`-sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[1,1]/n)),
                mean(dat$`Sweat Rate`+sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[1,1]/n))), byrow=TRUE, ncol=2))
(ISod=matrix(c(mean(dat$Sodium-sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[2,2]/n)),
                     mean(dat$Sodium+sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[2,2]/n))), byrow=TRUE, ncol=2))
(IPot=matrix(c(mean(dat$Potassium-sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[3,3]/n)),
                     mean(dat$Potassium+sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[3,3]/n))), byrow=TRUE, ncol=2))

(ISweatBon=matrix(c(mean(dat$`Sweat Rate`-abs(qt(alpha/(3),n-1))*sqrt(S[1,1]/n)),
                    mean(dat$`Sweat Rate`+abs(qt(alpha/(3),n-1))*sqrt(S[1,1]/n))), byrow=TRUE, ncol=2))
(ISodBon=matrix(c(mean(dat$Sodium-abs(qt(alpha/(3),n-1))*sqrt(S[2,2]/n)),
                  mean(dat$Sodium+abs(qt(alpha/(3),n-1))*sqrt(S[2,2]/n))), byrow=TRUE, ncol=2))
(IPotBon=matrix(c(mean(dat$Potassium-abs(qt(alpha/(3),n-1))*sqrt(S[3,3]/n)),
                  mean(dat$Potassium+abs(qt(alpha/(3),n-1))*sqrt(S[3,3]/n))), byrow=TRUE, ncol=2))

#5.18 ##############################################################
rm(list = ls())
alpha=.05
mu = matrix(c(500,50,3), byrow=TRUE, ncol = 1)
X <- read.table("T5-2.DAT", header=FALSE)
names(X) <- c("Social Science and History","Verbal","Science")
n <- nrow(X)
p <- ncol(X)
(xbar=matrix(c(mean(X[,1]),
               mean(X[,2]),
               mean(X[,3])), byrow=TRUE, ncol = 1))
S <- matrix(c(0,0,0,
              0,0,0,
              0,0,0),ncol=3)
for (i in 1:p){
  for (j in 1:p){
    S[i,j]=sum((X[,i]-xbar[i,])*(X[,j]-xbar[i,]))/(n-1)
  }
}
(S)
(Sinv=solve(S))
(T2=n*(t(xbar)-t(mu))%*%Sinv%*%(xbar-mu))

(Tdist <- qf(1-alpha,p,n-p)*(n-1)*p/(n-p))
lambda <- eigen(S)$values
(vecs <- eigen(S)$vectors)
Lengths <- matrix(c(sqrt(lambda[1])*sqrt(qf(1-alpha,p,n-p)*p*(n-1)/(n*(n-p))),
                    sqrt(lambda[2])*sqrt(qf(1-alpha,p,n-p)*p*(n-1)/(n*(n-p))),
                    sqrt(lambda[3])*sqrt(qf(1-alpha,p,n-p)*p*(n-1)/(n*(n-p)))),byrow=TRUE, ncol=1)
(Lengths)

qqnorm(X[,1])
png("QQSoc.png")
qqnorm(X[,1])
dev.off()

qqnorm(X[,2])
png("QQVerb.png")
qqnorm(X[,2])
dev.off()

qqnorm(X[,3])
png("QQSci.png")
qqnorm(X[,3])
dev.off()


plot(X$`Social Science and History`,X$Verbal,main="Social Vs Verbal", xlab="Social Science and History", ylab="Verbal")
png("SocialvHistory.png")
plot(X$`Social Science and History`,X$Verbal,main="Social Vs Verbal", xlab="Social Science and History", ylab="Verbal")
dev.off()

plot(X$`Social Science and History`,X$Science,main="Social Vs Science", xlab="Social Science and History", ylab="Science")
png("SocialvScience.png")
plot(X$`Social Science and History`,X$Verbal,main="Social Vs Science", xlab="Social Science and History", ylab="Science")
dev.off()

plot(X$Verbal,X$Science,main="Verbal vs Science", xlab="Verbal", ylab="Science")
png("VerbalvScience.png")
plot(X$Verbal,X$Science,main="Verbal vs Science", xlab="Verbal", ylab="Science")
dev.off()

#5.22 ##############################################################
rm(list = ls())
alpha=.05
X <- read.table("T5-13.dat", header=FALSE)
Xout <- X[-c(9, 21),]
names(X) <- c("Fuel","Repair","Capital")
xbar = matrix(c(mean(Xout$Fuel),mean(Xout$Repair),mean(Xout$Capital)), byrow=TRUE, ncol = 1)
n <- nrow(Xout)
p <- ncol(Xout)
S <- matrix(c(0,0,0,
              0,0,0,
              0,0,0),ncol=p)
for (i in 1:p){
  for (j in 1:p){
    S[i,j]=sum((Xout[,i]-xbar[i,])*(Xout[,j]-xbar[j,]))/(n-1)
  }
}
(S)
(Fstar<-qf(1-alpha,p,n-p))


qqnorm(X[,1])
png("QQFuel.png")
qqnorm(X[,1])
dev.off()

qqnorm(X[,2])
png("QQRepair.png")
qqnorm(X[,2])
dev.off()

qqnorm(X[,3])
png("QQCapital.png")
qqnorm(X[,3])
dev.off()

plot(X$Fuel,X$Repair,main="Fuel vs. Repair", xlab="Fuel", ylab="Repair")
png("FuelvRepair.png")
plot(X$Fuel,X$Repair,main="Fuel vs. Repair", xlab="Fuel", ylab="Repair")
dev.off()

plot(X$Fuel,X$Capital,main="Fuel vs. Capital", xlab="Fuel", ylab="Capital")
png("FuelvCapital.png")
plot(X$Fuel,X$Capital,main="Fuel vs. Capital", xlab="Fuel", ylab="Capital")
dev.off()

plot(X$Capital,X$Repair,main="Capital vs. Repair", xlab="Capital", ylab="Repair")
png("CapitalvRepair.png")
plot(X$Capital,X$Repair,main="Capital vs. Repair", xlab="Capital", ylab="Repair")
dev.off()





qqnorm(X[,1])
png("QQFuelOut.png")
qqnorm(Xout[,1])
dev.off()

qqnorm(X[,2])
png("QQRepairOut.png")
qqnorm(Xout[,2])
dev.off()

qqnorm(X[,3])
png("QQCapitalOut.png")
qqnorm(Xout[,3])
dev.off()

plot(Xout$Fuel,Xout$Repair,main="Fuel vs. Repair", xlab="Fuel", ylab="Repair")
png("FuelvRepairOut.png")
plot(Xout$Fuel,Xout$Repair,main="Fuel vs. Repair", xlab="Fuel", ylab="Repair")
dev.off()

plot(Xout$Fuel,Xout$Capital,main="Fuel vs. Capital", xlab="Fuel", ylab="Capital")
png("FuelvCapitalOut.png")
plot(Xout$Fuel,Xout$Capital,main="Fuel vs. Capital", xlab="Fuel", ylab="Capital")
dev.off()

plot(Xout$Capital,Xout$Repair,main="Capital vs. Repair", xlab="Capital", ylab="Repair")
png("CapitalvRepairOut.png")
plot(Xout$Capital,Xout$Repair,main="Capital vs. Repair", xlab="Capital", ylab="Repair")
dev.off()






(IFuel=matrix(c(mean(Xout$Fuel-sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[1,1]/n)),
                 mean(Xout$Fuel+sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[1,1]/n))), byrow=TRUE, ncol=2))
(IRepair=matrix(c(mean(Xout$Repair-sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[2,2]/n)),
               mean(Xout$Repair+sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[2,2]/n))), byrow=TRUE, ncol=2))
(ICapital=matrix(c(mean(Xout$Capital-sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[3,3]/n)),
               mean(Xout$Capital+sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[3,3]/n))), byrow=TRUE, ncol=2))

(IFuelBon=matrix(c(mean(Xout$Fuel-abs(qt(alpha/(2*p),n-1))*sqrt(S[1,1]/n)),
                    mean(Xout$Fuel+abs(qt(alpha/(2*p),n-1))*sqrt(S[1,1]/n))), byrow=TRUE, ncol=2))
(IRepairBon=matrix(c(mean(Xout$Repair-abs(qt(alpha/(2*p),n-1))*sqrt(S[2,2]/n)),
                  mean(Xout$Repair+abs(qt(alpha/(2*p),n-1))*sqrt(S[2,2]/n))), byrow=TRUE, ncol=2))
(ICapitalBon=matrix(c(mean(Xout$Capital-abs(qt(alpha/(2*p),n-1))*sqrt(S[3,3]/n)),
                  mean(Xout$Capital+abs(qt(alpha/(2*p),n-1))*sqrt(S[3,3]/n))), byrow=TRUE, ncol=2))

#5.23 #################
rm(list = ls())
alpha=.05
X <- read.table("T6-13.dat", header=FALSE)
names(X) <- c("MaxBreath","BasHeight","BasLength","NasHeight", "Time")
X<-cbind(X$MaxBreath,X$BasHeight,X$BasLength,X$NasHeight)
n <- nrow(X)
p <- ncol(X)
(x.bar <- t(matrix(1,ncol=n) %*% X)/n)
(D <- X - matrix(1,nrow=n) %*% t(x.bar))
(S <- (n-1)^(-1) * t(D)%*%D)
Sinv <- solve(S)

d <- c()
for (i in 1:n){
  d[i] <- (X[i,]-t(x.bar))%*%Sinv%*%t(X[i,]-t(x.bar))
}

# Example 4.13
chisq.quantiles <- qchisq(((1:n)-0.5)/n,df=p)
plot(y=sort(d),x=chisq.quantiles,pch=19, main="ChiSquare")
png("5chi.png")
plot(y=sort(d),x=chisq.quantiles,pch=19, main="ChiSquare")
dev.off()


qqnorm(X[,1])
png("QQMaxB.png")
qqnorm(X[,1])
dev.off()

qqnorm(X[,2])
png("QQBasH.png")
qqnorm(X[,2])
dev.off()

qqnorm(X[,3])
png("QQBasL.png")
qqnorm(X[,3])
dev.off()

qqnorm(X[,4])
png("QQNasH.png")
qqnorm(X[,4])
dev.off()


(Fstar<-qf(1-alpha,p,n-p))

(IMaxB=matrix(c(mean(X[,1]-sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[1,1]/n)),
                mean(X[,1]+sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[1,1]/n))), byrow=TRUE, ncol=2))
(IBasH=matrix(c(mean(X[,2]-sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[2,2]/n)),
                  mean(X[,2]+sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[2,2]/n))), byrow=TRUE, ncol=2))
(IBasL=matrix(c(mean(X[,3]-sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[3,3]/n)),
                   mean(X[,3]+sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[3,3]/n))), byrow=TRUE, ncol=2))
(INasH=matrix(c(mean(X[,4]-sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[3,3]/n)),
                mean(X[,4]+sqrt(Fstar*p*(n-1)/(n-p))*sqrt(S[3,3]/n))), byrow=TRUE, ncol=2))

(IMaxBBon=matrix(c(mean(X[,1]-abs(qt(alpha/(2*p),n-1))*sqrt(S[1,1]/n)),
                   mean(X[,1]+abs(qt(alpha/(2*p),n-1))*sqrt(S[1,1]/n))), byrow=TRUE, ncol=2))
(IBasHBon=matrix(c(mean(X[,2]-abs(qt(alpha/(2*p),n-1))*sqrt(S[2,2]/n)),
                     mean(X[,2]+abs(qt(alpha/(2*p),n-1))*sqrt(S[2,2]/n))), byrow=TRUE, ncol=2))
(IBasLBon=matrix(c(mean(X[,3]-abs(qt(alpha/(2*p),n-1))*sqrt(S[3,3]/n)),
                      mean(X[,3]+abs(qt(alpha/(2*p),n-1))*sqrt(S[3,3]/n))), byrow=TRUE, ncol=2))
(INasHBon=matrix(c(mean(X[,4]-abs(qt(alpha/(2*p),n-1))*sqrt(S[3,3]/n)),
                   mean(X[,4]+abs(qt(alpha/(2*p),n-1))*sqrt(S[3,3]/n))), byrow=TRUE, ncol=2))
