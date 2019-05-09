#setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Exam 3")

# Problem 1
rm(list=ls())
alpha <- .05
X <- read.table("Data/T6-8.dat", header = FALSE)
X <- as.matrix(X)
n <- nrow(X)
q <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D
Sinv <- solve(S)
C <- matrix(c(-1,-1,1,1,
              1,-1,1,-1,
              1,-1,-1,1), byrow = TRUE, ncol=4)
#Test Statistics
CX <- C%*%X.mean
CSC <- C%*%S%*%t(C)
CSCinv <- solve(CSC)
T2 <- n*t(CX)%*%CSCinv%*%CX
Tstar <- ((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1)#Critical Value

# Confidence Intervals
# Format Effect
C[1,]%*%X.mean - sqrt(((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))*sqrt(t(C[1,])%*%S%*%C[1,]/n)
C[1,]%*%X.mean + sqrt(((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))*sqrt(t(C[1,])%*%S%*%C[1,]/n)

# Parity Effect
C[2,]%*%X.mean - sqrt(((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))*sqrt(t(C[2,])%*%S%*%C[2,]/n)
C[2,]%*%X.mean + sqrt(((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))*sqrt(t(C[2,])%*%S%*%C[2,]/n)

# Interaction  Effect
C[3,]%*%X.mean - sqrt(((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))*sqrt(t(C[3,])%*%S%*%C[3,]/n)
C[3,]%*%X.mean + sqrt(((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))*sqrt(t(C[3,])%*%S%*%C[3,]/n)



#****************************************************************
#****************************************************************
#****************************************************************
# Problem 2
rm(list=ls())
X <- read.table("Data/T8-7.DAT", header = FALSE)
names(X) <- c("Family", "DistRd", "Cotton", "Maize", "Sorg", "Millet", "Bull", "Cattle", "Goats")
#First Scatter
plot(X$Family,X$DistRd)
png("famdistscatter.png")
plot(X$Family,X$DistRd, xlab = "Family", ylab= "Distance to Road", main = "Family vs DistRd")
dev.off()
#Using the scatterplot, we see that there are two points with DistRd>400 and one family point over 100.
#These are obvious outliers, so we parse for the relevant rows to remove.
bad1 <- which(X$DistRd>400)
bad1 <- c(bad1, which(X$Family > 100))
#Second Scatter
plot(X$DistRd, X$Cattle)
png("distcattlescatter.png")
plot(X$DistRd, X$Cattle, main = "DistRd vs Cattle", xlab = "Distance from Road", ylab = "Number of Cattle")
dev.off()
#Using the scatterplot, we see that there are points with DistRd>400 and one cattle point over 100.
#These are obvious outliers, so we parse for the relevant rows to remove.
bad2 <- which(X$DistRd>400)
bad2 <- c(bad2, which(X$Cattle > 100))
#Combine the outliers we found ignoring repeats and then clean the data:
(Outs <- unique(c(bad1,bad2)))
X <- X[-Outs,]
#### Now we Compute Correlation matrix R
X <- as.matrix(X)
n <- nrow(X)
q <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D
R <- cov2cor(S)
(l.r <- eigen(R)$values)
(e.r <- eigen(R)$vectors)
rownames(e.r) <- c("Family", "DistRd", "Cotton", "Maize", "Sorg", "Millet", "Bull", "Cattle", "Goats")
l.r/sum(l.r)
plot(l.r, main = "Scree Plot w/o Outliers", ylab = "EigenValues", xlab = "i")
lines(l.r)
png("farmscree.png")
plot(l.r, main = "Scree Plot w/o Outliers", ylab = "EigenValues", xlab = "i")
lines(l.r)
dev.off()
cumsum(l.r)/sum(l.r)
plot(cumsum(l.r)/sum(l.r), main = "Cumulative Variance", xlab="i", ylab = "Cumulative % Variance Explained")
lines(cumsum(l.r)/sum(l.r))

#********************************************************************
#********************************************************************
#********************************************************************
#Problem 3
rm(list=ls())
X <- read.csv("Data/creativeclass200711.csv", header = TRUE, sep=',')

#Get State -> Rural -> CCSHARE&BOHSHARE -> Numerics
WI <- X[which(X$State=="Wisconsin"),]
WI <- WI[which(WI$metro03==0),]
WI <- WI[,c(10,14)]
WI[,1] <- as.numeric(paste(WI$CCShare))
WI[,2] <- as.numeric(paste(WI$BohShare))
(WIbar <- matrix(c(mean(WI$CCShare), mean(WI$BohShare)), ncol=1))
WI<-as.matrix(WI)
nWI <- nrow(WI)
qWI <- ncol(WI)
D <- WI - matrix(1,nrow=nWI) %*% t(WIbar)
(SWI <- (nWI-1)^(-1) * t(D)%*%D)

MN <- X[which(X$State=="Minnesota"),]
MN <- MN[which(MN$metro03==0),]
MN <- MN[,c(10,14)]
MN[,1] <- as.numeric(paste(MN$CCShare))
MN[,2] <- as.numeric(paste(MN$BohShare))
(MNbar <- matrix(c(mean(MN$CCShare), mean(MN$BohShare)), ncol=1))
MN<-as.matrix(MN)
nMN <- nrow(MN)
qMN <- ncol(MN)
D <- MN - matrix(1,nrow=nMN) %*% t(MNbar)
(SMN <- (nMN-1)^(-1) * t(D)%*%D)

ND <- X[which(X$State=="North Dakota"),]
ND <- ND[which(ND$metro03==0),]
ND <- ND[,c(10,14)]
ND[,1] <- as.numeric(paste(ND$CCShare))
ND[,2] <- as.numeric(paste(ND$BohShare))
(NDbar <- matrix(c(mean(ND$CCShare), mean(ND$BohShare)), ncol=1))
ND<-as.matrix(ND)
nND <- nrow(ND)
qND <- ncol(ND)
D <- ND - matrix(1,nrow=nND) %*% t(NDbar)
(SND <- (nND-1)^(-1) * t(D)%*%D)

SD <- X[which(X$State=="South Dakota"),]
SD <- SD[which(SD$metro03==0),]
SD <- SD[,c(10,14)]
SD[,1] <- as.numeric(paste(SD$CCShare))
SD[,2] <- as.numeric(paste(SD$BohShare))
(SDbar <- matrix(c(mean(SD$CCShare), mean(SD$BohShare)), ncol=1))
SD<-as.matrix(SD)
nSD <- nrow(SD)
qSD <- ncol(SD)
D <- SD - matrix(1,nrow=nSD) %*% t(SDbar)
(SSD <- (nSD-1)^(-1) * t(D)%*%D)

MI <- X[which(X$State=="Michigan"),]
MI <- MI[which(MI$metro03==0),]
MI <- MI[,c(10,14)]
MI[,1] <- as.numeric(paste(MI$CCShare))
MI[,2] <- as.numeric(paste(MI$BohShare))
(MIbar <- matrix(c(mean(MI$CCShare), mean(MI$BohShare)), ncol=1))
MI<-as.matrix(MI)
nMI <- nrow(MI)
qMI <- ncol(MI)
D <- MI - matrix(1,nrow=nMI) %*% t(MIbar)
(SMI <- (nMI-1)^(-1) * t(D)%*%D)

#Compute Pooled S', pooled mean, and B matrix
(W <- (nWI-1)*SWI+(nMN-1)*SMN+(nND-1)*SND+(nSD-1)*SSD+(nMI-1)*SMI)
(xbar <- (1/(nWI+nMN+nND+nSD+nMI))*(nWI*WIbar+nMN*MNbar+nND*NDbar+nSD*SDbar+nMI*MIbar))
(B <- nWI * (WIbar-xbar)%*%t(WIbar-xbar) + nMN * (MNbar-xbar)%*%t(MNbar-xbar)+nND * (NDbar-xbar)%*%t(NDbar-xbar)+nSD * (SDbar-xbar)%*%t(SDbar-xbar)+nMI * (MIbar-xbar)%*%t(MIbar-xbar))
(Lambda <- det(W)/det(B+W))

n=nWI+nMN+nND+nSD+nMI #Total n
p=2 #Number of measure variables
g=5 #Number of states being compared
alpha=.01 #Test Level
#"Exact" Test
(test.statistic <- (n-p-2)/p * (1-sqrt(Lambda))/sqrt(Lambda))
(critical.value <- qf(1-alpha,2*p,2*(n-p-2)))
#"Large sample" test
(test.statistic.LS <- -(n-1-(p+g)/2)*log(Lambda))
(critical.value.LS <- qchisq(1-alpha,p*(g-1)))

#Intervals
tauWI <- (WIbar-xbar)
tauMN <- (MNbar-xbar)
tauND <- (NDbar-xbar)
tauSD <- (SDbar-xbar)
tauMI <- (MIbar-xbar)
alpha <- 0.05

#WI <-> MN
##CCShare
(tauWI[1,1] - tauMN[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nMN)*W[1,1]/(n-g))
(tauWI[1,1] - tauMN[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nMN)*W[1,1]/(n-g))
##BohShare
(tauWI[2,1] - tauMN[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nMN)*W[2,2]/(n-g))
(tauWI[2,1] - tauMN[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nMN)*W[2,2]/(n-g))

#WI <-> ND
##CCShare
(tauWI[1,1] - tauND[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nND)*W[1,1]/(n-g))
(tauWI[1,1] - tauND[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nND)*W[1,1]/(n-g))
##BohShare
(tauWI[2,1] - tauND[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nND)*W[2,2]/(n-g))
(tauWI[2,1] - tauND[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nND)*W[2,2]/(n-g))

#WI <-> SD
##CCShare
(tauWI[1,1] - tauSD[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nSD)*W[1,1]/(n-g))
(tauWI[1,1] - tauSD[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nSD)*W[1,1]/(n-g))
##BohShare
(tauWI[2,1] - tauSD[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nSD)*W[2,2]/(n-g))
(tauWI[2,1] - tauSD[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nSD)*W[2,2]/(n-g))

#WI <-> MI
##CCShare
(tauWI[1,1] - tauMI[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nMI)*W[1,1]/(n-g))
(tauWI[1,1] - tauMI[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nMI)*W[1,1]/(n-g))
##BohShare
(tauWI[2,1] - tauMI[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nMI)*W[2,2]/(n-g))
(tauWI[2,1] - tauMI[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nWI + 1/nMI)*W[2,2]/(n-g))

#ND <-> MN
##CCShare
(tauND[1,1] - tauMN[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nND + 1/nMN)*W[1,1]/(n-g))
(tauND[1,1] - tauMN[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nND + 1/nMN)*W[1,1]/(n-g))
##BohShare
(tauND[2,1] - tauMN[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nND + 1/nMN)*W[2,2]/(n-g))
(tauND[2,1] - tauMN[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nND + 1/nMN)*W[2,2]/(n-g))

#SD <-> MN
##CCShare
(tauSD[1,1] - tauMN[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nSD + 1/nMN)*W[1,1]/(n-g))
(tauSD[1,1] - tauMN[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nSD + 1/nMN)*W[1,1]/(n-g))
##BohShare
(tauSD[2,1] - tauMN[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nSD + 1/nMN)*W[2,2]/(n-g))
(tauSD[2,1] - tauMN[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nSD + 1/nMN)*W[2,2]/(n-g))

#MI <-> MN
##CCShare
(tauMI[1,1] - tauMN[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nMI + 1/nMN)*W[1,1]/(n-g))
(tauMI[1,1] - tauMN[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nMI + 1/nMN)*W[1,1]/(n-g))
##BohShare
(tauMI[2,1] - tauMN[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nMI + 1/nMN)*W[2,2]/(n-g))
(tauMI[2,1] - tauMN[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nMI + 1/nMN)*W[2,2]/(n-g))

#ND <-> SD
##CCShare
(tauND[1,1] - tauSD[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nND + 1/nSD)*W[1,1]/(n-g))
(tauND[1,1] - tauSD[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nND + 1/nSD)*W[1,1]/(n-g))
##BohShare
(tauND[2,1] - tauSD[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nND + 1/nSD)*W[2,2]/(n-g))
(tauND[2,1] - tauSD[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nND + 1/nSD)*W[2,2]/(n-g))

#ND <-> MI
##CCShare
(tauND[1,1] - tauMI[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nND + 1/nMI)*W[1,1]/(n-g))
(tauND[1,1] - tauMI[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nND + 1/nMI)*W[1,1]/(n-g))
##BohShare
(tauND[2,1] - tauMI[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nND + 1/nMI)*W[2,2]/(n-g))
(tauND[2,1] - tauMI[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nND + 1/nMI)*W[2,2]/(n-g))

#SD <-> MI
##CCShare
(tauSD[1,1] - tauMI[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nSD + 1/nMI)*W[1,1]/(n-g))
(tauSD[1,1] - tauMI[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nSD + 1/nMI)*W[1,1]/(n-g))
##BohShare
(tauSD[2,1] - tauMI[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/nSD + 1/nMI)*W[2,2]/(n-g))
(tauSD[2,1] - tauMI[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/nSD + 1/nMI)*W[2,2]/(n-g))

(Spooled <- (nWI-1)/(nWI+nMN+nND+nSD+nMI-g) * SWI + (nMN-1)/(nWI+nMN+nND+nSD+nMI-g) * SMN + (nND-1)/(nWI+nMN+nND+nSD+nMI-g) * SND +(nSD-1)/(nWI+nMN+nND+nSD+nMI-g) * SSD+(nMI-1)/(nWI+nMN+nND+nSD+nMI-g) * SMI)

u <- (1/(nWI-1) + 1/(nMN-1) + 1/(nND-1) + 1/(nSD-1) + 1/(nMI-1) - 1/(nWI + nMN + nND+nSD+nMI - g)) * ((2*p^2 + 3*p - 1)/(6*(p+1)*(g-1)))
(M <- (nWI + nMN + nND+nSD+nMI - g)*log(det(Spooled)) - ( (nWI-1)*log(det(SWI)) + (nMN-1)*log(det(SMN)) + (nND-1)*log(det(SND)) +(nSD-1)*log(det(SSD)) +(nMI-1)*log(det(SMI))) )

C <- (1-u)*M
v <- (1/2)*p*(p+1)*(g-1)
C.crit <- qchisq(1-alpha,v)
#This is not ideal for our test of equality of the covariance matrices, but we have larg sample sizes N so we proceed anyway


#******************************************************************
#******************************************************************
#******************************************************************
#Problem 4
rm(list=ls())
X <- read.csv("Data/VHA_Facility_Quality_and_Safety_Report_Hospital_Settings.csv", header = TRUE, sep=',')
X <- na.omit(X[9:16])
X <- as.matrix(X)
n <- nrow(X)
q <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D
R <- cov2cor(S)
(l <- eigen(S)$values)
(e <- eigen(S)$vectors)
plot(l)
lines(l)
png("VHAscree.png")
plot(l, main="VHA Data Scree Chart", ylab = "Eigenvalues")
lines(l)
dev.off()
cumsum(l)/sum(l)
P1 <- sqrt(l[1])*X%*%e[,1]
P2 <- sqrt(l[2])*X%*%e[,2]
P3 <- sqrt(l[3])*X%*%e[,3]
P <- matrix(c(P1,P2,P3), ncol=3)
