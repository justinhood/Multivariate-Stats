#setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Exam 1")

### 1 ################################################################
rm(list = ls())
psych <- read.table("./data/PSYCHPROFILE.DAT")
names(psych) <- c("Independence", "Support", "Benevolence", "Conformity", "Leadership", "Gender", "Socioeconomic")
X <- psych[,1:5]
xbar <- matrix(c(sum(X[,1])/length(X[,1]),
                 sum(X[,2])/length(X[,2]),
                 sum(X[,3])/length(X[,3]),
                 sum(X[,4])/length(X[,4]),
                 sum(X[,5])/length(X[,5])), ncol = 5)
(xbar)
(S <- matrix(rep(NA,5*5),ncol=5))
for(i in 1:5){
  for(j in 1:5){
    S[i,j]=sum((X[,i]-xbar[,i])*(X[,j]-xbar[,j]))/(length(X[,1])-1)
  }
}
(S)
R <- matrix(rep(NA,5*5),ncol=5)
for(i in 1:5){
  for(j in 1:5){
    R[i,j] = S[i,j]/(sqrt(S[i,i])*sqrt(S[j,j]))
  }
}
(R)
SN=((length(X[,1])-1)/length(X[,1]))*S
(SN)
RN <- matrix(rep(NA,5*5),ncol=5)
for(i in 1:5){
  for(j in 1:5){
    RN[i,j] = SN[i,j]/(sqrt(SN[i,i])*sqrt(SN[j,j]))
  }
}
(RN)

plot(X)
boxplot(X)
png("pairs.png")
plot(X)
dev.off()
png("box.png")
boxplot(X)
dev.off()



library(ggplot2)
qqnorm(X$Independence)
png("indepqq.png")
qqnorm(X$Independence)
dev.off()
qqnorm(X$Support)
qqnorm(X$Benevolence)
qqnorm(X$Conformity)
qqnorm(X$Leadership)
png("leaderqq.png")
qqnorm(X$Leadership)
dev.off()

#xo <- sort(X$Independence)
#n <- length(xo)
#problevels <- ((1:n)-0.5)/n
#q <- qnorm(problevels)
#plot(q,xo,pch=19)
#(rind=(sum((xo-mean(xo))*(q-mean(q))))/(sqrt(sum((xo-mean(xo))**2))*sqrt(sum((q-mean(q))**2))))

rQ <- c()
for(i in 1:5){
  xo <- sort(X[,i])
  n <- length(xo)
  problevels <- ((1:n)-0.5)/n
  q <- qnorm(problevels)
  rQ[i] <- sum((xo-mean(xo))*(q-mean(q)))/(sqrt(sum((xo-mean(xo))**2))*sqrt(sum((q-mean(q))**2)))
}
(rQ)

l.lambda <- function(lambda,n,x){
  if(lambda==0){
    x.lambda <- log(x)
  }
  else {
    x.lambda <- (x^lambda - 1)/lambda
  }
  return(
    (-n/2)*log(
      sum(
        (x.lambda-mean(x.lambda))^2)
    )+
      (lambda-1)*sum(log(x)))
}

lambdas <- seq(from=-1,to=1.50,by=.01)
l.lambdas <- c()
for (i in 1:length(lambdas)){
  l.lambdas[i] <- l.lambda(lambdas[i],n,X$Leadership)
}
plot(lambdas,l.lambdas)
lambdas[which(l.lambdas==max(l.lambdas))] # this is the transformation we should use

qqnorm(X$Leadership)
qqnorm((X$Leadership^0.38 - 1)/0.38)

LeadAdj <- (X$Leadership^0.38 - 1)/0.38

xo <- sort(LeadAdj)
n <- length(xo)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
sum((xo-mean(xo))*(q-mean(q)))/(sqrt(sum((xo-mean(xo))**2))*sqrt(sum((q-mean(q))**2)))

n <- length(X[,1])
(Sinv <- solve(S))
Xstar <- matrix(c(X$Independence,
                  X$Support,
                  X$Benevolence,
                  X$Conformity,
                  X$Leadership), byrow = TRUE, ncol=5)
d <- c()
for (i in 1:n){
  d[i] <- (Xstar[i,]-xbar)%*%Sinv%*%t(Xstar[i,]-xbar)
}
print(d)


chisq.quantiles <- qchisq(((1:n)-0.5)/n,df=5)
plot(y=sort(d),x=chisq.quantiles,pch=19)
png("psychchi.png")
plot(y=sort(d),x=chisq.quantiles,pch=19)
dev.off()


(cutoff=qchisq(0.50,5))
which(d<cutoff)
sum(d<cutoff)/n 

(xavg <- mean(xbar))
(Savg <- sum(S)/25)
### 2 #########################################
rm(list = ls())
sigma <- matrix(c(19254770, 9644284, 31090422,
                  9644284, 78716208, 32475410,
                  31090422, 32475410, 77679998), byrow=TRUE, ncol=3)
Astar <- matrix(c(19254770, 9644284,
                  9644284, 78716208), byrow=TRUE, ncol=2)
det(Astar)
det(sigma)
eigen(sigma)

## 3 #################################################
rm(list = ls())
A <- matrix(c(9, -2,
              -2, 6), byrow = TRUE, ncol=2)
det(A)
v <- eigen(A)
B <- v[1]$values[1]*v[2]$vectors[,1]%*%t(v[2]$vectors[,1])+v[1]$values[2]*v[2]$vectors[,2]%*%t(v[2]$vectors[,2])
A - B
## 4 #################################################
rm(list = ls())
A <- matrix(c(4, 4.001,
              4.001, 4.002), byrow=TRUE, ncol=2)
B <- matrix(c(4, 4.001,
              4.001, 4.002001), byrow=TRUE, ncol=2)
det(A)
det(B)
solve(A)
-3*solve(B)
### 5 #########################################################
rm(list = ls())
EQIDATA <- read.csv("./data/EQI_RESULTS_2013JULY22.CSV")[,1:10]
EQI <- EQIDATA[,5:10]
EQIbar <- matrix(c(sum(EQI[,1])/length(EQI[,1]),
                 sum(EQI[,2])/length(EQI[,2]),
                 sum(EQI[,3])/length(EQI[,3]),
                 sum(EQI[,4])/length(EQI[,4]),
                 sum(EQI[,5])/length(EQI[,5]),
                 sum(EQI[,6])/length(EQI[,6])), ncol = 6)
(EQIbar)
S <- matrix(rep(NA,6*6),ncol=6)
for(i in 1:6){
  for(j in 1:6){
    S[i,j]=sum((EQI[,i]-EQIbar[,i])*(EQI[,j]-EQIbar[,j]))/(length(EQI[,1])-1)
  }
}
(S)
R <- matrix(rep(NA,6*6),ncol=6)
for(i in 1:6){
  for(j in 1:6){
    R[i,j] = S[i,j]/(sqrt(S[i,i])*sqrt(S[j,j]))
  }
}
(R)
(det(S))
eigen(S)

EQIAL<-EQIDATA[EQIDATA$state=="AL",][,5:10]
EQICA<-EQIDATA[EQIDATA$state=="CA",][,5:10]
EQICT<-EQIDATA[EQIDATA$state=="CT",][,5:10]
EQIWI<-EQIDATA[EQIDATA$state=="WI",][,5:10]
dets <- c()
SAL <- matrix(rep(NA,5*5),ncol=5)
SCA <- matrix(rep(NA,5*5),ncol=5)
SCT <- matrix(rep(NA,5*5),ncol=5)
SWI <- matrix(rep(NA,5*5),ncol=5)

xALbar <- matrix(c(sum(EQIAL[,1])/length(EQIAL[,1]),
                   sum(EQIAL[,2])/length(EQIAL[,2]),
                   sum(EQIAL[,3])/length(EQIAL[,3]),
                   sum(EQIAL[,4])/length(EQIAL[,4]),
                   sum(EQIAL[,5])/length(EQIAL[,5])), ncol=5)
                   #sum(EQIAL[,6])/length(EQIAL[,6])), ncol = 6)
for(i in 1:5){
  for(j in 1:5){
    SAL[i,j]<- sum((EQIAL[,i]-xALbar[,i])*(EQIAL[,j]-xALbar[,j]))/(length(EQIAL[,1])-1)
  }
}
dets[1] <- det(SAL)

xCAbar <- matrix(c(sum(EQICA[,1])/length(EQICA[,1]),
                   sum(EQICA[,2])/length(EQICA[,2]),
                   sum(EQICA[,3])/length(EQICA[,3]),
                   sum(EQICA[,4])/length(EQICA[,4]),
                   sum(EQICA[,5])/length(EQICA[,5])), ncol=5)
                   #sum(EQICA[,6])/length(EQICA[,6])), ncol = 6)
for(i in 1:5){
  for(j in 1:5){
    SCA[i,j]<- sum((EQICA[,i]-xCAbar[,i])*(EQICA[,j]-xCAbar[,j]))/(length(EQICA[,1])-1)
  }
}
dets[2] <- det(SCA)

xCTbar <- matrix(c(sum(EQICT[,1])/length(EQICT[,1]),
                   sum(EQICT[,2])/length(EQICT[,2]),
                   sum(EQICT[,3])/length(EQICT[,3]),
                   sum(EQICT[,4])/length(EQICT[,4]),
                   sum(EQICT[,5])/length(EQICT[,5])), ncol=5)
                   #sum(EQICT[,6])/length(EQICT[,6])), ncol = 6)
for(i in 1:5){
  for(j in 1:5){
    SCT[i,j]<- sum((EQICT[,i]-xCTbar[,i])*(EQICT[,j]-xCTbar[,j]))/(length(EQICT[,1])-1)
  }
}
dets[3] <- det(SCT)

xWIbar <- matrix(c(sum(EQIWI[,1])/length(EQIWI[,1]),
                   sum(EQIWI[,2])/length(EQIWI[,2]),
                   sum(EQIWI[,3])/length(EQIWI[,3]),
                   sum(EQIWI[,4])/length(EQIWI[,4]),
                   sum(EQIWI[,5])/length(EQIWI[,5])), ncol=5)
                   #sum(EQIWI[,6])/length(EQIWI[,6])), ncol = 6)
for(i in 1:5){
  for(j in 1:5){
    SWI[i,j]<- sum((EQIWI[,i]-xWIbar[,i])*(EQIWI[,j]-xWIbar[,j]))/(length(EQIWI[,1])-1)
  }
}
dets[4] <- det(SWI)
(dets)



nevada <- EQIDATA[EQIDATA$state=="NV",]
orderedNV <- nevada[order(nevada$cat_rucc),]
condensed <- orderedNV[,c(2,5:10)]


nv <- c("Clark County", "Storey County", "Washoe County", "Carson City", "Douglas County", "Elko County", "Churchill County", "Humboldt County", "Lander County","Lyon County", "Mineral County", "Nye County", "White Pine County", "Esmeralda County", "Eureka County", "Lincoln County", "Pershing County")
f <- as.matrix(condensed[,2:7])
f <- t(f)
g <- data.frame(f)
names(g) <- nv
(C <- cor(g))
corrplot(C)
png("corr.png")
corrplot(C)
dev.off()

rQ <- c()
for(i in 1:6){
  xo <- sort(EQI[,i])
  n <- length(xo)
  problevels <- ((1:n)-0.5)/n
  q <- qnorm(problevels)
  rQ[i] <- sum((xo-mean(xo))*(q-mean(q)))/(sqrt(sum((xo-mean(xo))**2))*sqrt(sum((q-mean(q))**2)))
}
(rQ)
qqnorm(EQI$air_EQI_22July2013)
qqnorm(EQI$water_EQI_22July2013)
qqnorm(EQI$land_EQI_22July2013)
png("lowqq.png")
qqnorm(EQI$land_EQI_22July2013)
dev.off()
qqnorm(EQI$sociod_EQI_22July2013)
png("highqq.png")
qqnorm(EQI$sociod_EQI_22July2013)
dev.off()
qqnorm(EQI$built_EQI_22July2013)
qqnorm(EQI$EQI_22July2013)
cutoff <- 0.9953
which(rQ<cutoff)
