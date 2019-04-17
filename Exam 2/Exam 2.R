#setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Exam 2")

# 1
rm(list=ls())
X <- read.csv("./data/prodimp.csv",header = TRUE) #Get data



library(ggplot2) #Q plots 
qplot(X$level,X$productivity, main="qPlot of Productivity Data", xlab="Treatment", ylab="Productivity Improvement")
png("prodq.png")
qplot(X$level,X$productivity, main="qPlot of Productivity Data", xlab="Treatment", ylab="Productivity Improvement")
dev.off()

#make the mdoel
fit<- aov(X$productivity~X$level)
summary(fit)
#Factor and then recompute
X$level <- factor(X$level)
fit2<- aov(X$productivity~X$level)
summary(fit2)

#make the intervals
TukeyHSD(fit2, ordered = FALSE, conf.level = 0.95)
plot(TukeyHSD(fit2, ordered = FALSE, conf.level = 0.95))
png("prodtukey.png")
plot(TukeyHSD(fit2, ordered = FALSE, conf.level = 0.95))
dev.off()
#What are the means for comparison
(m1=mean(X$productivity[which(X$level==1)]))
(m2=mean(X$productivity[which(X$level==2)]))
(m3=mean(X$productivity[which(X$level==3)]))

#******************************************************************#
# 2
rm(list=ls())
X <- read.table("./data/T5-14.dat",header = FALSE)
X <- as.matrix(X)
first <- X[1:30,] #Just the good data

n <- nrow(first)
p <- ncol(first)

#Basic set of computations. Should probably just write a module for these
(xbar <- t(matrix(1,ncol=n) %*% first)/n)
(D <- first - matrix(1,nrow=n) %*% t(xbar))
(S <- (n-1)^(-1) * t(D)%*%D)
Sinv <- solve(S)
#D's
d <- c()
for (i in 1:length(X[,1])){
  d[i] <- (X[i,]-t(xbar))%*%Sinv%*%t(X[i,]-t(xbar))
}
#Chi Cutoffs
UCL95 <- qchisq(.95,df=p)
UCL99 <- qchisq(.99,df=p)
#For The In Spec Data Points
#################################################################
plot(x=1:n,y=d[1:n],ylim=c(0,max(UCL99,ceiling(max(d[1:n])))),main="In Spec T2", xlab="Sample",ylab="T2")
abline(h=UCL95,lty=3, col='blue')
abline(h=UCL99, lty=1, col='red')

png("inspec.png")
plot(x=1:n,y=d[1:n],ylim=c(0,max(UCL99,ceiling(max(d[1:n])))),main="In Spec T2",xlab="Sample",ylab="T2")
abline(h=UCL95,lty=3, col='blue')
abline(h=UCL99, lty=1, col='red')
dev.off()
####################################################################
#For the Out of Spec?
###################################################################
plot(x=1:length(d),y=d,ylim=c(0,max(UCL99,ceiling(max(d)))),main="All data T2", xlab="Sample",ylab="T2")
abline(h=UCL95,lty=3, col='blue')
abline(h=UCL99, lty=1, col='red')

png("outspec.png")
plot(x=1:length(d),y=d,ylim=c(0,max(UCL99,ceiling(max(d)))),main="All data T2", xlab="Sample",ylab="T2")
abline(h=UCL95,lty=3, col='blue')
abline(h=UCL99, lty=1, col='red')
dev.off()
#Find the bad data points
(which(d>UCL99))
#T test time
alpha=.05
mu0 <- matrix(c(0,0,0,0,0,0),ncol=1) #Zero vector
(T2 <- n*t(xbar-mu0)%*%Sinv%*%(xbar-mu0))
(Tcrit <- p*(n-1)*qf(1-alpha,p,n-p)/(n-p))
if(T2>Tcrit){
  print("Reject H0")
} else{
  print("Fail to Reject")
}

library(plotrix)
#Do the plots with all the data
dat <- X

# 1v2
X <- as.matrix(dat[,1:2])
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D

alpha <- 0.01
c2 <- qchisq(1-alpha,df=2)


angle <- atan(eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='X_2',xlab='X_1',xlim=c(-2,1),ylim=c(-1,.75))
points(X)
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)
png("12.png")
plot(0,pch='',ylab='X_2',xlab='X_1',xlim=c(-2,1),ylim=c(-1,.75))
points(X)
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)
dev.off()

#1v3
X <- as.matrix(dat[,c(1,3)])
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D

alpha <- 0.01
c2 <- qchisq(1-alpha,df=2)

angle <- atan(eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='X_3',xlab='X_1',xlim=c(-1.5,.5),ylim=c(-1.25,1.25))
points(X)
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)

# 1v4
X <- as.matrix(dat[,c(1,4)])
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D

alpha <- 0.01
c2 <- qchisq(1-alpha,df=2)

angle <- atan(eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='X_4',xlab='X_1',xlim=c(-1.25,.5),ylim=c(-1.25,1.25))
points(X)
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)

#1v5
X <- as.matrix(dat[,c(1,5)])
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D

alpha <- 0.01
c2 <- qchisq(1-alpha,df=2)

angle <- atan(eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='X_5',xlab='X_1',xlim=c(-1.25,.5),ylim=c(-1,2.25))
points(X)
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)


#1v6
X <- as.matrix(dat[,c(1,6)])
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D

alpha <- 0.01
c2 <- qchisq(1-alpha,df=2)

angle <- atan(eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='X_6',xlab='X_1',xlim=c(-1.25,.5),ylim=c(-1,1))
points(X)
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)

png("16.png")
plot(0,pch='',ylab='X_6',xlab='X_1',xlim=c(-1.25,.5),ylim=c(-1,1))
points(X)
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)
dev.off()


#2v6
X <- as.matrix(dat[,c(2,6)])
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D

alpha <- 0.01
c2 <- qchisq(1-alpha,df=2)

angle <- atan(eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='X_6',xlab='X_2',xlim=c(-1.25,.75),ylim=c(-.75,.85))
points(X)
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)


###################################################################
#3
rm(list=ls())
#Given calls
library(reshape2)
dat.fill <- read.csv("./data/BHO_Medication_Fill_Data__2010-2014.csv",header=TRUE)
tmp <- dcast(dat.fill, OMH.Region + Year + Quarter + Age.Group ~ Description.of.Metric)
dat <- tmp[,c(1:4,6,8,10)]
names(dat) <- c(names(dat)[1:4],"Mood30","Psychotrop30","Antipsych30")
head(dat)

# Look at dem histograms,
### Mood
hist(dat[,5],breaks=20,main="Mood30",xlab="% Refill")
png("moodhist.png")
hist(dat[,5],breaks=20,main="Mood30",xlab="% Refill")
dev.off()
### Psycho
hist(dat[,6],breaks=20,main="Psych30",xlab="% Refill")
png("psychhist.png")
hist(dat[,6],breaks=20,main="Psych30",xlab="% Refill")
dev.off()
## Anti
hist(dat[,7],breaks=20,main="Anti30",xlab="% Refill")
png("antihist.png")
hist(dat[,7],breaks=20,main="Anti30",xlab="% Refill")
dev.off()

#Q plot time
#QQ Mood
xmood <- sort(dat$Mood30)
n <- length(xmood)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xmood,pch=19)
png("moodq.png")
plot(q,xmood,pch=19)
dev.off()
#Compute rq statistic from text
q.bar <- mean(q)
x.bar <- mean(xmood)
critval95 <- 0.9913 # from table 4.2 p. 181
(rqmood <- sum((xmood-x.bar)*(q-q.bar))/(sqrt(sum((xmood-x.bar)^2))*sqrt(sum((q-q.bar)^2))))
if(rqmood>critval95){
  print("Its Normal Yo")
} else{
  print("Nice try bro")
}

## QQ Psycho
xpsych <- sort(dat$Psychotrop30)
n <- length(xpsych)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xpsych,pch=19)
png("psychq.png")
plot(q,xpsych,pch=19)
dev.off()
#Compute rq statistic from text
q.bar <- mean(q)
x.bar <- mean(xpsych)
critval95 <- 0.9913 # from table 4.2 p. 181
(rqpsych <- sum((xpsych-x.bar)*(q-q.bar))/(sqrt(sum((xpsych-x.bar)^2))*sqrt(sum((q-q.bar)^2))))
if(rqpsych>critval95){
  print("Its Normal Yo")
} else{
  print("Nice try bro")
}

## QQ Anti
xanti <- sort(dat$Antipsych30)
n <- length(xanti)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xanti,pch=19)
png("antiq.png")
plot(q,xanti,pch=19)
dev.off()
#Compute rq statistic from text
q.bar <- mean(q)
x.bar <- mean(xanti)
critval95 <- 0.9913 # from table 4.2 p. 181
(rqanti <- sum((xanti-x.bar)*(q-q.bar))/(sqrt(sum((xanti-x.bar)^2))*sqrt(sum((q-q.bar)^2))))
if(rqanti>critval95){
  print("Its Normal Yo")
} else{
  print("Nice try bro")
}

#Make the Chi Square plot now
n <- nrow(dat)
p <- ncol(dat)

reduced <- as.matrix(dat[,5:7])
(xbar <- t(matrix(1,ncol=n) %*% reduced)/n)
(D <- reduced - matrix(1,nrow=n) %*% t(xbar))
(S <- (n-1)^(-1) * t(D)%*%D)
Sinv <- solve(S)

d <- c()
for (i in 1:length(reduced[,1])){
  d[i] <- (reduced[i,]-t(xbar))%*%Sinv%*%t(reduced[i,]-t(xbar))
}

chisq.quantiles <- qchisq(((1:n)-0.5)/n,df=3)
f1<- summary(lm(sort(d)~chisq.quantiles))
plot(y=sort(d),x=chisq.quantiles,pch=19)
abline(f1$coefficients[1],f1$coefficients[2],col="red")
png("originalchi.png")
plot(y=sort(d),x=chisq.quantiles,pch=19)
abline(f1$coefficients[1],f1$coefficients[2],col="red")
dev.off()

#Looks like we have an outlier
which.max(d)
((dat[30,7]-mean(dat[,7]))/sqrt(var(dat[,7])))

#Because we cant see the data, lets just drop it from the analysis
newTest <- as.matrix(dat[,5:7])
newTest <- newTest[-30,]

#new Chi square
n <- nrow(newTest)
p <- ncol(newTest)

(xbar <- t(matrix(1,ncol=n) %*% newTest)/n)
(D <- newTest - matrix(1,nrow=n) %*% t(xbar))
(S <- (n-1)^(-1) * t(D)%*%D)
Sinv <- solve(S)

d <- c()
for (i in 1:length(newTest[,1])){
  d[i] <- (newTest[i,]-t(xbar))%*%Sinv%*%t(newTest[i,]-t(xbar))
}

chisq.quantiles <- qchisq(((1:n)-0.5)/n,df=3)
f2<-summary(lm(sort(d)~chisq.quantiles))
plot(y=sort(d),x=chisq.quantiles,pch=19)
abline(f2$coefficients[1],f2$coefficients[2],col="red")
png("newchi.png")
plot(y=sort(d),x=chisq.quantiles,pch=19)
abline(f2$coefficients[1],f2$coefficients[2],col="red")
dev.off()
#Better, lets check the qq stats

# NewQQ Mood
xmood <- sort(newTest[,1])
n <- length(xmood)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xmood,pch=19)
png("newmoodq.png")
plot(q,xmood,pch=19)
dev.off()
q.bar <- mean(q)
x.bar <- mean(xmood)
critval95 <- 0.9913 # from table 4.2 p. 181
(rqmood <- sum((xmood-x.bar)*(q-q.bar))/(sqrt(sum((xmood-x.bar)^2))*sqrt(sum((q-q.bar)^2))))
if(rqmood>critval95){
  print("Its Normal Yo")
} else{
  print("Nice try bro")
}
# New QQ Psycho
xpsych <- sort(newTest[,2])
n <- length(xpsych)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xpsych,pch=19)
png("newpsychq.png")
plot(q,xpsych,pch=19)
dev.off()
q.bar <- mean(q)
x.bar <- mean(xpsych)
critval95 <- 0.9913 # from table 4.2 p. 181
(rqpsych <- sum((xpsych-x.bar)*(q-q.bar))/(sqrt(sum((xpsych-x.bar)^2))*sqrt(sum((q-q.bar)^2))))
if(rqpsych>critval95){
  print("Its Normal Yo")
} else{
  print("Nice try bro")
}

# New QQ Anti
xanti <- sort(newTest[,3])
n <- length(xanti)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xanti,pch=19)
png("newantiq.png")
plot(q,xanti,pch=19)
dev.off()
q.bar <- mean(q)
x.bar <- mean(xanti)
critval95 <- 0.9913 # from table 4.2 p. 181
(rqanti <- sum((xanti-x.bar)*(q-q.bar))/(sqrt(sum((xanti-x.bar)^2))*sqrt(sum((q-q.bar)^2))))
if(rqanti>critval95){
  print("Its Normal Yo")
} else{
  print("Nice try bro")
}
#Psych and Anti still not quite normal. Look at power scaling

##Mood is good, alter Psych
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
n <- nrow(as.matrix(newTest[,2]))
lambdas <- seq(from=-1,to=10,by=.01)
l.lambdas <- c()
for (i in 1:length(lambdas)){
  l.lambdas[i] <- l.lambda(lambdas[i],n,newTest[,2])
}
plot(lambdas,l.lambdas)
(pstar=lambdas[which(l.lambdas==max(l.lambdas))]) # this is the transformation we should use

newPsych=(-1+newTest[,2]^pstar)/pstar
xpsych <- sort(newPsych)
n <- length(xpsych)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xpsych,pch=19)
png("scalepq.png")
plot(q,xpsych,pch=19)
dev.off()
q.bar <- mean(q)
x.bar <- mean(xpsych)
critval95 <- 0.9913 # from table 4.2 p. 181
(rqpsych2 <- sum((xpsych-x.bar)*(q-q.bar))/(sqrt(sum((xpsych-x.bar)^2))*sqrt(sum((q-q.bar)^2))))
if(rqpsych2>critval95){
  print("Its Normal Yo")
} else{
  print("Nice try bro")
}
##Now its normal. Good.

##Now that psych is good, do anti
n <- nrow(as.matrix(newTest[,3]))
lambdas <- seq(from=-1,to=10,by=.01)
l.lambdas <- c()
for (i in 1:length(lambdas)){
  l.lambdas[i] <- l.lambda(lambdas[i],n,newTest[,3])
}
plot(lambdas,l.lambdas)
(astar=lambdas[which(l.lambdas==max(l.lambdas))]) # this is the transformation we should use

newAnti = (-1+newTest[,3]^astar)/astar
xanti <- sort(newAnti)
n <- length(xanti)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xanti,pch=19)
png("scaleaq.png")
plot(q,xanti,pch=19)
dev.off()
q.bar <- mean(q)
x.bar <- mean(xanti)
critval95 <- 0.9913 # from table 4.2 p. 181
(rqanti2 <- sum((xanti-x.bar)*(q-q.bar))/(sqrt(sum((xanti-x.bar)^2))*sqrt(sum((q-q.bar)^2))))
if(rqanti2>critval95){
  print("Its Normal Yo")
} else{
  print("Nice try bro")
}
#Its not quite perfect, but is as good as we can get (log is much worse)
############ 
#Recombine our best attempts at normalcy
bestTry = matrix(c(newTest[,1],newPsych, newAnti),ncol=3)
#Make the chi-square on this data
n <- nrow(bestTry)
p <- ncol(bestTry)

(xbar <- t(matrix(1,ncol=n) %*% bestTry)/n)
(D <- bestTry - matrix(1,nrow=n) %*% t(xbar))
(S <- (n-1)^(-1) * t(D)%*%D)
Sinv <- solve(S)

d <- c()
for (i in 1:length(bestTry[,1])){
  d[i] <- (bestTry[i,]-t(xbar))%*%Sinv%*%t(bestTry[i,]-t(xbar))
}

chisq.quantiles <- qchisq(((1:n)-0.5)/n,df=3)
f3 <- summary(lm(sort(d)~chisq.quantiles))
plot(y=sort(d),x=chisq.quantiles,pch=19)
abline(f3$coefficients[1],f3$coefficients[2],col="red")

png("bestchi.png")
plot(y=sort(d),x=chisq.quantiles,pch=19)
abline(f3$coefficients[1],f3$coefficients[2],col="red")
dev.off()
#This seems to be the best in terms of linarity, so we continue.

## Now that the data is scaled again, we do our test.
n <- nrow(bestTry)
p <- ncol(bestTry)
alpha=.05
#Base null mu
mu0 <- matrix(c(75,75,75), ncol=1)
mu0[2] <-(-1+mu0[2]^pstar)/pstar #Transform to new scale
mu0[3] <-(-1+mu0[3]^astar)/astar # Transform to its new scale

#Do the test
(T2 <- n*t(xbar-mu0)%*%Sinv%*%(xbar-mu0))
(Tcrit <- p*(n-1)*qf(1-alpha,p,n-p)/(n-p))
if(T2>Tcrit){
  print("Reject H0")
} else{
  print("Fail to Reject")
}

################################################################
#4
rm(list=ls())
dat1 <- read.csv("./data/TTHM.csv",header=TRUE)
dat2 <- read.csv("./data/HAA5.csv",header=TRUE)
dat <- merge(dat1,dat2,by=names(dat1)[1],all=FALSE)
dat <- dat[,c(1:4,7)]
names(dat) <- c("Date","Year","Quarter","TTHM","HAA5")
head(dat)


#Histograms
hist(dat$TTHM)
png("thist.png")
hist(dat$TTHM)
dev.off()
hist(dat$HAA5)
png("hhist.png")
hist(dat$HAA5)
dev.off()
n <- nrow(dat)

#Xbar Charts
center <- mean(dat$TTHM)
UL <- center + 3*sd(dat$TTHM)
LL <- center - 3*sd(dat$TTHM)

plot(x=1:n,dat$TTHM,type="l",ylim=c(0,110),xlab="Observation Number",ylab="TTHM")
abline(h=center, col="green");abline(h=UL, col="red");abline(h=max(0,LL), col="red")
png("tthmxbar.png")
plot(x=1:n,dat$TTHM,type="l",ylim=c(0,110),xlab="Observation Number",ylab="TTHM")
abline(h=center, col="green");abline(h=UL, col="red");abline(h=max(0,LL), col="red")
dev.off()

center <- mean(dat$HAA5)
UL <- center + 3*sd(dat$HAA5)
LL <- center - 3*sd(dat$HAA5)

plot(x=1:n,dat$HAA5,type="l",ylim=c(0,110),xlab="Observation Number",ylab="HAA5")
abline(h=center, col="green");abline(h=UL, col="red");abline(h=max(0,LL), col="red")
png("haa5xbar.png")
plot(x=1:n,dat$HAA5,type="l",ylim=c(0,110),xlab="Observation Number",ylab="HAA5")
abline(h=center, col="green");abline(h=UL, col="red");abline(h=max(0,LL), col="red")
dev.off()

library(plotrix)
#Ellipse Chart
X <- as.matrix(dat[,4:5])
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D

alpha <- 0.01
c2 <- qchisq(1-alpha,df=2)


angle <- atan(eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='X_2',xlab='X_1',xlim=c(-5,100),ylim=c(-5,90))
points(X)
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)
png("waterelip.png")
plot(0,pch='',ylab='HAA5',xlab='TTHM',xlim=c(-5,100),ylim=c(-5,90))
points(X)
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)
dev.off()

Sinv <- solve(S)
d <- c()
for (i in 1:length(X[,1])){
  d[i] <- (X[i,]-t(X.mean))%*%Sinv%*%t(X[i,]-t(X.mean))
}
#Chi Cutoffs
UCL95 <- qchisq(.95,df=p)
UCL99 <- qchisq(.99,df=p)
#For The In Spec Data Points
#################################################################
plot(x=1:n,y=d[1:n],ylim=c(0,max(UCL99,ceiling(max(d[1:n])))),main="In Spec T2", xlab="Sample",ylab="T2")
abline(h=UCL95,lty=3, col='blue')
abline(h=UCL99, lty=1, col='red')

png("watert.png")
plot(x=1:n,y=d[1:n],ylim=c(0,max(UCL99,ceiling(max(d[1:n])))),main="In Spec T2",xlab="Sample",ylab="T2")
abline(h=UCL95,lty=3, col='blue')
abline(h=UCL99, lty=1, col='red')
dev.off()


n <- nrow(dat)
p <- 2
alpha=.05
#Base null mu
mu0 <- matrix(c(80,60), ncol=1)

#Do the test
(T2 <- n*t(X.mean-mu0)%*%Sinv%*%(X.mean-mu0))
(Tcrit <- p*(n-1)*qf(1-alpha,p,n-p)/(n-p))
if(T2>Tcrit){
  print("Reject H0")
} else{
  print("Fail to Reject")
}
