# setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Assignment 8")

## 8.2
rm(list=ls())
Sig=matrix(c(5,2,
             2,2), byrow = TRUE, ncol=2)
(Rho <- cov2cor(Sig))
(l.s <- eigen(Rho)$values)
(e.s <- eigen(Rho)$vectors)
(explained <- l.s/sum(l.s))

(l.s2 <- eigen(Sig)$values)
(e.s2 <- eigen(Sig)$vectors)
(explained2 <- l.s2/sum(l.s2))

(e.s[1,1]*sqrt(l.s[1]))
(e.s[2,1]*sqrt(l.s[1]))
(abs(e.s[1,2]*sqrt(l.s[2])))

## 8.6
rm(list=ls())
S <- matrix(c(7476.45, 303.62,
              303.62, 26.19),byrow=TRUE, ncol=2)
xbar <- matrix(c(155.60,
                 14.70), byrow=TRUE, ncol=1)

(l.s <- eigen(S)$values)
(e.s <- eigen(S)$vectors)
(explained <- l.s/sum(l.s))
library(plotrix)
c2 <- 1.4
angle <- atan(eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='X_2',xlab='X_1',xlim=c(45,300),ylim=c(5,25))
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=xbar[1],y=xbar[2],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)
png("ellip.png")
plot(0,pch='',ylab='X_2',xlab='X_1',xlim=c(45,300),ylim=c(5,25))
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=xbar[1],y=xbar[2],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)
dev.off

(Rho <- cov2cor(S))

## 8.7
(l.s <- eigen(Rho)$values)
(e.s <- eigen(Rho)$vectors)
(explained <- l.s/sum(l.s))

(e.s[1,1]*sqrt(l.s[1]))
(e.s[2,1]*sqrt(l.s[1]))

## 8.11
rm(list=ls())
X <- read.table("T8-5.DAT", header = FALSE)
X <- as.matrix(X)
X[,5] <- 10*X[,5]
n <- nrow(X)
p <- ncol(X)
(xbar <- t(matrix(1,ncol=n) %*% X)/n)
(D <- X - matrix(1,nrow=n) %*% t(xbar))
(S <- (n-1)^(-1) * t(D)%*%D)
Sinv <- solve(S)

(l.s <- eigen(S)$values)
(e.s <- eigen(S)$vectors)
(explained <- l.s/sum(l.s))
(cumulative <- explained)
for (i in 1:length(explained)){
  summer=0
  for (j in 1:i){
    summer=summer+explained[j]
  }
  cumulative[i]=summer
}
(cumulative)
(r11<-e.s[1,1]*sqrt(l.s[1])/sqrt(S[1,1]))
(r12<-e.s[2,1]*sqrt(l.s[1])/sqrt(S[2,2]))
(r13<-e.s[3,1]*sqrt(l.s[1])/sqrt(S[3,3]))
(r14<-e.s[4,1]*sqrt(l.s[1])/sqrt(S[4,4]))
(r15<-e.s[5,1]*sqrt(l.s[1])/sqrt(S[5,5]))

(r21<-e.s[1,2]*sqrt(l.s[2])/sqrt(S[1,1]))
(r22<-e.s[2,2]*sqrt(l.s[2])/sqrt(S[2,2]))
(r23<-e.s[3,2]*sqrt(l.s[2])/sqrt(S[3,3]))
(r24<-e.s[4,2]*sqrt(l.s[2])/sqrt(S[4,4]))
(r25<-e.s[5,2]*sqrt(l.s[2])/sqrt(S[5,5]))

## 8.14
rm(list=ls())
X <- read.table("T5-1.DAT", header = FALSE)
X <- as.matrix(X)
S <- matrix(c(2.879, 10.010, -1.810,
              10.010, 199.788, -5.640,
              -1.810, -5.640, 3.628), byrow = TRUE, ncol=3)

(l <- eigen(S)$values)
(e <- eigen(S)$vectors)

l/sum(l)
plot(l/sum(l))
png("eigprop.png")
plot(l/sum(l))
dev.off()

e[,1]*sqrt(l[1])/sqrt(diag(S))
e[,2]*sqrt(l[2])/sqrt(diag(S))
e[,3]*sqrt(l[3])/sqrt(diag(S))

qqnorm(X %*% e[,1])
png("1qq.png")
qqnorm(X %*% e[,1])
dev.off()

## 8.20
rm(list=ls())
X2 <- read.table("T8-6.DAT", header = FALSE)
X <-X2[,2:9]
X <- as.matrix(X)
n <- nrow(X)
p <- ncol(X)
(xbar <- t(matrix(1,ncol=n) %*% X)/n)
(D <- X - matrix(1,nrow=n) %*% t(xbar))
(S <- (n-1)^(-1) * t(D)%*%D)
Sinv <- solve(S)
(R <- cov2cor(S))
(l.s <- eigen(R)$values)
(e.s <- eigen(R)$vectors)
(explained <- l.s/sum(l.s))
(cumulative <- explained)
for (i in 1:length(explained)){
  summer=0
  for (j in 1:i){
    summer=summer+explained[j]
  }
  cumulative[i]=summer
}
(cumulative)
plot(l.s/sum(l.s), ylim = c(0,1))
lines(l.s/sum(l.s))
lines(cumulative, col="red")
png("cumvval.png")
plot(l.s/sum(l.s), ylim = c(0,1))
lines(l.s/sum(l.s))
lines(cumulative, col="red")
dev.off()

ranks <- X %*% -e.s[,1]
Y=data.frame(X2[,1],ranks)
(Y=Y[order(Y[,2]),])
### Womens data
X2 <- read.table("T1-9.dat", header = FALSE)
X <-X2[,2:8]
X <- as.matrix(X)
n <- nrow(X)
p <- ncol(X)
(xbar <- t(matrix(1,ncol=n) %*% X)/n)
(D <- X - matrix(1,nrow=n) %*% t(xbar))
(S <- (n-1)^(-1) * t(D)%*%D)
Sinv <- solve(S)
(R <- cov2cor(S))
(l.s <- eigen(R)$values)
(e.s <- eigen(R)$vectors)
(explained <- l.s/sum(l.s))
(cumulative <- explained)
for (i in 1:length(explained)){
  summer=0
  for (j in 1:i){
    summer=summer+explained[j]
  }
  cumulative[i]=summer
}
(cumulative)
plot(l.s/sum(l.s), ylim = c(0,1))
lines(l.s/sum(l.s))
lines(cumulative, col="red")


ranks <- X %*% -e.s[,1]
Y=data.frame(X2[,1],ranks)
(Y=Y[order(Y[,2]),])
## 8.21
rm(list=ls())
X2 <- read.table("T8-6.DAT", header = FALSE)
X <-X2[,2:9]
X[,1]=100/X[,1]
X[,2]=200/X[,2]
X[,3]=400/X[,3]
X[,4]=800/(X[,4]*60)
X[,5]=1500/(X[,5]*60)
X[,6]=5000/(X[,6]*60)
X[,7]=10000/(X[,7]*60)
X[,8]=42195/(X[,8]*60)

X <- as.matrix(X)
n <- nrow(X)
p <- ncol(X)
(xbar <- t(matrix(1,ncol=n) %*% X)/n)
(D <- X - matrix(1,nrow=n) %*% t(xbar))
(S <- (n-1)^(-1) * t(D)%*%D)
Sinv <- solve(S)
(l.s <- eigen(S)$values)
(e.s <- eigen(S)$vectors)
(explained <- l.s/sum(l.s))
(cumulative <- explained)
for (i in 1:length(explained)){
  summer=0
  for (j in 1:i){
    summer=summer+explained[j]
  }
  cumulative[i]=summer
}
(cumulative)
plot(l.s/sum(l.s), ylim = c(0,1))
lines(l.s/sum(l.s))
lines(cumulative, col="red")

ranks <- X %*% -e.s[,1]
Y=data.frame(X2[,1],ranks)
(Y=Y[order(Y[,2]),])
