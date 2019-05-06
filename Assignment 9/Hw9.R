#setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Assignment 9")
#9.1
rm(list=ls())
(p <- matrix(c(1.0, .63, .45,
              .63, 1.0, .35,
              .45, .35, 1.0), byrow=TRUE, ncol=3))
(Cov <- matrix(c(.19,0,0,
                0,.51,0,
                0,0,.75),byrow=TRUE, ncol=3))
(L <- matrix(c(.9,
              .7,
              .5),byrow=TRUE, ncol = 1))
(L%*%t(L)+Cov)
(err<- p-(L%*%t(L)+Cov))

#9.2
(h12 <- L[1]^2)
(h22 <- L[2]^2)
(h32 <- L[3]^2)

#9.3
lambda <- eigen(p)$values
vec <- eigen(p)$vectors
(L1 <- sqrt(lambda[1])*-vec[,1])
L1 <- as.matrix(L1)
(Psi <- p-(L1)%*%t(L1))
(lambda[1]/sum(lambda))

#9.4
(ptild = p-Cov)
lambdatild <- eigen(ptild)$values
vectild <- eigen(ptild)$vectors
(sqrt(lambdatild[1])*-vectild[,1])

#9.12
rm(list=ls())
S <- matrix(c(11.072, 8.019, 8.160,
              8.019, 6.417, 6.005,
              8.160,6.005,6.773),byrow = TRUE, ncol=3)
S<-10**-3*S
L <- matrix(c(.1022,
              .0752,
              .0765), byrow=TRUE, ncol = 1)
LL <- L%*%t(L)
n=24
Sn=S*(n-1)/n
Psi <- Sn-LL
Psi <- diag(Psi)
Psi <- diag(Psi)
(L^2)
sum(L^2)/sum(diag(Sn))
(Resid <- Sn-LL-Psi)

#9.13
p=3
m=1
(testStat <- (n-1-(2*p+4*m+5)/6)*log(det(LL+Psi)/det(Sn)))
#Test Condition m<(2*p+1-sqrt(8*p+1))/2
(2*p+1-sqrt(8*p+1))/2

#9.20
rm(list=ls())
X<-read.table("T1-5.dat", header = FALSE)
X<-data.frame(X[,1],X[,2],X[,5],X[,6])
X<-as.matrix(X)
n <- nrow(X)
p <- ncol(X)
(xbar <- t(matrix(1,ncol=n) %*% X)/n)
(D <- X - matrix(1,nrow=n) %*% t(xbar))
(S <- (n-1)^(-1) * t(D)%*%D)

lambda <- eigen(S)$values
vec <- eigen(S)$vectors
(sqrt(lambda[1])*-vec[,1])
((sqrt(lambda[1])*-vec[,1])^2)
(lambda[1])
(lambda[1]/sum(lambda))
(sqrt(lambda[2])*-vec[,2])
((sqrt(lambda[1])*-vec[,1])^2+(sqrt(lambda[2])*-vec[,2])^2)
(lambda[2]) 
((lambda[1]+lambda[2])/sum(lambda))
source("factanal2.R")
m1=factanal2(X, 1, rotation = "none")
m2=factanal2(X, 2, rotation = "none")
print(m2$loadings,cutoff=0)

Q = matrix(c((sqrt(lambda[1])*-vec[,1]),(sqrt(lambda[2])*-vec[,2])), ncol=2)
(Rot<- varimax(Q))
(A<-Q%*%Rot$rotmat)

(m22=factanal2(X, 2, rotation = "varimax"))

(R <- cov2cor(S))
lambdaR <- eigen(R)$values
vecR <- eigen(R)$vectors
(sqrt(lambdaR[1])*-vecR[,1])
((sqrt(lambdaR[1])*-vecR[,1])^2)
(lambdaR[1])
(lambdaR[1]/sum(lambdaR))
(sqrt(lambdaR[2])*-vecR[,2])
((sqrt(lambdaR[1])*-vecR[,1])^2+(sqrt(lambdaR[2])*-vecR[,2])^2)
(lambdaR[2]) 
((lambdaR[2])/sum(lambdaR))

m1R=factanal2(R, 2, rotation = "none")
