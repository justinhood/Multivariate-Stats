setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Assignment 3")

## 3.2 ###
rm(list = ls())
X <- matrix(c(3,4,
              6,-2,
              3,1),byrow = TRUE, ncol=2)
xbar <- matrix(c(sum(X[,1])/3,
                 sum(X[,2])/3), byrow = TRUE, ncol = 2)
plot(X[,1],X[,2], col="blue", main = "Points and X bar")
points(xbar[,1],xbar[,2], col="red")

png("32a.png")
plot(X[,1],X[,2], col="blue", main = "Points and X bar")
points(xbar[,1],xbar[,2], col="red")
dev.off()

d1=X[,1]-matrix(c(4,
                  4,
                  4), byrow=TRUE, ncol = 1)

d2=X[,2]-matrix(c(1,
                  1,
                  1), byrow=TRUE, ncol = 1)
library(rgl)
library(matlib)
library(MASS)
open3d()
vectors3d(X[,1],col="blue",lwd=2)
vectors3d(X[,2],col="blue",lwd=2)
vectors3d(t(d1),col="red",lwd=2)
vectors3d(t(d2),col="red",lwd=2)
axes3d() 

s11=(t(d1)%*%d1)/3
s22=(t(d2)%*%d2)/3
s12=(t(d1)%*%d2)/3

r12=s12/(sqrt(s11)*sqrt(s22))

## 3.5b ##
rm(list = ls())
(S <- matrix(rep(NA,2*2),ncol=2))
for(i in 1:2){
  for(j in 1:2){
    S[i,j]=sum((X[,i]-xbar[,i])*(X[,j]-xbar[,j]))/2
  }
}
(det(S))

## 3.9 ##
rm(list = ls())
X <- matrix(c(12,17,29,
              18,20,38,
              14,16,30,
              20,18,38,
              16,19,35), byrow=TRUE, ncol = 3)
xbar <- matrix(rep(NA,1*3),ncol=3)
for(i in 1:3){
  xbar[,i]=mean(X[,i])
}
xbar
Xadj=X
for(i in 1:5){
  for(j in 1:3){
    Xadj[i,j]=X[i,j]-xbar[,j]
  }
}
(Xadj)
a=matrix(c(1,1,-1),byrow=TRUE, ncol=3)
(Xadj%*%t(a))

(S <- matrix(rep(NA,3*3),ncol=3))
for(i in 1:3){
  for(j in 1:3){
    S[i,j]=sum((X[,i]-xbar[,i])*(X[,j]-xbar[,j]))/4
  }
}
(S)
(det(S))
(S%*%t(a))

## 3.14 ##
rm(list = ls())
X <- matrix(c(9,1,
              5,3,
              1,2), byrow=TRUE,ncol = 2)
b1=2*X[1,1]+3*X[1,2]
b2=2*X[2,1]+3*X[2,2]
b3=2*X[3,1]+3*X[3,2]
B=c(b1,b2,b3)
bmean=mean(B)
bvar=var(B)

c1=-1*X[1,1]+2*X[1,2]
c2=-1*X[2,1]+2*X[2,2]
c3=-1*X[3,1]+2*X[3,2]
C=c(c1,c2,c3)
cmean=mean(C)
cvar=var(C)

oldcov=cov(B,C)
##################
xbar=matrix(c(mean(X[,1]),mean(X[,2])),byrow = TRUE,ncol = 2)
Xadj=X
for (i in 1:3) {
  for (j in 1:2) {
    Xadj[i,j]=X[i,j]-xbar[,j]
  }
}
S=(t(Xadj)%*%Xadj)/2
b=matrix(c(2,3),byrow = TRUE,ncol=2)
c=matrix(c(-1,2),byrow = TRUE, ncol = 2)
(b%*%t(xbar))
(c%*%t(xbar))
(b%*%S%*%t(b))
(c%*%S%*%t(c))
(b%*%S%*%t(c))
