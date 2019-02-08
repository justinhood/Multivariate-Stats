setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Assignment 2")
## 2.1 ##
x <- matrix(c(5,1,3),ncol=1)
y <- matrix(c(-1,3,1),ncol=1)
Lx <- sqrt(t(x) %*% x)
Ly <- sqrt(t(y) %*% y)
ctheta <- (t(x)%*%y)/(Lx*Ly)
theta <- acos(ctheta)
deg=theta*180/pi
(prj <- (t(x)%*%y)/(Ly**2))
prj <- prj[1,1] * y


## 2.5 ##
rm(list = ls())
library(MASS)
Q <- matrix(c(5/13, 12/13,
              -12/13, 5/13),byrow = TRUE, ncol=2)
fractions(Q)
(t(Q)%*%Q)
(Q%*%t(Q))
t(Q)%*%Q == Q%*%t(Q)

## 2.6 ##
rm(list = ls())
A <- matrix(c(9, -2,
              -2, 6), byrow = TRUE, ncol=2)
A == t(A)
eigen(A)$values

## 2.8 ##
rm(list = ls())
A <- matrix(c(1, 2,
              2, -2), byrow = TRUE, ncol=2)
(lambda <- eigen(A)$values)
(vecs <- eigen(A)$vectors)
(spec <- lambda[1]*vecs[,1]%*%t(vecs[,1])+
  lambda[2]*vecs[,2]%*%t(vecs[,2]))
A == round(spec)

## 2.9 ##
(inver=solve(A))
fractions(inver)
(inver %*% A)
(inverlambda <- eigen(inver)$values)
fractions(inverlambda)
(invere <- eigen(inver)$vectors)
(vecs)

## 2.15 ##
rm(list = ls())
A <- matrix(c(3, -1,
              -1, 3), byrow = TRUE, ncol=2)
(eigen(A)$values)

## 2.18 ##
rm(list = ls())
A <- matrix(c(4, -sqrt(2),
              -sqrt(2), 3), byrow = TRUE, ncol=2)
(eigen(A)$values)
(eigen(A)$vectors)

## 2.21 ##
rm(list = ls())
A <- matrix(c(1,1,
              2,-2,
              2,2), byrow = TRUE, ncol=2)
B <- t(A)%*%A
(eigen(B)$values)
(eigen(B)$vectors)

(C <- A%*%t(A))
(eigen(C)$values)
(eigen(C)$vectors)

## 2.25 ##
rm(list = ls())
sig <- matrix(c(25, -2, 4,
              -2,4,1,
              4,1,9), byrow = TRUE, ncol=3)
(vhalf <- matrix(rep(NA,3*3),ncol=3))
for(i in 1:3){
  for(j in 1:3){
    if(i == j){
    vhalf[i,j] = sqrt(sig[i,i])
    } else{
      vhalf[i,j] = 0
    }
  }
}
(vhalf)
(vinv=solve(vhalf))
rho <- matrix(rep(NA,3*3),ncol=3)
for(i in 1:3){
  for(j in 1:3){
    rho[i,j]=sig[i,j]/(sqrt(sig[i,i])*sqrt(sig[j,j]))
  }
}
(rho)
fractions(rho)
(vhalf %*% rho %*% vhalf)
(sig == round(vhalf %*% rho %*% vhalf))
