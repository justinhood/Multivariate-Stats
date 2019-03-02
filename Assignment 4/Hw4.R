setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Assignment 4")
### 4.18 ###
X <- matrix(c(3,6,
              4,4,
              5,7,
              4,7), byrow=TRUE, ncol = 2)

xbar <- matrix(c(mean(X[,1]),
                      mean(X[,2])), byrow=TRUE, ncol = 1)

(S <- matrix(c(0,0,
               0,0),ncol=2))

for (i in 1:4) {
  S=S+(X[i,]-xbar)%*%t(X[i,]-xbar)
}
(S)

## 4.25 ##
rm(list = ls())
X1 <- matrix(c(108.28,152.36,95.04,65.45,62.97,263.99,265.19,285.06,92.01,165.68),ncol=1)
X2 <- matrix(c(17.05,16.59,10.91,14.14,9.52,25.33,18.54,15.73,8.10,11.13),ncol=1)
X3 <- matrix(c(1484.1,750.33,766.42,1110.46,1031.29,195.26,193.83,191.11,1175.16,211.15),ncol=1)
X <- cbind(X1,X2,X3)

n <- nrow(X)
p <- ncol(X)
(x.bar <- t(matrix(1,ncol=n) %*% X)/n)
(D <- X - matrix(1,nrow=n) %*% t(x.bar))
(S <- (n-1)^(-1) * t(D)%*%D)
(Sinv <- solve(S))

d <- c()
for (i in 1:n){
  d[i] <- (X[i,]-t(x.bar))%*%Sinv%*%t(X[i,]-t(x.bar))
}
print(d)

chisq.quantiles <- c(.3518,.7978,1.2125,1.6416,2.1095,2.6430,3.2831,4.1083,5.317,7.8147)
plot(y=sort(d),x=chisq.quantiles,pch=19)
png("425.png")
plot(y=sort(d),x=chisq.quantiles,pch=19)
dev.off()

### 4.25 ###
rm(list = ls())
X <- scan("T1-5.dat")
X <- matrix(X,ncol = 7, byrow=TRUE)
Xadj <- X[,5:6]
X
Xadj

n <- nrow(Xadj)
p <- ncol(Xadj)
(x.bar <- t(matrix(1,ncol=n) %*% Xadj)/n)
(D <- Xadj - matrix(1,nrow=n) %*% t(x.bar))
(S <- (n-1)^(-1) * t(D)%*%D)
(Sinv <- solve(S))

d <- c()
for (i in 1:n){
  d[i] <- (Xadj[i,]-t(x.bar))%*%Sinv%*%t(Xadj[i,]-t(x.bar))
}
print(d)
(cutoff <- qchisq(0.50,2))
which(d<cutoff)
sum(d<cutoff)/n # 40% of data fall within 50% contour

chisq.quantiles <- qchisq(((1:n)-0.5)/n,df=2)
plot(y=sort(d),x=chisq.quantiles,pch=19)
png("429.png")
plot(y=sort(d),x=chisq.quantiles,pch=19)
dev.off()
