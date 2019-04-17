#setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Assignment 6")

###6.4
rm(list=ls())

dat <- read.table("T6-1.dat")
names(dat) <- c("BOD.com","SS.com","BOD.sta","SS.sta")
altered <- log(dat)

alpha <- 0.05
n <- nrow(dat)
p <- 2

dj1 <- altered$BOD.com - altered$BOD.sta
dj2 <- altered$SS.com - altered$SS.sta
d <- matrix(c(dj1,dj2),ncol=2,nrow=n)
(d.bar <- matrix(c(mean(dj1),mean(dj2)),ncol=1))
(Sd <- var(d))
Sinv<-solve(Sd)

(T2 <- n*t(d.bar)%*%solve(Sd)%*%d.bar)

(T2.crit <- p*(n-1)/(n-p)*qf(1-alpha,p,n-p))

ifelse(T2>T2.crit,"Reject H0","Fail to Reject H0")

#Simultaneous CI
d.bar[1,1] - sqrt(p*(n-1)/(n-p)*qf(1-alpha,p,n-p))*sqrt(Sd[1,1]/n)
d.bar[1,1] + sqrt(p*(n-1)/(n-p)*qf(1-alpha,p,n-p))*sqrt(Sd[1,1]/n)

d.bar[2,1] - sqrt(p*(n-1)/(n-p)*qf(1-alpha,p,n-p))*sqrt(Sd[2,2]/n)
d.bar[2,1] + sqrt(p*(n-1)/(n-p)*qf(1-alpha,p,n-p))*sqrt(Sd[2,2]/n)


#### Bonferroni
tstar=abs(qt(alpha/(2*p),n-1))
(lower=d.bar[1,1]-tstar*sqrt(Sd[1,1]/n))
(upper=d.bar[1,1]+tstar*sqrt(Sd[1,1]/n))

(lower=d.bar[2,1]-tstar*sqrt(Sd[2,2]/n))
(upper=d.bar[2,1]+tstar*sqrt(Sd[2,2]/n))

d12=t(d[1,]-d.bar[1,1])%*%solve(Sd)%*%(d[1,]-d.bar[1,1])

dist <- c()
for (i in 1:n){
  dist[i] <- (d[i,]-t(d.bar))%*%Sinv%*%t(d[i,]-t(d.bar))
}


chisq.quantiles <- qchisq(((1:n)-0.5)/n,df=2)

plot(y=sort(dist),x=chisq.quantiles,pch=19)
png("64chi.png")
plot(y=sort(dist),x=chisq.quantiles,pch=19)
dev.off()


### 6.5
rm(list=ls())
alpha=0.05
x.bar <- matrix(c(46.1,
                  57.3,
                  50.4),byrow=TRUE, ncol=1)
n=40
S <- matrix(c(101.3, 63.0, 71.0,
              63.0, 80.2, 55.6,
              71.0, 55.6, 97.4), byrow = TRUE, ncol=3)
C <- matrix(c(1,-1,0,
              0,1,-1), byrow = TRUE, ncol=3)
inv <- solve(C%*%S%*%t(C))
T2 <- n*t(C%*%x.bar)%*%inv%*%(C%*%x.bar)
(T2.crit <- (nrow(S)-1)*(n-1)/(n-nrow(S)+1)*qf(1-alpha,nrow(S)-1,n-nrow(S)+1))

t(C[1,])%*%x.bar-sqrt((n-1)*(nrow(S)-1)*qf(1-alpha,nrow(S)-1,n-nrow(S)+1)/(n+1-nrow(S)))*sqrt(t(C[1,])%*%S%*%C[1,]/n)
t(C[1,])%*%x.bar+sqrt((n-1)*(nrow(S)-1)*qf(1-alpha,nrow(S)-1,n-nrow(S)+1)/(n+1-nrow(S)))*sqrt(t(C[1,])%*%S%*%C[1,]/n)

t(C[2,])%*%x.bar-sqrt((n-1)*(nrow(S)-1)*qf(1-alpha,nrow(S)-1,n-nrow(S)+1)/(n+1-nrow(S)))*sqrt(t(C[2,])%*%S%*%C[2,]/n)
t(C[2,])%*%x.bar+sqrt((n-1)*(nrow(S)-1)*qf(1-alpha,nrow(S)-1,n-nrow(S)+1)/(n+1-nrow(S)))*sqrt(t(C[2,])%*%S%*%C[2,]/n)

C <- matrix(c(1,-1,0,
              0,1,-1,
              -1,0,1), byrow = TRUE, ncol=3)

t(C[3,])%*%x.bar-sqrt((n-1)*(nrow(S)-1)*qf(1-alpha,nrow(S)-1,n-nrow(S)+1)/(n+1-nrow(S)))*sqrt(t(C[3,])%*%S%*%C[3,]/n)
t(C[3,])%*%x.bar+sqrt((n-1)*(nrow(S)-1)*qf(1-alpha,nrow(S)-1,n-nrow(S)+1)/(n+1-nrow(S)))*sqrt(t(C[3,])%*%S%*%C[3,]/n)

### 6.16
rm(list=ls())
alpha=0.05
X <- read.table("T4-3.DAT")
X <- X[,1:4]
(x.bar <- matrix(c(mean(X[,1]),mean(X[,2]),mean(X[,3]),mean(X[,4])),ncol=1))
n=nrow(X)
(S <- matrix(rep(NA,4*4),ncol=4))
for(i in 1:4){
  for(j in 1:4){
    S[i,j]=sum((X[,i]-x.bar[i,])*(X[,j]-x.bar[j,]))/(n-1)
  }
}
(S)
C <- matrix(c(1,-1,0,0,
              0,1,-1,0,
              0,0,1,-1), byrow = TRUE, ncol=4)
inv <- solve(C%*%S%*%t(C))
T2 <- n*t(C%*%x.bar)%*%inv%*%(C%*%x.bar)
(T2.crit <- (nrow(S)-1)*(n-1)/(n-nrow(S)+1)*qf(1-alpha,nrow(S)-1,n-nrow(S)+1))

C <- matrix(c(1,-1,0,0,
              0,1,-1,0,
              0,0,1,-1,
              1,1,-1,-1), byrow = TRUE, ncol=4)

t(C[4,])%*%x.bar-sqrt((n-1)*(nrow(S)-1)*qf(1-alpha,nrow(S)-1,n-nrow(S)+1)/(n+1-nrow(S)))*sqrt(t(C[4,])%*%S%*%C[4,]/n)
t(C[4,])%*%x.bar+sqrt((n-1)*(nrow(S)-1)*qf(1-alpha,nrow(S)-1,n-nrow(S)+1)/(n+1-nrow(S)))*sqrt(t(C[4,])%*%S%*%C[4,]/n)

### 6.24
rm(list=ls())
alpha=0.05
X <- read.table("T6-13.dat")
names(X) <- c("MaxBreath","BasHeight","BasLength","NasHeight","Time")
(a <- table(X[,5]))
n1 <- 30
n2 <- 30
n3 <- 30
n <- n1 + n2 + n3
p <- 4
g <- 3

xbar1 <- matrix(c(mean(X[1:30,1]),mean(X[1:30,2]),mean(X[1:30,3]),mean(X[1:30,4])), byrow = TRUE, ncol=1)
xbar2 <- matrix(c(mean(X[31:60,1]),mean(X[31:60,2]),mean(X[31:60,3]),mean(X[31:60,4])), byrow = TRUE, ncol=1)
xbar3 <- matrix(c(mean(X[61:90,1]),mean(X[61:90,2]),mean(X[61:90,3]),mean(X[61:90,4])), byrow = TRUE, ncol=1)

t1=X[1:30,]
t2=X[31:60,]
t3=X[61:90,]

(S1 <- matrix(rep(NA,4*4),ncol=4))
for(i in 1:4){
  for(j in 1:4){
    S1[i,j]=sum((t1[,i]-xbar1[i,])*(t1[,j]-xbar1[j,]))/(30-1)
  }
}
S1

(S2 <- matrix(rep(NA,4*4),ncol=4))
for(i in 1:4){
  for(j in 1:4){
    S2[i,j]=sum((t2[,i]-xbar2[i,])*(t2[,j]-xbar2[j,]))/(30-1)
  }
}
S2

(S3 <- matrix(rep(NA,4*4),ncol=4))
for(i in 1:4){
  for(j in 1:4){
    S3[i,j]=sum((t3[,i]-xbar3[i,])*(t3[,j]-xbar3[j,]))/(30-1)
  }
}
S3

(W <- (n1-1)*S1 + (n2-1)*S2 + (n3-1)*S3)
(xbar <- (1/(n1+n2+n3))*(n1*xbar1 + n2*xbar2 + n3*xbar3))
(B <- n1 * (xbar1-xbar)%*%t(xbar1-xbar) + n2 * (xbar2-xbar)%*%t(xbar2-xbar) + n3 * (xbar3-xbar)%*%t(xbar3-xbar))

(Lambda <- det(W)/det(B+W))

(test.statistic <- (n-p-2)/p * (1-sqrt(Lambda))/sqrt(Lambda))
(critical.value <- qf(1-alpha,2*p,2*(n-p-2)))

tau1 <- (xbar1-xbar)
tau2 <- (xbar2-xbar)
tau3 <- (xbar3-xbar)

# Max Breath
(tau1[1,1] - tau2[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[1,1]/(n-g))
(tau1[1,1] - tau2[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[1,1]/(n-g))

(tau1[1,1] - tau3[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[1,1]/(n-g))
(tau1[1,1] - tau3[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[1,1]/(n-g))

(tau2[1,1] - tau3[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[1,1]/(n-g))
(tau2[1,1] - tau3[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[1,1]/(n-g))

#BasHeight
(tau1[2,1] - tau2[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[2,2]/(n-g))
(tau1[2,1] - tau2[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[2,2]/(n-g))

(tau1[2,1] - tau3[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[2,2]/(n-g))
(tau1[2,1] - tau3[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[2,2]/(n-g))

(tau2[2,1] - tau3[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[2,2]/(n-g))
(tau2[2,1] - tau3[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[2,2]/(n-g))

#BasLength
(tau1[3,1] - tau2[3,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[3,3]/(n-g))
(tau1[3,1] - tau2[3,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[3,3]/(n-g))

(tau1[3,1] - tau3[3,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[3,3]/(n-g))
(tau1[3,1] - tau3[3,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[3,3]/(n-g))

(tau2[3,1] - tau3[3,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[3,3]/(n-g))
(tau2[3,1] - tau3[3,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[3,3]/(n-g))

#NasHeight
(tau1[4,1] - tau2[4,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[4,4]/(n-g))
(tau1[4,1] - tau2[4,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[4,4]/(n-g))

(tau1[4,1] - tau3[4,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[4,4]/(n-g))
(tau1[4,1] - tau3[4,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[4,4]/(n-g))

(tau2[4,1] - tau3[4,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[4,4]/(n-g))
(tau2[4,1] - tau3[4,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[4,4]/(n-g))

png("624pair.png")
pairs(X[,1:4])
dev.off()

#### 6.31
rm(list=ls())
alpha=.05
X <- read.table("T6-17.dat")
names(X) <- c("Location", "Variety", "x1","x2","x3")

dat2 <- X
dat2[,1] <- as.factor(dat2[,1])
dat2[,2] <- as.factor(dat2[,2])

fit <- manova(cbind(x1,x2,x3) ~ Location * Variety, data=dat2)
(fit.summary <- summary(fit)) 
fit.summary$SS
anova(fit)


xo <- sort(fit$residuals[,1])
n <- length(xo)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xo,pch=19)
png("631r1.png")
plot(q,xo,pch=19)
dev.off()

xo <- sort(fit$residuals[,2])
n <- length(xo)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xo,pch=19)
png("631r2.png")
plot(q,xo,pch=19)
dev.off()

xo <- sort(fit$residuals[,3])
n <- length(xo)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xo,pch=19)
png("631r3.png")
plot(q,xo,pch=19)
dev.off()


fit2 <-aov(x1 ~ Location+Variety+Location*Variety, data=dat2)
(fit2.summary <- summary(fit2)) 
anova(fit2)

fit3 <-aov(x2 ~ Location+Variety+Location*Variety, data=dat2)
(fit3.summary <- summary(fit3)) 
anova(fit3)

fit4 <-aov(x3 ~ Location+Variety+Location*Variety, data=dat2)
(fit4.summary <- summary(fit4)) 
anova(fit4)

loc2=X[X[,1]==2,]
t1=loc2[1:2,]
t2=loc2[3:4,]
t3=loc2[5:6,]
n1 <- 2
n2 <- 2
n3 <- 2
n <- n1 + n2 + n3
p <- 3
g <- 1

xbar5 <- matrix(c(mean(loc2[1:2,3]),mean(loc2[1:2,4]),mean(loc2[1:2,5])), byrow = TRUE, ncol=1)
xbar6 <- matrix(c(mean(loc2[3:4,3]),mean(loc2[3:4,4]),mean(loc2[3:4,5])), byrow = TRUE, ncol=1)
xbar8 <- matrix(c(mean(loc2[5:6,3]),mean(loc2[5:6,4]),mean(loc2[5:6,5])), byrow = TRUE, ncol=1)

(S1 <- matrix(rep(NA,3*3),ncol=3))
for(i in 1:3){
  for(j in 1:3){
    S1[i,j]=sum((t1[,i]-xbar5[i,])*(t1[,j]-xbar5[j,]))/(2-1)
  }
}
S1

(S2 <- matrix(rep(NA,3*3),ncol=3))
for(i in 1:3){
  for(j in 1:3){
    S2[i,j]=sum((t2[,i]-xbar6[i,])*(t2[,j]-xbar6[j,]))/(2-1)
  }
}
S2

(S3 <- matrix(rep(NA,3*3),ncol=3))
for(i in 1:3){
  for(j in 1:3){
    S3[i,j]=sum((t3[,i]-xbar8[i,])*(t3[,j]-xbar8[j,]))/(2-1)
  }
}
S3

(W <- (n1-1)*S1 + (n2-1)*S2 + (n3-1)*S3)
(xbar <- (1/(n1+n2+n3))*(n1*xbar5 + n2*xbar6 + n3*xbar8))
(B <- n1 * (xbar5-xbar)%*%t(xbar5-xbar) + n2 * (xbar6-xbar)%*%t(xbar6-xbar) + n3 * (xbar8-xbar)%*%t(xbar8-xbar))

(Lambda <- det(W)/det(B+W))

(test.statistic <- (n-p-2)/p * (1-sqrt(Lambda))/sqrt(Lambda))
(critical.value <- qf(1-alpha,2*p,2*(n-p-2)))

tau1 <- (xbar5-xbar)
tau2 <- (xbar6-xbar)
tau3 <- (xbar8-xbar)

(tau1[1,1] - tau2[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[1,1]/(n-g))
(tau1[1,1] - tau2[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[1,1]/(n-g))

(tau1[1,1] - tau3[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[1,1]/(n-g))
(tau1[1,1] - tau3[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[1,1]/(n-g))

(tau2[1,1] - tau3[1,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[1,1]/(n-g))
(tau2[1,1] - tau3[1,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[1,1]/(n-g))



(tau1[2,1] - tau2[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[2,2]/(n-g))
(tau1[2,1] - tau2[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[2,2]/(n-g))

(tau1[2,1] - tau3[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[2,2]/(n-g))
(tau1[2,1] - tau3[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[2,2]/(n-g))

(tau2[2,1] - tau3[2,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[2,2]/(n-g))
(tau2[2,1] - tau3[2,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[2,2]/(n-g))



(tau1[3,1] - tau2[3,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[3,3]/(n-g))
(tau1[3,1] - tau2[3,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[3,3]/(n-g))

(tau1[3,1] - tau3[3,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[3,3]/(n-g))
(tau1[3,1] - tau3[3,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[3,3]/(n-g))

(tau2[3,1] - tau3[3,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[3,3]/(n-g))
(tau2[3,1] - tau3[3,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n2)*W[3,3]/(n-g))


### 6.32
rm(list=ls())
alpha=.05
X <- matrix(c(10.35, 25.93, 'SS', 'p',
              13.41, 38.63, 'JL', 'p',
              7.78, 25.15, 'LP', 'p',
              10.4, 24.25, 'SS', 'm',
              17.78, 41.45, 'JL', 'm',
              10.4, 29.20, 'LP', 'm'), byrow=TRUE, ncol=4)
X <- data.frame(X)
names(X) <- c("560CM", "720CM", "Species","Nutrient")

dat2 <- X
dat2[,3] <- as.factor(dat2[,3])
dat2[,4] <- as.factor(dat2[,4])


fit <- manova(cbind('560CM', '720CM') ~ Species * Nutrient, data=dat2)
(fit.summary <- summary(fit)) 
fit.summary$SS
anova(fit)