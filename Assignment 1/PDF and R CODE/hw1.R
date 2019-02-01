setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Assignment 1")
#####1.2######
x1=c(1,2,3,3,4,5,6,8,9,11)
x2=c(18.95,19.00,17.95,15.54,14.00,12.95,8.94,7.49,6.00,3.99)
library(ggplot2)
library(ggExtra)
dat=data.frame(x1,x2)
names(dat)<- c("Age","Cost")
png("12scatter.png")
ggMarginal(ggplot(dat, aes(Age, Cost))+geom_point() + theme_gray() + ggtitle("Used Car Cost vs Age") + xlab("Age") + ylab("Cost (in Thousands)"), type="histogram", fill="steelblue", col="darkblue")
dev.off()
(mean(dat$Age))
mean(dat$Cost)
s11=sum((dat$Age-mean(dat$Age))**2)/10
s22=sum((dat$Cost-mean(dat$Cost))**2)/10
s12=sum((dat$Age-mean(dat$Age))*(dat$Cost-mean(dat$Cost)))/10
r12=s12/(sqrt(s11)*sqrt(s22))

#### 1.4 ####
dat <- read.table("P1-4.DAT", header = FALSE)
names(dat)<-c("x1","x2","x3")
png("14scatter.png")
ggMarginal(ggplot(dat, aes(x1, x2))+geom_point() + theme_gray() + ggtitle("x_1 vs x_2, Sales vs Profits") + xlab("Sales (billions)") + ylab("Profits (billions)"), type="histogram", fill="steelblue", col="darkblue")
dev.off()
mean(dat$x1)
mean(dat$x2)
s11=sum((dat$x1-mean(dat$x1))**2)/length(dat$x1)
s22=sum((dat$x2-mean(dat$x2))**2)/length(dat$x2)
s12=sum((dat$x1-mean(dat$x1))*(dat$x2-mean(dat$x2)))/10
r12=s12/(sqrt(s11)*sqrt(s22))

### 1.5 ####
png("1523scatter.png")
ggMarginal(ggplot(dat, aes(x2, x3))+geom_point() + theme_gray() + ggtitle("x_2 vs x_3, Profits vs Assets") + xlab("Profits (billions)") + ylab("Assets (billions)"), type="histogram", fill="steelblue", col="darkblue")
dev.off()

png("1513scatter.png")
ggMarginal(ggplot(dat, aes(x1, x3))+geom_point() + theme_gray() + ggtitle("x_1 vs x_3, Sales vs Assets") + xlab("Sales (billions)") + ylab("Assets (billions)"), type="histogram", fill="steelblue", col="darkblue")
dev.off()

mean(dat$x1)
mean(dat$x2)
mean(dat$x3)

s33=sum((dat$x3-mean(dat$x3))**2)/length(dat$x3)
s13=sum((dat$x1-mean(dat$x1))*(dat$x3-mean(dat$x3)))/10
s23=sum((dat$x2-mean(dat$x2))*(dat$x3-mean(dat$x3)))/10
r13=s13/(sqrt(s11)*sqrt(s33))
r23=s23/(sqrt(s22)*sqrt(s33))

### 1.17 ####
dat <- read.table("T1-9.dat", header = FALSE)
names(dat) <- c("Country","100m","200m","400m","800m(min)","1500m(min)","3000m(min)", "Marathon(min)")
dat=dat[,2:8]
summary(dat)

SN <- matrix(rep(NA,7*7),ncol=7)
for(i in 1:7){
  for(j in 1:7){
    SN[i,j] = (1/54)*sum((dat[,i]-mean(dat[,i]))*(dat[,j]-mean(dat[,j])))
  }
}
print(SN)

RN <- matrix(rep(NA,7*7),ncol=7)
for(i in 1:7){
  for(j in 1:7){
    RN[i,j] = SN[i,j]/(sqrt(SN[i,i])*sqrt(SN[j,j]))
  }
}
print(RN)

### 1.18 ###
speed=dat
speed$`800m(min)`=speed$`800m(min)`*60
speed$`1500m(min)`=speed$`1500m(min)`*60
speed$`3000m(min)`=speed$`3000m(min)`*60
speed$`Marathon(min)`=speed$`Marathon(min)`*60

speed$`100m`=100/speed$`100m`
speed$`200m`=200/speed$`200m`
speed$`400m`=400/speed$`400m`
speed$`800m(min)`=800/speed$`800m(min)`
speed$`1500m(min)`=1500/speed$`1500m(min)`
speed$`3000m(min)`=3000/speed$`3000m(min)`
speed$`Marathon(min)`=42195/speed$`Marathon(min)`

summary(speed)
SSpeed <- matrix(rep(NA,7*7),ncol=7)
for(i in 1:7){
  for(j in 1:7){
    SSpeed[i,j] = (1/54)*sum((speed[,i]-mean(speed[,i]))*(speed[,j]-mean(speed[,j])))
  }
}
print(SSpeed)

RSpeed <- matrix(rep(NA,7*7),ncol=7)
for(i in 1:7){
  for(j in 1:7){
    RSpeed[i,j] = SSpeed[i,j]/(sqrt(SSpeed[i,i])*sqrt(SSpeed[j,j]))
  }
}
print(RSpeed)
