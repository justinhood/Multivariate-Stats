setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Assignment 5")

### 1 ###
rm(list=ls())
set.seed(1) # Lock the random
EPA <- read.csv("FuelEconomy.csv", header = TRUE)
economy <- EPA$Mileage # No scaling Needed

(x.bar <- mean(economy))
(s <- sd(economy))
n <- length(economy)

#Plot the data
png("1a.png")
plot(EPA, main="MPG by Car")
abline(h=c(x.bar,x.bar-s,x.bar+s), col=c("red", "blue", "blue"), lty=c(1,2,2), lwd=c(1,2,2))
dev.off()


# T-interval
alpha <- 0.05 #1-0.95=.05
#         mean +-          t*              Scaler
(lowert = x.bar - qt(1-alpha/2,df=(n-1))*s/sqrt(n))
(uppert = x.bar + qt(1-alpha/2,df=(n-1))*s/sqrt(n))


#Bootstrap
#Run 1000 straps
n.boot <- 1000
x.bars <- c()
#First hundred for dot plotting ease of parsing
for (i in 1:100){
  x.bars[i] <- mean(sample(economy,size=n,replace=TRUE))
}

library(ggplot2)
xbardf <- data.frame(bootmeans=x.bars)
ggplot(xbardf,aes(x=bootmeans)) + geom_dotplot()

#Save This plot
png("econdot.png")
ggplot(xbardf,aes(x=bootmeans)) + geom_dotplot()
dev.off()

#Finish the sampling
for (i in 101:n.boot){
  x.bars[i] <- mean(sample(economy,size=n,replace=TRUE))
}

#Easy Histogram Command
hist(x.bars)

#Save PNG
png("econhist.png")
hist(x.bars)
dev.off()

#CI construction by quantiles
quantile(x.bars,0.025);quantile(x.bars,0.975) # 95% CI
quantile(x.bars,0.05);quantile(x.bars,0.95) # 90% CI



### 2 ###############################################################
rm(list=ls())
dat <- read.csv("EG12-15EYES.csv",header = TRUE)
eyes <- data.frame(dat$Group, dat$Score)
names(eyes)<-c("EyeColor", "Score")

#Plots
library(ggplot2)
qplot(eyes$EyeColor,eyes$Score)
#PNG
png("2qplot.png")
qplot(eyes$EyeColor,eyes$Score)
dev.off()

plot(Score~EyeColor,data=eyes)
#PNG
png("2stem.png")
plot(Score~EyeColor,data=eyes)
dev.off()


fit <- aov(Score~EyeColor,data=eyes)
anova(fit)

#Tukey Intervals
TukeyHSD(fit, ordered = FALSE, conf.level = 0.95)
plot(TukeyHSD(fit, ordered = FALSE, conf.level = 0.95))
#png
png("tukey.png")
plot(TukeyHSD(fit, ordered = FALSE, conf.level = 0.95))
dev.off()


## 3 ###################################################
rm(list=ls())
soft <- read.table("PremiumDistribution.txt")
names(soft) <- c("LapseTime", "Agent", "Transaction")

#Plots
library(ggplot2)
qplot(soft$Agent, soft$LapseTime)
#PNG
png("3qplot.png")
qplot(soft$Agent, soft$LapseTime)
dev.off()

boxplot(LapseTime~Agent,data=soft)
#PNG
png("3stem.png")
boxplot(LapseTime~Agent,data=soft)
dev.off()


soft$Agent <- factor(soft$Agent)
fit <- aov(LapseTime~Agent,data=soft)
anova(fit)
plot(LapseTime~Agent,data=soft)

#png
png("3box.png")
plot(LapseTime~Agent,data=soft)
dev.off()

#CI's
alpha <- 0.05
TukeyHSD(fit, ordered = FALSE, conf.level = 1-alpha)
plot(TukeyHSD(fit, ordered = FALSE, conf.level = 1-alpha))
#png
png("3tukey.png")
plot(TukeyHSD(fit, ordered = FALSE, conf.level = 1-alpha))
dev.off()


m <- choose(5,2) # Common approach is evenly dividing the alpha level among each confidence interval to be made
t.test(soft$LapseTime[which(soft$Agent==1)],
       soft$LapseTime[which(soft$Agent==2)],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(soft$LapseTime[which(soft$Agent==1)],
       soft$LapseTime[which(soft$Agent==3)],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(soft$LapseTime[which(soft$Agent==1)],
       soft$LapseTime[which(soft$Agent==4)],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(soft$LapseTime[which(soft$Agent==1)],
       soft$LapseTime[which(soft$Agent==5)],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(soft$LapseTime[which(soft$Agent==2)],
       soft$LapseTime[which(soft$Agent==3)],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(soft$LapseTime[which(soft$Agent==2)],
       soft$LapseTime[which(soft$Agent==4)],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(soft$LapseTime[which(soft$Agent==2)],
       soft$LapseTime[which(soft$Agent==5)],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(soft$LapseTime[which(soft$Agent==3)],
       soft$LapseTime[which(soft$Agent==4)],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(soft$LapseTime[which(soft$Agent==3)],
       soft$LapseTime[which(soft$Agent==5)],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(soft$LapseTime[which(soft$Agent==4)],
       soft$LapseTime[which(soft$Agent==5)],
       conf.level=1-alpha/m)$conf.int[1:2]

## 4 #################################################
rm(list=ls())
set.seed(1) # Lock the random
soft <- read.table("PremiumDistribution.txt")
names(soft) <- c("LapseTime", "Agent", "Transaction")
n<- length(soft[,1])
soft$Agent <- factor(soft$Agent)
#Bootstrap
#Run 1000 straps
n.boot <- 1000
R2s <- c()
#First hundred for dot plotting ease of parsing
for (i in 1:100){
  #Get Row numbers w/repeat
  X=sample(nrow(soft),size=n,replace=TRUE)
  #Null data frame
  dat = data.frame()
  #Extract by row into data frame
  for (j in 1:length(X)){
    dat[j,1]=soft[X[j],1]
    dat[j,2]=soft[X[j],2]
    dat[j,3]=soft[X[j],3]
  }
  names(dat) <- c("LapseTime", "Agent", "Transaction")
  #Run ANOVA on sample frame
  fit <- aov(LapseTime~Agent,data=dat)
  #compute R^2
  R2s[i] <- 1-anova(fit)$`Sum Sq`[2]/sum(anova(fit)$`Sum Sq`)
}
#Plot R2's
library(ggplot2)
R2 <- data.frame(bootmeans=R2s)
ggplot(R2,aes(x=bootmeans)) + geom_dotplot()

#Save This plot
png("R2dot.png")
ggplot(R2,aes(x=bootmeans)) + geom_dotplot()
dev.off()

#Finish the sampling
for (i in 101:n.boot){
  #Get Row numbers w/repeat
  X=sample(nrow(soft),size=n,replace=TRUE)
  #Null data frame
  dat = data.frame()
  #Extract by row into data frame
  for (j in 1:length(X)){
    dat[j,1]=soft[X[j],1]
    dat[j,2]=soft[X[j],2]
    dat[j,3]=soft[X[j],3]
  }
  names(dat) <- c("LapseTime", "Agent", "Transaction")
  #Run ANOVA on sample frame
  fit <- aov(LapseTime~Agent,data=dat)
  #compute R^2
  R2s[i] <- 1-anova(fit)$`Sum Sq`[2]/sum(anova(fit)$`Sum Sq`)
}

#Easy Histogram Command
hist(R2s)

#Save PNG
png("R2hist.png")
hist(R2s)
dev.off()

#CI construction by quantiles
quantile(R2s,0.025);quantile(R2s,0.975) # 95% CI
quantile(R2s,0.05);quantile(R2s,0.95) # 90% CI

## 5 ##############################################################
rm(list=ls())
head(OrchardSprays)
str(OrchardSprays)
boxplot(decrease~treatment, data=OrchardSprays, xlab="treatment", ylab="decrease mean (effect size)")
png("5box.png")
boxplot(decrease~treatment, data=OrchardSprays, xlab="treatment", ylab="decrease mean (effect size)")
dev.off()
treatbar <- c()
treatbar[1] <- mean(OrchardSprays$decrease[OrchardSprays$treatment == "A"])
treatbar[2] <- mean(OrchardSprays$decrease[OrchardSprays$treatment == "B"])
treatbar[3] <- mean(OrchardSprays$decrease[OrchardSprays$treatment == "C"])
treatbar[4] <- mean(OrchardSprays$decrease[OrchardSprays$treatment == "D"])
treatbar[5] <- mean(OrchardSprays$decrease[OrchardSprays$treatment == "E"])
treatbar[6] <- mean(OrchardSprays$decrease[OrchardSprays$treatment == "F"])
treatbar[7] <- mean(OrchardSprays$decrease[OrchardSprays$treatment == "G"])
treatbar[8] <- mean(OrchardSprays$decrease[OrchardSprays$treatment == "H"])

treatdev <- c()
treatdev[1] <- sd(OrchardSprays$decrease[OrchardSprays$treatment == "A"])
treatdev[2] <- sd(OrchardSprays$decrease[OrchardSprays$treatment == "B"])
treatdev[3] <- sd(OrchardSprays$decrease[OrchardSprays$treatment == "C"])
treatdev[4] <- sd(OrchardSprays$decrease[OrchardSprays$treatment == "D"])
treatdev[5] <- sd(OrchardSprays$decrease[OrchardSprays$treatment == "E"])
treatdev[6] <- sd(OrchardSprays$decrease[OrchardSprays$treatment == "F"])
treatdev[7] <- sd(OrchardSprays$decrease[OrchardSprays$treatment == "G"])
treatdev[8] <- sd(OrchardSprays$decrease[OrchardSprays$treatment == "H"])


#### Kruskal-Wallis Test####
rm(list=ls())

dat <- data.frame(OrchardSprays,ranks=rank(OrchardSprays$decrease))
n <- nrow(dat)
ni <- 8 # for each

Ri <- by(dat$ranks,
         INDICES=OrchardSprays$treatment,
         FUN=sum)
(H <-  - 3*(n+1) + 12/(n*(n+1)) * sum(Ri^2/ni))
(Hcrit <- qchisq(1-0.05,4-1))

# Using built-in functions
kruskal.test(x=OrchardSprays$decrease,g=OrchardSprays$treatment)

alpha <- .05
m <- choose(8,2) # Common approach is evenly dividing the alpha level among each confidence interval to be made
t.test(OrchardSprays$decrease[which(OrchardSprays$treatment=="A")],
       OrchardSprays$decrease[which(OrchardSprays$treatment=="B")],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(OrchardSprays$decrease[which(OrchardSprays$treatment=="A")],
       OrchardSprays$decrease[which(OrchardSprays$treatment=="C")],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(OrchardSprays$decrease[which(OrchardSprays$treatment=="A")],
       OrchardSprays$decrease[which(OrchardSprays$treatment=="D")],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(OrchardSprays$decrease[which(OrchardSprays$treatment=="B")],
       OrchardSprays$decrease[which(OrchardSprays$treatment=="C")],
       conf.level=1-alpha/m)$conf.int[1:2]
t.test(OrchardSprays$decrease[which(OrchardSprays$treatment=="B")],
       OrchardSprays$decrease[which(OrchardSprays$treatment=="D")],
       conf.level=1-alpha/m)$conf.int[1:2]
