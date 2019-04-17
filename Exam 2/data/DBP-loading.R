rm(list=ls())
dat1 <- read.csv("TTHM.csv",header=TRUE)
dat2 <- read.csv("HAA5.csv",header=TRUE)
dat <- merge(dat1,dat2,by=names(dat1)[1],all=FALSE)
dat <- dat[,c(1:4,7)]
names(dat) <- c("Date","Year","Quarter","TTHM","HAA5")
head(dat)
