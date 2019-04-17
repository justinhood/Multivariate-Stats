#install.packages("reshape2")
library(reshape2)
dat.fill <- read.csv("BHO_Medication_Fill_Data__2010-2014.csv",header=TRUE)
tmp <- dcast(dat.fill, OMH.Region + Year + Quarter + Age.Group ~ Description.of.Metric)
dat <- tmp[,c(1:4,6,8,10)]
names(dat) <- c(names(dat)[1:4],"Mood30","Psychotrop30","Antipsych30")
head(dat)