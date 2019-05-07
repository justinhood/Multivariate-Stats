#setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Project")


#Make my n=2 dimenstional data/hypershpere image
library(plotrix)
p1=data.frame(x=double(),y=double()) #Points
names(p1)=c("X", "Y")
while(dim(p1)[1]<200){ #Get 200 points for the first bubble by making circle
  x=runif(1, min=-1, max=1)
  y=runif(1, min=-1, max=1)
  if((x**2+y**2)<=1){
    p1[nrow(p1)+1,]=c(x,y)
  }
}

while(dim(p1)[1]<400){ #200 more but with a shift to a new center now
  x=runif(1, min=-1, max=1)
  y=runif(1, min=-1, max=1)
  if((x**2+y**2)<=1){
    p1[nrow(p1)+1,]=c(x+3,y+3)
  }
}
plot(p1)
png("circles.png")
plot(p1)
dev.off()


####### Example 12.4
library(sem)
rm(list=ls())
X <- readMoments("T12-3.DAT",diag=TRUE) #Language Linakage data
colnames(X) <- c("E","N", "Da", "Du", "G","Fr", "Sp", "I", "P", "H", "Fi")
rownames(X) <- c("E","N", "Da", "Du", "G","Fr", "Sp", "I", "P", "H", "Fi")

D <- 10-X #Convert to distance
C <- hclust(as.dist(D), method = "average") #Do an average distance cluster analysis
plot(C) #Dendro Plot
png("LanguageDendo.png")
plot(C)
dev.off()

### Example 12.5
rm(list = ls())
D <- matrix(c(0, 9, 3, 6,11,
              9,0, 7 , 5, 10,
              3,7,0, 9,2,
              6,5,9,0,8,
              11,10,2,8,0), byrow=TRUE, ncol=5) #Complete linkage example
C <- hclust(as.dist(D), method = "complete") #Desired Result
plot(C)
#This plot shows that the first cluster will be (35), then (24), then (124),
#and finally (12345)
min(D[D>0]) #Find the smallest value in distance matrix
Dnew=D
#Corresponds to (35) cluster. Remove relevant rows, and add new row,
Dnew <- Dnew[-c(3,5),-c(3,5)]
d351 <- max(D[3,1],D[5,1])
d352 <- max(D[3,2],D[5,2])
d354 <- max(D[3,4],D[5,4])
Dnew <- cbind(c(d351,d352,d354),Dnew)
Dnew <- rbind(c(0,d351,d352,d354),Dnew)
D <- Dnew
min(D[D>0])#Find the smallest value in distance matrix
#Corresponds to (24) cluster. Remove relevant rows, and add new row,
Dnew <- Dnew[-c(3,4),-c(3,4)]
d2435 <- max(D[1,3],D[1,4])
d241 <- max(D[3,2],D[4,2])
Dnew <- cbind(Dnew[,1],c(0,d241),Dnew[,2])
Dnew <- rbind(Dnew[1,],c(d2435,0,d241),Dnew[2,])
min(D[D>0])
#Corresponds to (124) cluster. Remove relevant rows, and add new row,
Dnew <- Dnew[-c(2,3),-c(2,3)]
d24135 <- max(D[1,3],D[1,2])
Dnew <- cbind(c(0,d24135),c(d24135,0))
D<- Dnew
min(D[D>0])#Find the smallest value in distance matrix
#Corresponds to the (12345) cluster, we are done.
#We see here that our results matched perfectly the complete linkage method in R


### Example 12.7
rm(list=ls())
library(sem)
X <- readMoments("T12-5.DAT",diag=TRUE) #Correlations in public utility data
C <- hclust(as.dist(X), method = "complete")
plot(C)



#######################################################################
#######################################################################
############# Some Worked Problems ####################################
#######################################################################
#######################################################################
#12 .1
'''
Consider the binary values:
{1 if from South, 0 else}
{1 if elected first term, 0 else}
{1 if Democrat, 0 else}
{1 if Prior Experience, 0 else}
{1 if former Vice, 0 else}
Then, consider the following,
'''
rm(list = ls())
X <- matrix(c(0,1,0,0,0,
              1,1,1,0,0,
              0,0,0,1,1,
              0,1,0,1,1,
              1,0,1,1,1,
              0,1,1,1,0), byrow=TRUE, ncol=5)
#Reagan Carter
a=1
b=3
sim1 <- matrix(c(0,0,
                 0,0), byrow=TRUE, ncol=2)
for(i in 1:5){
  if(X[a,i]==1 && X[b,i]==1){
    sim1[1,1]<- sim1[1,1]+1
  }
  if(X[a,i]==1 && X[b,i]==0){
    sim1[1,2]<- sim1[1,2]+1
  }
  if(X[a,i]==0 && X[b,i]==1){
    sim1[2,1]<- sim1[2,1]+1
  }
  if(X[a,i]==0 && X[b,i]==0){
    sim1[2,2]<- sim1[2,2]+1
  }
}
(simCoef1 <- (sim1[1,1]+sim1[2,2])/sum(sim1))

## 12.7 Stock Data
rm(list=ls())
X <- matrix(c(1,.63,.51,.12,.16,
              .63,1,.57,.32,.21,
              .51,.57,1,.18,.15,
              .12,.32,.18,1,.68,
              .16,.21,.15,.68,1), byrow = TRUE, ncol=5)
row.names(X) <- c("JP Morgan", "Citibank","Wells Fargo", "Royal DutchShell","ExxonMobil")
C <- hclust(as.dist(X), method = "complete")
plot(C)
#We see based on this analysis that JP Morgan and Royal DutchShell are closely related, followed by Wells fargo and ExxonMobil
C <- hclust(as.dist(X), method = "average")
plot(C)
#We see that there is a similar initial breakdown in this dendrogram, but uwing the average method actually causes the two pairs to merge before adding Citibank now

## 12.15
#This one will be k-means clustering
rm(list=ls())
X <- read.table("T11-9.dat")
X<-X[-c(1,2)]
k1 <- kmeans(X[,1:5],2)
print(k1)
k2 <- kmeans(X[,1:5],3)
print(k2)
k3 <- kmeans(X[,1:5],4)
print(k3)
hist(k1$cluster)
hist(k2$cluster)
hist(k3$cluster)
#We see that the data clusters into roughly even groups, and that each successive addtion of a group is likely just a natural partition of an existing group

# 12.26 Maloi Farm
rm(list=ls())
X <- read.table("T8-7.DAT", header=FALSE)
X <- X[-c(25,34,69,72),]#,-c(1,7:9)]
D <- matrix(0, nrow= 72, ncol=72)
for(i in 1:72){
  for(j in i:72){
    d <- sqrt(sum((X[i,]-X[j,])^2))
    D[i,j] <- d
    
  }
}
D<-t(D)
D <- as.dist(D)
avg <- hclust(D, method="average")
plot(avg)
ward <- hclust(D, method = "ward.D")
plot(ward)
## Based on the Ward method, we see tat there are around 7 base clusters of farms based on Euclidean Distance
library(cluster)
k1 <- kmeans(X,5)
s=silhouette(k1$cl,D)
plot(s)
clusplot(X,k1$cluster)

k2 <- kmeans(X,6)
s2=silhouette(k2$cl,D)
plot(s2)
clusplot(X,k2$cluster)

#We see that the data is well wrapped within the 5 or 6 clusters that we are using.