#setwd("~/Desktop/PSM/Spring 2019/Multivariate-Stats/Project")
library(plotrix)
p1=data.frame(x=double(),y=double())
names(p1)=c("X", "Y")
while(dim(p1)[1]<200){
  x=runif(1, min=-1, max=1)
  y=runif(1, min=-1, max=1)
  if((x**2+y**2)<=1){
    p1[nrow(p1)+1,]=c(x,y)
  }
}
while(dim(p1)[1]<400){
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


####### 12.4
library(sem)
rm(list=ls())
X <- readMoments("T12-3.DAT",diag=TRUE)
colnames(X) <- c("E","N", "Da", "Du", "G","Fr", "Sp", "I", "P", "H", "Fi")
rownames(X) <- c("E","N", "Da", "Du", "G","Fr", "Sp", "I", "P", "H", "Fi")
'''
for(i in 1:dim(X)[1]){
  for(j in i:dim(X)[2]){
    X[i,j]=X[j,i]
  }
}'''

D <- 10-X
C <- hclust(as.dist(D), method = "average")
plot(C)
png("LanguageDendo.png")
plot(C)
dev.off()

### 12.5
rm(list = ls())
D <- matrix(c(0, 9, 3, 6,11,
              9,0, 7 , 5, 10,
              3,7,0, 9,2,
              6,5,9,0,8,
              11,10,2,8,0), byrow=TRUE, ncol=5)
C <- hclust(as.dist(D), method = "complete")
plot(C)
min(D[D>0])
Dnew=D
#Corresponds to (35) cluster. Remove relevant rows, and add new row,
Dnew <- Dnew[-c(3,5),-c(3,5)]
d351 <- max(D[3,1],D[5,1])
d352 <- max(D[3,2],D[5,2])
d354 <- max(D[3,4],D[5,4])
Dnew <- cbind(c(d351,d352,d354),Dnew)
Dnew <- rbind(c(0,d351,d352,d354),Dnew)
D <- Dnew
min(D[D>0])
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
min(D[D>0])
#Corresponds to the (12345) cluster, we are done.
### 12.7
rm(list=ls())
X <- readMoments("T12-5.DAT",diag=TRUE)

C <- hclust(as.dist(X), method = "complete")
plot(C)
