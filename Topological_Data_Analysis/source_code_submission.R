#Author: Lawrence Adu-Gyamfi
# No: 1484610
# Date: 30/05/1990

library(rgl)
library(data.table)
library(TDA)
library(rgl)
library(plot3D)

# set working directory (path to points to be included in file_path below)
file_path <- "##"
setwd(file_path)

# load file
points1 <- fread("points1.csv", data.table=F)
points2 <- fread("points2.csv", data.table=F)
points3 <- fread("points3.csv", data.table=F)

# quick plots
plot(points1, col=rainbow(1000))
plot(points2)
plot(points3)

# #3-D plots
plot3d(points1,xlab="x-axis",ylab="y-axis", zlab="z-axis",main="Points1 3-D", col=rainbow(1000), top=F)
plot3d(points2,xlab="x-axis",ylab="y-axis", zlab="z-axis",main="Points2 3-D", top=F)
plot3d(points3,xlab="x-axis",ylab="y-axis", zlab="z-axis",main="Points3 3-D")

# #3-D plots
par(mfrow3d(1,3))
scatter3D(points1$x1, points1$x2,points1$x3,xlab="x-axis",ylab="y-axis", 
      zlab="z-axis", main="Points1 3-D", theta=45)
scatter3D(x=points2$x1, y=points2$x2, z=points2$x3,main="Points2 3-D", 
          xlab="x-axis",ylab="y-axis", zlab="z-axis", theta=45, d=5)
scatter3D(x=points3$x1, y=points3$x2, z=points3$x3,main="Points3 3-D", theta=45, d=2)


# set limits of axis
Xlim <- c(-3,3)
Ylim <- c(-2,2)
Zlim <- c(-2,2)
by=0.10

Xseq <- seq(from=Xlim[1], to=Xlim[2], by=by)
Yseq <- seq(from=Ylim[1], to=Ylim[2], by=by)
Zseq <- seq(from=Zlim[1], to=Zlim[2], by=by)

# create grid
Grid <- expand.grid(Xseq, Yseq, Zseq)

#calculate 95% confidence band using bootstrap
band_1 <- bootstrapBand(X = points1, FUN=kde, Grid =Grid, 
                        B = 100,parallel = T, alpha = 0.05, h = 0.3)
band_2 <- bootstrapBand(X = points2, FUN=kde, Grid =Grid, 
                        B = 100,parallel = T, alpha = 0.05, h = 0.3)
band_3 <- bootstrapBand(X = points3, FUN=kde, Grid =Grid, 
                        B = 100,parallel = T, alpha = 0.05, h = 0.3)


# Calculate the persistence homology and make the plots

#using KDE

par(mfrow=c(1,3))
DiagKDE_1 <- gridDiag(X=points1[,1:3], FUN=kde,h=0.3,sublevel=F, maxdimension = 2, lim=cbind(Xlim,Ylim,Zlim), by=by, library="Dionysus",printProgress=F)
plot(x=DiagKDE_1[["diagram"]], main="KDE Diagram point1", barcode = F, rotated = F,
     band=2*band_1[["width"]])

#point 2
DiagKDE_2 <- gridDiag(X=points2[,1:3], FUN=kde,h=0.3,sublevel=F, maxdimension = 2, lim=cbind(Xlim,Ylim,Zlim), by=by, library="Dionysus",printProgress=F)
plot(x=DiagKDE_2[["diagram"]], main="KDE Diagram point 2", barcode = F, rotated = F,
     band =2*band_2[["width"]])

#points 3
DiagKDE_3 <- gridDiag(X=points3[,1:3], FUN=kde,h=0.3,sublevel=F, maxdimension = 2, lim=cbind(Xlim,Ylim,Zlim), by=by, library="Dionysus",printProgress=F)
plot(x=DiagKDE_3[["diagram"]], main="KDE Diagram point 3", barcode = F, rotated = F,
     band = 2*band_3[["width"]])

# Using VIETORIS RIPS

par(mfrow=c(1,3))
DiagVR_1 <- ripsDiag(X=points1[,1:3],maxdimension=2, maxscale=1.5,dist='euclidean', library="GUDHI",printProgress=F)
plot(x=DiagVR_1[["diagram"]], main="Vietoris-Rips Diagram Points 1", barcode = F)

DiagVR_2 <- ripsDiag(X=points2[,1:3],maxdimension=2, maxscale=1,dist='euclidean', library="GUDHI",printProgress=F)
plot(x=DiagVR_2[["diagram"]], main="Vietoris-Rips Diagram - Points 2", barcode = F)

DiagVR_3 <- ripsDiag(X=points3[,1:3],maxdimension=2, maxscale=1,dist='euclidean', library="GUDHI",printProgress=F)
plot(x=DiagVR_3[["diagram"]], main="Vietoris-Rips Diagram Points 3", barcode = F)
