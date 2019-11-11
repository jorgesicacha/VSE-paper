library(sp)
library(rgeos)
library(INLA)
library(dplyr)
library(raster)
library(pbapply)
library(reshape)
library(tiff)

## Initial data manipulation ##

datagbif <- readRDS("moosegbif.rds")
newcrs <- CRS("+proj=robin +datum=WGS84 +units=km")
norwaybu <-raster::getData("GADM",country="NOR",level=1)
norwaybu <- spTransform(norwaybu,newcrs)
norway <- SpatialPolygons(list(norwaybu@polygons[[6]])) ## Hedmark polygon
proj4string(norway) <- CRS(proj4string(norwaybu))
datagbif.no <- datagbif
datagbif.no <- datagbif.no[which(datagbif.no$basisOfRecord=="HUMAN_OBSERVATION"),] #Only human observation records
dates <- paste(datagbif.no$day,datagbif.no$month,datagbif.no$year,sep="/")
dates <- as.Date(dates,format = "%d/%m/%Y")
datagbif.no <- datagbif.no[which(dates > "2000-01-01"),]
coordsgbif.no <- datagbif.no[,c("lon","lat")]
coordsgbif.no<- distinct(coordsgbif.no)
coordsgbif.no <- coordsgbif.no[complete.cases(coordsgbif.no),]
sppointscoords.no <- SpatialPoints(coordsgbif.no,proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
sppointscoords.no <- spTransform(sppointscoords.no,newcrs)
sppoints.gbif <- gIntersection(norway,sppointscoords.no)


## Making the models ##

source('spde-book-functions.R') ##Taken from Krainski book
premesh<- inla.sp2segment(norway)
mesh <- inla.mesh.2d(boundary = premesh, loc = sppoints.gbif, offset = c(10, -0.2),
                     max.edge = c(7,50),cutoff=4,min.angle = 30)

nv <- mesh$n
meshpoints <- mesh$loc[,1:2]
par(mar = c(0, 0, 0, 0))
plot(mesh, asp = 1, main = '')
points(sppoints.gbif, col = 2, pch = 19)

#### Defining PC priors ####
spde <- inla.spde2.pcmatern(mesh, 
                            prior.sigma = c(1, 0.05), 
                            prior.range = c(15, 0.05))


#### Constructing the dual mesh #### Partly taken from Krainski book
dmesh <- book.mesh.dual(mesh)

#### Converting the domain polygon into a SpatialPolygons class ####
domain.polys <- norway@polygons
domainSP0 <- SpatialPolygons(domain.polys)

w1 <- pblapply(1:length(dmesh), function(i) {
  gIntersects(dmesh[i, ], domainSP0)})
w1.1 <- do.call(rbind,w1)
table(w1.1)
w1.1.1 <- which(w1.1)
length(w1.1.1)
## Let's focus on the largest polygon
areas.no <- c()
for(i in 1:length(norway@polygons[[1]]@Polygons)){
  areas.no[i]<- norway@polygons[[1]]@Polygons[[i]]@area
  print(i)
}
summary(areas.no)
which.max(areas.no)
sum(areas.no)
domain.polys1 <- norway@polygons[[1]]@Polygons[[which.max(areas.no)]]
domainSP1 <- SpatialPolygons(list(Polygons(list(domain.polys1),1)))

w2 <- pblapply(1:length(dmesh), function(i) {
  gIntersects(dmesh[i, ], domainSP1)})
w2.1 <- do.call(rbind,w2)
table(w2.1)
w2.1.1 <- which(w2.1)
length(w2.1.1)
w2.2 <- pblapply(w2.1.1, function(i) {
  gArea(gIntersection(dmesh[i, ], domainSP1))})
w2.2.1 <- do.call(rbind,w2.2) ##Needed areas

no.int <- setdiff(1:length(dmesh),w1.1.1)
df.val0 <- data.frame(poly=no.int,area=0)
int.1 <- w2.1.1
df.val1 <- data.frame(poly=int.1,area=w2.2.1)
df.val <- rbind(df.val0,df.val1)
melt.areas <- melt(df.val,id="poly")
areas.final <- cast(melt.areas,poly~variable,sum)
sum(areas.final$area) ##For validation
sum(areas.no)
w <- areas.final$area
#### Summary of these weights ####
sum(w)
table(w > 0)

##################################
par(mar = c(2, 2, 1, 1), mgp = 2:0)
plot(mesh$loc, asp = 1,  pch = 19, xlab = '', ylab = '',col=(w==0)+1) 
plot(dmesh, add = TRUE)
points(sppoints.gbif, col = "blue", pch = 19,cex=0.1)
##################################

### Getting the covariates ###

ext <- extent(c(min(mesh$loc[,1])-10,max(mesh$loc[,1])+10,min(mesh$loc[,2])-10,max(mesh$loc[,2])+10))
ext2_0 <- spTransform(SpatialPoints(matrix(ext,ncol=2,byrow=F),proj4string =CRS(proj4string(norway)) ),CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
ext2 <- extent(c(min(ext2_0@coords[,1]),max(ext2_0@coords[,1]),min(ext2_0@coords[,2]),max(ext2_0@coords[,2])))


## TRI ##

tri <- raster("current_30arcsec_tri.tif")
tri <- projectRaster(tri,crs=newcrs,method = "ngb")
tri_hed <- crop(tri,ext)
tri_hed1 <- mask(tri_hed, norway)
plot(tri_hed1,col=viridis(256))

## RAD ##
for(i in 1:12){
  if(i<10){assign(paste("rad",i,sep=""),raster(paste("wc2.0_30s_srad_0",i,".tif",sep="")))}
  else{assign(paste("rad",i,sep=""),raster(paste("wc2.0_30s_srad_",i,".tif",sep="")))}
  assign(paste("rad",i,sep=""),crop(get(paste("rad",i,sep="")),ext2))
  assign(paste("rad",i,sep=""),projectRaster(get(paste("rad",i,sep="")),crs=newcrs,method = "ngb"))
}

rad_hed <- mean(rad1,rad2,rad3,rad4,rad5,rad6,rad7,rad8,rad9,rad10,rad11,rad12)
rad_hed1 <- mask(rad_hed,norway)
plot(rad_hed1,col=viridis(256))

#### Distance ####
dist_hed <- raster("dist_hed.tif")
dist_hed1 <- mask(dist_hed,norway)
plot(dist_hed1, col=viridis(64, direction=-1))

## Points for prediction ##
a <- min(norway@polygons[[1]]@Polygons[[which.max(areas.no)]]@coords[,1])
b <- max(norway@polygons[[1]]@Polygons[[which.max(areas.no)]]@coords[,1])
c <- min(norway@polygons[[1]]@Polygons[[which.max(areas.no)]]@coords[,2])
d <- max(norway@polygons[[1]]@Polygons[[which.max(areas.no)]]@coords[,2])
x.pred <- seq(a-10,b+10,length.out = 1000)
y.pred <- seq(c-10,d+10,length.out = 1000)
locs.pred <- expand.grid(x.pred,y.pred)
locs.pred<- gIntersection(norway,SpatialPoints(locs.pred,proj4string = CRS(proj4string(norway))))


covpoints.mesh=matrix(nrow=mesh$n,ncol=5)
covpoints.obs=matrix(nrow=nrow(sppoints.gbif@coords),ncol=13)
covpoints.pred = matrix(nrow=nrow(locs.pred@coords),ncol=13)
covpoints.mesh[,1:2] <- meshpoints
covpoints.obs[,1:2] <- sppoints.gbif@coords
covpoints.pred[,1:2] <- locs.pred@coords

for(j in 3:ncol(covpoints.mesh)){
  if(j==3){
    covpoints.mesh[,j] <- extract(tri_hed,covpoints.mesh[,1:2])
    covpoints.obs[,j] <- extract(tri_hed,covpoints.obs[,1:2])
    covpoints.pred[,j] <- extract(tri_hed,covpoints.pred[,1:2])
  }
  if(j==4){
    covpoints.mesh[,j] <- extract(rad_hed,covpoints.mesh[,1:2])
    covpoints.obs[,j] <- extract(rad_hed,covpoints.obs[,1:2])
    covpoints.pred[,j] <- extract(rad_hed,covpoints.pred[,1:2])}
  if(j==5){
    covpoints.mesh[,j] <- extract(dist_hed,covpoints.mesh[,1:2])
    covpoints.obs[,j] <- extract(dist_hed,covpoints.obs[,1:2])
    covpoints.pred[,j] <- extract(dist_hed,covpoints.pred[,1:2])}
}

mean.tri <- mean(covpoints.pred[,3]);sd.tri <- sd(covpoints.pred[,3])
mean.rad <- mean(covpoints.pred[,4]);sd.rad <- sd(covpoints.pred[,4])

covpoints.mesh[,3] <- (covpoints.mesh[,3]-mean.tri)/sd.tri
covpoints.obs[,3] <- (covpoints.obs[,3]-mean.tri)/sd.tri
covpoints.pred[,3] <- (covpoints.pred[,3]-mean.tri)/sd.tri
covpoints.mesh[,4] <- (covpoints.mesh[,4]-mean.rad)/sd.rad
covpoints.obs[,4] <- (covpoints.obs[,4]-mean.rad)/sd.rad
covpoints.pred[,4] <- (covpoints.pred[,4]-mean.rad)/sd.rad

dist.obs <- -0.5*(covpoints.obs[,5])^2
dist.pred <- -0.5*(covpoints.pred[,5])^2
dist.mesh <- -0.5*(covpoints.mesh[,5])^2

n.c <- nrow(sppoints.gbif@coords)
y.pp.c <- rep(0:1, c(nv, n.c))
e.pp.c <- c(w, rep(0, n.c))

#### Defining the projection matrix ####
imat <- Diagonal(nv, rep(1, nv))
lmat.c <- inla.spde.make.A(mesh, sppoints.gbif@coords)
A.pp.c <- rbind(imat, lmat.c)

#### Defining the stack ####
stk.pp.c <- inla.stack(
  data = list(y = y.pp.c, e = e.pp.c), 
  A = list(1, A.pp.c),
  effects = list(list(b0 = 1,TRI=c(covpoints.mesh[,3],covpoints.obs[,3]),RAD=c(covpoints.mesh[,4],covpoints.obs[,4]),DIST=c(dist.mesh,dist.obs)), list(i = 1:nv)),
  tag = 'pp.c')
############################

#### Fitting the LGCP with INLA ... and summarizing the results ####
naive.model <- model2 <-  inla(y ~ 0 + b0  + RAD + TRI  +  f(i, model = spde),
                               family = 'poisson', data = inla.stack.data(stk.pp.c),
                               control.predictor = list(A = inla.stack.A(stk.pp.c)),
                               control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE),verbose=T,
                               control.inla = list(strategy = "gaussian", int.strategy = "eb"),E = inla.stack.data(stk.pp.c)$e)

bprior = list(theta = list(prior = "normal", param=c(1,0.05)))

vse.model <- model3 <- inla(y ~ 0 + b0  + RAD + TRI + + f(i, model = spde)+ f(DIST,model="clinear",range=c(0,Inf),hyper = bprior),
                            family = 'poisson', data = inla.stack.data(stk.pp.c),
                            control.predictor = list(A = inla.stack.A(stk.pp.c)),
                            control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE),verbose=T,
                            control.inla = list(strategy = "gaussian", int.strategy = "eb"),E = inla.stack.data(stk.pp.c)$e)

###########################
## Prediction in Hedmark ##
###########################

lmat.pred <- inla.spde.make.A(mesh, as.matrix(locs.pred@coords))
A.pred <- lmat.pred
stk.pred <- inla.stack(
  data = list(y = NA,e=rep(0,nrow(locs.pred@coords))),
  A = list(1,A.pred), 
  effects = list(list(b0 = 1,TRI=covpoints.pred[,3],RAD=covpoints.pred[,4],DIST=c(dist.pred)), list(i = 1:nv)),
  tag = 'pred.mu')
stk.full <- inla.stack(stk.pp.c,stk.pred)

####################
## Naïve Approach ##
####################
formu.naive <- model2$.args$formula
modforpred.naive <- inla(formu.naive,
                   family = 'poisson', data = inla.stack.data(stk.full),
                   control.mode = list(theta = model2$mode$theta,restart=FALSE),
                   control.predictor = list(A = inla.stack.A(stk.full),compute=TRUE),
                   control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE),verbose=TRUE,E = inla.stack.data(stk.full)$e,
                   control.inla = list(strategy = "gaussian", int.strategy = "eb"))


ind.pred <- inla.stack.index(stk.full, 'pred.mu')$data
df.pred <- modforpred.naive$summary.fitted.values[ind.pred,]
df.pred$medianpred <- (df.pred$`0.5quant`)
plotdfnew.naive=data.frame(xcoord=locs.pred@coords[,1],ycoord=locs.pred@coords[,2],median=df.pred$medianpred,mean=df.pred$mean,sd=df.pred$sd)
plotdfnew1.naive <- SpatialPointsDataFrame(plotdfnew.naive[,c("xcoord","ycoord")],data=plotdfnew.naive[,-(1:2)],proj4string = newcrs)
r <- raster(plotdfnew1.naive)
r1<-disaggregate(r, fact=res(r)/c(0.5,0.5))
rastxrespvar.naive <- rasterize(plotdfnew1.naive@coords,r1,plotdfnew1.naive$median, fun=mean,na.rm=T)
rastsd.naive <- rasterize(plotdfnew1.naive@coords,r1,plotdfnew1.naive$sd, fun=mean,na.rm=T)
par(mfrow=c(1,1))
plot(rastxrespvar.naive,col=viridis(64),axes=F,box=F)
plot(rastsd.naive,col=viridis(64),axes=F,box=F)


##################
## VSE Approach ##
##################
formu.vse <- model3$.args$formula
modforpred.vse <- inla(formu.vse,
                   family = 'poisson', data = inla.stack.data(stk.full),
                   control.mode = list(theta = model3$mode$theta,restart=FALSE),
                   control.predictor = list(A = inla.stack.A(stk.full),compute=TRUE),
                   control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE),verbose=TRUE,E = inla.stack.data(stk.full)$e,
                   control.inla = list(strategy = "gaussian", int.strategy = "eb"))


ind.pred <- inla.stack.index(stk.full, 'pred.mu')$data
df.pred <- modforpred.vse$summary.fitted.values[ind.pred,]
df.pred$medianpred <- (df.pred$`0.5quant`)
plotdfnew.vse=data.frame(xcoord=locs.pred@coords[,1],ycoord=locs.pred@coords[,2],median=df.pred$medianpred,mean=df.pred$mean,sd=df.pred$sd)
plotdfnew1.vse <- SpatialPointsDataFrame(plotdfnew.vse[,c("xcoord","ycoord")],data=plotdfnew.vse[,-(1:2)],proj4string = newcrs)

df.pred$DIST <- dist.pred
q.s <- modforpred.vse$summary.random$DIST[which(modforpred.vse$summary.random$DIST$ID %in% dist.pred),]
q.s1 <- c()
for(i in 1:nrow(df.pred)){q.s1[i] <- q.s$`0.5quant`[which(q.s$ID == df.pred$DIST[i])]}
df.pred$medianpred1 <- df.pred$`0.5quant` - q.s1

########################################
## Uncertainty for lambda(s)q(s)/q(s) ##
########################################

## Run in server ##
marg.pred <- modforpred$marginals.fitted.values[ind.pred]
marg.q.pred <- modforpred$marginals.random$DIST[which(modforpred$summary.random$DIST$ID %in% dist.pred)] 
B <- 10000
samp.lambda <-  matrix(nrow=length(marg.pred),ncol=B)
samp.q <- matrix(nrow=length(marg.pred),ncol=B)
for(i in 1:length(marg.pred)){
  if(max(marg.pred[[i]])==Inf){samp.lambda[i,]<-rep(modforpred$summary.fitted.values[ind.pred,][i,4],B)}
  else{a <- inla.rmarginal(marg.pred[[i]],n=B);samp.lambda[i,]<-a}
  if(max(marg.q.pred[[which(q.s$ID == df.pred$DIST[i])]])==Inf){samp.q[i,]<-rep(modforpred$summary.random$DIST[which(q.s$ID == df.pred$DIST[i]),5],B) }
  else{samp.q[i,] <- inla.rmarginal(marg.q.pred[[which(q.s$ID == df.pred$DIST[i])]],n=B)}
}

samp.lambdaoq <- samp.lambda-samp.q
postmed <- apply(samp.lambdaoq,1,median)
postq1 <-  apply(samp.lambdaoq,1,quantile,0.25)
postq3 <- apply(samp.lambdaoq,1,quantile,0.75)
postsd <- apply(samp.lambdaoq,1,sd)

## Results from server ##
vseunc <- readRDS("plotdfnew_vse.rds")
plotdfnew.vse <- vseunc
r <- raster(plotdfnew1.vse)
r1<-disaggregate(r, fact=res(r)/c(0.5,0.5))
rastxrespmed1.vse <- rasterize(cbind(plotdfnew.vse$xcoord,plotdfnew.vse$ycoord),r1,plotdfnew.vse$postmed, fun=mean)
rastxrespsd.vse <- rasterize(cbind(plotdfnew.vse$xcoord,plotdfnew.vse$ycoord),r1,plotdfnew.vse$sd, fun=max)
par(mfrow=c(1,1))
plot(rastxrespmed1.vse,col=viridis(64))
plot(rastxrespsd.vse,col=viridis(64))


plot(rastdiff <- -rastxrespvar.naive+rastxrespmed1.vse,axes=F,box=F,col=inferno(64, direction = 1))
plot(rastsddiff <- -rastsd.naive+rastxrespsd,axes=F,box=F,col=inferno(64, direction = 1))

par(mfrow=c(1,1))
transdist <- list()
qsum <- matrix(nrow=length(model3$marginals.random$DIST),ncol=3)
for(i in 1:length(model3$marginals.random$DIST)){
  transdist[[i]] <- inla.tmarginal(function(x) exp(x),model3$marginals.random$DIST[[i]]);
  qsum[i,] <- c(inla.zmarginal(transdist[[i]],silent = T)[[3]],inla.zmarginal(transdist[[i]],silent = T)[[5]],inla.zmarginal(transdist[[i]],silent = T)[[7]]) }

qsum.df <- data.frame(qsum)
names(qsum.df) <- c("q0.025","q0.5","q0.975")
qsum.df$dist <- sqrt(-2*model3$summary.random$DIST[,1])

par(mar=c(5,3,3,3))
plot((sqrt(-2*model3$summary.random$DIST[,1])),exp(model3$summary.random$DIST[,5]),type="l",xlab="Distance (Km)",ylab="q(s)")
lines((sqrt(-2*model3$summary.random$DIST[,1])),exp(model3$summary.random$DIST[,4]),lty=2,col="blue")
lines((sqrt(-2*model3$summary.random$DIST[,1])),exp(model3$summary.random$DIST[,6]),lty=2,col="blue")

### Focus on Zones ##
Zone1c <- SpatialPoints(data.frame(x=875,y=6535),proj4string = newcrs)
z1a <- c(861.7529,6399.844)
z1b <- c(903.3762,6399.844)
z1c <- c(861.7529,6372.493)
z1d <- c(903.3762,6372.493)
Z1 <- SpatialPolygons(list(Polygons(list(Polygon(rbind(z1a,z1b,z1d,z1c))),ID="1")),proj4string = newcrs)
Zone2c <- SpatialPoints(data.frame(x=910,y=6360),proj4string = newcrs)
z2a <-c(770.9316,6586.453)
z2b <- c(814.7563,6586.453)
z2c <- c(770.9316,6559.866)
z2d <- c(814.7563,6559.866)
Z2 <- SpatialPolygons(list(Polygons(list(Polygon(rbind(z2a,z2b,z2d,z2c))),ID="1")),proj4string = newcrs)
plot(Z2,add=T)
points(Zone2c)
sppz1 <- gIntersection(sppoints.gbif,Z1)
sppz2 <- gIntersection(sppoints.gbif,Z2)


par(mfrow=c(2,2),mar=c(1,0,0,0))
plot(crop(rastxrespvar.naive,Z1),zlim=c(-6,-2),axes=F,box=F,legend=F,col=viridis(64))
points(sppz1,cex=0.5,pch=19)
plot(crop(rastxrespmed1.vse,Z1),zlim=c(-6,-2),axes=F,box=F,col=viridis(64))
points(sppz1,cex=0.5,pch=19)
plot(crop(rastsd.naive,Z1),zlim=c(0.3,0.75),axes=F,box=F,legend=F,col=viridis(64))
points(sppz1,cex=0.5,pch=19)
plot(crop(rastxrespsd.vse,Z1),zlim=c(0.3,0.75),axes=F,box=F,col=viridis(64))
points(sppz1,cex=0.5,pch=19)

par(mfrow=c(2,2))
plot(crop(rastxrespvar.naive,Z2),zlim=c(-6.5,-3.5),axes=F,box=F,legend=F,col=viridis(64))
points(sppz2,cex=0.5,pch=19)
plot(crop(rastxrespmed1.vse,Z2),zlim=c(-6.5,-3.5),axes=F,box=F,col=viridis(64))
points(sppz2,cex=0.5,pch=19)
plot(crop(rastsd.naive,Z2),zlim=c(0.4,0.85),axes=F,box=F,legend=F,col=viridis(64))
points(sppz2,cex=0.5,pch=19)
plot(crop(rastxrespsd.vse,Z2),zlim=c(0.4,0.85),axes=F,box=F,col=viridis(64))
points(sppz2,cex=0.5,pch=19)
