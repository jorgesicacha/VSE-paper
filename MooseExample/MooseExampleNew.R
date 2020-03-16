options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")

library(sp)
library(rgeos)
library(INLA)
library(dplyr)
library(raster)
library(pbapply)
library(reshape)
library(tiff)
library(ggplot2)
library(gridExtra)

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

source('spde-book-functions.R')
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
x.pred <- seq(a-10,b+10,length.out = 500)
y.pred <- seq(c-10,d+10,length.out =500)
locs.pred <- expand.grid(x.pred,y.pred)
locs.pred<- gIntersection(norway,SpatialPoints(locs.pred,proj4string = CRS(proj4string(norway))))
points.pred <- locs.pred@coords

covpoints.mesh=matrix(nrow=mesh$n,ncol=5)
covpoints.obs=matrix(nrow=nrow(sppoints.gbif@coords),ncol=5)
covpoints.pred = matrix(nrow=nrow(locs.pred@coords),ncol=5)
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

vse.model <- model3 <- inla(y ~ 0 + b0  + RAD + TRI  + f(i, model = spde)+ f(DIST,model="clinear",range=c(0,Inf),hyper = bprior),
                            family = 'poisson', data = inla.stack.data(stk.pp.c),
                            control.predictor = list(A = inla.stack.A(stk.pp.c)),
                            control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE),verbose=T,
                            control.inla = list(strategy = "gaussian", int.strategy = "eb"),E = inla.stack.data(stk.pp.c)$e)


library(PCDSpline)
Ibasis <- Ispline(c(covpoints.obs[,5],covpoints.mesh[,5],covpoints.pred[,5]), order=3, knots=seq(0,8,1))

coefmat.obs <- -t(Ibasis[c(1,2,3,4,5,6,7,8,9,10),1:length(covpoints.obs[,5])])
coefmat.mesh <- -t(Ibasis[c(1,2,3,4,5,6,7,8,9,10),length(covpoints.obs[,5])+(1:length(covpoints.mesh[,5]))])
coefmat.pred <- -t(Ibasis[c(1,2,3,4,5,6,7,8,9,10),length(covpoints.obs[,5])+length(covpoints.mesh[,5])+(1:length(covpoints.pred[,5]))])

coefmat.all <- rbind(coefmat.obs,coefmat.mesh,coefmat.pred)

stk.pp.c.1 <- inla.stack(
  data = list(y = y.pp.c, e = e.pp.c), 
  A = list(1, A.pp.c),
  effects = list(list(b0 = 1,TRI=c(covpoints.mesh[,3],covpoints.obs[,3]),RAD=c(covpoints.mesh[,4],covpoints.obs[,4]),
                      dist1=c(coefmat.mesh[,2],coefmat.obs[,2]),dist2=c(coefmat.mesh[,3],coefmat.obs[,3]),dist3=c(coefmat.mesh[,5],coefmat.obs[,5])),
                 list(i = 1:nv)),
  tag = 'pp.c') 

bprior = list(theta = list(prior = "normal", param=c(1,0.05)))

vse.model.spl <- model4 <- inla(y ~ 0 + b0  + TRI + RAD +   f(i, model = spde)+ f(dist1,model="clinear",range=c(0,Inf),hyper=bprior)+
                f(dist2,model="clinear",range=c(0,Inf),hyper=bprior)+f(dist3,model="clinear",range=c(0,Inf),hyper = bprior),
              family = 'poisson', data = inla.stack.data(stk.pp.c.1),
              control.predictor = list(A = inla.stack.A(stk.pp.c.1)),
              control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE),verbose = T,
              control.inla = list(strategy = "gaussian", int.strategy = "eb"),
              E = inla.stack.data(stk.pp.c)$e)
  

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
stk.full1 <- inla.stack(stk.pp.c.1,stk.pred)

#################
## Naïve Model ##
#################
formu.naive <- model2$.args$formula
modforpred.naive <- inla(formu.naive,
                         family = 'poisson', data = inla.stack.data(stk.full),
                         control.mode = list(theta = model2$mode$theta,restart=FALSE),
                         control.predictor = list(A = inla.stack.A(stk.full),compute=TRUE),
                         control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE),verbose=TRUE,E = inla.stack.data(stk.full)$e,
                         control.inla = list(strategy = "gaussian", int.strategy = "eb"))


ind.pred <- inla.stack.index(stk.full, 'pred.mu')$data
df.pred <- data.frame(points.pred,modforpred.naive$summary.fitted.values[ind.pred,])
m1 <- lapply(modforpred.naive$marginals.fitted.values[ind.pred], function(x) inla.tmarginal(exp,x))
m2 <- lapply(m1, function(x) inla.qmarginal(0.5,x))
rm(m1)
names(df.pred)[c(1:2,5:7)] <- c("x","y","quant0025","quant05","quant0975")
df.pred$int <- do.call("rbind",m2)

pred.naive <- ggplot()+
  geom_point(data=df.pred,mapping=aes(x,y,color=quant05),size=2.5)

sd.naive <- ggplot()+
  geom_point(data=df.pred,mapping=aes(x,y,color=sd),size=2.5)

grid.arrange(pred.naive,sd.naive,nrow=1)

saveRDS(df.pred,"dfpred500NewNew.rds")

###############
## VSE Model ##
###############

formu.vse <- model3$.args$formula
modforpred.vse <- inla(formu.vse,
                       family = 'poisson', data = inla.stack.data(stk.full),
                       control.mode = list(theta = model3$mode$theta,restart=FALSE),
                       control.predictor = list(A = inla.stack.A(stk.full),compute=TRUE),
                       control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE),verbose=TRUE,E = inla.stack.data(stk.full)$e,
                       control.inla = list(strategy = "gaussian", int.strategy = "eb"))


ind.pred <- inla.stack.index(stk.full, 'pred.mu')$data
df.pred.vse <- data.frame(points.pred,modforpred.vse$summary.fitted.values[ind.pred,])
m1.vse <- lapply(modforpred.vse$marginals.fitted.values[ind.pred], function(x) inla.tmarginal(exp,x))
m2.vse <- lapply(m1.vse, function(x) inla.qmarginal(0.5,x))
rm(m1.vse)
names(df.pred.vse)[c(1:2,5:7)] <- c("x","y","quant0025","quant05","quant0975")
df.pred.vse$int <- do.call("rbind",m2.vse)

ind.q <- c()
for(i in 1:length(dist.pred)){
  ind.q[i] <- which(modforpred.vse$summary.random$DIST$ID == dist.pred[i])
  print(i)
}

q.s <- modforpred.vse$summary.random$DIST[ind.q,]
q.s.marg <- modforpred.vse$marginals.random$DIST[ind.q]
B <- 10000
s5 <- list()
s6 <- list()

for(i in 1:length(dist.pred)){
  s1 <- inla.tmarginal(exp,marginal = modforpred.vse$marginals.fitted.values[ind.pred][[i]])
  s2 <- inla.tmarginal(exp,marginal = q.s.marg[[i]])
  s3 <- inla.rmarginal(s1,n=B)
  s3.1 <- inla.rmarginal(marginal =modforpred.vse$marginals.fitted.values[ind.pred][[i]],n=B)
  s4 <- inla.rmarginal(s2,n=B)
  s4.1 <- inla.rmarginal(marginal=q.s.marg[[i]],n=B)
  s5[[i]] <- (s3/s4)
  s6[[i]] <- (s3.1-s4.1)
  print(i)
}

pred.median.vse <- do.call("rbind",lapply(s5,median))
pred.sd.vse <- do.call("rbind",lapply(s5,sd))
pred.log.median.vse <- do.call("rbind",lapply(s6,median))
pred.log.sd.vse <- do.call("rbind",lapply(s6,sd))

df.pred.vse2 <- data.frame(points.pred,modforpred.vse$summary.fitted.values[ind.pred,])
df.pred.vse2$logint <- pred.log.median.vse
df.pred.vse2$int <- pred.median.vse
df.pred.vse2$logsd <- pred.log.sd.vse
df.pred.vse2$sd <- pred.sd.vse

pred.vse <- ggplot()+
  geom_point(data=df.pred.vse2,mapping = aes(x,y,color=logint),size=2.5)

sd.vse <- ggplot()+
  geom_point(data=df.pred.vse2,mapping=aes(x,y,color=logsd),size=2.5)

grid.arrange(pred.vse,sd.vse,nrow=1)

##############################################
## Plotting the relationship between d and q ##
##############################################
s2 <- matrix(nrow=length(dist.pred),ncol=3)
for(i in 1:length(dist.pred)){
  s2_0 <- inla.tmarginal(exp,marginal = q.s.marg[[i]])
  s2[i,] <- inla.qmarginal(p=c(0.025,0.5,0.975),marginal = s2_0)
  print(i)
}

dfplotq.vse <- data.frame(dist=covpoints.pred[,5],q0025=s2[,1],q05=s2[,2],q975=s2[,3])

#####################
### VSE-SPL Model ###
#####################

  stk.pred <- inla.stack(
    data = list(y = NA,e=rep(0,nrow(locs.pred@coords))),
    A = list(1,A.pred), 
    effects = list(list(b0 = 1,TRI=covpoints.pred[,3],RAD=covpoints.pred[,4],dist1=coefmat.pred[,2],dist2=coefmat.pred[,3],dist3=coefmat.pred[,5]), list(i = 1:nv)),
    tag = 'pred.mu')
  stk.full <- inla.stack(stk.pp.c.1,stk.pred)
  
  formu.vse.spl <- vse.model.spl$.args$formula
  modforpred.vse.spl <- inla(formu.vse.spl,
                             family = 'poisson', data = inla.stack.data(stk.full),
                             control.mode = list(theta = vse.model.spl$mode$theta,restart=FALSE),
                             control.predictor = list(A = inla.stack.A(stk.full),compute=TRUE),
                             control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE),verbose=T,E = inla.stack.data(stk.full)$e,
                             control.inla = list(strategy = "gaussian", int.strategy = "eb"))
  
  
  ind.pred <- inla.stack.index(stk.full, 'pred.mu')$data
  df.pred.vse.spl <- data.frame(points.pred,modforpred.vse.spl$summary.fitted.values[ind.pred,])


m1.vse.spl <- lapply(modforpred.vse.spl$marginals.fitted.values[ind.pred], function(x) inla.tmarginal(exp,x))
m2.vse.spl <- lapply(m1.vse.spl, function(x) inla.qmarginal(0.5,x))
rm(m1.vse.spl)
names(df.pred.vse.spl)[c(1:2,5:7)] <- c("x","y","quant0025","quant05","quant0975")
df.pred.vse.spl$dist <- rup(covpoints.pred[,5])
df.pred.vse.spl$int <- do.call("rbind",m2.vse.spl)
df.pred.vse.spl$areacell <- cell.pred.areas
df.pred.vse.spl$count <- df.pred.vse.spl$int * df.pred.vse.spl$areacell
sum(df.pred.vse.spl$count)
summdfpred.vse.spl <- group_by(df.pred.vse.spl,dist)
summdfpred.vse.spl <- summarise(summdfpred.vse.spl,avgint = sum(count))

ggplot()+
  geom_bar(data=summdfpred.vse.spl,mapping = aes(x=dist,y=avgint/sum(avgint)),stat="identity",alpha=0.3,col="red")+
  geom_histogram(data=df.just.all,aes(x=dist.obs,y=stat(width*density)),breaks=seq(0,8,0.1),alpha=0.3,col="blue")+
  theme_bw()

ind.q.spl1 <- c()
ind.q.spl2 <- c()
ind.q.spl3 <- c()
combins1 <- matrix(c(2,3,5),ncol=3,byrow=T)
x=1
for(i in 1:length(dist.pred)){
  ind.q.spl1[i] <- which(modforpred.vse.spl$summary.random$dist1$ID == coefmat.pred[i,combins1[x,1]])
  ind.q.spl2[i] <- which(modforpred.vse.spl$summary.random$dist2$ID == coefmat.pred[i,combins1[x,2]])
  ind.q.spl3[i] <- which(modforpred.vse.spl$summary.random$dist3$ID == coefmat.pred[i,combins1[x,3]])
  print(i)
} 

q.s.spl1 <- modforpred.vse.spl$summary.random$dist1[ind.q.spl1,]
q.s.spl2 <- modforpred.vse.spl$summary.random$dist2[ind.q.spl2,]
q.s.spl3 <- modforpred.vse.spl$summary.random$dist3[ind.q.spl3,]
q.s.marg.1 <- modforpred.vse.spl$marginals.random$dist1[ind.q.spl1]
q.s.marg.2 <- modforpred.vse.spl$marginals.random$dist2[ind.q.spl2]
q.s.marg.3 <- modforpred.vse.spl$marginals.random$dist3[ind.q.spl3]
B <- 10000
s5<- list()
s6 <- list()
for(i in 1:length(dist.pred)){
  
  
  s1 <- inla.tmarginal(exp,marginal = modforpred.vse.spl$marginals.fitted.values[ind.pred][[i]])
  s2.1 <-inla.tmarginal(exp,marginal = q.s.marg.1[[i]])
  s2.2 <-inla.tmarginal(exp,marginal = q.s.marg.2[[i]])
  s2.3 <-inla.tmarginal(exp,marginal = q.s.marg.3[[i]])
  s3 <- inla.rmarginal(s1,n=B)
  s3.1 <- inla.rmarginal(n=B,marginal = modforpred.vse.spl$marginals.fitted.values[ind.pred][[i]])
  s4.1 <- inla.rmarginal(s2.1,n=B)
  s4.1.1 <- inla.rmarginal(n=B,marginal = q.s.marg.1[[i]])
  s4.2 <- inla.rmarginal(s2.2,n=B)
  s4.2.1 <- inla.rmarginal(n=B,marginal = q.s.marg.2[[i]])
  s4.3 <- inla.rmarginal(s2.3,n=B)
  s4.3.1 <- inla.rmarginal(n=B,marginal = q.s.marg.3[[i]])
  s4 <- s4.1*s4.2*s4.3
  s4_1 <- s4.1.1 + s4.2.1 + s4.3.1
  s5[[i]] <- s3/s4
  s6[[i]] <- s3.1 - s4_1
  #}
  print(i)
}

pred.median.vse.spl <- do.call("rbind",lapply(s5,median))
pred.sd.vse.spl <- do.call("rbind",lapply(s5,sd))
pred.log.median.vse.spl <- do.call("rbind",lapply(s6,median))
pred.log.sd.vse.spl <- do.call("rbind",lapply(s6,sd))



df.pred.vse2.spl <- data.frame(points.pred,modforpred.vse.spl$summary.fitted.values[ind.pred,])
df.pred.vse2.spl$logint <- pred.log.median.vse.spl
df.pred.vse2.spl$int <- pred.median.vse.spl
df.pred.vse2.spl$logsd <- pred.log.sd.vse.spl
df.pred.vse2.spl$sd <- pred.sd.vse.spl

pred.vse.spl <- ggplot()+
  geom_point(data=df.pred.vse2.spl,mapping = aes(x,y,color=logint),size=2.5)

sd.vse.spl <- ggplot()+
  geom_point(data=df.pred.vse2.spl,mapping=aes(x,y,color=logsd),size=2.5)

grid.arrange(pred.vse.spl,sd.vse.spl,nrow=1)


###############################################
## Plotting the relationship between d and q ##
###############################################
s2spl <- matrix(nrow=length(dist.pred),ncol=3)
for(i in 1:length(dist.pred)){
  s2_0_1 <- inla.tmarginal(exp,marginal = q.s.marg.1[[i]])
  s2_0_2 <- inla.tmarginal(exp,marginal = q.s.marg.2[[i]])
  s2_0_3 <- inla.tmarginal(exp,marginal = q.s.marg.3[[i]])
  
  s4.1 <- inla.rmarginal(s2_0_1,n=B)
  s4.2 <- inla.rmarginal(s2_0_2,n=B)
  s4.3 <- inla.rmarginal(s2_0_3,n=B)
  s4 <- s4.1*s4.2*s4.3
  
  s2spl[i,] <- quantile(x = s4,probs=c(0.025,0.5,0.975))
  print(i)
}

dfplotq.vsespl <- data.frame(dist=covpoints.pred[,5],q0025=s2spl[,1],q05=s2spl[,2],q975=s2spl[,3])


colors <- c("VSE" = "blue", "EVSE" = "red")
ggplot()+
  geom_line(data=dfplotq.vse,mapping = aes(dist,q05,col="VSE"))+
  geom_line(data=dfplotq.vsespl,mapping = aes(dist,q05,col="EVSE"))+
  theme_bw()+
  labs(y = "q(s)",x="d(s)",color="Model")+
  scale_color_manual(values = colors)

