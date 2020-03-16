library(sp)
library(rgeos)
library(INLA)
library(dplyr)
library(raster)
library(pbapply)
library(reshape)
library(tiff)
library(spatstat)
library(maptools)
library(gdata)
library(ggplot2)
library(gridExtra)
library(PCDSpline)
library(foreach)
library(doParallel)


setwd("/Users/jorgespa/Documents/Research/Work for paper 2") ## Remove later ##

## Polygon of Hedmark ##
newcrs <- CRS("+proj=robin +datum=WGS84 +units=km")
norwaybu <-raster::getData("GADM",country="NOR",level=1)
norwaybu <- spTransform(norwaybu,newcrs)
norway <- SpatialPolygons(list(norwaybu@polygons[[6]]))
proj4string(norway) <- CRS(proj4string(norwaybu))


## This is a grid, useful for prediction ##
a <- min(norway@polygons[[1]]@Polygons[[1]]@coords[,1])
b <- max(norway@polygons[[1]]@Polygons[[1]]@coords[,1])
c <- min(norway@polygons[[1]]@Polygons[[1]]@coords[,2])
d <- max(norway@polygons[[1]]@Polygons[[1]]@coords[,2])
x.pred <- seq(a-10,b+10,length.out = 100)
y.pred <- seq(c-10,d+10,length.out = 100)
locs.pred <- expand.grid(x.pred,y.pred)
locs.pred<- gIntersection(norway,SpatialPoints(locs.pred,proj4string = CRS(proj4string(norway))))

rast.dist <- raster("dist_hed.tif")

## Simulate a LGCP on Hedmark ##
win <- as.owin(norway)
x0 <- unique(locs.pred@coords[,1])
x0 <- x0[order(x0)]
y0 <- unique(locs.pred@coords[,2])
y0 <- y0[order(y0)]

ext <- extent(c(634.5525,1031.97,6249.143,6671.252))
ext2_0 <- spTransform(SpatialPoints(matrix(ext,ncol=2,byrow=F),proj4string =CRS(proj4string(norway)) ),CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
ext2 <- extent(c(min(ext2_0@coords[,1]),max(ext2_0@coords[,1]),min(ext2_0@coords[,2]),max(ext2_0@coords[,2])))
for(i in 1:12){
  if(i<10){assign(paste("rad",i,sep=""),raster(paste("wc2.0_30s_srad_0",i,".tif",sep="")))}
  else{assign(paste("rad",i,sep=""),raster(paste("wc2.0_30s_srad_",i,".tif",sep="")))}
  assign(paste("rad",i,sep=""),crop(get(paste("rad",i,sep="")),ext2))
  assign(paste("rad",i,sep=""),projectRaster(get(paste("rad",i,sep="")),crs=newcrs))
}

rad_hed <- mean(rad1,rad2,rad3,rad4,rad5,rad6,rad7,rad8,rad9,rad10,rad11,rad12)
rad_hed1 <- (rad_hed-8893.427)/218.3044 #Same mean/sd as in the example.
rad_hed1 <- crop(rad_hed1,norway)
rad_hed1 <- mask(rad_hed1,norway)

simulation <- function(pars,thinpar,zetap,BB,d1,scenario,basis.idxs = c(1,2,3)){
  
  if(scenario==1){beta0 <- pars[1];beta1 <- pars[2];sigma2x <- pars[3];range <- pars[4];zeta <- thinpar;zetap <- zetap;B=100}
  if(scenario==2){beta0 <- pars[1];beta1 <- pars[2];sigma2x <- pars[3];range <- pars[4];zeta <- thinpar;zetap <- zetap;d1=d1;p0=exp(-0.5*zeta*(d1^2));B=100}
  lg.s.c <- list()
  results <- list()
  naive <- list()
  vse <- list()
  vse.spl <- list()
  results0 <- list()
  naive0 <- list()
  vse0 <- list()
  vse.spl0  <-list()
  results0_0 <- list()
  naive0_0 <- list()
  vse0_0 <- list()
  vse.spl0_0 <- list()
  sampsize <- matrix(nrow=BB,ncol=length(zeta))
  set.seed(201019)
  for(k in 1:BB){
    lg.s.c[[k]] <- rLGCP('matern',beta0+beta1*as.im.RasterLayer(rad_hed1), var = sigma2x, scale = range / sqrt(8), 
                         nu = 1)
    
    ObsPoints_0 <- SpatialPointsDataFrame(data.frame(x=lg.s.c[[k]]$x,y=lg.s.c[[k]]$y),data= data.frame(ID=1:lg.s.c[[k]]$n),proj4string = newcrs)
    dist2road.obs0 <- extract(rast.dist,ObsPoints_0)
    dist.obs0 <- -0.5*(dist2road.obs0)^2
    
    
    ## Thin the PP based on the detection function
    
    thinlgcp.s <- function(pp,p.s){
      samp <- c()
      for(i in 1:nrow(pp)){samp[i] <- rbinom(1,1,p.s[i])}
      return(list(pp=pp[which(samp==1),],samp=samp))
    }
    
    for(m in 1:length(zeta)){
      if(scenario==1){thinprob.d <- exp(zeta[m]*dist.obs0)}
      if(scenario==2){thinprob.d <- ifelse(dist2road.obs0<=d1,exp(-0.5*zeta[m]*(dist2road.obs0^2)),exp(-0.5*zeta[m]*(d1^2)))}
      
      xy.c <- data.frame(x=lg.s.c[[k]]$x,y=lg.s.c[[k]]$y)
      tpp1.p.d <- thinlgcp.s(xy.c,thinprob.d)
      samp <- tpp1.p.d$samp
      tpp1.p.d <- tpp1.p.d$pp
      sampsize[k,m] <- nrow(tpp1.p.d)
      
      ### Making the mesh ###
      premesh<- inla.sp2segment(norway)
      mesh <- inla.mesh.2d(loc = tpp1.p.d, boundary = norway, max.edge = c(10,50),offset = c(1000, -0.2),cutoff=4,min.angle = 30)
      nv <- mesh$n
      meshpoints <- mesh$loc[,1:2]
      
      #### Constructing the dual mesh. Partially taken from Krainski et al 2018 ####
      source('spde-book-functions.R')
      dmesh <- book.mesh.dual(mesh)
      
      domain.polys <- norway@polygons
      domainSP0 <- SpatialPolygons(domain.polys)
      
      w1 <- pblapply(1:length(dmesh), function(i) {
        gIntersects(dmesh[i, ], domainSP0)})
      w1.1 <- do.call(rbind,w1)
      table(w1.1)
      w1.1.1 <- which(w1.1)
      length(w1.1.1)
      
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
      w2.2.1 <- do.call(rbind,w2.2) 
      
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
      
      points.mesh <- data.frame(x=meshpoints[,1],y=meshpoints[,2])
      points.obs <- tpp1.p.d
      
      MeshPoints <- SpatialPointsDataFrame(points.mesh ,data= data.frame(ID=1:nrow(points.mesh)),proj4string = newcrs)
      
      ## Distances to road ##
      dist2road.mesh <- extract(rast.dist,points.mesh)
      dist2road.obs <- extract(rast.dist,points.obs)
      dist.obs <- -0.5*(dist2road.obs)^2
      dist.mesh <- -0.5*(dist2road.mesh)^2
      
      ## Naïve model ##
      
      n.c <- nrow(tpp1.p.d)
      y.pp.c <- rep(0:1, c(nv, n.c))
      e.pp.c <- c(w, rep(0, n.c))
      
      #### Projection matrix ####
      imat <- Diagonal(nv, rep(1, nv))
      lmat.c <- inla.spde.make.A(mesh, as.matrix(tpp1.p.d))
      A.pp.c <- rbind(imat, lmat.c)
      
      cov.mesh <- extract(rad_hed1,points.mesh)
      cov.obs <- extract(rad_hed1,points.obs)
      
      #### Stack ####
      stk.pp.c <- inla.stack(
        data = list(y = y.pp.c, e = e.pp.c), 
        A = list(1, A.pp.c),
        effects = list(list(b0 = 1,cov=c(cov.mesh,cov.obs),dist=c(dist.mesh,dist.obs)), list(i = 1:nv)),
        tag = 'pp.c')
      
      spde <- inla.spde2.pcmatern(mesh, 
                                  prior.sigma = c(3, 0.05), 
                                  prior.range = c(10, 0.1))
      
      naive.model <- inla(y ~ 0 + b0  + cov + f(i, model = spde),
                          family = 'poisson', data = inla.stack.data(stk.pp.c), 
                          control.predictor = list(A = inla.stack.A(stk.pp.c)), 
                          control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE),verbose = F,
                          control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                          E = inla.stack.data(stk.pp.c)$e)
      
      
      naive.model$summary.fixed
      naive.model$summary.hyperpar
      
       
      #### VSE model ####
      
      bprior = list(theta = list(prior = "normal", param=c(zetap[m],0.05)))
      vse.model <- inla(y ~ 0 + b0  + cov +  f(i, model = spde)+ f(dist,model="clinear",range=c(0,Inf),hyper = bprior),
                        family = 'poisson', data = inla.stack.data(stk.pp.c), 
                        control.predictor = list(A = inla.stack.A(stk.pp.c)), 
                        control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE),verbose = F,
                        control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                        E = inla.stack.data(stk.pp.c)$e)
      
      
      #### EVSE Model ####
      
      Ibasis <- Ispline(c(dist2road.obs,dist2road.mesh), order=3, knots=seq(0,8,1))
      
      coefmat.obs <- -t(Ibasis[c(1,2,3,4,5,6,7,8,9,10),1:length(dist2road.obs)])
      coefmat.mesh <- -t(Ibasis[c(1,2,3,4,5,6,7,8,9,10),length(dist2road.obs)+(1:length(dist2road.mesh))])
      
      coefmat.all <- rbind(coefmat.obs,coefmat.mesh)
      
      stk.pp.c.1 <- inla.stack(
        data = list(y = y.pp.c, e = e.pp.c), 
        A = list(1, A.pp.c),
        effects = list(list(b0 = 1,cov=c(cov.mesh,cov.obs),dist1=c(coefmat.mesh[,basis.idxs[1]],coefmat.obs[,basis.idxs[1]]),dist2=c(coefmat.mesh[,basis.idxs[2]],coefmat.obs[,basis.idxs[2]]),dist3=c(coefmat.mesh[,basis.idxs[3]],coefmat.obs[,basis.idxs[3]])
        ),
        list(i = 1:nv)),
        tag = 'pp.c') 
      
      
      bprior = list(theta = list(prior = "normal", param=c(1,0.05)))
      
      vse.model.spl <- inla(y ~ 0 + b0  + cov +  f(i, model = spde)+ f(dist1,model="clinear",range=c(0,Inf),hyper=bprior)+
                              f(dist2,model="clinear",range=c(0,Inf),hyper=bprior)+f(dist3,model="clinear",range=c(0,Inf),hyper = bprior),
                            family = 'poisson', data = inla.stack.data(stk.pp.c.1),
                            control.predictor = list(A = inla.stack.A(stk.pp.c.1)),
                            control.compute = list(dic=TRUE,cpo=TRUE,waic=TRUE),verbose = F,
                            control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                            E = inla.stack.data(stk.pp.c)$e)
      
      vse.model.spl$summary.fixed
      vse.model.spl$summary.hyperpar
      
      results0[[m]] <- data.frame(real.bo=beta0,naive.b0=naive.model$summary.fixed$`0.5quant`[1],vse.b0=vse.model$summary.fixed$`0.5quant`[1],vse.spl.b0=vse.model.spl$summary.fixed$`0.5quant`[1],real.b1=beta1,naive.b1=naive.model$summary.fixed$`0.5quant`[2],vse.b1=vse.model$summary.fixed$`0.5quant`[2],vse.spl.b1=vse.model.spl$summary.fixed$`0.5quant`[2],real.range=range,naive.range=naive.model$summary.hyperpar$`0.5quant`[1],vse.range=vse.model$summary.hyperpar$`0.5quant`[1],vse.spl.range=vse.model.spl$summary.hyperpar$`0.5quant`[1],real.sd=sqrt(sigma2x),naive.sd=naive.model$summary.hyperpar$`0.5quant`[2],vse.sd=vse.model$summary.hyperpar$`0.5quant`[2],vse.spl.sd=vse.model.spl$summary.hyperpar$`0.5quant`[2],naive.dic=naive.model$dic$dic,vsl.dic=vse.model$dic$dic,vsl.spl.dic=vse.model.spl$dic$dic,naive.waic=naive.model$waic$waic,vse.waic=vse.model$waic$waic,vse.spl.waic=vse.model.spl$waic$waic,naive.cpo=-sum(log(naive.model$cpo$cpo)),vse.cpo=-sum(log(vse.model$cpo$cpo)),vse.spl.cpo=-sum(log(vse.model.spl$cpo$cpo))) 
      naive0[[m]] <- naive.model
      vse0[[m]] <- vse.model
      vse.spl0[[m]] <- vse.model.spl
      print(c(k,m))
    }
    results[[k]] <- results0
    naive[[k]] <- naive0
    vse[[k]] <- vse0
    vse.spl[[k]] <- vse.spl0
    
  }
  return(list(results=results,nv=naive,vse=vse,sampsize=sampsize))
}

pars1 <- data.frame(beta0=-4.25,beta1=0.82,s2=0.7,range=34)
pars <- as.numeric(pars1)
basis.comb <- matrix(c(2,4,6),ncol=3,byrow = T)
simu_results1 <- simulation(pars=pars,thinpar=c(0,1,8,16),zetap = rep(1,4),BB=1,scenario=1,basis.idxs=basis.comb[1,])

options <-  matrix(c(2,4,6,0.5),ncol=4,byrow = T)
simu_results2 <- simulation(pars=pars,thinpar=c(0,1,8,16),zetap = rep(1,4),d1=options[1,4],BB=1,scenario=2,basis.idxs=options[1,1:3])
