###Case5
rm(list=ls())

require(spatstat)
require(expm)
library(ggmap);library(maps);library(data.table);library(fossil)
library(maptools);library(raster);library(spatialEco);library(fields);library(MASS)
library(MCMCpack);library(cowplot);library(spdep);library(phyclust)
library(maptools);library(rgdal);library(sf);library(matrixStats);library(netmeta)

library(MCMCglmm)
library(matrixStats)
library(CARBayes)



#Dahl
getDahl <- function(MFMfit, burn)
{
  ################################################################
  
  ## Input: MFMfit = the result from CDMFM_new ##
  ##        burn = the number of burn-in iterations ##
  
  ## Output: 
  ##         zout = estimated clustering configuration, a n by 1 vector##
  ##         Qout = estimated probability matrix, a k by k matrix ##
  
  #################################################################
  iters <- MFMfit$Iterates[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x[[1]]
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  DahlAns
}

#MFM
CDMFM_new1 <- function(data, niterations, alpha, kappa, GAMMA, initNClusters,VN,V,X,Mpart)
{
  ################################################################
  
  ## Model: n_{i}|z,lambda \sim Possion(X%*%beta_{z_i}) ##
  ##        beta_{r} \sim MLGamma(0,IV,alpha,kappa), r = 1,...,k ##
  ##        P(z_i = j) = \pi_j, j = 1,...,k ##
  ##        \pi \sim Dirichlet_k(GAMMA,...,GAMMA) ##
  ##        k-1 \sim possion(1) ##
  
  ################################################################
  
  ## Input: data = the vector of responses ##
  ##        niterations = the total number of iterations in MFM-MLG ##
  ##        alpha, kappa = hyperparameters (shape, rate) for the prior on elements in beta in Multivariate Log Gamma distribution ##
  ##        GAMMA = the parameter in Dirichlet distribution that controls the relative size of clusters ##
  ##        initNClusters = the initial number of clusters ##
  
  ## Output: 
  ##         zout = clustering configuration, a n by 1 vector##
  ##         phiout = possion parameters, a p by k vector ##
  
  #################################################################
  n = length(data)
  #gamma <- GAMMA
  p=dim(X)[2]
  #V=sqrtm(sig2.V*diag(p))
  IV=solve(V)
  
  #initialization of clustering configuration
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  
  #phi<-rgamma(initNClusters, shape = alpha, rate = beta)
  phi=matrix(0,p,initNClusters)
  for(j in 1:initNClusters){
    phi[,j]=MLG(matrix(0,p,1),V,alpha,kappa) 
  }
  History <- vector("list", niterations)
  
  ##start Gibb's sampling
  for (iter in 1:niterations)
  {
    ## update z ##
    clusterSizes = table(as.factor(clusterAssign))
    nClusters = length(clusterSizes)
    for (i in 1:n)
    { #determine whether ith component is a singleton 
      cur.cluster.i = clusterAssign[i]
      if (clusterSizes[clusterAssign[i]] > 1){
        # not a singleton, have |C|+1 choices
        c.counts.noi = clusterSizes  #c.counts.noi corresponds to |C|
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          #if(nClusters>1){
          (GAMMA+c.counts.noi[x])*dpois(data[i],exp(X[i,]%*%phi[,x]))#}else{
          #(GAMMA+c.counts.noi[x])*dpois(data[i],exp(X[i,]%*%as.matrix(phi)))
          #}
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        clusterProbs[nClusters+1]<-Mpart[[i]]$result*GAMMA*exp(VN[nClusters+1]-VN[nClusters])
        
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i > nClusters)
        {
          #phinew = rep(0,nClusters+1)
          #phinew[1:nClusters] = phi
          phinew = matrix(0,p,nClusters+1)
          phinew[,1:nClusters] = phi
          phinew[,nClusters+1] = MLG(matrix(0,p,1),V,alpha,kappa)
          phi = phinew
          clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
          nClusters <- length(clusterSizes)} else
          {phi = phi
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes)}
      } else {
        # a singleton, have |C| choices
        c.counts.noi = clusterSizes
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1 - GAMMA# can offset the gamma adding later
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          #if(nClusters > 1){
          (GAMMA+c.counts.noi[x])*dpois(data[i],exp(X[i,]%*%phi[,x]))#}else{
          #(GAMMA+c.counts.noi[x])*dpois(data[i],exp(X[i,]%*%as.matrix(phi)))
          #}
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        clusterProbs[nClusters+1]<-Mpart[[i]]$result*GAMMA*exp(VN[nClusters+1]-VN[nClusters])
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        # remove the empty cluster
        if (cluster.i > nClusters)
        {      clusterAssign[i] <- cur.cluster.i #put the new cluster in the place of the only singleten one
        clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
        } else
        {      
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i, clusterAssign-1, clusterAssign) # to delete the previous group index
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes) 
          #phi = phi[-cur.cluster.i]}
          phi = phi[,-cur.cluster.i]}
      }
      #print(i)
    }
    # end for loop over subjects i
    ## update phi ##
    #phi = rep(0, nClusters)
    phi = matrix(0,p,nClusters)
    for (r in 1:nClusters){
      # #set up for cMLG
      # Hphi=rbind(IV,X) #could precompute
      # alpha_phi=rbind(alpha*matrix(1,p,1),sum(data[clusterAssign == r])*matrix(1,n,1)) 
      # kappa_phi=rbind(kappa*matrix(1,p,1),sum(clusterAssign == r)*matrix(1,n,1))       
      # phi[,r]=solve(t(Hphi)%*%Hphi)%*%t(Hphi)%*%log(rgamma((p+n),alpha_phi,rate=kappa_phi)) #+1e-4
      # #phi[,r]=solve(t(Hphi)%*%Hphi)%*%t(Hphi)%*%MLG(matrix(0,(n+p),1),diag(n+p),alpha_phi,kappa_phi) #+1e-4
      
      Hphi=rbind(IV,X[which(clusterAssign == r),])
      alpha_phi=rbind(alpha*matrix(1,p,1),as.matrix(data[clusterAssign == r])) 
      kappa_phi=rbind(kappa*matrix(1,p,1),matrix(1,sum(clusterAssign == r),1))       
      phi[,r]=solve(t(Hphi)%*%Hphi)%*%t(Hphi)%*%log(rgamma((p+sum(clusterAssign == r)),alpha_phi,rate=kappa_phi)) #+1e-4
      
    }
    History[[iter]] <- list(zout = clusterAssign,phiout = phi)
    cat(" iteration:", iter,"\n",clusterAssign,"\n")
  }# for loop over iterations
  
  list(Iterates = History)
}

#Precompute Mpart
CNC=function(data,sig2.V,sig2.V1,alpha,kappa,alpha1,kappa1,X){
  library(cubature)
  X=as.matrix(X)
  X=t(X)
  nn=dim(X)[1]
  #nn=1
  p=dim(X)[2]
  V=sqrtm(sig2.V*diag(p))
  #COEF=1/sqrt(det(V%*%t(V)))*(kappa^alpha/gamma(alpha))^p
  logCOEF= -1/2*log(det(V%*%t(V)))+p*(alpha*log(kappa)-lgamma(alpha))
  COEF=exp(logCOEF)
  
  IV=solve(V)
  H=rbind(IV,X)
  
  VHinv=solve(sqrtm(sig2.V1*diag(nn+p)))      #sig2.V1: could change the value
  
  Q2=VHinv[,-(1:p)]*sig2.V1
  
  #numeric constant (using numeric integral)
  V1=solve(cbind(H,Q2))
  #if(data <=0){
  #f=function(x){1/det(V1%*%t(V1))*(kappa1^alpha1/gamma(alpha1))^p*1/gamma(data)*exp(alpha1*matrix(1,1,(nn+p))%*%solve(V1)%*%c(x[1],x[2],rep(0,nn))-matrix(kappa1,1,(nn+p))%*%exp(solve(V1)%*%c(x[1],x[2],rep(0,nn))))}
  #NC=adaptIntegrate(f, lowerLimit = rep(-Inf,p), upperLimit = rep(Inf,p))$integral
  
  #M1=det(cbind(H,Q2))*(kappa1^alpha1/gamma(alpha1))^p*1/gamma(data)/NC
  #result=COEF/M1
  #}else{
  #gammadata=gamma(170)
  #f=function(x){1/det(V1%*%t(V1))*(kappa1^alpha1/gamma(alpha1))^p*1/gamma(data)*exp(alpha1*matrix(1,1,(nn+p))%*%solve(V1)%*%c(x[1],x[2],rep(0,nn))-matrix(kappa1,1,(nn+p))%*%exp(solve(V1)%*%c(x[1],x[2],rep(0,nn))))}
  
  #In the "f", x[1],...x[p] ; if p=2, c(x[1],x[2],rep(0,nn));;; if p=1, c(x[1],rep(0,nn))
  f=function(x){exp(alpha1*matrix(1,1,(nn+p))%*%solve(V1)%*%c(x[1],x[2],rep(0,nn))-matrix(kappa1,1,(nn+p))%*%exp(solve(V1)%*%c(x[1],x[2],rep(0,nn))))}
  NC=adaptIntegrate(f, lowerLimit = rep(-Inf,p), upperLimit = rep(Inf,p))$integral
  logConstant= -log(det(V1%*%t(V1)))+p*(alpha1*log(kappa1)-lgamma(alpha1))-lgamma(data)
  NC=logConstant+log(NC) 
  
  logM1=log(det(cbind(H,Q2)))+p*(alpha1*log(kappa1)-lgamma(alpha1))-lgamma(data)-NC
  result=COEF/exp(logM1)
  #}
  return(list(result=result,V=V,alpha=alpha,kappa=kappa))
}

#MLG
MLG=function(C,V,shape, rate){
  dim=length(C)
  omega=matrix(0,dim,1)
  omega=log(rgamma(dim,shape=shape,rate=rate))
  #omega=log(1/rate)+log(rgamma(dim,shape=shape,rate=1)) #or omega=log(scale)+log(rgamma(dim,shape=shape,scale=1))
  q=C+V%*%omega
  return(q)
}

#cut
cut = function(d,c)
{as.numeric(d<=c)}

#MRF_MFM
CDMFM_MLG <- function(data,distance,neighbour, niterations, alpha, kappa, GAMMA, initNClusters,VN,V,X,Mpart,tune)
{
  ################################################################
  
  ## Model: n_{i}|z,lambda \sim Possion(X%*%beta_{z_i}) ##
  ##        beta_{r} \sim MLGamma(0,IV,alpha,kappa), r = 1,...,k ##
  ##        P(z_i = j) = \pi_j, j = 1,...,k ##
  ##        \pi \sim Dirichlet_k(GAMMA,...,GAMMA) ##
  ##        k-1 \sim possion(1) ##
  
  ################################################################
  
  ## Input: data = the vector of responses ##
  ##        niterations = the total number of iterations in MFM-MLG ##
  ##        alpha, kappa = hyperparameters (shape, rate) for the prior on elements in beta in Multivariate Log Gamma distribution ##
  ##        GAMMA = the parameter in Dirichlet distribution that controls the relative size of clusters ##
  ##        initNClusters = the initial number of clusters ##
  
  ## Output: 
  ##         zout = clustering configuration, a n by 1 vector##
  ##         phiout = possion parameters, a p by k vector ##
  
  #################################################################
  n = length(data)
  #gamma <- GAMMA
  p=dim(X)[2]
  #V=sqrtm(sig2.V*diag(p))
  IV=solve(V)
  
  #initialization of clustering configuration
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  
  #phi<-rgamma(initNClusters, shape = alpha, rate = beta)
  phi=matrix(0,p,initNClusters)
  for(j in 1:initNClusters){
    phi[,j]=MLG(matrix(0,p,1),V,alpha,kappa) 
  }
  History <- vector("list", niterations)
  
  ##start Gibb's sampling
  for (iter in 1:niterations)
  {
    ## update z ##
    clusterSizes = table(as.factor(clusterAssign))
    nClusters = length(clusterSizes)
    for (i in 1:n)
    { #determine whether ith component is a singleton 
      cur.cluster.i = clusterAssign[i]
      if (clusterSizes[clusterAssign[i]] > 1){
        # not a singleton, have |C|+1 choices
        c.counts.noi = clusterSizes  #c.counts.noi corresponds to |C|
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          #if(nClusters>1){
          #(GAMMA+c.counts.noi[x])*dpois(data[i],exp(X[i,]%*%phi[,x]))#}else{
          #(GAMMA+c.counts.noi[x])*dpois(data[i],exp(X[i,]%*%as.matrix(phi)))
          #}
          #weight = sum(cut(distance[i,clusterAssign_temp==x],ll))
          cost = exp(tune*sum(cut(distance[i,clusterAssign_temp==x],neighbour)))
          (GAMMA+c.counts.noi[x])*cost*dpois(data[i],exp(X[i,]%*%phi[,x]))
          
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        clusterProbs[nClusters+1]<-Mpart[[i]]$result*GAMMA*exp(VN[nClusters+1]-VN[nClusters])
        
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i > nClusters)
        {
          #phinew = rep(0,nClusters+1)
          #phinew[1:nClusters] = phi
          phinew = matrix(0,p,nClusters+1)
          phinew[,1:nClusters] = phi
          phinew[,nClusters+1] = MLG(matrix(0,p,1),V,alpha,kappa)
          phi = phinew
          clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
          nClusters <- length(clusterSizes)} else
          {phi = phi
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes)}
      } else {
        # a singleton, have |C| choices
        c.counts.noi = clusterSizes
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1 - GAMMA# can offset the gamma adding later
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          #if(nClusters > 1){
          #(GAMMA+c.counts.noi[x])*dpois(data[i],exp(X[i,]%*%phi[,x]))#}else{
          #(GAMMA+c.counts.noi[x])*dpois(data[i],exp(X[i,]%*%as.matrix(phi)))
          #}
          #weight = sum(cut(distance[i,clusterAssign_temp==x],ll))
          cost = exp(tune*sum(cut(distance[i,clusterAssign_temp==x],neighbour)))
          (GAMMA+c.counts.noi[x])*cost*dpois(data[i],exp(X[i,]%*%phi[,x]))
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        clusterProbs[nClusters+1]<-Mpart[[i]]$result*GAMMA*exp(VN[nClusters+1]-VN[nClusters])
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        # remove the empty cluster
        if (cluster.i > nClusters)
        {      clusterAssign[i] <- cur.cluster.i #put the new cluster in the place of the only singleten one
        clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
        } else
        {      
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i, clusterAssign-1, clusterAssign) # to delete the previous group index
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes) 
          #phi = phi[-cur.cluster.i]}
          phi = phi[,-cur.cluster.i]}
      }
      #print(i)
    }
    # end for loop over subjects i
    ## update phi ##
    #phi = rep(0, nClusters)
    phi = matrix(0,p,nClusters)
    for (r in 1:nClusters){
      # #set up for cMLG
      # Hphi=rbind(IV,X) #could precompute
      # alpha_phi=rbind(alpha*matrix(1,p,1),sum(data[clusterAssign == r])*matrix(1,n,1)) 
      # kappa_phi=rbind(kappa*matrix(1,p,1),sum(clusterAssign == r)*matrix(1,n,1))       
      # phi[,r]=solve(t(Hphi)%*%Hphi)%*%t(Hphi)%*%log(rgamma((p+n),alpha_phi,rate=kappa_phi)) #+1e-4
      # #phi[,r]=solve(t(Hphi)%*%Hphi)%*%t(Hphi)%*%MLG(matrix(0,(n+p),1),diag(n+p),alpha_phi,kappa_phi) #+1e-4
      
      Hphi=rbind(IV,X[which(clusterAssign == r),])
      alpha_phi=rbind(alpha*matrix(1,p,1),as.matrix(data[clusterAssign == r])) 
      kappa_phi=rbind(kappa*matrix(1,p,1),matrix(1,sum(clusterAssign == r),1))       
      phi[,r]=solve(t(Hphi)%*%Hphi)%*%t(Hphi)%*%log(rgamma((p+sum(clusterAssign == r)),alpha_phi,rate=kappa_phi)) #+1e-4
      
    }
    History[[iter]] <- list(zout = clusterAssign,phiout = phi)
    cat(" iteration:", iter,"\n",clusterAssign,"\n")
  }# for loop over iterations
  
  list(Iterates = History)
}

######Draw Map#####
centroids <- readRDS("GAcentroids.rds")
centroids <- as.data.frame(centroids)

a <- 1:159
map.county <- map_data('county')
counties   <- unique(map.county[, 5:6])
louiscounties <- counties[counties$region == "georgia", ]
age_map <- data.frame(
  state_names = louiscounties$region,
  county_names = louiscounties$subregion,
  Region = factor(a)
)
age_map <- data.table(age_map)
setkey(age_map, state_names, county_names)

map.county <-
  data.table(map.county[map.county$region == "georgia", ])
setkey(map.county, region, subregion)

map.df <- map.county[age_map]

for (i in 1:nrow(centroids)){
  if (centroids$x[i] - 2 * centroids$y[i] < -150) {
    map.df[map.df$subregion == row.names(centroids)[i], "Region"] <- "1"
  } else if (centroids$x[i] + centroids$y[i] > -51) {
    map.df[map.df$subregion == row.names(centroids)[i], "Region"] <- "2"
  } else {
    map.df[map.df$subregion == row.names(centroids)[i], "Region"] <- "3"
  }
}


map.df$Region <- factor(map.df$Region)
map.df$Cluster <- map.df$Region

asm=NA
for (i in 1:nrow(centroids)) {
  if ((centroids$x[i] - 2 * centroids$y[i] < -152) | (centroids$x[i] - 2 * centroids$y[i] > -147.5)) {
    asm[i] <- 1
  } else {
    asm[i] <- 2
  }
}

p=2
betaMat <- t(matrix(nrow = 159, ncol = p, byrow = TRUE)) #n:159, p=2

for (i in 1:159) {
  ## cluster 1
  betaMat[,asm == 1] <- rep(1.5,p)
  ## cluster 2
  betaMat[,asm == 2] <- rep(1,p)
}

n=dim(betaMat)[2]
#Precompute VN
GAMMA = 1; lambda = 1; N = n 
VN<-0
tmax = n+10
for (t in 1:tmax)
{
  rr = log(0)
  for (k in t:259)
  {
    b = sum(log((k-t+1):k))-sum(log((k*GAMMA):(k*GAMMA+N-1))) + dpois(k-1, lambda, log = TRUE)  #lambda or exp(X%*%betaMat)
    m = max(b,rr)
    rr = log(exp(rr-m) + exp(b-m)) + m
  }
  VN[t] = rr
}

###Setting
n.iter=3000
burn_in=1000
times=1
Totalrep=50
Cluster_MFM=matrix(0,Totalrep,1)
LPML_MFM=matrix(0,Totalrep,1)
RI.MFM=matrix(0,n.iter,Totalrep)

Cluster_MRF025=matrix(0,Totalrep,1)
LPML_MRF025=matrix(0,Totalrep,1)
RI.MRF025=matrix(0,n.iter,Totalrep)
beta025=list()

Cluster_MRF05=matrix(0,Totalrep,1)
LPML_MRF05=matrix(0,Totalrep,1)
RI.MRF05=matrix(0,n.iter,Totalrep)
beta05=list()

Cluster_MRF075=matrix(0,Totalrep,1)
LPML_MRF075=matrix(0,Totalrep,1)
RI.MRF075=matrix(0,n.iter,Totalrep)
beta075=list()

Cluster_MRF1=matrix(0,Totalrep,1)
LPML_MRF1=matrix(0,Totalrep,1)
RI.MRF1=matrix(0,n.iter,Totalrep)
beta1=list()

LPML_LGP=matrix(0,Totalrep,1)
LPML_CAR=matrix(0,Totalrep,1)


###
Countyshape<-readShapePoly('Counties_Georgia.shp')
A<-poly2nb(Countyshape)
adjacency<-nb2mat(A, zero.policy=TRUE)
distance = netdistance(adjacency)
diag(distance) = 10000 

X=array(runif(n*p*Totalrep,1,2), dim=c(n,p,Totalrep))
sim.data <- array(0,dim=c(n,Totalrep))

fit_MFM=list()
fita_MRF025=list()
fita_MRF05=list()
fita_MRF075=list()
fita_MRF1=list()


result_MFM=list()
Result_MFMMLG025=list()
Result_MFMMLG05=list()
Result_MFMMLG075=list()
Result_MFMMLG1=list()

###Run MCMC; Contain MFM and MRF
for(j in 1:Totalrep){
  
  ###Generate Data
  X[,,j]=matrix(runif(n*p,1,2),n,p)
  
  while (TRUE) {
    for(i in 1:n){
      sim.data[i,j]=rpois(1,exp(X[i,,j]%*%betaMat[,i]))
    }
    check=length(which(sim.data[,j]==0))
    if(check==0){
      break
    }
  }
  
  Mpart <- vector("list", n)
  for(i in 1:n){
    #Mpart[i]=CNC(data=sim.data[i],sig2.V=10,sig2.V1=10,alpha=1,kappa=10,alpha1=1,kappa1=10,X[i,])
    Mpart[[i]]=CNC(data=sim.data[i,j],sig2.V=10,sig2.V1=10,alpha=1,kappa=1,alpha1=1,kappa1=1,X=X[i,,j])
    if ((i%%10) == 0){print(i)}
  }
  
  ####MFM####
  fit_MFM[[j]]=CDMFM_new1(data=sim.data[,j], niterations=n.iter, alpha=Mpart[[1]]$alpha, kappa=Mpart[[1]]$kappa, 
                          GAMMA=GAMMA, initNClusters=4,VN=VN,V=Mpart[[1]]$V,X=X[,,j],Mpart = Mpart)
  result_MFM[[j]] =getDahl(fit_MFM[[j]],burn = burn_in)[[1]]
  MCR_MFM = sapply(1:times, function(x) {length(unique(result_MFM[[j]]))})
  Cluster_MFM[j]=MCR_MFM
  
  ###Calculate CPO
  CPOinv=matrix(0,n,(n.iter-(burn_in)))
  CPO=rep(0,n)
  for(i in 1:n){
    for(s in (burn_in+1):n.iter){
      CPOinv[i,(s-burn_in)]= -dpois(sim.data[i,j],lambda = exp(X[i,,j]%*%fit_MFM[[j]][[1]][[s]]$phiout[,fit_MFM[[j]][[1]][[s]]$zout[i]]),log=TRUE)
    }
    #CPO[i]<-1/mean(CPOinv[i,])
    CPO[i]<- log((n.iter-(burn_in)))-logSumExp(CPOinv[i,])    #log scale
  }
  LPML_MFM[j]=sum(CPO)
  
  for(q in 1:n.iter){
    RI.MFM[q,j]=rand.index(asm,fit_MFM[[j]][[1]][[q]]$zout)
  }
  
  ####MRF 0.25####
  fita_MRF025[[j]]=CDMFM_MLG(data = sim.data[,j], distance = distance , neighbour=4 , niterations =n.iter,
                             alpha = 1, kappa = 1, GAMMA = 1, initNClusters = 4,VN = VN, V = Mpart[[1]]$V,
                             X = X[,,j],Mpart = Mpart, tune = 0.25)
  Result_MFMMLG025[[j]]=getDahl(fita_MRF025[[j]],burn = burn_in)[[1]]
  beta025[[j]]=getDahl(fita_MRF025[[j]],burn = burn_in)[[2]]
  MCR_MFMMLG = sapply(1:times, function(x) {length(unique(Result_MFMMLG025[[j]]))})
  Cluster_MRF025[j]=MCR_MFMMLG
  
  ###Calculate CPO
  CPOinv=matrix(0,n,(n.iter-(burn_in)))
  CPO=rep(0,n)
  for(i in 1:n){
    for(s in (burn_in+1):n.iter){
      CPOinv[i,(s-burn_in)]= -dpois(sim.data[i,j],lambda = exp(X[i,,j]%*%fita_MRF025[[j]][[1]][[s]]$phiout[,fita_MRF025[[j]][[1]][[s]]$zout[i]]),log=TRUE)
    }
    CPO[i]<- log((n.iter-(burn_in)))-logSumExp(CPOinv[i,])      #log scale
  }
  LPML_MRF025[j]=sum(CPO)
  
  #Save the RI for each iteration over replicate
  for(q in 1:n.iter){
    RI.MRF025[q,j]=rand.index(asm,fita_MRF025[[j]][[1]][[q]]$zout)
  }
  
  ####MRF 0.5####
  fita_MRF05[[j]]=CDMFM_MLG(data = sim.data[,j], distance = distance , neighbour=4 , niterations =n.iter,
                            alpha = 1, kappa = 1, GAMMA = 1, initNClusters = 4,VN = VN, V = Mpart[[1]]$V,
                            X = X[,,j],Mpart = Mpart, tune = 0.5)
  Result_MFMMLG05[[j]]=getDahl(fita_MRF05[[j]],burn = burn_in)[[1]]
  beta05[[j]]=getDahl(fita_MRF05[[j]],burn = burn_in)[[2]]
  MCR_MFMMLG = sapply(1:times, function(x) {length(unique(Result_MFMMLG05[[j]]))})
  Cluster_MRF05[j]=MCR_MFMMLG
  
  ###Calculate CPO
  CPOinv=matrix(0,n,(n.iter-(burn_in)))
  CPO=rep(0,n)
  for(i in 1:n){
    for(s in (burn_in+1):n.iter){
      CPOinv[i,(s-burn_in)]= -dpois(sim.data[i,j],lambda = exp(X[i,,j]%*%fita_MRF05[[j]][[1]][[s]]$phiout[,fita_MRF05[[j]][[1]][[s]]$zout[i]]),log=TRUE)
    }
    CPO[i]<- log((n.iter-(burn_in)))-logSumExp(CPOinv[i,])      #log scale
  }
  LPML_MRF05[j]=sum(CPO)
  
  #Save the RI for each iteration over replicate
  for(q in 1:n.iter){
    RI.MRF05[q,j]=rand.index(asm,fita_MRF05[[j]][[1]][[q]]$zout)
  }
  
  
  ####MRF 0.75####
  fita_MRF075[[j]]=CDMFM_MLG(data = sim.data[,j], distance = distance , neighbour=4 , niterations =n.iter,
                             alpha = 1, kappa = 1, GAMMA = 1, initNClusters = 4,VN = VN, V = Mpart[[1]]$V,
                             X = X[,,j],Mpart = Mpart, tune = 0.75)
  Result_MFMMLG075[[j]]=getDahl(fita_MRF075[[j]],burn = burn_in)[[1]]
  beta075[[j]]=getDahl(fita_MRF075[[j]],burn = burn_in)[[2]]
  MCR_MFMMLG = sapply(1:times, function(x) {length(unique(Result_MFMMLG075[[j]]))})
  Cluster_MRF075[j]=MCR_MFMMLG
  
  ###Calculate CPO
  CPOinv=matrix(0,n,(n.iter-(burn_in)))
  CPO=rep(0,n)
  for(i in 1:n){
    for(s in (burn_in+1):n.iter){
      CPOinv[i,(s-burn_in)]= -dpois(sim.data[i,j],lambda = exp(X[i,,j]%*%fita_MRF075[[j]][[1]][[s]]$phiout[,fita_MRF075[[j]][[1]][[s]]$zout[i]]),log=TRUE)
    }
    CPO[i]<- log((n.iter-(burn_in)))-logSumExp(CPOinv[i,])      #log scale
  }
  LPML_MRF075[j]=sum(CPO)
  
  #Save the RI for each iteration over replicate
  for(q in 1:n.iter){
    RI.MRF075[q,j]=rand.index(asm,fita_MRF075[[j]][[1]][[q]]$zout)
  }
  
  ####MRF 1####
  fita_MRF1[[j]]=CDMFM_MLG(data = sim.data[,j], distance = distance , neighbour=4 , niterations =n.iter,
                           alpha = 1, kappa = 1, GAMMA = 1, initNClusters = 4,VN = VN, V = Mpart[[1]]$V,
                           X = X[,,j],Mpart = Mpart, tune = 1)
  Result_MFMMLG1[[j]]=getDahl(fita_MRF1[[j]],burn = burn_in)[[1]]
  beta1[[j]]=getDahl(fita_MRF1[[j]],burn = burn_in)[[2]]
  MCR_MFMMLG = sapply(1:times, function(x) {length(unique(Result_MFMMLG1[[j]]))})
  Cluster_MRF1[j]=MCR_MFMMLG
  
  ###Calculate CPO
  CPOinv=matrix(0,n,(n.iter-(burn_in)))
  CPO=rep(0,n)
  for(i in 1:n){
    for(s in (burn_in+1):n.iter){
      CPOinv[i,(s-burn_in)]= -dpois(sim.data[i,j],lambda = exp(X[i,,j]%*%fita_MRF1[[j]][[1]][[s]]$phiout[,fita_MRF1[[j]][[1]][[s]]$zout[i]]),log=TRUE)
    }
    CPO[i]<- log((n.iter-(burn_in)))-logSumExp(CPOinv[i,])      #log scale
  }
  LPML_MRF1[j]=sum(CPO)
  
  #Save the RI for each iteration over replicate
  for(q in 1:n.iter){
    RI.MRF1[q,j]=rand.index(asm,fita_MRF1[[j]][[1]][[q]]$zout)
  }
  
  #LGP
  #organize data
  y.LGP = sim.data[,j]
  y.LGP = data.frame(y.LGP)
  
  obs.LGP = X[,,j]
  #obs.LGP = X
  dat.LGP = cbind(y.LGP)
  dat2.LGP = data.frame(dat.LGP)
  
  
  #prior specifications
  PBV.yfixed <- diag(1) * 1
  prior.m2a.5.1 <- list(B = list(mu = rep(0, 1),V = PBV.yfixed))
  
  
  #fit Poisson LGP using pacakge MCMCglmm
  model1.LGP <- MCMCglmm(y.LGP ~ 1, family = "poisson"
                         ,data=dat2.LGP,pl = TRUE,pr=TRUE,prior = prior.m2a.5.1,singular.ok=TRUE
                         ,nitt=n.iter,burnin=burn_in,thin=1)
  
  lams=t(model1.LGP$Liab)
  ###Calculate CPO
  CPOinv.LGP=matrix(0,n,(n.iter-(burn_in)))
  CPO.LGP=rep(0,n)
  for(i in 1:n){
    for(s in (burn_in+1):n.iter){
      CPOinv.LGP[i,(s-burn_in)]= -dpois(sim.data[i,j],lambda = exp(lams[i,(s-burn_in)]),log=TRUE)
    }
    CPO.LGP[i]<- log((n.iter-(burn_in)))-logSumExp(CPOinv.LGP[i,])      #log scale
  }
  LPML_LGP[j]=sum(CPO.LGP)
  
  #CAR
  formula.CAR <- sim.data[,j]~X[,,j]
  W.CAR=1*exp(rdist(centroids))
  model.CAR <- S.CARbym(formula=formula.CAR, family="poisson",W=W.CAR, burnin=burn_in, n.sample=n.iter)
  LPML_CAR[j]=model.CAR$modelfit[5]
  
  print(j)
}





save.image("~/MRF/Case5withEst.RData")