library(igraph)
library(intergraph)
library(network)
library(glmnet)

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
Totalrep=50
fused_error = rep(0,Totalrep)
n=dim(betaMat)[2]
p=2
betaMat <- t(matrix(nrow = 159, ncol = p, byrow = TRUE)) #n:159, p=2
Countyshape<-readShapePoly('Counties_Georgia.shp')
A<-poly2nb(Countyshape)
adjacency<-nb2mat(A, zero.policy=TRUE)
distance = netdistance(adjacency)
diag(distance) = 10000 

X=array(runif(n*p*Totalrep,1,2), dim=c(n,p,Totalrep))
sim.data <- array(0,dim=c(n,Totalrep))


for (i in 1:159) {
  ## cluster 1
  betaMat[,asm == 1] <- rep(1.5,p)
  ## cluster 2
  betaMat[,asm == 2] <- rep(1,p)
}

X=array(runif(n*p*Totalrep,1,2), dim=c(n,p,Totalrep))
sim.data <- array(0,dim=c(n,Totalrep))

Georgia_graph =  graph_from_adjacency_matrix(adjacency!=0,mode ="undirected")
Georgia_tree = mst(Georgia_graph)


D=as_adjacency_matrix(Georgia_tree)

diag(D)=-1

D= as.matrix(cbind(D,D))

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
  
  Xnew = cbind(diag(X[,1,j]),diag(X[,2,j]))
  
  flasso = genlasso(y = sim.data[,j], X= Xnew, D=D)
  
  fused_error[j]=min(apply(flasso$beta -c(betaMat[1,],betaMat[2,]),2, function(x) sqrt(sum(x^2)/n/2)))
}

mean(fused_error)