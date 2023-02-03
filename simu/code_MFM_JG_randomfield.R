require(spatstat)
require(ngspatial)
require(netmeta)
cut = function(d,c)
{as.numeric(d<=c)}
#fit2 = CDMFM_new1(datasimple, niterations = 500, alpha = 1, beta = 1, GAMMA = 1, LAMBDA = 1, initNClusters = 5,VN = VN)

## Dahl's method to summarize the samples from the MCMC
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

## diagonal on distance matrix is very large value
## function for Collapsed sampler
## ll is the parameter in logistic decay function
## neighhour is parameter in random field
CDMFM_new1 <- function(data, distance,ll,neighbour,niterations, alpha, beta, GAMMA, LAMBDA, initNClusters,VN)
{
  ## Model: n_{i}|z,lambda \sim Possion(lambda_{z_i}) ##
  ##        phi_{r} \sim Gamma(alpha,beta), r = 1,...,k ##
  ##        P(z_i = j) = \pi_j, j = 1,...,k ##
  ##        \pi \sim Dirichlet_k(GAMMA,...,GAMMA) ##
  ##        k-1 \sim possion(1) ##
  
  
  
  ################################################################
  
  ## Input: data = the vector of responses ##
  ##        niterations = the total number of iterations in MFM-SBM ##
  ##        alpha, beta = hyperparameters (shape, rate) for the prior on elements in lambda in Gamma distribution ##
  ##        GAMMA = the parameter in Dirichlet distribution that controls the relative size of clusters ##
  ##        LAMBDA = the parameter for Poisson distrition ##
  ##        initNClusters = the initial number of clusters ##
  
  
  ## Output: 
  ##         zout = clustering configuration, a n^2 by 1 vector##
  ##         phiout = possion parameters, a k by 1 vector ##
  
  #################################################################
  n = length(data)
  #precomputation for prespecified coefficient VN
  lambda <- LAMBDA
  gamma <- GAMMA
  N=n ## n is the number of oberservations
  # VN<-0
  # tmax = n+10
  # for (t in 1:tmax)
  # {
  #   r = log(0)
  #   for (k in t:500)
  #   {
  #     b = sum(log((k-t+1):k))-sum(log((k*gamma):(k*gamma+N-1))) + dpois(k-1, lambda, log = TRUE)
  #     m = max(b,r)
  #     r = log(exp(r-m) + exp(b-m)) + m
  #   }
  #   VN[t] = r
  # }
  #initialization of clustering configuration
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  
  phi<-rgamma(initNClusters, shape = alpha, rate = beta)
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
        #c.counts.noi = clusterSizes  #c.counts.noi corresponds to |C|
        #c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          #weight = sum(exp(-distance[i,clusterAssign_temp==x]))
          weight = sum(cut(distance[i,clusterAssign_temp==x],ll))
          cost = exp(0.5*sum(cut(distance[i,clusterAssign_temp==x],neighbour)))
          (GAMMA+weight)*cost*dpois(data[i],phi[x])
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        clusterProbs[nClusters+1]<-GAMMA*beta^alpha/gamma(alpha)*gamma(alpha+data[i])/(beta+1)^(data[i]+alpha)/factorial(data[i])*exp(VN[nClusters+1]-VN[nClusters])
        
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i > nClusters)
        {
          phinew = rep(0,nClusters+1)
          phinew[1:nClusters] = phi
          phinew[nClusters+1] = rgamma(1, shape = alpha, rate = beta)
          phi = phinew
          clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
          nClusters <- length(clusterSizes)} else
          {phi = phi
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes)}
      } else {
        # a singleton, have |C| choices
        #c.counts.noi = clusterSizes
        #c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1 - GAMMA# can offset the gamma adding later
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          #weight = sum(exp(-distance[i,clusterAssign_temp==x]))
          weight = sum(cut(distance[i,clusterAssign_temp==x],ll))
          cost = exp(0.5*sum(cut(distance[i,clusterAssign_temp==x],neighbour)))
          (GAMMA+weight)*cost*dpois(data[i],phi[x])
        })
        clusterProbs[clusterAssign[i]] = 0
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        clusterProbs[nClusters+1]<-GAMMA*beta^alpha/gamma(alpha)*gamma(alpha+data[i])/(beta+1)^(data[i]+alpha)/factorial(data[i])*exp(VN[nClusters]-VN[nClusters-1])
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
          phi = phi[-cur.cluster.i]}
      }
    }
    # end for loop over subjects i
    ## update phi ##
    phi = rep(0, nClusters)
    AA = rep(0,nClusters)
    NN = rep(0,nClusters)
    for (r in 1:nClusters){
      phi[r] = rgamma(1,alpha + sum(data[clusterAssign == r]), beta + sum(clusterAssign == r))
    }
    History[[iter]] <- list(zout = clusterAssign,phiout = phi)
    cat(" iteration:", iter,"\n",clusterAssign,"\n")
  }# for loop over iterations
  
  list(Iterates = History)
}


gamma = 1; lambda = 1; N = 100
VN<-0
tmax = 100+10
for (t in 1:tmax)
{
  r = log(0)
  for (k in t:500)
  {
    b = sum(log((k-t+1):k))-sum(log((k*gamma):(k*gamma+N-1))) + dpois(k-1, lambda, log = TRUE)
    m = max(b,r)
    r = log(exp(r-m) + exp(b-m)) + m
  }
  VN[t] = r
}

aaa = adjacency.matrix(10)

distance = netdistance(aaa)
diag(distance) = 10000
#distance = (distance-1)/10

result = CDMFM_new1(data = as.vector(count), distance = distance, ll = 150,neighbour = 8, niterations = 2000, alpha = 1, beta = 1, GAMMA = 1, LAMBDA = 1, initNClusters = 1,VN = VN)

# require(ngspatial)