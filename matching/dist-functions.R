smahal <- function(z,X){
  X <- as.matrix(X)
  n <- dim(X)[1]
  rownames(X) <- 1:n
  k <- dim(X)[2]
  m <- sum(z)
  for(j in 1:k){
    X[,j] <- rank(X[,j])
  }
  cv <- cov(X)
  vuntied <- var(1:n)
  rat <- sqrt(vuntied/diag(cv))
  cv <- diag(rat)%*%cv%*%diag(rat)
  out <- matrix(NA,m,n-m)
  Xc <- X[z==0,]
  Xt <- X[z==1,]
  rownames(out) <- rownames(X)[z==1]
  colnames(out) <- rownames(X)[z==0]
  icov <- ginv(cv)
  for(i in 1:m){
    out[i,] <- mahalanobis(Xc,Xt[i,],icov,inverted=T)
  }
  out
}

# Function for adding a propensity score caliper to a distance matrix dmat
# calipersd is the caliper in terms of standard deviation of the logit propensity score

addcaliper <- function(dmat,z,logitp,calipersd=.2,penalty=1000){
  sd.logitp <- sd(logitp)
  adif <- abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif <- (adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat <- dmat+adif*penalty
  dmat
}
# Add penalty for not matching on a covariate
addpenalty.binarycov <- function(dmat,x,z,penalty=2){
  adif <- abs(outer(x[z==1],x[z==0],"-"))
  dmat <- dmat+adif*penalty
  dmat
}

stratified.matched.set <- function(match.vec.list, z){
  subject.match.order <- c()
  stratum.numeric <- c()
  #treated <- treated[rownames(X.mat.full)]
  for(i in 1:length(match.vec.list)){
    match.vec.tmp <- match.vec.list[[i]]
    subject.match.order <- c(subject.match.order, as.numeric(names(match.vec.tmp)))
    stratum.short.tmp <- substr(match.vec.tmp, start = 3, stop = 10)
    stratum.numeric <- c(stratum.numeric, (as.numeric(stratum.short.tmp) + i*10000)) # makes sure we don't repeat names for the matched sets across strata
  }
  sort.unique.stratum <- sort(unique(stratum.numeric))
  stratum.myindex.matchvecorder <- rep(0, length(stratum.numeric))
  for(i in 1:length(sort.unique.stratum)){
    stratum.myindex.matchvecorder[stratum.numeric == sort.unique.stratum[i]] <- i
  }
  
  #stratum.myindex <- rep(0, length(treated))
  stratum.myindex <- rep(0,  length(z))
  stratum.max <- 0
  stratum.myindex[subject.match.order] <- stratum.myindex.matchvecorder
  stratum.max <- max(stratum.myindex, na.rm = TRUE)
  
  #names(stratum.myindex) <- names(treated)
  names(stratum.myindex) <- names(z)
  stratum.myindex[stratum.myindex == 0] <- NA
  stratum.myindex
}