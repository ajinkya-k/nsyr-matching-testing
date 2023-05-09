library(MASS)
library(optmatch)
library(Hmisc)

# Rank based Mahalanobis distance
# Function for computing
# rank based Mahalanobis distance.  Prevents an outlier from
# inflating the variance for a variable, thereby decreasing its importance.
# Also, the variances are not permitted to decrease as ties
# become more common, so that, for example, it is not more important
# to match on a rare binary variable than on a common binary variable
# z is a vector, length(z)=n, with z=1 for treated, z=0 for control
# X is a matrix with n rows containing variables in the distance

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


# Function to find best VR match


# Function to extract matched sets
# matched.set <- function(match.vec, X.mat.full, treated){

# [SKD]: 20 November 2018 ... Add a line at the end to replace 0's with NA

matched.set <- function(match.vec, z){
  subject.match.order <- c()
  stratum.numeric <- c()
  #treated <- treated[rownames(X.mat.full)]

  subject.match.order <- as.numeric(names(match.vec))
  stratum.short.tmp <- substr(match.vec, start = 3, stop = 10)
  stratum.numeric <- as.numeric(stratum.short.tmp)
  
  sort.unique.stratum <- sort(unique(stratum.numeric))
  stratum.myindex.matchvecorder <- rep(0, length(stratum.numeric))
  for(i in 1:length(sort.unique.stratum)){
    # if ( sum(stratum.numeric == sort.unique.stratum[i], na.rm = T) > 5) print(i) 
    stratum.myindex.matchvecorder[stratum.numeric == sort.unique.stratum[i]] <- i
  }
  #stratum.myindex <- rep(0, length(treated))
  #stratum.myindex <- rep(0, length(z))
  stratum.myindex <- rep(NA, length(z))
  stratum.max <- 0
  stratum.myindex[subject.match.order] <- stratum.myindex.matchvecorder
  stratum.max <- max(stratum.myindex, na.rm = TRUE)
  names(stratum.myindex) <- names(z)
  stratum.myindex
}


# Update: 27 November 2018. unmatched individuals will be saved as NA not 0
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



# Let us build a data frame containing all of the covariate data and also the matched set


# [SKD]: November 20, 2018: we will allow matched.sets to include NA values to signify that someone is unmatched
# 

stratified.df <- function(cov.data, out.data, matched.sets, treated){
  
  #n.strata <- length(unique(matched.sets[!is.na(matched.set)])) # number of actual matched sets
  if(any(is.na(matched.sets))){
    stratum.list <- c(unique(matched.sets[!is.na(matched.sets)]), NA)
    n.strata <- length(stratum.list) - 1
  } else{
    stratum.list <- unique(matched.sets)
    n.strata <- length(stratum.list)
  }
  #stratum.list <- c(unique(matched.sets[!is.na(matched.sets)]), NA) # does not include the NA's
  
  #n.subjects <- sum(matched.sets != 0)
  n.subjects <- nrow(cov.data)
  
  if("block" %in% colnames(cov.data)){
    
    strat.cov <- data.frame("treated" = rep(NA, times = n.subjects),
                            "block" = rep(NA, times = n.subjects),
                            "stratum" = rep(NA, times = n.subjects),
                            "weight" = rep(NA, times = n.subjects))
    strat.out <- data.frame("treated" = rep(NA, times = n.subjects),
                            "block" = rep(NA, times = n.subjects),
                            "stratum" = rep(NA, times = n.subjects),
                            "weight" = rep(NA, times = n.subjects))    
  } else{
    
    strat.cov <- data.frame("treated" = rep(NA, times = n.subjects),
                            "stratum" = rep(NA, times = n.subjects),
                            "weight" = rep(NA, times = n.subjects))
    strat.out <- data.frame("treated" = rep(NA, times = n.subjects),
                            "stratum" = rep(NA, times = n.subjects),
                            "weight" = rep(NA, times = n.subjects))
  }


  #rownames(strat.cov) <- names(matched.sets[matched.sets != 0])
  #rownames(strat.out) <- names(matched.sets[matched.sets != 0])
  rownames(strat.cov) <- names(matched.sets)
  rownames(strat.out) <- names(matched.sets)
  
  cov.names <- colnames(cov.data)[!colnames(cov.data) %in% c("cohort", "block", "treated")]
  out.names <- colnames(out.data)[!colnames(out.data) %in% c("cohort", "block", "treated")]
  
  strat.cov[,cov.names] <- NA
  strat.out[,out.names] <- NA
  #print("starting loop")
  #for(i in unik.stratum.list){
    #print(paste("i = ", i))
  # this loops over only the valid matched sets (i.e. not the unmatched individuals)
  for(i in 1:n.strata){
    subjects.tmp <- names(which(matched.sets == stratum.list[i]))
    treated.tmp <- names(which(treated[subjects.tmp] == 1))
    control.tmp <- names(which(treated[subjects.tmp] == 0))
    n.stratum <- length(subjects.tmp)
    n.treated <- length(treated.tmp)

    strat.cov[treated.tmp, "treated"] <- 1
    strat.cov[control.tmp, "treated"] <- 0
    if("block" %in% colnames(cov.data)) strat.cov[subjects.tmp, "block"] <- cov.data[subjects.tmp,"block"]
    strat.cov[subjects.tmp,"stratum"] <- i
    
    strat.out[treated.tmp, "treated"] <- 1
    strat.out[control.tmp, "treated"] <- 0
    if("block" %in% colnames(cov.data)) strat.out[subjects.tmp, "block"] <- cov.data[subjects.tmp,"block"]
    strat.out[subjects.tmp, "stratum"] <- i
    
    # weighting to treated population
    strat.cov[treated.tmp, "weight"] <- 1
    strat.cov[control.tmp, "weight"] <- n.treated/(n.stratum - n.treated)
    
    strat.out[treated.tmp, "weight"] <- 1
    strat.out[control.tmp, "weight"] <- n.treated/(n.stratum - n.treated)
 
    
    #print(paste(" added weights"))
    
    strat.cov[subjects.tmp, cov.names] <- cov.data[subjects.tmp, cov.names]
    strat.out[subjects.tmp, out.names] <- out.data[subjects.tmp, out.names]
  }
  
  # Now we need to add in all of the individuals whose matched.set is na
  if(any(is.na(matched.sets))){
    subjects.tmp <- names(which(is.na(matched.sets)))
    treated.tmp <- names(which(treated[subjects.tmp] == 1))
    control.tmp <- names(which(treated[subjects.tmp] == 0))
    n.stratum <- length(subjects.tmp)
    n.treated <- length(treated.tmp)
    
    strat.cov[treated.tmp, "treated"] <- 1
    strat.cov[control.tmp, "treated"] <- 0
    if("block" %in% colnames(cov.data)) strat.cov[subjects.tmp, "block"] <- cov.data[subjects.tmp,"block"]
    #strat.cov[subjects.tmp, "stratum"] <- ifelse(i != 0, i, NA)
    strat.cov[subjects.tmp,"stratum"] <- NA
    
    strat.out[treated.tmp, "treated"] <- 1
    strat.out[control.tmp, "treated"] <- 0
    if("block" %in% colnames(cov.data)) strat.out[subjects.tmp, "block"] <- cov.data[subjects.tmp,"block"]
    strat.out[subjects.tmp, "stratum"] <- NA
    
    # weighting to treated population
    strat.cov[subjects.tmp, "weight"] <- 0
    strat.out[subjects.tmp, "weight"] <- 0

    strat.cov[subjects.tmp, cov.names] <- cov.data[subjects.tmp, cov.names]
    strat.out[subjects.tmp, out.names] <- out.data[subjects.tmp, out.names]
  }

  strat.cov <- strat.cov[order(strat.cov$stratum),]
  strat.out <- strat.out[order(strat.out$stratum),]
  results <- list("cov" = strat.cov, "out" = strat.out)
  #results <- list("cov" = strat.cov)
  return(results)
}

# 
num.balance.table <- function(match.obj){
  cov.data <- match.obj$cov
  #out.data <- match.obj$out
  cov.names <- colnames(cov.data)[-(1:3)] # ignore the first three columns
  n.cov <- length(cov.names)
  # For now pretend everything is continuous
  match.table <- data.frame("Treated.Mean" = NA, "After.Treated.Mean" = NA, "Before.Control.Mean" = NA, "After.Control.Mean" = NA, "Before.Std.diff" = NA, "After.Std.diff" = NA)
  match.table.missing <- data.frame("Treated.Mean" = NA, "Before.Control.Mean" = NA, "After.Control.Mean" = NA, "Before.Std.diff" = NA, "After.Std.diff" = NA)
  raw.treated.index <- which(cov.data[,"treated"] == 1)
  raw.control.index <- which(cov.data[,"treated"] == 0)

  after.treated.index <- which(cov.data[,"treated"] == 1 & !is.na(cov.data[,"stratum"]))
  after.control.index <- which(cov.data[,"treated"] == 0 & !is.na(cov.data[,"stratum"]))
  for(v in cov.names){
    combined.sd <- sqrt(0.5*var(cov.data[raw.treated.index,v], na.rm = TRUE) + 0.5*var(cov.data[raw.control.index,v], na.rm = TRUE))
    treated.mean <- mean(cov.data[raw.treated.index,v], na.rm = TRUE)
    after.treated.mean <- wtd.mean(cov.data[after.treated.index,v], weights = cov.data[after.treated.index,"weight"])
    before.control.mean <- mean(cov.data[raw.control.index,v], na.rm = TRUE)
    after.control.mean <- wtd.mean(cov.data[after.control.index,v], weights = cov.data[after.control.index,"weight"])
    std.diff.before <- (treated.mean - before.control.mean)/combined.sd
    std.diff.after <- (after.treated.mean - after.control.mean)/combined.sd
    
    if(grepl(".missing", v)){
      #TODO: add after.treated.mean in the mising table
      match.table.missing[v,] <- c(treated.mean, before.control.mean, after.control.mean, std.diff.before, std.diff.after)
    } else{
      match.table[v,] <- c(treated.mean, after.treated.mean, before.control.mean, after.control.mean, std.diff.before, std.diff.after)
    }
  }
  match.table <- match.table[-1,]
  match.table.missing <- match.table.missing[-1,]
  return(list(match.table = match.table, match.table.missing = match.table.missing))
}

