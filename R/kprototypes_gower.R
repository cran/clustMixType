kproto_gower <- function(x, k, lambda = NULL, iter.max = 100, na.rm = "yes", keep.data = TRUE, verbose = TRUE){
  # # enable input of tibbles #...done within kproto()
  # if(is_tibble(x) == TRUE){x <- as.data.frame(x)}
  
  # initial error checks
  if(!is.data.frame(x)) stop("x should be a data frame!")
  if(ncol(x) < 2) stop("For clustering x should contain at least two variables!")
  #  if(iter.max < 1 | nstart < 1) stop("iter.max and nstart must not be specified < 1!") # no parameter within kproto_gower()
  if(!is.null(lambda)){
    if(any(lambda < 0)) stop("lambda must be specified >= 0!")
    if(!any(lambda > 0)) stop("lambda must be specified > 0 for at least two variables!")
    if(length(lambda) != ncol(x)) {
      warning("For gower distance length(lambda) must match the # of variables. lambda will be ignored.")
      lambda <- NULL
      }
    }
  # check for numeric and factor variables
  numvars <- sapply(x, is.numeric)
  anynum <- any(numvars)
  ordvars <- sapply(x, is.ordered)
  anyord <- any(ordvars)
  catvars <- sapply(x, is.factor) & !ordvars
  anyfact <- any(catvars)
  
  # initialize lookup table to store standardized ranks / ranges of numeric variables for predict.kproto()
  lookup <- list() 
  
  # vector of ranges for normalization  
  if(any(numvars)) rgnums <- sapply(x[, numvars, drop = FALSE], function(z) diff(range(z, na.rm = TRUE)))
  if(any(ordvars)){
    xord   <- x[, ordvars, drop = FALSE] # store original variables 
    # ...and replace ordered variables by their ranks
    for(jord in which(ordvars)) x[,jord] <- rank(x[,jord], na.last = "keep")
    rgords <- sapply(x[, ordvars, drop = FALSE], function(z) diff(range(z, na.rm = TRUE)))
  }

  # initialize prototypes
  if(!is.data.frame(k)){
    if (length(k) == 1){
      if(as.integer(k) != k){k <- as.integer(k); warning(paste("k has been set to", k,"!"))}
      if(sum(complete.cases(x)) < k) stop("Data frame has less complete observations than clusters!")
      ids <- sample(row.names(x[complete.cases(x),]), k)
      protos <- x[ids,]
    }
    if (length(k) > 1){
      if(nrow(x) < length(k)) stop("Data frame has less observations than clusters!")
      ids <- k
      k <- length(ids)
      if(length(unique(ids)) != length(ids)) stop("If k is specified as a vector it should contain different indices!")
      if(any(ids<1)|any(ids>nrow(x))) stop("If k is specified as a vector all elements must be valid indices of x!")
      #check for integer
      protos <- x[ids,]
      if(any(!complete.cases(protos))) stop("Choose initial prototypes without missing values!")
    }
    rm(ids)
  }
  if(is.data.frame(k)){
    if(nrow(x) < nrow(k)) stop("Data frame has less observations than clusters!")
    if(length(names(k)) != length(names(x))) stop("k and x have different numbers of columns!")
    if(any(names(k) != names(x))) stop("k and x have different column names!")
    if(anynum) {if( any(sapply(k, is.numeric) != numvars)) stop("Numeric variables of k and x do not match!")}
    if(anyfact) {if( any(sapply(k, is.factor) != catvars)) stop("Factor variables of k and x do not match!")}

    if(anyord){
      # ...and replace ordered variables for initial prototypes by their ranks according to 
      for(jord in 1:ncol(xord)){
        levsj <- unique(xord[,jord])
        rgsj  <- unique(x[,which(ordvars)[jord]])
        names(rgsj) <- levsj
        if(!all(k[,which(ordvars)[jord]] %in% levsj)) stop(paste("Some levels of initial prototypes do not match levels in variable", names(k)[which(ordvars)[jord]],"!")) 
        k[,which(ordvars)[jord]] <- as.numeric(rgsj[k[,which(ordvars)[jord]]])
      }
    }

    protos <- k
    if(any(!complete.cases(protos))) stop("Prototypes with missing values. Choose initial prototypes without missing values!")
    k <- nrow(protos)
  }
  if(k < 1) stop("Number of clusters k must not be smaller than 1!")
  
  if(length(lambda) > 0){
    if(length(lambda) != sum(c(numvars,catvars,ordvars))) {
      warning("For gower distance if lambda is specified, its length must be the sum of numeric and factor variables in the data frame! Otherwise it will be ignored.")
      lambda <- NULL
    }
  }
  
  # initialize clusters
  clusters  <- numeric(nrow(x)) 
  tot.dists <- NULL
  moved   <- NULL
  iter <- 1
  
  # check for any equal prototypes and reduce cluster number in case of occurence
  if(k > 1){
    keep.protos <- rep(TRUE,k)
    for(l in 1:(k-1)){
      for(m in (l+1):k){
        d1 <- d2 <- d3 <- 0
        if(any(numvars)) d1 <- sum((protos[l,numvars, drop = FALSE]-protos[m,numvars, drop = FALSE])^2) # euclidean for numerics
        if(any(catvars)) d2 <- sum(protos[l,catvars, drop = FALSE] != protos[m,catvars, drop = FALSE]) # compare levels 
        if(any(ordvars)) d3 <- sum((protos[l,ordvars, drop = FALSE]-protos[m,ordvars, drop = FALSE])^2) # euclidean for ranks of ordinals
        if((d1+d2+d3) == 0) keep.protos[m] <- FALSE 
      }
    }
    if(!all(keep.protos)){
      protos <- protos[keep.protos,]
      k <- sum(keep.protos)
      if(verbose) message("Equal prototypes merged. Cluster number reduced to:", k, "\n\n")      
    }
  }

  # special case only one cluster
  if(k == 1){clusters <- rep(1, nrow(x)); size  <- table(clusters); iter <- iter.max} # REM: named vector size is needed later...
  
  # start iterations for standard case (i.e. k > 1)
  while(iter < iter.max){
    
    # compute distances 
    nrows <- nrow(x)
    dists <- matrix(NA, nrow=nrows, ncol = k)
    for(i in 1:k){
      
      # in case of no numeric / factor / ordinal variables set:
      d1 <- d2 <- d3 <- rep(0, nrows)
      
      if(any(numvars)){
        d1 <- abs(x[, numvars, drop = FALSE] - matrix(rep(as.numeric(protos[i, numvars, drop = FALSE]), nrows), nrow=nrows, byrow=T))
        for(jnum in 1:ncol(d1)) d1[,jnum] <- d1[,jnum] / rgnums[jnum] 
        d1[is.na(d1)] <- 0
        if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
        if(is.null(lambda)) d1 <- rowSums(d1)
      }
      
      if(any(catvars)){
        d2 <- sapply(which(catvars), function(j) return(x[,j] != rep(protos[i,j], nrows)) )
        d2[is.na(d2)] <- FALSE
        if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
        if(is.null(lambda)) d2 <- rowSums(d2)
      }

      if(any(ordvars)){
        d3 <- abs(x[, ordvars, drop = FALSE] - matrix(rep(as.numeric(protos[i, ordvars, drop = FALSE]), nrows), nrow=nrows, byrow=T))
        for(jord in 1:ncol(d3)) d3[,jord] <- d3[,jord] / rgords[jord] 
        d3[is.na(d3)] <- 0
        if(length(lambda) > 1) d3 <- as.matrix(d3) %*% lambda[ordvars]
        if(is.null(lambda)) d3 <- rowSums(d3)
      }
      
      dists[,i] <- d1 + d2 + d3
    }
    
    # assign clusters 
    old.clusters  <- clusters
    # clusters      <- apply(dists, 1, function(z) which.min(z))
    clusters      <- apply(dists, 1, function(z) {a <- which(z == min(z)); if (length(a)>1) a <- sample(a,1); return(a)}) # sample in case of multiple minima
    size          <- table(clusters)  
    min.dists     <- apply(cbind(clusters, dists), 1, function(z) z[z[1]+1])
    within        <- as.numeric(by(min.dists, clusters, sum))
    tot.within    <- sum(within)
    # prevent from empty classes
    #tot.within    <- numeric(k)
    #totw.list     <- by(min.dists, clusters, sum) 
    #tot.within[names(totw.list)] <- as.numeric(totw.list)
    
    # ...check for empty clusters and eventually reduce number of prototypes    
    if (length(size) < k){
      k <- length(size)
      protos <- protos[1:length(size),]  
      if(verbose) cat("Empty clusters occur. Cluster number reduced to:", k, "\n\n")
    }
    
    # trace
    tot.dists <- c(tot.dists, sum(tot.within))      
    moved <- c(moved, sum(clusters != old.clusters))
   
    # compute new prototypes
    remids <- as.integer(names(size))
    for(i in remids){
      some_vals <- sapply(x[clusters == i, , drop = FALSE], function(z) !all(is.na(z))) # only update variables if not all values are NA
      if(any(some_vals & numvars)){
        protos[which(remids == i), some_vals & numvars] <- sapply(x[clusters == i, some_vals & numvars, drop = FALSE], median, na.rm = TRUE)
      }
      if(any(some_vals & catvars)){
        protos[which(remids == i), some_vals & catvars] <- sapply(x[clusters == i, some_vals & catvars, drop = FALSE], function(z) levels(z)[which.max(table(z))])
      }
      if(any(some_vals & ordvars)){
        protos[which(remids == i), some_vals & ordvars] <- sapply(x[clusters == i, some_vals & ordvars, drop = FALSE], median, na.rm = TRUE)
      }
    }
    
    if(k == 1){clusters <- rep(1, length(clusters)); size <- table(clusters); iter <- iter.max; break}
    
    # check for any equal prototypes and reduce cluster number in case of occurence
    if(iter == (iter.max-1)){ # REM: for last iteration equal prototypes are allowed. otherwise less prototypes than assigned clusters.
      keep.protos <- rep(TRUE,k)
      for(l in 1:(k-1)){
        for(m in (l+1):k){        
          d1 <- d2 <- d3 <- 0
          if(any(numvars)) d1 <- sum((protos[l,numvars, drop = FALSE]-protos[m,numvars, drop = FALSE])^2) # euclidean for numerics
          if(any(catvars)) d2 <- sum(protos[l,catvars, drop = FALSE] != protos[m,catvars, drop = FALSE]) # compare levels for categorics 
          if(any(ordvars)) d3 <- sum((protos[l,ordvars, drop = FALSE]-protos[m,ordvars, drop = FALSE])^2) # euclidean for ranks of ordinals
          if((d1+d2+d3) == 0) keep.protos[m] <- FALSE 
        }
      }
      if(!all(keep.protos)){
        protos <- protos[keep.protos,]
        k <- sum(keep.protos)
        if(verbose) cat("Equal prototypes merged. Cluster number reduced to:", k, "\n\n")      
      }
    }

    # add stopping rules
    if(moved[length(moved)] ==  0) break
    
    if(k == 1){clusters <- rep(1, length(clusters)); size <- table(clusters); iter <- iter.max; break}
    
    #cat("iter", iter, "moved", moved[length(moved)], "tot.dists",tot.dists[length(tot.dists)],"\n" )      
    iter <- iter+1
  }

  ### Final update of prototypes and dists
  if(iter == iter.max){ # otherwise there have been no moves anymore and prototypes correspond to cluster assignments 
    # compute new prototypes
    remids <- as.integer(names(size))
    for(i in remids){
      some_vals <- sapply(x[clusters == i, , drop = FALSE], function(z) !all(is.na(z))) # only update variables if not all values are NA
      if(any(some_vals & numvars)){
        protos[which(remids == i), some_vals & numvars] <- sapply(x[clusters == i, some_vals & numvars, drop = FALSE], median, na.rm = TRUE)
      }
      if(any(some_vals & catvars)){
        protos[which(remids == i), some_vals & catvars] <- sapply(x[clusters == i, some_vals & catvars, drop = FALSE], function(z) levels(z)[which.max(table(z))])
      }
      if(any(some_vals & ordvars)){
        protos[which(remids == i), some_vals & ordvars] <- sapply(x[clusters == i, some_vals & ordvars, drop = FALSE], median, na.rm = TRUE)
      }
    }
    
    # compute distances 
    nrows <- nrow(x)
    dists <- matrix(NA, nrow=nrows, ncol = k)
    for(i in 1:k){
      
      # in case of no numeric / factor / ordinal variables set:
      d1 <- d2 <- d3 <- rep(0, nrows)
      
      if(any(numvars)){
        d1 <- abs(x[, numvars, drop = FALSE] - matrix(rep(as.numeric(protos[i, numvars, drop = FALSE]), nrows), nrow=nrows, byrow=T))
        for(jnum in 1:ncol(d1)) d1[,jnum] <- d1[,jnum] / rgnums[jnum] 
        d1[is.na(d1)] <- 0
        if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
        if(is.null(lambda)) d1 <- rowSums(d1)
      }
      
      if(any(catvars)){
        d2 <- sapply(which(catvars), function(j) return(x[,j] != rep(protos[i,j], nrows)) )
        d2[is.na(d2)] <- FALSE
        if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
        if(is.null(lambda)) d2 <- rowSums(d2)
      }
      
      if(any(ordvars)){
        d3 <- abs(x[, ordvars, drop = FALSE] - matrix(rep(as.numeric(protos[i, ordvars, drop = FALSE]), nrows), nrow=nrows, byrow=T))
        for(jord in 1:ncol(d3)) d3[,jord] <- d3[,jord] / rgords[jord] 
        d3[is.na(d3)] <- 0
        if(length(lambda) > 1) d3 <- as.matrix(d3) %*% lambda[ordvars]
        if(is.null(lambda)) d3 <- rowSums(d3)
      }
      
      dists[,i] <- d1 + d2 + d3
    }

        
    size          <- table(clusters)  
    min.dists     <- apply(cbind(clusters, dists), 1, function(z) z[z[1]+1])
    within        <- as.numeric(by(min.dists, clusters, sum))
    tot.within    <- sum(within)
  }
  
  # replace ranks of ordinal variables from prototypes by the level with the closest rank
  if(any(ordvars)){
    protos.ord <- protos[,ordvars, drop = FALSE]
    for(jord in 1:ncol(protos.ord)){
      protos[,which(ordvars)[jord]] <- sapply(protos.ord[,jord], function(z) xord[which.min(abs(rank(xord[,jord]) - z)), jord])
    }
    for(jord in which(ordvars)) protos[,jord] <- factor(protos[,jord], levels = levels(protos[,jord]), ordered = TRUE)
    
    # store standardizied ranks in lookup table for predict.kproto()
    for(jord in 1:ncol(xord)){
      df <- data.frame(unique(xord[,jord]),
                       unique(x[,which(ordvars)[jord], drop = FALSE] / rgords[jord])
                       )
      names(df) <- c("level", "value")
      lookup[[names(xord)[jord]]] <- df
    }  
  }
  # store ranges of numeric variables for predict.kproto
  if(any(numvars)) lookup$num_ranges <- rgnums

  if(na.rm == "no"){
    allNAs <- apply(x,1,function(z) all(is.na(z))) # note: not passed from kproto!  
    if(sum(allNAs) > 0){
      clusters[allNAs] <- NA
      dists[allNAs,] <- NA
    }
  }

  names(clusters) <- row.names(dists) <- row.names(x)
  rownames(protos) <- NULL
  # create result: 
  res <- list(cluster = clusters,  
              centers = protos, 
              lambda = lambda, 
              size = size,
              withinss = within,
              tot.withinss = tot.within,   
              dists = dists, 
              iter = iter,
              stdization = lookup,
              trace = list(tot.dists = tot.dists, moved = moved))
  
  # if(keep.data) res$data = x  # ...done by kproto() 
  # class(res) <- "kproto"
  return(res)
}
