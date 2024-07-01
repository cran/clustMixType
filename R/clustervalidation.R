


calc.dist <- function(type, lambda, data1, data2 = NULL, all.dists = FALSE, data_plus = NULL){
  
  if(type == "huang"){# calculation with standard distance of k-prototypes
    
    # if(is.null(lambda)){
    #   lambda <- lambdaest(data, num.method = 1, fac.method = 1, outtype = "numeric", verbose = TRUE)
    # }
    
    if(!is.null(data2)){#contains the input two datasets?
      
      #tbd: check if structure of data1 and data2 are equal
      data <- rbind(data1, data2)
      
      if(nrow(data1) == 1 & nrow(data2) == 1){#calculation the distance between two observations (due to saving time)
        
        #numerical variables
        d1 <- (data[1, data_plus$numvars, drop = FALSE]-data[2, data_plus$numvars, drop = FALSE])^2
        d1[is.na(d1)] <- 0
        if(length(lambda) == 1){
          d1 <- sum(d1)
        }else{ #length(lambda) > 1
          d1 <- as.matrix(d1, nrow = 1) %*% matrix(lambda[data_plus$numvars], ncol = 1)
        }
        
        #categorical variables
        d2 <- sapply(which(data_plus$catvars), function(j) return(data[1,j] != data[2,j]))
        d2[is.na(d2)] <- FALSE
        if(length(lambda) == 1){
          d2 <- lambda * sum(d2)
        }else{#length(lambda) > 1
          d2 <- matrix(d2, nrow = 1) %*% matrix(lambda[data_plus$catvars], ncol = 1)
        }
        
        return(d1 + d2)
      }else{#calculate distances between objects of two datasets
        
        #numerical variables
        d1 <- (data1[rep(1:nrow(data1), each=nrow(data2)), data_plus$numvars, drop = FALSE] -
                 data2[rep(1:nrow(data2), times=nrow(data1)), data_plus$numvars, drop = FALSE])^2
        d1[is.na(d1)] <- 0
        if(length(lambda) == 1){
          d1 <- rowSums(d1)
        }else{ #length(lambda) > 1
          d1 <- as.matrix(d1) %*% lambda[data_plus$numvars]
        }
        
        #categorical variables
        d2 <- data1[rep(1:nrow(data1), each=nrow(data2)), data_plus$catvars, drop = FALSE] != 
          data2[rep(1:nrow(data2), times = nrow(data1)), data_plus$catvars, drop = FALSE]
        d2[is.na(d2)] <- FALSE
        if(length(lambda) == 1){
          d2 <- lambda * rowSums(d2)
        }else{ #length(lambda) > 1
          d2 <- as.matrix(d2) %*% lambda[data_plus$catvars]
        }
        
        return(sum(d1 + d2))
      }
    }else{ #calculate distances between the objects of one dataset
      
      d1 <- (data1[rep(1:(nrow(data1)-1),times=(nrow(data1)-1):1, each=1), data_plus$numvars, drop = FALSE] -
               data1[unlist(lapply(2:nrow(data1), seq, to=nrow(data1))), data_plus$numvars, drop = FALSE])^2
      d1[is.na(d1)] <- 0
      if(length(lambda) == 1) d1 <- rowSums(d1)
      if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[data_plus$numvars]
      
      d2 <- data1[rep(1:(nrow(data1)-1),times=(nrow(data1)-1):1, each=1), data_plus$catvars, drop = FALSE] != 
        data1[unlist(lapply(2:nrow(data1), seq, to=nrow(data1))), data_plus$catvars, drop = FALSE]
      d2[is.na(d2)] <- FALSE
      if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
      if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[data_plus$catvars]
      
      if(all.dists == TRUE){
        m <- matrix(rep(0,nrow(data1)^2), ncol = nrow(data1))
        
        m[lower.tri(m)] <- d1+d2
        m <- t(m)
        m[lower.tri(m)] <- d1+d2
        return(m)
      }else{
        return(sum(d1 + d2))
      }
    }
  }else{ #calculation with gower distance
    
    if(!is.null(data2)){#contains the input two datasets?
      
      #tbd: check if structure of data1 and data2 are equal
      #data <- rbind(data1, data2)
      
      if(nrow(data1) == 1 & nrow(data2) == 1){#calculation the distance between two observations (due to saving time)
        
        d1 <- d2 <- d3 <- 0
        
        if(any(data_plus$numvars)){
          d1 <- abs(data1[1, data_plus$numvars, drop = FALSE] - data2[1, data_plus$numvars, drop = FALSE])
          for(jnum in 1:ncol(d1)) d1[,jnum] <- d1[,jnum] / data_plus$rgnums[jnum] 
          d1[is.na(d1)] <- 0
          if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[data_plus$numvars]
          if(is.null(lambda)) d1 <- rowSums(d1)
        }
        
        if(any(data_plus$catvars)){
          d2 <- data1[1, data_plus$catvars, drop = FALSE] != data2[1, data_plus$catvars, drop = FALSE]
          d2[is.na(d2)] <- FALSE
          if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[data_plus$catvars]
          if(is.null(lambda)) d2 <- rowSums(d2)
        }
        
        if(any(data_plus$ordvars)){
          d3 <- abs(data1[1, data_plus$ordvars, drop = FALSE] - data2[1, data_plus$ordvars, drop = FALSE])
          for(jord in 1:ncol(d3)) d3[,jord] <- d3[,jord] / data_plus$rgords[jord] 
          d3[is.na(d3)] <- 0
          if(length(lambda) > 1) d3 <- as.matrix(d3) %*% lambda[data_plus$ordvars]
          if(is.null(lambda)) d3 <- rowSums(d3)
        }
        
        return(d1 + d2 + d3)
        
      }else{#calculate distances between objects of two datasets
        
        # in case of no numeric / factor / ordinal variables set:
        d1 <- d2 <- d3 <- rep(0, nrow(data1)*nrow(data2))
        
        if(any(data_plus$numvars)){
          d1 <- abs(data1[rep(1:nrow(data1), each=nrow(data2)), data_plus$numvars, drop = FALSE] - 
                      data2[rep(1:nrow(data2), times=nrow(data1)), data_plus$numvars, drop = FALSE])
          for(jnum in 1:ncol(d1)) d1[,jnum] <- d1[,jnum] / data_plus$rgnums[jnum] 
          #time efficient alternative for codeline before:
          # d1 <- d1 / rep(data_plus$rgnums, rep(nrow(d1), ncol(d1)))
          d1[is.na(d1)] <- 0
          if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[data_plus$numvars]
          if(is.null(lambda)) d1 <- rowSums(d1)
        }
        
        if(any(data_plus$catvars)){
          d2 <- data1[rep(1:nrow(data1), each=nrow(data2)), data_plus$catvars, drop = FALSE] != 
            data2[rep(1:nrow(data2), times = nrow(data1)), data_plus$catvars, drop = FALSE]
          d2[is.na(d2)] <- FALSE
          if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[data_plus$catvars]
          if(is.null(lambda)) d2 <- rowSums(d2)
        }
        
        if(any(data_plus$ordvars)){
          d3 <- abs(data1[rep(1:nrow(data1), each=nrow(data2)), data_plus$ordvars, drop = FALSE] - 
                      data2[rep(1:nrow(data2), times=nrow(data1)), data_plus$ordvars, drop = FALSE])
          for(jord in 1:ncol(d3)) d3[,jord] <- d3[,jord] / data_plus$rgords[jord] 
          #time efficient alternative for codeline before:
          # d3 <- d3 / rep(data_plus$rgords, rep(nrow(d3), ncol(d3)))
          d3[is.na(d3)] <- 0
          if(length(lambda) > 1) d3 <- as.matrix(d3) %*% lambda[data_plus$ordvars]
          if(is.null(lambda)) d3 <- rowSums(d3)
        }
        
        return(sum(d1 + d2 + d3))
      }
    }else{ #calculate distances between the objects of one dataset
      
      d1 <- d2 <- d3 <- rep(0,length(rep(1:(nrow(data1)-1),times=(nrow(data1)-1):1, each=1)))
      
      if(any(data_plus$numvars)){
        d1 <- abs(data1[rep(1:(nrow(data1)-1),times=(nrow(data1)-1):1, each=1), data_plus$numvars, drop = FALSE] - 
                    data1[unlist(lapply(2:nrow(data1), seq, to=nrow(data1))), data_plus$numvars, drop = FALSE])
        for(jnum in 1:ncol(d1)) d1[,jnum] <- d1[,jnum] / data_plus$rgnums[jnum] 
        #time efficient alternative for codeline before:
        # d1 <- d1 / rep(rgnums, rep(nrow(d1), ncol(d1)))
        d1[is.na(d1)] <- 0
        if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[data_plus$numvars]
        if(is.null(lambda)) d1 <- rowSums(d1)
      }
      
      if(any(data_plus$catvars)){
        d2 <- data1[rep(1:(nrow(data1)-1),times=(nrow(data1)-1):1, each=1), data_plus$catvars, drop = FALSE] != 
          data1[unlist(lapply(2:nrow(data1), seq, to=nrow(data1))), data_plus$catvars, drop = FALSE]
        d2[is.na(d2)] <- FALSE
        if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[data_plus$catvars]
        if(is.null(lambda)) d2 <- rowSums(d2)
      }
      
      if(any(data_plus$ordvars)){
        d3 <- abs(data1[rep(1:(nrow(data1)-1),times=(nrow(data1)-1):1, each=1), data_plus$ordvars, drop = FALSE] - 
                    data1[unlist(lapply(2:nrow(data1), seq, to=nrow(data1))), data_plus$ordvars, drop = FALSE])
        for(jord in 1:ncol(d3)) d3[,jord] <- d3[,jord] / data_plus$rgords[jord] 
        d3[is.na(d3)] <- 0
        if(length(lambda) > 1) d3 <- as.matrix(d3) %*% lambda[data_plus$ordvars]
        if(is.null(lambda)) d3 <- rowSums(d3)
      }
      
      if(all.dists == TRUE){ #return distance matrix
        m <- matrix(rep(0,nrow(data1)^2), ncol = nrow(data1))
        
        m[lower.tri(m)] <- d1+d2+d3
        m <- t(m)
        m[lower.tri(m)] <- d1+d2+d3
        return(m)
      }else{ #return the sum of all dists between objects
        return(sum(d1 + d2 + d3))
      }
      
    }
  }
}


create.S_w <- function(object, data_plus){
  S_w <- 0
  for(k in seq_along(object$size)){
    x_k <- object$data[which(object$cluster==k),]
    if(object$size[k]>1){
      S_w <- S_w + calc.dist(lambda = object$lambda, type = object$type, 
                             data1 = x_k, 
                             data_plus = data_plus)
    }
  }
  return(S_w)
}


create.S_b <- function(object, data_plus){
  S_b <- 0
  for(k in 1:(length(object$size)-1)){
    x_k <- object$data[which(object$cluster==k),]
    for(l in (k+1):length(object$size)){
      S_b <- S_b + calc.dist(type = object$type, lambda = object$lambda,
                             data1 = x_k, 
                             data2 = object$data[which(object$cluster==l),], 
                             data_plus = data_plus)
    }
  }
  return(S_b)
}


create.N_w <- function(object){
  sum((object$size*(object$size-1))/2)
}




cindex_kproto <- function(object = NULL, type = NULL, data = NULL, k = NULL, S_sort = NULL, kp_obj = "optimal", lambda = NULL, data_plus = NULL, verbose = FALSE, ...){
  
  if(!is.null(object)){# it follows: determine the index value of the given cluster partition
    
    if(is.null(S_sort)){
      n <- length(object$cluster)
      S_all <- matrix(numeric(), ncol = n, nrow = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          S_all[i,j] <- calc.dist(type = object$type, lambda = object$lambda,
                                  data1 = object$data[i,], data2 = object$data[j,],
                                  data_plus = data_plus)
          
        }
      }
      S_sort <- sort(S_all)
    }
    
    S_w <- create.S_w(object, data_plus = data_plus)
    
    N_w <- create.N_w(object)
    S_min <- sum(head(S_sort,n = N_w))
    S_max <- sum(tail(S_sort,n = N_w))
    
    if(S_min != S_max){
      index <- (S_w - S_min)/(S_max - S_min)
    }
    
    return(index)
  }else{# determine index-optimal cluster partition
    
    n <- nrow(data)
    S_all <- matrix(numeric(), ncol=n, nrow=n)
    for(i in 1:(n - 1)){
      for(j in (i + 1):n){
        S_all[i,j] <- calc.dist(type = type, lambda = lambda, 
                                data1 = data[i,], data2 = data[j,],
                                data_plus = data_plus)
      }
    }
    S_sort <- sort(S_all)
    
    #calculate all kproto objects for k
    if(type == "gower"){
      # ...replace ranks with original ordered values
      data[,which(data_plus$ordvars)] <- data_plus$data_ord
    }
    object <- kproto(x = data, k = k[1], type = type, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
    if(type == "gower"){
      # ...replace ranks with original ordered values
      for(jord in which(data_plus$ordvars)) object$data[,jord] <- rank(data[,jord])
    }
    trace_kp <- list(list("index" = cindex_kproto(object = object, S_sort = S_sort, data_plus = data_plus), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      if(type == "gower"){
        # ...replace ranks with original ordered values
        data[,which(data_plus$ordvars)] <- data_plus$data_ord
      }
      object <- kproto(x = data, k = q, type = type, keep.data = TRUE, lambda = lambda, verbose = verbose)#, ...)
      if(type == "gower"){
        # ...replace ranks with original ordered values
        for(jord in which(data_plus$ordvars)) object$data[,jord] <- rank(data[,jord])
      }
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = cindex_kproto(object = object, S_sort = S_sort, data_plus = data_plus), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a cluster partition with same number of cluster but different validation index
        index_value <- cindex_kproto(object = object, S_sort = S_sort, data_plus = data_plus)
        if(!(index_value %in% unlist(lapply(trace_kp, `[[`, 1))[which(unlist(lapply(trace_kp, `[[`, 2)) == length(object$size))])){
          trace_kp <- c(trace_kp, list(list("index" = index_value, "k" = length(object$size), "object" = object)))
        }
      }
    }
    
    # save all index values for comparison and as summary for output
    indices <- unlist(lapply(trace_kp, `[[`, 1))
    names(indices) <- lapply(trace_kp, `[[`, 2)
    
    # find the optimal k, if it is ambiguously: sample
    k_m <- which(indices == min(indices))
    if(length(k_m)>1){k_m <- sample(k_m,1)}
    k_opt <- as.integer(names(indices[k_m]))
    index_opt <- indices[k_m]
    
    output <- switch(kp_obj,
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp, "type" = type),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object, "type" = type))
    return(output)
  }
}


dunn_kproto <- function(object = NULL, type = NULL, data = NULL, k = NULL, kp_obj = "optimal", lambda = NULL, data_plus = NULL, verbose = FALSE, ...){
  
  cond <- FALSE
  if(!is.null(object)){
    if(!is.null(object$lambda)){
      cond <- all(!as.logical(data_plus$numvars*object$lambda))
    }
  }else{
    if(!is.null(lambda)){
      cond <- all(!as.logical(data_plus$numvars*lambda))
    }
  }
  if(cond) message("As a result of the choice of lambda: No numeric variables in x! Index calculation might result in NA...\n")
  
  if(!is.null(object)){
    
    k <- length(object$size)
    
    #determine d(C_i,C_j)
    min_CiCj <- matrix(numeric(k*k), ncol = k, nrow = k)
    for(i in 1:(k-1)){
      xi <- object$data[which(object$cluster == i),]
      for(j in (i+1):k){
        xj <- object$data[which(object$cluster == j),]
        min_ij <- calc.dist(type = object$type, lambda = object$lambda, 
                            data1 = xi[1,], data2 = xj[1,],
                            data_plus = data_plus)
        for(l in 1:object$size[i]){
          for(m in 1:object$size[j]){
            min_neu <- calc.dist(type = object$type, lambda = object$lambda, 
                                 data1 = xi[l,], data2 = xj[m,],
                                 data_plus = data_plus)
            if(min_neu < min_ij){
              min_ij <- min_neu
            }
          }
        }
        min_CiCj[i,j] <- min_ij
      }
    }
    
    if(length(min_CiCj[min_CiCj > 0]) > 0){
      numerator <- min(min_CiCj[min_CiCj > 0])
    }else{
      return(NA)
    }
    
    
    #determine diam(C_k)
    max_diam <- numeric(k)
    for(p in 1:k){
      xi <- object$data[which(object$cluster == p),]
      if(object$size[p] > 1){
        max_ij <- calc.dist(type = object$type, lambda = object$lambda, 
                            data1 = xi[1,], data2 = xi[2,],
                            data_plus = data_plus)
        for(l in 1:(object$size[p]-1)){
          for(m in (l+1):(object$size[p])){
            max_neu <- calc.dist(type = object$type, lambda = object$lambda, 
                                 data1 = xi[l,], data2 = xi[m,],
                                 data_plus = data_plus)
            if(max_neu > max_ij){
              max_ij <- max_neu
            }
          }
        }
        max_diam[p] <- max_ij
      }
    }
    denominator <- max(max_diam)
    
    if(is.finite(numerator/denominator)){
      return(numerator/denominator)
    }else{
      return(NA)
    }
    
    
  }else{
    n <- nrow(data)
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], type = type, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
    trace_kp <- list(list("index" = dunn_kproto(object = object, data_plus = data_plus), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, type = type, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = dunn_kproto(object = object, data_plus = data_plus), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a clusterpartition with same number of cluster but different validation index
        index_value <- dunn_kproto(object = object, data_plus = data_plus)
        if(!(index_value %in% unlist(lapply(trace_kp, `[[`, 1))[which(unlist(lapply(trace_kp, `[[`, 2)) == length(object$size))])){
          trace_kp <- c(trace_kp, list(list("index" = index_value, "k" = length(object$size), "object" = object)))
        }
      }
    }
    
    # save all index values for comparison and as summary for output
    indices <- unlist(lapply(trace_kp, `[[`, 1))
    names(indices) <- lapply(trace_kp, `[[`, 2)
    
    # find the optimal k, if it is ambiguously: sample
    if(all(is.na(indices))){return(NA)} # returning NA if dunn index couldn't be calculated (for all k), otherwise choose partition with highest index value
    k_m <- which(indices == max(indices, na.rm = TRUE))
    if(length(k_m)>1){k_m <- sample(k_m, 1)}
    k_opt <- as.integer(names(indices[k_m]))
    index_opt <- indices[k_m]
    
    output <- switch(kp_obj,
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp, "type" = type),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object, "type" = type))
    return(output)
    
  }
}


gamma_kproto <- function(object = NULL, type = NULL, data = NULL, k = NULL, dists = NULL, kp_obj = "optimal", lambda = NULL, data_plus = NULL, verbose = FALSE, ...){
  
  if(!is.null(object)){
    
    if(is.null(dists)){
      n <- nrow(object$data)
      dists <- matrix(numeric(), nrow = n, ncol = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- calc.dist(type = object$type, lambda = object$lambda, 
                                  data1 = object$data[i,], data2 = object$data[j,],
                                  data_plus = data_plus)
        }
      }
    }
    
    s_plus <- 0
    s_minus <- 0
    dist_within <- 0
    dist_between <- 0
    for(k in as.numeric(names(object$size))){
      dist_within <- c(dist_within, na.omit(as.vector(dists[as.numeric(which(object$cluster == k)), as.numeric(which(object$cluster == k))])))
      dist_between <- c(dist_between, na.omit(as.vector(dists[as.numeric(which(object$cluster == k)), as.numeric(which(object$cluster != k))])))
    }
    dist_within <- dist_within[-1]
    for(m in 1:length(dist_within)){
      s_plus <- s_plus + sum(dist_within[m] < dist_between)
      s_minus <- s_minus + sum(dist_within[m] > dist_between)
    }
    index <- (s_plus - s_minus)/(s_plus + s_minus)
    
    return(index)
  }else{
    n <- nrow(data)
    
    if(is.null(dists)){
      dists <- matrix(numeric(), nrow = n, ncol = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- calc.dist(type = type, lambda = lambda, 
                                  data1 = data[i,], data2 = data[j,],
                                  data_plus = data_plus)
        }
      }
    }
    
    index <- numeric(n)
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], type = type, keep.data = TRUE,  lambda = lambda, verbose = verbose, ...)
    trace_kp <- list(list("index" = gamma_kproto(object = object, dists = dists, data_plus = data_plus), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, type = type, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = gamma_kproto(object = object, dists = dists, data_plus = data_plus), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a clusterpartition with same number of cluster but different validation index
        index_value <- gamma_kproto(object = object, dists = dists, data_plus = data_plus)
        if(!(index_value %in% unlist(lapply(trace_kp, `[[`, 1))[which(unlist(lapply(trace_kp, `[[`, 2)) == length(object$size))])){
          trace_kp <- c(trace_kp, list(list("index" = index_value, "k" = length(object$size), "object" = object)))
        }
      }
    }
    
    # save all index values for comparison and as summary for output
    indices <- unlist(lapply(trace_kp, `[[`, 1))
    names(indices) <- lapply(trace_kp, `[[`, 2)
    
    # find the optimal k, if it is ambiguously: sample
    k_m <- which(indices == max(indices))
    if(length(k_m)>1){k_m <- sample(k_m,1)}
    k_opt <- as.integer(names(indices[k_m]))
    index_opt <- indices[k_m]
    
    output <- switch(kp_obj,
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp, "type" = type),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object, "type" = type))
    return(output)
    
  }
}


gplus_kproto <- function(object = NULL, type = NULL, data = NULL, k = NULL, dists = NULL, kp_obj = "optimal", lambda = NULL, data_plus = NULL, verbose = FALSE, ...){
  
  if(!is.null(object)){
    
    n <- nrow(object$data)
    if(is.null(dists)){
      dists <- matrix(numeric(), nrow=n, ncol=n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- calc.dist(type = object$type, lambda = object$lambda, 
                                  data1 = object$data[i,], data2 = object$data[j,],
                                  data_plus = data_plus)
        }
      }
    }
    
    s_minus <- 0
    dist_within <- 0
    dist_between <- 0
    for(k in as.numeric(names(object$size))){
      dist_within <- c(dist_within, na.omit(as.vector(dists[as.numeric(which(object$cluster == k)),as.numeric(which(object$cluster == k))])))
      dist_between <- c(dist_between, na.omit(as.vector(dists[as.numeric(which(object$cluster == k)),as.numeric(which(object$cluster != k))])))
    }
    dist_within <- dist_within[-1]
    for(m in 1:length(dist_within)){
      s_minus <- s_minus + sum(dist_within[m] > dist_between)
    }
    N_t <- n*(n - 1)/2
    
    index <- (2 * s_minus)/(N_t * (N_t - 1))
    
    return(index)
  }else{
    n <- nrow(data)
    
    if(is.null(dists)){
      dists <- matrix(numeric(), nrow = n, ncol = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- calc.dist(type = type, lambda = lambda, 
                                  data1 = data[i,], data2 = data[j,],
                                  data_plus = data_plus)
        }
      }
    }
    
    if(is.null(k)){k <- 2:sqrt(n)}
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], type = type, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
    trace_kp <- list(list("index" = gplus_kproto(object = object, dists = dists, data_plus = data_plus), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, type = type, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = gplus_kproto(object = object, dists = dists, data_plus = data_plus), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a clusterpartition with same number of cluster but different validation index
        index_value <- gplus_kproto(object = object, dists = dists, data_plus = data_plus)
        if(!(index_value %in% unlist(lapply(trace_kp, `[[`, 1))[which(unlist(lapply(trace_kp, `[[`, 2)) == length(object$size))])){
          trace_kp <- c(trace_kp, list(list("index" = index_value, "k" = length(object$size), "object" = object)))
        }
      }
    }
    
    # save all index values for comparison and as summary for output
    indices <- unlist(lapply(trace_kp, `[[`, 1))
    names(indices) <- lapply(trace_kp, `[[`, 2)
    
    # find the optimal k, if it is ambiguously: sample
    k_m <- which(indices == min(indices))
    if(length(k_m)>1){k_m <- sample(k_m,1)}
    k_opt <- as.integer(names(indices[k_m]))
    index_opt <- indices[k_m]
    
    output <- switch(kp_obj,
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp, "type" = type),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object, "type" = type))
    return(output)
  }
}


mcclain_kproto <- function(object = NULL, type = NULL, data = NULL, k = NULL, kp_obj = "optimal", lambda = NULL, data_plus = NULL, verbose = FALSE, ...){
  
  if(!is.null(object)){
    n <- nrow(object$data)
    S_w_kpa <- create.S_w(object, data_plus = data_plus)
    S_b_kpa <- create.S_b(object, data_plus = data_plus)
    N_w <- create.N_w(object)
    N_t <- n*(n - 1)/2
    N_b <- N_t-N_w
    
    index <- (S_w_kpa/N_w)/(S_b_kpa/N_b)
    
    return(index)
  }else{
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], type = type, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
    trace_kp <- list(list("index" = mcclain_kproto(object = object, data_plus = data_plus), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, type = type, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = mcclain_kproto(object = object, data_plus = data_plus), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a cluster partition with same number of cluster but different validation index
        index_value <- mcclain_kproto(object = object, data_plus = data_plus)
        if(!(index_value %in% unlist(lapply(trace_kp, `[[`, 1))[which(unlist(lapply(trace_kp, `[[`, 2)) == length(object$size))])){
          trace_kp <- c(trace_kp, list(list("index" = index_value, "k" = length(object$size), "object" = object)))
        }
      }
    }
    
    # save all index values for comparison and as summary for output
    indices <- unlist(lapply(trace_kp, `[[`, 1))
    names(indices) <- lapply(trace_kp, `[[`, 2)
    
    # find the optimal k, if it is ambiguously: sample
    k_m <- which(indices == min(indices))
    if(length(k_m)>1){k_m <- sample(k_m,1)}
    k_opt <- as.integer(names(indices[k_m]))
    index_opt <- indices[k_m]
    
    output <- switch(kp_obj,
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp, "type" = type),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object, "type" = type))
    return(output)
  }
}


ptbiserial_kproto <- function(object = NULL, type = NULL, data = NULL, k = NULL, s_d = NULL, kp_obj = "optimal", lambda = NULL, data_plus = NULL, verbose = FALSE, ...){
  
  if(!is.null(object)){
    
    n <- nrow(object$data)
    if(is.null(s_d)){
      S_all <- matrix(numeric(), ncol = n, nrow = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          S_all[i,j] <- calc.dist(type = object$type, lambda = object$lambda, 
                                  data1 = object$data[i,], data2 = object$data[j,],
                                  data_plus = data_plus)
        }
      }
      s_d <- sd(S_all, na.rm = TRUE)
    }
    
    S_w_kpa <- create.S_w(object, data_plus = data_plus)
    S_b_kpa <- create.S_b(object, data_plus = data_plus)
    N_w <- create.N_w(object)
    N_t <- n*(n-1)/2
    N_b <- N_t-N_w
    
    index <- ((S_b_kpa/N_b) - (S_w_kpa/N_w)) * sqrt(N_w * N_b/N_t^2) / s_d
    
    return(index)
  }else{
    
    n <- nrow(data)
    if(is.null(s_d)){
      S_all <- matrix(numeric(), ncol = n, nrow = n)
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          S_all[i,j] <- calc.dist(type = type, lambda = lambda,
                                  data1 = data[i,], data2 = data[j,],
                                  data_plus = data_plus)
        }
      }
      s_d <- sd(S_all, na.rm = TRUE)
    }
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], type = type, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
    trace_kp <- list(list("index" = ptbiserial_kproto(object = object, s_d = s_d, data_plus = data_plus), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, type = type, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = ptbiserial_kproto(object = object, s_d = s_d, data_plus = data_plus), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a cluster partition with same number of cluster but different validation index
        index_value <- ptbiserial_kproto(object = object, s_d = s_d, data_plus = data_plus)
        if(!(index_value %in% unlist(lapply(trace_kp, `[[`, 1))[which(unlist(lapply(trace_kp, `[[`, 2)) == length(object$size))])){
          trace_kp <- c(trace_kp, list(list("index" = index_value, "k" = length(object$size), "object" = object)))
        }
      }
    }
    
    # save all index values for comparison and as summary for output
    indices <- unlist(lapply(trace_kp, `[[`, 1))
    names(indices) <- lapply(trace_kp, `[[`, 2)
    
    # find the optimal k, if it is ambiguously: sample
    k_m <- which(indices == max(indices))
    if(length(k_m)>1){k_m <- sample(k_m,1)}
    k_opt <- as.integer(names(indices[k_m]))
    index_opt <- indices[k_m]
    
    output <- switch(kp_obj,
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp, "type" = type),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object, "type" = type))
    return(output)
  }
}


silhouette_kproto <- function(object = NULL, type = NULL, data = NULL, k = NULL, kp_obj = "optimal", lambda = NULL, data_plus = NULL, verbose = FALSE, ...){
  
  if(!is.null(object)){
    
    dists <- calc.dist(type = object$type, lambda = object$lambda, 
                       data1 = object$data, all.dists = TRUE,
                       data_plus = data_plus)
    
    cluster_dists <- matrix(numeric(nrow(object$data)*length(table(object$cluster))), 
                            nrow = nrow(object$data), ncol = length(table(object$cluster)))
    for(i in 1:length(table(object$cluster))){
      if(!(length(which(object$cluster == i)) == 1)){
        cluster_dists[,i] <- rowMeans(dists[,which(object$cluster == i)])
      }else{
        cluster_dists[,i] <- dists[,which(object$cluster == i)]
      }
    }
    
    #determine ai, bi and si
    a <- numeric(nrow(object$data))
    b <- numeric(nrow(object$data))
    s <- numeric(nrow(object$data))
    for(i in 1:nrow(object$data)){
      if(is.na(object$cluster[i])){
        # special case: no cluster assignment cluster[i] (usually resulting of all variables NA in x[i])
        s[i] <- NA
      }else{
        a[i] <- cluster_dists[i, object$cluster[i]]
        b[i] <- min(cluster_dists[i, -object$cluster[i]])
        
        if(max(a[i], b[i], na.rm = TRUE) == 0){
          # special case: a[i]=0 and b[i]=0 => s = 0, since x[i] lies equally far away (distance = 0) from both the clusters
          s[i] <- 0
        }else{
          s[i] <- (b[i] - a[i])/max(a[i],b[i], na.rm = TRUE)
        }
      }
      
    }
    if(any(table(object$cluster) == 1)){
      for(i in which(object$cluster %in% as.integer(which(table(object$cluster) == 1)))){
        s[i] <- 0
      }
      cat(length(which(object$cluster %in% as.integer(which(table(object$cluster) == 1))))," cluster with only one observation\n")
    }
    
    index <- mean(s, na.rm = TRUE)
    
    return(index)
  }else{
    #n <- nrow(data)
    
    #calculate all kproto objects for k
    object <- kproto(x = data, type = type, k = k[1], keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
    trace_kp <- list(list("index" = silhouette_kproto(object = object, data_plus = data_plus), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, type = type, k = q, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = silhouette_kproto(object = object, data_plus = data_plus), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a cluster partition with same number of cluster but different validation index
        index_value <- silhouette_kproto(object = object, data_plus = data_plus)
        if(!(index_value %in% unlist(lapply(trace_kp, `[[`, 1))[which(unlist(lapply(trace_kp, `[[`, 2)) == length(object$size))])){
          trace_kp <- c(trace_kp, list(list("index" = index_value, "k" = length(object$size), "object" = object)))
        }
      }
    }
    
    # save all index values for comparison and as summary for output
    indices <- unlist(lapply(trace_kp, `[[`, 1))
    names(indices) <- lapply(trace_kp, `[[`, 2)
    
    # find the optimal k, if it is ambiguously: sample
    k_m <- which(indices == max(indices, na.rm = TRUE))
    if(length(k_m)>1){k_m <- sample(k_m,1)}
    k_opt <- as.integer(names(indices[k_m]))
    index_opt <- indices[k_m]
    
    output <- switch(kp_obj,
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp, "type" = type),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object, "type" = type))
    return(output)
  }
}


tau_kproto <- function(object = NULL, type = NULL, data = NULL, k = NULL, dists = NULL, kp_obj = "optimal", lambda = NULL, data_plus = NULL, verbose = FALSE, ...){
  
  if(!is.null(object)){
    
    if(is.null(dists)){
      n <- nrow(object$data)
      dists <- matrix(numeric(), nrow = n, ncol = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- calc.dist(type = object$type, lambda = object$lambda, 
                                  data1 = object$data[i,], data2 = object$data[j,],
                                  data_plus = data_plus)
        }
      }
    }
    
    s_plus <- 0
    s_minus <- 0
    dist_within <- 0
    dist_between <- 0
    for(k in as.numeric(names(object$size))){
      dist_within <- c(dist_within, na.omit(as.vector(dists[as.numeric(which(object$cluster == k)),as.numeric(which(object$cluster == k))])))
      dist_between <- c(dist_between, na.omit(as.vector(dists[as.numeric(which(object$cluster == k)),as.numeric(which(object$cluster != k))])))
    }
    dist_within <- dist_within[-1]
    for(m in 1:length(dist_within)){
      s_plus <- s_plus + sum(dist_within[m] < dist_between)
      s_minus <- s_minus + sum(dist_within[m] > dist_between)
    }
    
    n <- nrow(object$data)
    N_t <- n * (n - 1)/2
    M <- combn(1:n, m = 2)
    M <- rbind(M, apply(X = M, MARGIN = 2, function(x){
      if(any(is.na(c(object$cluster[x[1]], object$cluster[x[2]])))){return(-1)}else{
        if(object$cluster[x[1]] == object$cluster[x[2]]){return(1)}else{return(-1)}
      }
    }))
    t <- 0
    for(i in 1:(ncol(M)-1)){
      t <- t + length(which(M[3,i]*M[3,-(1:i)] > 0))
    }
    
    index <- (s_plus - s_minus)/sqrt((N_t * (N_t - 1) * 0.5 - t) * N_t * (N_t-1) * 0.5)
    
    return(index)
  }else{
    n <- nrow(data)
    
    if(is.null(dists)){
      dists <- matrix(numeric(), nrow = n, ncol = n)
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          dists[i,j] <- calc.dist(type = type, lambda = lambda, 
                                  data1 = data[i,], data2 = data[j,],
                                  data_plus = data_plus)
        }
      }
    }
    
    if(is.null(k)){k <- 2:sqrt(n)}
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], type = type, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
    trace_kp <- list(list("index" = tau_kproto(object = object, dists = dists, data_plus = data_plus), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, type = type, keep.data = TRUE, lambda = lambda, verbose = verbose, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = tau_kproto(object = object, dists = dists, data_plus = data_plus), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a clusterpartition with same number of cluster but different validation index
        index_value <- tau_kproto(object = object, dists = dists, data_plus = data_plus)
        if(!(index_value %in% unlist(lapply(trace_kp, `[[`, 1))[which(unlist(lapply(trace_kp, `[[`, 2)) == length(object$size))])){
          trace_kp <- c(trace_kp, list(list("index" = index_value, "k" = length(object$size), "object" = object)))
        }
      }
    }
    
    # save all index values for comparison and as summary for output
    indices <- unlist(lapply(trace_kp, `[[`, 1))
    names(indices) <- lapply(trace_kp, `[[`, 2)
    
    # find the optimal k, if it is ambiguously: sample
    k_m <- which(indices == max(indices))
    if(length(k_m)>1){k_m <- sample(k_m,1)}
    k_opt <- as.integer(names(indices[k_m]))
    index_opt <- indices[k_m]
    
    output <- switch(kp_obj,
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp, "type" = type),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object, "type" = type))
    return(output)
  }
}






#' @title Validating k Prototypes Clustering
#'
#' @description Calculating the preferred validation index for a k-Prototypes clustering with k clusters or computing the optimal number of clusters based on the choosen index for k-Prototype clustering. Possible validation indices are: \code{cindex}, \code{dunn}, \code{gamma}, \code{gplus}, \code{mcclain}, \code{ptbiserial}, \code{silhouette} and \code{tau}.
#' 
#' @param method Character specifying the validation index: \code{cindex}, \code{dunn}, \code{gamma}, \code{gplus}, \code{mcclain}, \code{ptbiserial}, \code{silhouette} (default) or \code{tau}.
#' @param object Object of class \code{kproto} resulting from a call with \code{kproto(..., keep.data=TRUE)}.
#' @param data Original data; only required if \code{object == NULL} and neglected if \code{object != NULL}.
#' @param type Character, to specify the distance for clustering; either \code{"huang"} or \code{"gower"}.
#' @param k Vector specifying the search range for optimum number of clusters; if \code{NULL} the range will set as \code{2:sqrt(n)}. Only required if \code{object == NULL} and neglected if \code{object != NULL}.
#' @param lambda Factor to trade off between Euclidean distance of numeric variables and simple matching coefficient between categorical variables.
#' @param kp_obj character either "optimal" or "all": Output of the index-optimal clustering (kp_obj == "optimal") or all computed cluster partitions (kp_obj == "all"); only required if \code{object != NULL}.
#' @param verbose Logical, whether additional information about process should be printed.
#' @param ... Further arguments passed to \code{\link[clustMixType]{kproto}}, like:
#'   \itemize{
#'     \item \code{nstart}: If > 1 repetitive computations of \code{kproto} with random initializations are computed.
#'     \item \code{na.rm}: Character, either \code{"yes"} to strip \code{NA} values for complete case analysis, \code{"no"} to keep and ignore \code{NA} values, \code{"imp.internal"} to impute the \code{NAs} within the algorithm or \code{"imp.onestep"} to apply the algorithm ignoring the \code{NAs} and impute them after the partition is determined.
#'   }
#' 
#' @details More information about the implemented validation indices:
#'   \itemize{
#'     \item \code{cindex} \deqn{Cindex = \frac{S_w-S_{min}}{S_{max}-S_{min}}} \cr
#' For \eqn{S_{min}} and \eqn{S_{max}} it is necessary to calculate the distances between all pairs of points in the entire data set (\eqn{\frac{n(n-1)}{2}}). 
#' \eqn{S_{min}} is the sum of the "total number of pairs of objects belonging to the same cluster" smallest distances and 
#' \eqn{S_{max}} is the sum of the "total number of pairs of objects belonging to the same cluster" largest distances. \eqn{S_w} is the sum of the within-cluster distances. \cr
#' The minimum value of the index is used to indicate the optimal number of clusters.
#' 
#'     \item \code{dunn} \deqn{Dunn = \frac{\min_{1 \leq i < j \leq q} d(C_i, C_j)}{\max_{1 \leq k \leq q} diam(C_k)}} \cr
#' The following applies: The dissimilarity between the two clusters \eqn{C_i} and \eqn{C_j} is defined as \eqn{d(C_i, C_j)=\min_{x \in C_i, y \in C_j} d(x,y)} and
#' the diameter of a cluster is defined as \eqn{diam(C_k)=\max_{x,y \in C} d(x,y)}. \cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#' 
#'     \item \code{gamma} \deqn{Gamma = \frac{s(+)-s(-)}{s(+)+s(-)}} \cr 
#' Comparisons are made between all within-cluster dissimilarities and all between-cluster dissimilarities. 
#' \eqn{s(+)} is the number of concordant comparisons and \eqn{s(-)} is the number of discordant comparisons.
#' A comparison is named concordant (resp. discordant) if a within-cluster dissimilarity is strictly less (resp. strictly greater) than a between-cluster dissimilarity.\cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#' 
#'     \item \code{gplus} \deqn{Gplus = \frac{2 \cdot s(-)}{\frac{n(n-1)}{2} \cdot (\frac{n(n-1)}{2}-1)}} \cr 
#' Comparisons are made between all within-cluster dissimilarities and all between-cluster dissimilarities. 
#' \eqn{s(-)} is the number of discordant comparisons and a comparison is named discordant if a within-cluster 
#' dissimilarity is strictly greater than a between-cluster dissimilarity. \cr
#' The minimum value of the index is used to indicate the optimal number of clusters.
#' 
#'     \item \code{mcclain} \deqn{McClain = \frac{\bar{S}_w}{\bar{S}_b}} \cr 
#' \eqn{\bar{S}_w} is the sum of within-cluster distances divided by the number of within-cluster distances and 
#' \eqn{\bar{S}_b} is the sum of between-cluster distances divided by the number of between-cluster distances.\cr
#' The minimum value of the index is used to indicate the optimal number of clusters.
#' 
#'     \item\code{ptbiserial} \deqn{Ptbiserial = \frac{(\bar{S}_b-\bar{S}_w) \cdot (\frac{N_w \cdot N_b}{N_t^2})^{0.5}}{s_d}} \cr 
#' \eqn{\bar{S}_w} is the sum of within-cluster distances divided by the number of within-cluster distances and 
#' \eqn{\bar{S}_b} is the sum of between-cluster distances divided by the number of between-cluster distances.\cr
#' \eqn{N_t} is the total number of pairs of objects in the data, \eqn{N_w} is the total number of pairs of 
#' objects belonging to the same cluster and \eqn{N_b} is the total number of pairs of objects belonging to different clusters.
#' \eqn{s_d} is the standard deviation of all distances.\cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#' 
#'     \item \code{silhouette} \deqn{Silhouette = \frac{1}{n} \sum_{i=1}^n \frac{b(i)-a(i)}{max(a(i),b(i))}} \cr 
#' \eqn{a(i)} is the average dissimilarity of the i\emph{th} object to all other objects of the same/own cluster.
#' \eqn{b(i)=min(d(i,C))}, where \eqn{d(i,C)} is the average dissimilarity of the i\emph{th} object to all the other clusters except the own/same cluster.\cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#'     
#'     \item \code{tau} \deqn{Tau = \frac{s(+) - s(-)}{((\frac{N_t(N_t-1)}{2}-t)\frac{N_t(N_t-1)}{2})^{0.5}}} \cr 
#' Comparisons are made between all within-cluster dissimilarities and all between-cluster dissimilarities. 
#' \eqn{s(+)} is the number of concordant comparisons and \eqn{s(-)} is the number of discordant comparisons.
#' A comparison is named concordant (resp. discordant) if a within-cluster dissimilarity is strictly less 
#' (resp. strictly greater) than a between-cluster dissimilarity.\cr
#' \eqn{N_t} is the total number of distances \eqn{\frac{n(n-1)}{2}} and \eqn{t} is the number of comparisons 
#' of two pairs of objects where both pairs represent within-cluster comparisons or both pairs are between-cluster
#' comparisons. \cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#'    
#'   }
#' 
#' @return For computing the optimal number of clusters based on the choosen validation index for k-Prototype clustering the output contains:
#' @return \item{k_opt}{optimal number of clusters (sampled in case of ambiguity)}
#' @return \item{index_opt}{index value of the index optimal clustering}
#' @return \item{indices}{calculated indices for \eqn{k=2,...,k_{max}}}
#' @return \item{kp_obj}{if(kp_obj == "optimal") the kproto object of the index optimal clustering and if(kp_obj == "all") all kproto which were calculated}
#' @return For computing the index-value for a given k-Prototype clustering the output contains:
#' @return \item{index}{calculated index-value}
#'
#' @author Rabea Aschenbruck
#' 
#' @references \itemize{
#'     \item Aschenbruck, R., Szepannek, G. (2020): 
#'     Cluster Validation for Mixed-Type Data. 
#'     \emph{Archives of Data Science, Series A, Vol 6, Issue 1}.
#'     \doi{10.5445/KSP/1000098011/02}.
#'     
#'     \item Charrad, M., Ghazzali, N., Boiteau, V., Niknafs, A. (2014): 
#'     NbClust: An R Package for Determining the Relevant Number of Clusters in a Data Set. 
#'     \emph{Journal of Statistical Software, Vol 61, Issue 6}.
#'     \doi{10.18637/jss.v061.i06}.
#'   }
#'
#' @examples
#' \dontrun{
#' # generate toy data with factors and numerics
#' n   <- 10
#' prb <- 0.99
#' muk <- 2.5 
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' 
#' # calculate optimal number of cluster, index values and clusterpartition with Silhouette-index
#' val <- validation_kproto(method = "silhouette", data = x, k = 3:5, nstart = 5)
#' 
#' 
#' # apply k-prototypes
#' kpres <- kproto(x, 4, keep.data = TRUE)
#' 
#' # calculate cindex-value for the given clusterpartition
#' cindex_value <- validation_kproto(method = "cindex", object = kpres)
#' }
#' 
#' @rdname validation_kproto
#' 
#' @importFrom stats na.omit
#' @importFrom utils combn
#' @importFrom utils head
#' @importFrom utils tail
#' 
#' 
#' @export
validation_kproto <- function(method = "silhouette", object = NULL, data = NULL, type = "huang", k = NULL, lambda = NULL, kp_obj = "optimal", verbose = FALSE, ...){

  if(is.null(method)) stop("validation methode must be choosen!")
  if(!(method %in% c("cindex", "dunn", "gamma", "gplus", "mcclain", "ptbiserial", "silhouette", "tau"))) stop("choose one of these methods: cindex, dunn, gamma, gplus, mcclain, ptbiserial, silhouette, tau")
  
  if(is.null(data) & is.null(object)) stop("'data' or 'object' must be given!")
  
  if(!is.null(object) & !inherits(object, "kproto")) stop("object must be type of 'kproto'")
  if(!is.null(object) & is.null(object$data)) stop("kproto_object should have the original data included (kproto(..., keep.data = TRUE))")
  if(!is.null(object)){
    
    if(object$type == "gower"){
      data_plus <- list()
      
      data_plus$numvars <- sapply(object$data, is.numeric)
      data_plus$ordvars <- sapply(object$data, is.ordered)
      data_plus$catvars <- sapply(object$data, is.factor) & !data_plus$ordvars
      
      # ranges for normalization and transformed data for gower's distance calculation
      if(any(data_plus$numvars)) data_plus$rgnums <- sapply(object$data[, data_plus$numvars, drop = FALSE], function(z) diff(range(z)))
      if(any(data_plus$ordvars)){
        data_plus$data_ord <- object$data[, data_plus$ordvars, drop = FALSE] # store original variables 
        # ...and replace ordered variables by their ranks
        for(jord in which(data_plus$ordvars)) object$data[,jord] <- rank(object$data[,jord])
        data_plus$rgords <- sapply(object$data[, data_plus$ordvars, drop = FALSE], function(z) diff(range(z)))
      }
    }
    if(object$type == "huang"){
      data_plus <- list()
      data_plus$numvars <- sapply(object$data, is.numeric)
      data_plus$catvars <- sapply(object$data, is.factor)
    }
  }
  
  if(is.null(object)){
    if(!is.data.frame(data)) stop("data should be a data frame!")
    if(ncol(data) < 2) stop("for clustering data should contain at least two variables!")
    if(nrow(data) < 4) stop("for clustering data should contain at least four objects!")
    
    if(type %in% c("Huang", "standard")) type <- "huang"
    if(type == "Gower") type <- "gower"
    if(!type %in% c("huang", "gower")) stop("Argument type must be either 'huang' or 'gower'!")
    
    if(type == "gower"){
      type <- "gower"
      data_plus <- list()
      
      data_plus$numvars <- sapply(data, is.numeric)
      data_plus$ordvars <- sapply(data, is.ordered)
      data_plus$catvars <- sapply(data, is.factor) & !data_plus$ordvars
      
      # ranges for normalization and transformed data for gower's distance calculation
      if(any(data_plus$numvars)) data_plus$rgnums <- sapply(data[, data_plus$numvars, drop = FALSE], function(z) diff(range(z)))
      if(any(data_plus$ordvars)){
        data_plus$data_ord <- data[, data_plus$ordvars, drop = FALSE] # store original variables 
        # ...and replace ordered variables by their ranks
        for(jord in which(data_plus$ordvars)) data[,jord] <- rank(data[,jord])
        data_plus$rgords <- sapply(data[, data_plus$ordvars, drop = FALSE], function(z) diff(range(z)))
      }
      
      if(length(lambda) > 0 & length(lambda) != sum(c(data_plus$numvars, data_plus$catvars, data_plus$ordvars))) {
        warning("For gower distance if lambda is specified, its length must be the sum of numeric and factor variables in the data frame!")
        lambda <- NULL
      }
    }
    
    if(type == "huang"){
      data_plus <- list()
      
      data_plus$numvars <- sapply(data, is.numeric)
      data_plus$catvars <- sapply(data, is.factor)
      
      if(!is.null(lambda)){
        if(length(lambda) > 1){
          if(length(lambda) != sum(c(data_plus$numvars, data_plus$catvars))) stop("If lambda is a vector, its length should be the sum of numeric and factor variables in the data frame!")
        }else{
          if(length(lambda) == 1) {if(lambda == 0) stop("lambda has to be a value != 0. For automatic calculation use lambda = NULL (default setting)!")}
        }
      }else{
        lambda <- lambdaest(x = data, num.method = 1, fac.method = 1, outtype = "numeric", verbose = FALSE)
        if(verbose) cat("Estimated lambda:", lambda, "\n\n")
      }
    }
    
    if(!is.null(k) & length(k) == 1){
      stop("k should be the search range for optimum number of clusters, e.g. c(2:sqrt(n))")
    }
    if(length(k) > 1){
      if(nrow(data) < max(k)) stop("Data frame has less observations than clusters!")
      if(any(k <= 1) | any(k >= nrow(data))) stop("Elements of k must be greater than 1 and strictly less than n!")
      if(all(!as.integer(k)==k)) stop("Elements of k must be type of integer")
    }
    if(1 %in% k){
      if(method %in% c("cindex", "mcclain", "ptbiserial")) stop(paste("Index calculation not possible for k = 1. As the", method, "index is based on the relation of within-cluster distances and between-cluster distances, choose k > 1."))
      if(method %in% c("gamma", "gplus", "tau")) stop(paste("Index calculation not possible for k = 1. Since the", method, "index insists on comparing within-cluster distances and between-cluster distances, choose k > 1."))
      if(method == "dunn") stop("Index calculation not possible for k = 1. Calculation of the dunn index requires the determination of the distance between two clusters, so you have to set k > 1.")
      if(method == "silhouette") stop("Index calculation not possible for k = 1. The Silhouette index rates the average within-cluster distance for the own cluster and for the best alternative cluster for every object, so you have to set k > 1.")
    }
    if(length(object$size) == 1){
      if(method %in% c("cindex", "mcclain", "ptbiserial")) stop(paste("Index calculation not possible for kproto-objects with only one cluster. As the", method, "index is based on the relation of within-cluster distances and between-cluster distances a kproto-object which has been computed with k > 1 is required."))
      if(method %in% c("gamma", "gplus", "tau")) stop(paste("Index calculation not possible for kproto-objects with only one cluster. Since the", method, "index insists on comparing within-cluster distances and between-cluster distances a kproto-object which has been computed with k > 1 is required."))
      if(method == "dunn") stop("Index calculation not possible for kproto-objects with only one cluster. Calculation of the dunn index requires the determination of the distance between two clusters. You have to insert a kproto-object which has been computed with k > 1.")
      if(method == "silhouette") stop("Index calculation not possible for kproto-objects with only one cluster. The Silhouette index rates the average within-cluster distance for the own cluster and for the best alternative cluster for every object. You have to insert a kproto-object which has been computed with k > 1.")
    }
    if(is.null(k)){
      k <- 2:sqrt(nrow(data))
    }
  }
  
  if(!(kp_obj %in% c("optimal","all"))) stop("kp_obj must either be 'optimal' or 'all'!")
  
  output <- switch(method,
                   "cindex" = cindex_kproto(object = object, type = type, data = data, k = k, kp_obj = kp_obj, lambda = lambda, data_plus = data_plus, ...),
                   "dunn" = dunn_kproto(object = object, type = type, data = data, k = k, kp_obj = kp_obj, lambda = lambda, data_plus = data_plus, ...),
                   "gamma" = gamma_kproto(object = object, type = type, data = data, k = k, kp_obj = kp_obj, lambda = lambda, data_plus = data_plus, ...),
                   "gplus" = gplus_kproto(object = object, type = type, data = data, k = k, kp_obj = kp_obj, lambda = lambda, data_plus = data_plus, ...),
                   "mcclain" = mcclain_kproto(object = object, type = type, data = data, k = k, kp_obj = kp_obj, lambda = lambda, data_plus = data_plus, ...),
                   "ptbiserial" = ptbiserial_kproto(object = object, type = type, data = data, k = k, kp_obj = kp_obj, lambda = lambda, data_plus = data_plus, ...),
                   "silhouette" = silhouette_kproto(object = object, type = type, data = data, k = k, kp_obj = kp_obj, lambda = lambda, data_plus = data_plus, ...),
                   "tau" = tau_kproto(object = object, type = type, data = data, k = k, kp_obj = kp_obj, lambda = lambda, data_plus = data_plus, ...))
  
  return(output)
}











