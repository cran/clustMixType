


create.S_w_kpa <- function(object){
  numvars <- sapply(object$data, is.numeric)
  catvars <- sapply(object$data, is.factor)
  S_w <- 0
  for(k in seq_along(object$size)){
    x_k <- object$data[which(object$cluster==k),]
    if(object$size[k]>1){

      d1 <- (x_k[rep(1:(object$size[k]-1),times=(object$size[k]-1):1, each=1), numvars, drop = FALSE] -
               x_k[unlist(lapply(2:object$size[k], seq, to=object$size[k])), numvars, drop = FALSE])^2
      if(length(object$lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
      if(length(object$lambda) > 1) d1 <- as.matrix(d1) %*% object$lambda[numvars]
      
      d2 <- x_k[rep(1:(object$size[k]-1),times=(object$size[k]-1):1, each=1),
                catvars, drop = FALSE] != x_k[unlist(lapply(2:object$size[k], seq, to=object$size[k])), catvars, drop = FALSE]
      d2[is.na(d2)] <- FALSE
      if(length(object$lambda) == 1) d2 <- object$lambda * rowSums(d2)
      if(length(object$lambda) > 1) d2 <- as.matrix(d2) %*% object$lambda[catvars]
      
      S_w <- S_w + sum(d1 + d2)
    }
  }
  return(S_w)
}


create.S_b_kpa <- function(object){
  numvars <- sapply(object$data, is.numeric)
  catvars <- sapply(object$data, is.factor)
  S_b <- 0
  for(k in 1:(length(object$size)-1)){
    x_k <- object$data[which(object$cluster==k),]
    for(l in (k+1):length(object$size)){
      x_l <- object$data[which(object$cluster==l),]
      
      d1 <- (x_k[rep(1:object$size[k], each=object$size[l]), numvars, drop = FALSE] -
               x_l[rep(1:object$size[l], times=object$size[k]), numvars, drop = FALSE])^2
      if(length(object$lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
      if(length(object$lambda) > 1) d1 <- as.matrix(d1) %*% object$lambda[numvars]

      d2 <- x_k[rep(1:object$size[k], each=object$size[l]),
                catvars, drop = FALSE] != x_l[rep(1:object$size[l], times = object$size[k]), catvars, drop = FALSE]
      d2[is.na(d2)] <- FALSE
      if(length(object$lambda) == 1) d2 <- object$lambda * rowSums(d2)
      if(length(object$lambda) > 1) d2 <- as.matrix(d2) %*% object$lambda[catvars]

      S_b <- S_b + sum(d1 + d2)
    }
  }
  return(S_b)
}


create.dist_kpa <- function(lambda = NULL, data1, data2){
  data <- rbind(data1, data2)
  numvars <- sapply(data, is.numeric)
  catvars <- sapply(data, is.factor)
  
  d1 <- (data[1,numvars, drop = FALSE]-data[2,numvars, drop = FALSE])^2
  if(length(lambda) == 1) d1 <- sum(d1)
  if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
  
  d2 <- sapply(which(catvars), function(j) return(data[1,j] != data[2,j]))
  if(length(lambda) == 1) d2 <- lambda * sum(d2)
  if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
  
  return(d1 + d2)
}


create.N_w <- function(object){
  sum((object$size*(object$size-1))/2)
}




cindex_kproto <- function(object = NULL, data = NULL, k = NULL, S_sort = NULL, kp_obj = "optimal", ...){
  
  if(!is.null(object)){
    if(is.null(object$data)) stop("object should have the original data included (kproto(..., keep.data = TRUE))")
    
    if(is.null(S_sort)){
      n <- length(object$cluster)
      S_all <- matrix(numeric(), ncol = n, nrow = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          S_all[i,j] <- create.dist_kpa(lambda = object$lambda, data1 = object$data[i,], data2 = object$data[j,])
        }
      }
      S_sort <- sort(S_all)
    }
    
    S_w_kpa <- create.S_w_kpa(object)
    
    N_w <- create.N_w(object)
    S_min <- sum(head(S_sort,n = N_w))
    S_max <- sum(tail(S_sort,n = N_w))
    
    if(S_min != S_max){
      index <- (S_w_kpa - S_min)/(S_max - S_min)
    }
    
    return(index)
  }else{
    n <- nrow(data)
    p <- ncol(data)
    
    numvars <- sapply(data, is.numeric)
    anynum <- any(numvars)
    catvars <- sapply(data, is.factor)
    anyfact <- any(catvars)
    vnum <- mean(sapply(data[,numvars, drop = FALSE], var, na.rm = TRUE))
    vcat <- mean(sapply(data[,catvars, drop = FALSE], function(z) return(1-sum((table(z)/sum(!is.na(z)))^2))))
    lambda <- vnum/vcat
    
    S_all <- matrix(numeric(), ncol=n, nrow=n)
    for(i in 1:(n - 1)){
      for(j in (i + 1):n){
        S_all[i,j] <- create.dist_kpa(lambda = lambda, data1 = data[i,], data2 = data[j,])
      }
    }
    S_sort <- sort(S_all)
    
    if(is.null(k)){k <- 2:sqrt(n)}
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], keep.data = TRUE, ...)
    trace_kp <- list(list("index" = cindex_kproto(object = object, S_sort = S_sort), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, keep.data = TRUE, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = cindex_kproto(object = object, S_sort = S_sort), 
                                        "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a clusterpartition with same number of cluster but different validation index
        index_value <- cindex_kproto(object = object, S_sort = S_sort)
        if(trace_kp[[which(unlist(lapply(trace_kp, `[[`, 2)) == 2)]]$index != index_value){
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
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object))
    return(output)
  }
}


dunn_kproto <- function(object = NULL, data = NULL, k = NULL, kp_obj = "optimal", ...){
  
  if(!is.null(object)){
    if(is.null(object$data)) stop("object should have the original data included (kproto(..., keep.data = TRUE))")
    
    k <- length(object$size)
    
    #determine d(C_i,C_j)
    min_CiCj <- matrix(numeric(k*k),ncol = k,nrow = k)
    for(i in 1:(k-1)){
      xi <- object$data[which(object$cluster == i),]
      for(j in (i+1):k){
        xj <- object$data[which(object$cluster == j),]
        min_ij <- create.dist_kpa(object$lambda, data1 = xi[1,], data2 = xj[1,])
        for(l in 1:object$size[i]){
          for(m in 1:object$size[j]){
            min_neu <- create.dist_kpa(object$lambda, data1 = xi[l,], data2 = xj[m,])
            if(min_neu < min_ij){
              min_ij <- min_neu
            }
          }
        }
        min_CiCj[i,j] <- min_ij
      }
    }
    Zaehler <- min(min_CiCj[min_CiCj > 0])
    
    #determine diam(C_k)
    max_diam <- numeric(k)
    for(p in 1:k){
      xi <- object$data[which(object$cluster == p),]
      if(object$size[p] > 1){
        max_ij <- create.dist_kpa(object$lambda, data1 = xi[1,], data2 = xi[2,])
        for(l in 1:(object$size[p]-1)){
          for(m in (l+1):(object$size[p])){
            max_neu <- create.dist_kpa(object$lambda, data1 = xi[l,], data2 = xi[m,])
            if(max_neu > max_ij){
              max_ij <- max_neu
            }
          }
        }
        max_diam[p] <- max_ij
      }
    }
    Nenner <- max(max_diam)
    
    if(is.finite(Zaehler/Nenner)){
      return(Zaehler/Nenner)
    }else{
      return(NA)
    }
    
    
  }else{
    n <- nrow(data) 
    
    if(is.null(k)){k <- 2:sqrt(n)}
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], keep.data = TRUE, ...)
    trace_kp <- list(list("index" = dunn_kproto(object = object), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, keep.data = TRUE, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = dunn_kproto(object = object), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a clusterpartition with same number of cluster but different validation index
        index_value <- dunn_kproto(object = object)
        if(trace_kp[[which(unlist(lapply(trace_kp, `[[`, 2)) == 2)]]$index != index_value){
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
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object))
    return(output)

  }
}


gamma_kproto <- function(object = NULL, data = NULL, k = NULL, dists = NULL, kp_obj = "optimal", ...){
  
  if(!is.null(object)){
    if(is.null(object$data)) stop("object should have the original data included (kproto(..., keep.data = TRUE))")
    
    if(is.null(dists)){
      n <- nrow(object$data)
      dists <- matrix(numeric(), nrow = n, ncol = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- create.dist_kpa(object$lambda, data1 = object$data[i,], data2 = object$data[j,])
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
    
    numvars <- sapply(data, is.numeric)
    anynum <- any(numvars)
    catvars <- sapply(data, is.factor)
    anyfact <- any(catvars)
    vnum <- mean(sapply(data[,numvars, drop = FALSE], var, na.rm = TRUE))
    vcat <- mean(sapply(data[,catvars, drop = FALSE], function(z) return(1 - sum((table(z)/sum(!is.na(z)))^2))))
    lambda <- vnum/vcat
    
    if(is.null(dists)){
      dists <- matrix(numeric(), nrow = n, ncol = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- create.dist_kpa(lambda, data1 = data[i,], data2 = data[j,])
        }
      }
    }
    
    index <- numeric(n)
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], keep.data = TRUE, ...)
    trace_kp <- list(list("index" = gamma_kproto(object = object, dists = dists), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, keep.data = TRUE, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = gamma_kproto(object = object, dists = dists), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a clusterpartition with same number of cluster but different validation index
        index_value <- gamma_kproto(object = object, dists = dists)
        if(trace_kp[[which(unlist(lapply(trace_kp, `[[`, 2)) == 2)]]$index != index_value){
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
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object))
    return(output)

  }
}


gplus_kproto <- function(object = NULL, data = NULL, k = NULL, dists = NULL, kp_obj = "optimal", ...){
  
  if(!is.null(object)){
    if(is.null(object$data)) stop("object should have the original data included (kproto(..., keep.data = TRUE))")
    
    n <- nrow(object$data)
    if(is.null(dists)){
      dists <- matrix(numeric(), nrow=n, ncol=n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- create.dist_kpa(object$lambda, data1 = object$data[i,], data2 = object$data[j,])
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
    
    numvars <- sapply(data, is.numeric)
    anynum <- any(numvars)
    catvars <- sapply(data, is.factor)
    anyfact <- any(catvars)
    vnum <- mean(sapply(data[,numvars, drop = FALSE], var, na.rm = TRUE))
    vcat <- mean(sapply(data[,catvars, drop = FALSE], function(z) return(1 - sum((table(z)/sum(!is.na(z)))^2))))
    lambda <- vnum/vcat
    
    if(is.null(dists)){
      dists <- matrix(numeric(), nrow = n, ncol = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- create.dist_kpa(lambda, data1 = data[i,], data2 = data[j,])
        }
      }
    }
    
    if(is.null(k)){k <- 2:sqrt(n)}
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], keep.data = TRUE, ...)
    trace_kp <- list(list("index" = gplus_kproto(object = object, dists = dists), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, keep.data = TRUE, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = gplus_kproto(object = object, dists = dists), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a clusterpartition with same number of cluster but different validation index
        index_value <- gplus_kproto(object = object, dists = dists)
        if(trace_kp[[which(unlist(lapply(trace_kp, `[[`, 2)) == 2)]]$index != index_value){
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
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object))
    return(output)
  }
}


mcclain_kproto <- function(object = NULL, data = NULL, k = NULL, kp_obj = "optimal", ...){
  
  if(!is.null(object)){
    if(is.null(object$data)) stop("kproto_object should have the original data included (kproto(..., keep.data = TRUE))")
    
    n <- nrow(object$data)
    S_w_kpa <- create.S_w_kpa(object)
    S_b_kpa <- create.S_b_kpa(object)
    N_w <- create.N_w(object)
    N_t <- n*(n - 1)/2
    N_b <- N_t-N_w
    
    index <- (S_w_kpa/N_w)/(S_b_kpa/N_b)
    
    return(index)
  }else{
    n <- nrow(data) 
    
    if(is.null(k)){k <- 2:sqrt(n)}
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], keep.data = TRUE, ...)
    trace_kp <- list(list("index" = mcclain_kproto(object = object), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, keep.data = TRUE, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = mcclain_kproto(object = object), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a clusterpartition with same number of cluster but different validation index
        index_value <- mcclain_kproto(object = object)
        if(trace_kp[[which(unlist(lapply(trace_kp, `[[`, 2)) == 2)]]$index != index_value){
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
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object))
    return(output)
  }
}


ptbiserial_kproto <- function(object = NULL, data = NULL, k = NULL, s_d = NULL, kp_obj = "optimal", ...){
  
  if(!is.null(object)){
    if(is.null(object$data)) stop("object should have the original data included (kproto(...,keep.data=TRUE))")
    
    n <- nrow(object$data)
    if(is.null(s_d)){
      S_all <- matrix(numeric(), ncol = n, nrow = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          S_all[i,j] <- create.dist_kpa(lambda = object$lambda, data1 = object$data[i,], data2 = object$data[j,])
        }
      }
      s_d <- sd(S_all, na.rm = TRUE)
    }
    
    S_w_kpa <- create.S_w_kpa(object)
    S_b_kpa <- create.S_b_kpa(object)
    N_w <- create.N_w(object)
    N_t <- n*(n-1)/2
    N_b <- N_t-N_w
    
    index <- ((S_b_kpa/N_b) - (S_w_kpa/N_w)) * sqrt(N_w * N_b/N_t^2) / s_d
    
    return(index)
  }else{
    n <- nrow(data)
    
    numvars <- sapply(data, is.numeric)
    anynum <- any(numvars)
    catvars <- sapply(data, is.factor)
    anyfact <- any(catvars)
    vnum <- mean(sapply(data[,numvars, drop = FALSE], var, na.rm = TRUE))
    vcat <- mean(sapply(data[,catvars, drop = FALSE], function(z) return(1-sum((table(z)/sum(!is.na(z)))^2))))
    lambda <- vnum/vcat
    
    if(is.null(s_d)){
      S_all <- matrix(numeric(), ncol = n, nrow = n)
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          S_all[i,j] <- create.dist_kpa(lambda = lambda, data1 = data[i,], data2 = data[j,])
        }
      }
      s_d <- sd(S_all, na.rm = TRUE)
    }
    
    if(is.null(k)){k <- 2:sqrt(n)}
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], keep.data = TRUE, ...)
    trace_kp <- list(list("index" = ptbiserial_kproto(object = object, s_d = s_d), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, keep.data = TRUE, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = ptbiserial_kproto(object = object, s_d = s_d), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a clusterpartition with same number of cluster but different validation index
        index_value <- ptbiserial_kproto(object = object, s_d = s_d)
        if(trace_kp[[which(unlist(lapply(trace_kp, `[[`, 2)) == 2)]]$index != index_value){
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
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object))
    return(output)
  }
}


silhouette_kproto <- function(object = NULL, data = NULL, k = NULL, kp_obj = "optimal", ...){
  
  if(!is.null(object)){
    if(is.null(object$data)) stop("object should have the original data included (kproto(...,keep.data=TRUE))")
    
    n <- nrow(object$data)
    x <- object$data
    cluster <- object$cluster
    k <- length(table(cluster))
    lambda <- object$lambda
    numvars <- sapply(x, is.numeric)
    catvars <- sapply(x, is.factor)
    
    
    protos <- x
    nrows <- nrow(x)
    dists <- matrix(NA, nrow = nrows, ncol = nrows)
    for(i in 1:nrows){
      #distances of the numeric variables
      d1 <- (x[,numvars, drop = FALSE] - matrix(rep(as.numeric(protos[i, numvars, drop = FALSE]), nrows), nrow = nrows, byrow = TRUE))^2
      if(length(lambda) == 1) d1 <- rowSums(d1)
      if(length(lambda) > 1) d1 <- d1 %*% lambda[numvars]
      
      #distances of the categorical variances
      d2 <- sapply(which(catvars), function(j) return(x[,j] != rep(protos[i,j], nrows)))
      if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
      if(length(lambda) > 1) d2 <- d2 %*% lambda[catvars]
      
      dists[,i] <- d1 + d2
    }
    
    cluster_dists <- matrix(numeric(nrows*k), nrow = nrows, ncol = k)
    for(i in 1:k){
      if(!(length(which(cluster == i)) == 1)){
        cluster_dists[,i] <- rowMeans(dists[,which(cluster == i)])
      }else{
        cluster_dists[,i] <- dists[,which(cluster == i)]
      }
    }
    
    #determine ai, bi and si
    a <- numeric(nrows)
    b <- numeric(nrows)
    s <- numeric(nrows)
    for(i in 1:nrows){
      a[i] <- cluster_dists[i,cluster[i]]
      b[i] <- min(cluster_dists[i,-cluster[i]])
      s[i] <- (b[i] - a[i])/max(a[i],b[i])
    }
    if(any(table(cluster) == 1)){
      for(i in which(cluster %in% as.integer(which(table(cluster) == 1)))){
        s[i] <- 0
      }
      cat(length(which(cluster %in% as.integer(which(table(cluster) == 1))))," Cluster mit nur einem Element\n")
    }
    
    index <- mean(s)
    
    return(index)
  }else{
    n <- nrow(data) 
    
    if(is.null(k)){k <- 2:sqrt(n)}
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], keep.data = TRUE, ...)
    trace_kp <- list(list("index" = silhouette_kproto(object = object), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, keep.data = TRUE, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = silhouette_kproto(object = object), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a clusterpartition with same number of cluster but different validation index
        index_value <- silhouette_kproto(object = object)
        if(trace_kp[[which(unlist(lapply(trace_kp, `[[`, 2)) == 2)]]$index != index_value){
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
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object))
    return(output)
  }
}


tau_kproto <- function(object = NULL, data = NULL, k = NULL, dists = NULL, kp_obj = "optimal", ...){
  
  if(!is.null(object)){
    if(is.null(object$data)) stop("object should have the original data included (kproto(..., keep.data = TRUE))")
    
    if(is.null(dists)){
      n <- nrow(object$data)
      dists <- matrix(numeric(), nrow = n, ncol = n)
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          dists[i,j] <- create.dist_kpa(lambda = object$lambda, data1 = object$data[i,], data2 = object$data[j,])
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
    M <- rbind(M,apply(X = M, MARGIN = 2, function(x) if(object$cluster[x[1]] == object$cluster[x[2]]){return(1)}else{return(-1)}))
    t <- 0
    for(i in 1:(ncol(M)-1)){
      t <- t + length(which(M[3,i]*M[3,-(1:i)] > 0))
    }
    
    index <- (s_plus - s_minus)/sqrt((N_t * (N_t - 1) * 0.5 - t) * N_t * (N_t-1) * 0.5)
    
    return(index)
  }else{
    n <- nrow(data)
    
    numvars <- sapply(data, is.numeric)
    anynum <- any(numvars)
    catvars <- sapply(data, is.factor)
    anyfact <- any(catvars)
    vnum <- mean(sapply(data[, numvars, drop = FALSE], var, na.rm = TRUE))
    vcat <- mean(sapply(data[, catvars, drop = FALSE], function(z) return(1 - sum((table(z)/sum(!is.na(z)))^2))))
    lambda <- vnum/vcat
    
    if(is.null(dists)){
      dists <- matrix(numeric(), nrow = n, ncol = n)
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          dists[i,j] <- create.dist_kpa(lambda = lambda, data1 = data[i,], data2 = data[j,])
        }
      }
    }
    
    if(is.null(k)){k <- 2:sqrt(n)}
    
    #calculate all kproto objects for k
    object <- kproto(x = data, k = k[1], keep.data = TRUE, ...)
    trace_kp <- list(list("index" = tau_kproto(object = object, dists = dists), 
                          "k" = length(object$size), "object" = object))
    for(q in k[-1]){
      object <- kproto(x = data, k = q, keep.data = TRUE, ...)
      #save kproto object, if there isn't an object for this number of cluster
      if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
        trace_kp <- c(trace_kp, list(list("index" = tau_kproto(object = object, dists = dists), 
                                          "k" = length(object$size), "object" = object)))
      }else{
        # save kproto object if there is a clusterpartition with same number of cluster but different validation index
        index_value <- tau_kproto(object = object, dists = dists)
        if(trace_kp[[which(unlist(lapply(trace_kp, `[[`, 2)) == 2)]]$index != index_value){
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
                     "all" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp),
                     "optimal" = list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices, "kp_obj" = trace_kp[[k_m]]$object))
    return(output)
  }
}






#' @title Validating k Prototypes Clustering
#'
#' @description Calculating the prefered validation index for a k-Prototypes clustering with k clusters or computing the optimal number of clusters based on the choosen index for k-Prototype clustering. Possible validation indices are: \code{cindex}, \code{dunn}, \code{gamma}, \code{gplus}, \code{mcclain}, \code{ptbiserial}, \code{silhouette} and \code{tau}.
#' 
#' @param method character specifying the validation index: \code{cindex}, \code{dunn}, \code{gamma}, \code{gplus}, \code{mcclain}, \code{ptbiserial}, \code{silhouette} or \code{tau}.
#' @param object Object of class \code{kproto} resulting from a call with \code{kproto(..., keep.data=TRUE)}
#' @param data Original data; only required if \code{object == NULL} and neglected if \code{object != NULL}.
#' @param k Vector specifying the search range for optimum number of clusters; if \code{NULL} the range will set as \code{2:sqrt(n)}. Only required if \code{object == NULL} and neglected if \code{object != NULL}.
#' @param kp_obj character either "optimal" or "all": Output of the index-optimal clustering (kp_obj == "optimal") or all computed clusterpartitions (kp_obj == "all"); only required if \code{object != NULL}.
#' @param ... Further arguments passed to \code{\link[clustMixType]{kproto}}, like:
#'   \itemize{
#'     \item \code{nstart}: If > 1 repetetive computations of \code{kproto} with random initializations are computed.
#'     \item \code{lambda}: Factor to trade off between Euclidean distance of numeric variables and simple matching coefficient between categorical variables.
#'     \item \code{verbose}: Logical whether information about the cluster procedure should be given. Caution: If \code{verbose=FALSE}, the reduction of the number of clusters is not mentioned.
#'   }
#' 
#' @details More information about the implemented validation indices:
#'   \itemize{
#'     \item{\code{cindex}} {\deqn{Cindex = \frac{S_w-S_{min}}{S_{max}-S_{min}}} \cr
#' For \eqn{S_{min}} and \eqn{S_{max}} it is nessesary to calculate the distances between all pairs of points in the entire data set (\eqn{\frac{n(n-1)}{2}}). 
#' \eqn{S_{min}} is the sum of the "total number of pairs of objects belonging to the same cluster" smallest distances and 
#' \eqn{S_{max}} is the sum of the "total number of pairs of objects belonging to the same cluster" largest distances. \eqn{S_w} is the sum of the within-cluster distances. \cr
#' The minimum value of the index is used to indicate the optimal number of clusters.
#'     }
#'     \item{\code{dunn}}{\deqn{Dunn = \frac{\min_{1 \leq i < j \leq q} d(C_i, C_j)}{\max_{1 \leq k \leq q} diam(C_k)}} \cr
#' The following applies: The dissimilarity between the two clusters \eqn{C_i} and \eqn{C_j} is defined as \eqn{d(C_i, C_j)=\min_{x \in C_i, y \in C_j} d(x,y)} and
#' the diameter of a cluster is defined as \eqn{diam(C_k)=\max_{x,y \in C} d(x,y)}. \cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#'     }
#'     \item{\code{gamma}}{\deqn{Gamma = \frac{s(+)-s(-)}{s(+)+s(-)}} \cr 
#' Comparisons are made between all within-cluster dissimilarities and all between-cluster dissimilarities. 
#' \eqn{s(+)} is the number of concordant comparisons and \eqn{s(-)} is the number of discordant comparisons.
#' A comparison is named concordant (resp. discordant) if a within-cluster dissimilarity is strictly less (resp. strictly greater) than a between-cluster dissimilarity.\cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#'     }
#'     \item{\code{gplus}}{\deqn{Gplus = \frac{2 \cdot s(-)}{\frac{n(n-1)}{2} \cdot (\frac{n(n-1)}{2}-1)}} \cr 
#' Comparisons are made between all within-cluster dissimilarities and all between-cluster dissimilarities. 
#' \eqn{s(-)} is the number of discordant comparisons and a comparison is named discordant if a within-cluster 
#' dissimilarity is strictly greater than a between-cluster dissimilarity. \cr
#' The minimum value of the index is used to indicate the optimal number of clusters.
#'     }
#'     \item{\code{mcclain}}{\deqn{McClain = \frac{\bar{S}_w}{\bar{S}_b}} \cr 
#' \eqn{\bar{S}_w} is the sum of within-cluster distances divided by the number of within-cluster distances and 
#' \eqn{\bar{S}_b} is the sum of between-cluster distances divided by the number of between-cluster distances.\cr
#' The minimum value of the index is used to indicate the optimal number of clusters.
#'     }
#'     \item{\code{ptbiserial}}{\deqn{Ptbiserial = \frac{(\bar{S}_b-\bar{S}_w) \cdot (\frac{N_w \cdot N_b}{N_t^2})^{0.5}}{s_d}} \cr 
#' \eqn{\bar{S}_w} is the sum of within-cluster distances divided by the number of within-cluster distances and 
#' \eqn{\bar{S}_b} is the sum of between-cluster distances divided by the number of between-cluster distances.\cr
#' \eqn{N_t} is the total number of pairs of objects in the data, \eqn{N_w} is the total number of pairs of 
#' objects belonging to the samecluster and \eqn{N_b} is the total number of pairs of objects belonging to different clusters.
#' \eqn{s_d} is the standard deviation of all distances.\cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#'     }
#'     \item{\code{silhouette}}{\deqn{Silhouette = \frac{1}{n} \sum_{i=1}^n \frac{b(i)-a(i)}{max(a(i),b(i))}} \cr 
#' \eqn{a(i)} is the average dissimilarity of the i\emph{th} object to all other objects of the same/own cluster.
#' \eqn{b(i)=min(d(i,C))}, where \eqn{d(i,C)} is the average dissimilarity of the i\emph{th} object to all the other clusters except the own/same cluster.\cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#'     }
#'     \item{\code{tau}}{\deqn{Tau = \frac{s(+) - s(-)}{((\frac{N_t(N_t-1)}{2}-t)\frac{N_t(N_t-1)}{2})^{0.5}}} \cr 
#' Comparisons are made between all within-cluster dissimilarities and all between-cluster dissimilarities. 
#' \eqn{s(+)} is the number of concordant comparisons and \eqn{s(-)} is the number of discordant comparisons.
#' A comparison is named concordant (resp. discordant) if a within-cluster dissimilarity is strictly less 
#' (resp. strictly greater) than a between-cluster dissimilarity.\cr
#' \eqn{N_t} is the total number of distances \eqn{\frac{n(n-1)}{2}} and \eqn{t} is the number of comparisons 
#' of two pairs of objects where both pairs represent within-cluster comparisons or both pairs are between-cluster
#' comparisons. \cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#'     }
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
#'     \item Charrad, M., Ghazzali, N., Boiteau, V., Niknafs, A. (2014): 
#'     NbClust: An R Package for Determining the Relevant Number of Clusters in a Data Set. 
#'     \emph{Journal of Statistical Software, Vol 61, Issue 6}.
#'     \href{http://www.jstatsoft.org/v61/i06/}{\emph{www.jstatsoft.org}}.
#'   }
#'
#' @examples
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
validation_kproto <- function(method = NULL, object = NULL, data = NULL, k = NULL, kp_obj = "optimal", ...){
  
  if(is.null(method)) stop("validation methode must be choosen!")
  if(!(method %in% c("cindex", "dunn", "gamma", "gplus", "mcclain", "ptbiserial", "silhouette", "tau"))) stop("choose one of these methods: cindex, dunn, gamma, gplus, mcclain, ptbiserial, silhouette, tau")
  
  if(is.null(data) && is.null(object)) stop("data or object muss be given!")
  
  if(!is.null(data) && !is.data.frame(data)) stop("data should be a data frame!")
  if(!is.null(data) && ncol(data) < 2) stop("for clustering data should contain at least two variables!")
  if(!is.null(data) && nrow(data) < 4) stop("for clustering data should contain at least four objects!")
  
  if(!is.null(object) && class(object) != "kproto") stop("object must be type of 'kproto'")
  
  if(!is.null(k) && length(k) == 1){
    stop("k should be the search range for optimum number of clusters, e.g. c(2:sqrt(n))")
  }
  if(length(k) > 1){
    if(nrow(data) < max(k)) stop("Data frame has less observations than clusters!")
    if(any(k < 1) | any(k > nrow(data))) stop("Elements of k must be greater than 1 and strictly less than n!")
    if(all(!as.integer(k)==k)) stop("Elements of k must be type of integer")
  }
  
  if(!(kp_obj %in% c("optimal","all"))) stop("kp_obj must either be optimal or all!")
  
  output <- switch(method,
                   "cindex" = cindex_kproto(object = object, data = data, k = k, kp_obj = kp_obj, ...),
                   "dunn" = dunn_kproto(object = object, data = data, k = k, kp_obj = kp_obj, ...),
                   "gamma" = gamma_kproto(object = object, data = data, k = k, kp_obj = kp_obj, ...),
                   "gplus" = gplus_kproto(object = object, data = data, k = k, kp_obj = kp_obj, ...),
                   "mcclain" = mcclain_kproto(object = object, data = data, k = k, kp_obj = kp_obj, ...),
                   "ptbiserial" = ptbiserial_kproto(object = object, data = data, k = k, kp_obj = kp_obj, ...),
                   "silhouette" = silhouette_kproto(object = object, data = data, k = k, kp_obj = kp_obj, ...),
                   "tau" = tau_kproto(object = object, data = data, k = k, kp_obj = kp_obj, ...))
  
  return(output)
}











