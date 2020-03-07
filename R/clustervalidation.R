





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





#' @title Validating k Prototypes Clustering: Cindex
#'
#' @description Calculating the Cindex for a k-Prototypes clustering with k clusters or computing the optimal number of clusters based on the Cindex for k-Prototype clustering.
#' 
#' @details \deqn{Cindex = \frac{S_w-S_{min}}{S_{max}-S_{min}}} \cr
#' For \eqn{S_{min}} and \eqn{S_{max}} it is nessesary to calculate the distances between all pairs of points in the entire data set (\eqn{\frac{n(n-1)}{2}}). 
#' \eqn{S_{min}} is the sum of the "total number of pairs of objects belonging to the same cluster" smallest distances and 
#' \eqn{S_{max}} is the sum of the "total number of pairs of objects belonging to the same cluster" largest distances. \eqn{S_w} is the sum of the within-cluster distances. \cr
#' The minimum value of the index is used to indicate the optimal number of clusters.
#' 
#' 
#' @param object Object of class \code{kproto} resulting from a call with \code{kproto(..., keep.data=TRUE)}
#' @param data Original data; only required if \code{object == NULL}.
#' @param k Vector specifying the search range for optimum number of clusters; if \code{NULL} the range will set as \code{2:sqrt(n)}. Only required if \code{object == NULL}.
#' @param S_sort for internal purposes only
#' @param ... Further arguments passed to \code{\link[clustMixType]{kproto}}, like:
#'   \itemize{
#'     \item \code{nstart}: If > 1 repetetive computations of \code{kproto} with random initializations are computed.
#'     \item \code{lambda}: Factor to trade off between Euclidean distance of numeric variables and simple matching coefficient between categorical variables.
#'     \item \code{verbose}: Logical whether information about the cluster procedure should be given. Caution: If \code{verbose=FALSE}, the reduction of the number of clusters is not mentioned.
#'   }
#'
#' @return For computing the optimal number of clusters based on the Cindex for k-Prototype clustering the output contains:
#' @return \item{k_opt}{optimal number of clusters}
#' @return \item{indices}{calculated indices for \eqn{k=2,...,k_{max}}}
#' @return For computing the Cindex-value for a given k-Prototype clustering the output contains:
#' @return \item{index}{calculated index-value}
#'
#' @references \itemize{
#'     \item Charrad, M., Ghazzali, N., Boiteau, V., Niknafs, A. (2014): 
#'     NbClust: An R Package for Determining the Relevant Number of Clusters in a Data Set. 
#'     \emph{Journal of Statistical Software, Vol 61, Issue 6}.
#'     \href{http://www.jstatsoft.org/v61/i06/}{\emph{www.jstatsoft.org}}.
#'   }
#'
#' @seealso Other clustervalidation indices: \code{\link[clustMixType]{dunn_kproto}},
#' \code{\link[clustMixType]{gamma_kproto}}, \code{\link[clustMixType]{gplus_kproto}},
#' \code{\link[clustMixType]{mcclain_kproto}}, \code{\link[clustMixType]{ptbiserial_kproto}},
#' \code{\link[clustMixType]{silhouette_kproto}}, \code{\link[clustMixType]{tau_kproto}}
#'
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 10
#' prb <- 0.99
#' muk <- 2.5 
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' 
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' 
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' 
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' # apply k prototyps
#' kpres <- kproto(x, 4, keep.data=TRUE)
#' 
#' # calculate cindex-value
#' cindex_value <- cindex_kproto(object = kpres)
#' 
#' # calculate optimal number of cluster
#' k_opt <- cindex_kproto(data = x, k = 3:5, nstart = 5, verbose=FALSE)
#' 
#' @author Rabea Aschenbruck
#' 
#' @rdname cindex_kproto
#' @importFrom utils head
#' @importFrom utils tail
#' 
#' @export
cindex_kproto <- function(object = NULL, data = NULL, k = NULL, S_sort = NULL, ...){
  
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
    
    indices <- rep(NA,length(k))
    names(indices) <- k
    for(q in k){
      object <- kproto(x = data, k = q, keep.data = TRUE, ...)
      if(is.na(indices[length(object$size) == names(indices)])){
        indices[length(object$size)==names(indices)] <- cindex_kproto(object = object, S_sort = S_sort)
      }
    }
    
    k_opt <- as.numeric(names(sort(indices)[1]))
    
    return(list("k_opt" = k_opt, "indices" = indices))
  }
}





#' @title Validating k Prototypes Clustering: Dunn index
#'
#' @description Calculating the Dunn index for a k-Prototypes clustering with k clusters or 
#' computing the optimal number of clusters based on the Dunn index for k-Prototype clustering.
#' 
#' @details \deqn{Dunn = \frac{\min_{1 \leq i < j \leq q} d(C_i, C_j)}{\max_{1 \leq k \leq q} diam(C_k)}} \cr
#' The following applies: The dissimilarity between the two clusters \eqn{C_i} and \eqn{C_j} is defined as \eqn{d(C_i, C_j)=\min_{x \in C_i, y \in C_j} d(x,y)} and
#' the diameter of a cluster is defined as \eqn{diam(C_k)=\max_{x,y \in C} d(x,y)}. \cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#' 
#' @param object Object of class \code{kproto} resulting from a call with \code{kproto(..., keep.data=TRUE)}
#' @param data Original data; only required if \code{object == NULL}.
#' @param k Vector specifying the search range for optimum number of clusters; if \code{NULL} the range will set as \code{2:sqrt(n)}. Only required if \code{object == NULL}.
#' @param ... Further arguments passed to \code{\link[clustMixType]{kproto}}, like:
#'   \itemize{
#'     \item \code{nstart}: If > 1 repetetive computations of \code{kproto} with random initializations are computed.
#'     \item \code{lambda}: Factor to trade off between Euclidean distance of numeric variables and simple matching coefficient between categorical variables.
#'     \item \code{verbose}: Logical whether information about the cluster procedure should be given. Caution: If \code{verbose=FALSE}, the reduction of the number of clusters is not mentioned.
#'   }
#'
#' @return For computing the optimal number of clusters based on the Dunn index for k-Prototype clustering the output contains:
#' @return \item{k_opt}{optimal number of clusters}
#' @return \item{indices}{calculated indices for \eqn{k=2,...,k_{max}}}
#' @return For computing the Dunn index-value for a given k-Prototype clustering the output contains:
#' @return \item{index}{calculated index-value}
#'
#' @references \itemize{
#'     \item Charrad, M., Ghazzali, N., Boiteau, V., Niknafs, A. (2014): 
#'     NbClust: An R Package for Determining the Relevant Number of Clusters in a Data Set. 
#'     \href{http://www.jstatsoft.org/v61/i06/}{\emph{Journal of Statistical Software, Vol 61, Issue 6}}.
#'   }
#'
#' @seealso Other clustervalidation indices: \code{\link[clustMixType]{dunn_kproto}},
#' \code{\link[clustMixType]{gamma_kproto}}, \code{\link[clustMixType]{gplus_kproto}},
#' \code{\link[clustMixType]{mcclain_kproto}}, \code{\link[clustMixType]{ptbiserial_kproto}},
#' \code{\link[clustMixType]{silhouette_kproto}}, \code{\link[clustMixType]{tau_kproto}}
#'
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 10
#' prb <- 0.99
#' muk <- 2.5 
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' 
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' 
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' 
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' # apply k prototyps
#' kpres <- kproto(x, 4, keep.data=TRUE)
#' 
#' # calculate index-value
#' dunn_value <- dunn_kproto(object = kpres)
#' 
#' \dontrun{
#' # calculate optimal number of cluster
#' k_opt <- dunn_kproto(data = x, k = 3:5, nstart = 5, verbose = FALSE)
#' }
#' 
#' 
#' @author Rabea Aschenbruck
#' 
#' @rdname dunn_kproto
#' 
#' @export
dunn_kproto <- function(object = NULL, data = NULL, k = NULL, ...){
  
  if(is.null(data) && is.null(object)) stop("data or object muss be given!")
  
  if(!is.null(data) && !is.data.frame(data)) stop("data should be a data frame!")
  if(!is.null(data) && ncol(data) < 2) stop("For clustering data should contain at least two variables!")
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
    
    indices <- rep(NA,length(k))
    names(indices) <- k
    for(q in k){
      object <- kproto(x = data, k = q, ...)
      if(length(object$size)>1){
        if(is.na(indices[length(object$size) == names(indices)])){
          indices[length(object$size) == names(indices)] <- dunn_kproto(object = object)
        }
      }
    }
    
    k_opt <- as.numeric(names(sort(indices, decreasing = TRUE)[1]))
    
    return(list("k_opt" = k_opt, "indices" = indices))
  }
}




#' @title Validating k Prototypes Clustering: Gamma index
#'
#' @description Calculating the Gamma index for a k-Prototypes clustering with k clusters or 
#' computing the optimal number of clusters based on the Gamma index for k-Prototype clustering.
#' 
#' @details \deqn{Gamma = \frac{s(+)-s(-)}{s(+)+s(-)}} \cr 
#' Comparisons are made between all within-cluster dissimilarities and all between-cluster dissimilarities. 
#' \eqn{s(+)} is the number of concordant comparisons and \eqn{s(-)} is the number of discordant comparisons.
#' A comparison is named concordant (resp. discordant) if a within-cluster dissimilarity is strictly less (resp. strictly greater) than a between-cluster dissimilarity.\cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#' 
#' @param object Object of class \code{kproto} resulting from a call with \code{kproto(..., keep.data=TRUE)}
#' @param data Original data; only required if \code{object == NULL}.
#' @param k Vector specifying the search range for optimum number of clusters; if \code{NULL} the range will set as \code{2:sqrt(n)}. Only required if \code{object == NULL}.
#' @param dists for internal purposes only
#' @param ... Further arguments passed to \code{\link[clustMixType]{kproto}}, like:
#'   \itemize{
#'     \item \code{nstart}: If > 1 repetetive computations of \code{kproto} with random initializations are computed.
#'     \item \code{lambda}: Factor to trade off between Euclidean distance of numeric variables and simple matching coefficient between categorical variables.
#'     \item \code{verbose}: Logical whether information about the cluster procedure should be given. Caution: If \code{verbose=FALSE}, the reduction of the number of clusters is not mentioned.
#'   }
#'
#' @return For computing the optimal number of clusters based on the Gamma index for k-Prototype clustering the output contains:
#' @return \item{k_opt}{optimal number of clusters}
#' @return \item{indices}{calculated indices for \eqn{k=2,...,k_{max}}}
#' @return For computing the Gamma index-value for a given k-Prototype clustering the output contains:
#' @return \item{index}{calculated index-value}
#'
#' @references \itemize{
#'     \item Charrad, M., Ghazzali, N., Boiteau, V., Niknafs, A. (2014): 
#'     NbClust: An R Package for Determining the Relevant Number of Clusters in a Data Set. 
#'     \href{http://www.jstatsoft.org/v61/i06/}{\emph{Journal of Statistical Software, Vol 61, Issue 6}}.
#'   }
#'
#' @seealso Other clustervalidation indices: \code{\link[clustMixType]{dunn_kproto}},
#' \code{\link[clustMixType]{dunn_kproto}}, \code{\link[clustMixType]{gplus_kproto}},
#' \code{\link[clustMixType]{mcclain_kproto}}, \code{\link[clustMixType]{ptbiserial_kproto}},
#' \code{\link[clustMixType]{silhouette_kproto}}, \code{\link[clustMixType]{tau_kproto}}
#' 
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 10
#' prb <- 0.99
#' muk <- 2.5 
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' 
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' 
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' 
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' # apply k prototyps
#' kpres <- kproto(x, 4, keep.data=TRUE)
#' 
#' # calculate index-value
#' gamma_value <- gamma_kproto(object = kpres)
#' 
#' # calculate optimal number of cluster
#' k_opt <- gamma_kproto(data = x, k = 3:5, nstart = 5, verbose = FALSE)
#' 
#' @author Rabea Aschenbruck
#' 
#' @rdname gamma_kproto
#' @importFrom stats na.omit
#' 
#' @export
gamma_kproto <- function(object = NULL, data = NULL, k = NULL, dists = NULL, ...){
  
  if(is.null(data) && is.null(object)) stop("data or object muss be given!")
  
  if(!is.null(data) && !is.data.frame(data)) stop("data should be a data frame!")
  if(!is.null(data) && ncol(data) < 2) stop("For clustering data should contain at least two variables!")
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
    
    if(is.null(k)){k <- 2:sqrt(n)}
    
    indices <- rep(NA,length(k))
    names(indices) <- k
    for(q in k){
      object <- kproto(x = data, k = q, ...)
      if(is.na(indices[length(object$size) == names(indices)])){
        indices[length(object$size) == names(indices)] <- gamma_kproto(object = object, dists = dists)
      }
    }
    
    k_opt <- as.numeric(names(sort(indices, decreasing = TRUE)[1]))
  
    return(list("k_opt" = k_opt, "indices" = indices))
  }
}




#' @title Validating k Prototypes Clustering: Gplus index
#'
#' @description Calculating the Gplus index for a k-Prototypes clustering with k clusters or 
#' computing the optimal number of clusters based on the Gplus index for k-Prototype clustering.
#' 
#' @details \deqn{Gplus = \frac{2 \cdot s(-)}{\frac{n(n-1)}{2} \cdot (\frac{n(n-1)}{2}-1)}} \cr 
#' Comparisons are made between all within-cluster dissimilarities and all between-cluster dissimilarities. 
#' \eqn{s(-)} is the number of discordant comparisons and a comparison is named discordant if a within-cluster 
#' dissimilarity is strictly greater than a between-cluster dissimilarity. \cr
#' The minimum value of the index is used to indicate the optimal number of clusters.
#' 
#' @param object Object of class \code{kproto} resulting from a call with \code{kproto(..., keep.data=TRUE)}
#' @param data Original data; only required if \code{object == NULL}.
#' @param k Vector specifying the search range for optimum number of clusters; if \code{NULL} the range will set as \code{2:sqrt(n)}. Only required if \code{object == NULL}.
#' @param dists for internal purposes only
#' @param ... Further arguments passed to \code{\link[clustMixType]{kproto}}, like:
#'   \itemize{
#'     \item \code{nstart}: If > 1 repetetive computations of \code{kproto} with random initializations are computed.
#'     \item \code{lambda}: Factor to trade off between Euclidean distance of numeric variables and simple matching coefficient between categorical variables.
#'     \item \code{verbose}: Logical whether information about the cluster procedure should be given. Caution: If \code{verbose=FALSE}, the reduction of the number of clusters is not mentioned.
#'   }
#'
#' @return For computing the optimal number of clusters based on the Gplus index for k-Prototype clustering the output contains:
#' @return \item{k_opt}{optimal number of clusters}
#' @return \item{indices}{calculated indices for \eqn{k=2,...,k_{max}}}
#' @return For computing the Gplus index-value for a given k-Prototype clustering the output contains:
#' @return \item{index}{calculated index-value}
#'
#' @references \itemize{
#'     \item Charrad, M., Ghazzali, N., Boiteau, V., Niknafs, A. (2014): 
#'     NbClust: An R Package for Determining the Relevant Number of Clusters in a Data Set. 
#'     \href{http://www.jstatsoft.org/v61/i06/}{\emph{Journal of Statistical Software, Vol 61, Issue 6}}.
#'   }
#'
#' @seealso Other clustervalidation indices: \code{\link[clustMixType]{dunn_kproto}},
#' \code{\link[clustMixType]{dunn_kproto}}, \code{\link[clustMixType]{gamma_kproto}},
#' \code{\link[clustMixType]{mcclain_kproto}}, \code{\link[clustMixType]{ptbiserial_kproto}},
#' \code{\link[clustMixType]{silhouette_kproto}}, \code{\link[clustMixType]{tau_kproto}}
#' 
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 10
#' prb <- 0.99
#' muk <- 2.5
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' 
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' 
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' 
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' # apply k prototyps
#' kpres <- kproto(x, 4, keep.data=TRUE)
#' 
#' # calculate index-value
#' gplus_value <- gplus_kproto(object = kpres)
#' 
#' # calculate optimal number of cluster
#' k_opt <- gplus_kproto(data = x, k = 3:5, nstart = 5, verbose = FALSE)
#' 
#' @author Rabea Aschenbruck
#' 
#' @rdname gplus_kproto
#' @importFrom stats na.omit
#' 
#' @export
gplus_kproto <- function(object = NULL, data = NULL, k = NULL, dists = NULL, ...){
  
  if(is.null(data) && is.null(object)) stop("data or object muss be given!")
  
  if(!is.null(data) && !is.data.frame(data)) stop("data should be a data frame!")
  if(!is.null(data) && ncol(data) < 2) stop("For clustering data should contain at least two variables!")
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
    
    indices <- rep(NA,length(k))
    names(indices) <- k
    for(q in k){
      object <- kproto(x = data, k = q, ...)
      if(is.na(indices[length(object$size) == names(indices)])){
        indices[length(object$size) == names(indices)] <- gplus_kproto(object = object, dists = dists)
      }
    }
    
    k_opt <- as.numeric(names(sort(indices)[1]))
    
    return(list("k_opt" = k_opt, "indices" = indices))
  }
}





#' @title Validating k Prototypes Clustering: McClain index
#'
#' @description Calculating the McClain index for a k-Prototypes clustering with k clusters or 
#' computing the optimal number of clusters based on the McClain index for k-Prototype clustering.
#' 
#' @details \deqn{McClain = \frac{\bar{S}_w}{\bar{S}_b}} \cr 
#' \eqn{\bar{S}_w} is the sum of within-cluster distances divided by the number of within-cluster distances and 
#' \eqn{\bar{S}_b} is the sum of between-cluster distances divided by the number of between-cluster distances.\cr
#' The minimum value of the index is used to indicate the optimal number of clusters.
#' 
#' @param object Object of class \code{kproto} resulting from a call with \code{kproto(..., keep.data=TRUE)}
#' @param data Original data; only required if \code{object == NULL}.
#' @param k Vector specifying the search range for optimum number of clusters; if \code{NULL} the range will set as \code{2:sqrt(n)}. Only required if \code{object == NULL}.
#' @param ... Further arguments passed to \code{\link[clustMixType]{kproto}}, like:
#'   \itemize{
#'     \item \code{nstart}: If > 1 repetetive computations of \code{kproto} with random initializations are computed.
#'     \item \code{lambda}: Factor to trade off between Euclidean distance of numeric variables and simple matching coefficient between categorical variables.
#'     \item \code{verbose}: Logical whether information about the cluster procedure should be given. Caution: If \code{verbose=FALSE}, the reduction of the number of clusters is not mentioned.
#'   }
#'
#' @return For computing the optimal number of clusters based on the McClain index for k-Prototype clustering the output contains:
#' @return \item{k_opt}{optimal number of clusters}
#' @return \item{indices}{calculated indices for \eqn{k=2,...,k_{max}}}
#' @return For computing the McClain index-value for a given k-Prototype clustering the output contains:
#' @return \item{index}{calculated index-value}
#'
#' @references \itemize{
#'     \item Charrad, M., Ghazzali, N., Boiteau, V., Niknafs, A. (2014): 
#'     NbClust: An R Package for Determining the Relevant Number of Clusters in a Data Set. 
#'     \href{http://www.jstatsoft.org/v61/i06/}{\emph{Journal of Statistical Software, Vol 61, Issue 6}}.
#'   }
#'
#' @seealso Other clustervalidation indices: \code{\link[clustMixType]{dunn_kproto}},
#' \code{\link[clustMixType]{dunn_kproto}}, \code{\link[clustMixType]{gamma_kproto}},
#' \code{\link[clustMixType]{gplus_kproto}}, \code{\link[clustMixType]{ptbiserial_kproto}},
#' \code{\link[clustMixType]{silhouette_kproto}}, \code{\link[clustMixType]{tau_kproto}}
#'
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 10
#' prb <- 0.99
#' muk <- 2.5
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' 
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' 
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' 
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' # apply k prototyps
#' kpres <- kproto(x, 4, keep.data=TRUE)
#' 
#' # calculate index-value
#' mcclain_value <- mcclain_kproto(object = kpres)
#' 
#' # calculate optimal number of cluster
#' k_opt <- mcclain_kproto(data = x, k = 3:5, nstart = 5, verbose = FALSE)
#' 
#' @author Rabea Aschenbruck
#' 
#' @rdname mcclain_kproto
#' 
#' @export
mcclain_kproto <- function(object = NULL, data = NULL, k = NULL, ...){
  
  if(is.null(data) && is.null(object)) stop("data or object muss be given!")
  
  if(!is.null(data) && !is.data.frame(data)) stop("data should be a data frame!")
  if(!is.null(data) && ncol(data) < 2) stop("For clustering data should contain at least two variables!")
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
    
    indices <- rep(NA,length(k))
    names(indices) <- k
    for(q in k){
      object <- kproto(x = data, k = q, ...)
      if(is.na(indices[length(object$size) == names(indices)])){
        indices[length(object$size) == names(indices)] <- mcclain_kproto(object = object)
      }
    }
    
    k_opt <- as.numeric(names(sort(indices)[1]))
    
    return(list("k_opt" = k_opt, "indices" = indices))
  }
}




#' @title Validating k Prototypes Clustering: Ptbiserial index
#'
#' @description Calculating the Ptbiserial index for a k-Prototypes clustering with k clusters or 
#' computing the optimal number of clusters based on the Ptbiserial index for k-Prototype clustering.
#' 
#' @details \deqn{Ptbiserial = \frac{(\bar{S}_b-\bar{S}_w) \cdot (\frac{N_w \cdot N_b}{N_t^2})^{0.5}}{s_d}} \cr 
#' \eqn{\bar{S}_w} is the sum of within-cluster distances divided by the number of within-cluster distances and 
#' \eqn{\bar{S}_b} is the sum of between-cluster distances divided by the number of between-cluster distances.\cr
#' \eqn{N_t} is the total number of pairs of objects in the data, \eqn{N_w} is the total number of pairs of 
#' objects belonging to the samecluster and \eqn{N_b} is the total number of pairs of objects belonging to different clusters.
#' \eqn{s_d} is the standard deviation of all distances.\cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#' 
#' @param object Object of class \code{kproto} resulting from a call with \code{kproto(..., keep.data=TRUE)}
#' @param data Original data; only required if \code{object == NULL}.
#' @param k Vector specifying the search range for optimum number of clusters; if \code{NULL} the range will set as \code{2:sqrt(n)}. Only required if \code{object == NULL}.
#' @param s_d for internal purposes only
#' @param ... Further arguments passed to \code{\link[clustMixType]{kproto}}, like:
#'   \itemize{
#'     \item \code{nstart}: If > 1 repetetive computations of \code{kproto} with random initializations are computed.
#'     \item \code{lambda}: Factor to trade off between Euclidean distance of numeric variables and simple matching coefficient between categorical variables.
#'     \item \code{verbose}: Logical whether information about the cluster procedure should be given. Caution: If \code{verbose=FALSE}, the reduction of the number of clusters is not mentioned.
#'   }
#'
#' @return For computing the optimal number of clusters based on the Ptbiserial index for k-Prototype clustering the output contains:
#' @return \item{k_opt}{optimal number of clusters}
#' @return \item{indices}{calculated indices for \eqn{k=2,...,k_{max}}}
#' @return For computing the Ptbiserial index-value for a given k-Prototype clustering the output contains:
#' @return \item{index}{calculated index-value}
#'
#' @references \itemize{
#'     \item Charrad, M., Ghazzali, N., Boiteau, V., Niknafs, A. (2014): 
#'     NbClust: An R Package for Determining the Relevant Number of Clusters in a Data Set. 
#'     \href{http://www.jstatsoft.org/v61/i06/}{\emph{Journal of Statistical Software, Vol 61, Issue 6}}.
#'   }
#'
#' @seealso Other clustervalidation indices: \code{\link[clustMixType]{dunn_kproto}},
#' \code{\link[clustMixType]{dunn_kproto}}, \code{\link[clustMixType]{gamma_kproto}},
#' \code{\link[clustMixType]{gplus_kproto}}, \code{\link[clustMixType]{mcclain_kproto}},
#' \code{\link[clustMixType]{silhouette_kproto}}, \code{\link[clustMixType]{tau_kproto}}
#' 
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 10
#' prb <- 0.99
#' muk <- 2.5
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' 
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' 
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' 
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' # apply k prototyps
#' kpres <- kproto(x, 4, keep.data=TRUE)
#' 
#' # calculate index-value
#' Ptbiserial_value <- ptbiserial_kproto(object = kpres)
#' 
#' # calculate optimal number of cluster
#' k_opt <- ptbiserial_kproto(data = x, k = 3:5, nstart = 5, verbose = FALSE)
#' 
#' @author Rabea Aschenbruck
#' 
#' @rdname ptbiserial_kproto
#' 
#' @export
ptbiserial_kproto <- function(object = NULL, data = NULL, k = NULL, s_d = NULL, ...){
  
  if(is.null(data) && is.null(object)) stop("data or object muss be given!")
  
  if(!is.null(data) && !is.data.frame(data)) stop("data should be a data frame!")
  if(!is.null(data) && ncol(data) < 2) stop("For clustering data should contain at least two variables!")
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
    
    indices <- rep(NA,length(k))
    names(indices) <- k
    for(q in k){
      object <- kproto(x = data, k=q, ...)
      if(is.na(indices[length(object$size) == names(indices)])){
        indices[length(object$size) == names(indices)] <- ptbiserial_kproto(object = object, s_d = s_d)
      }
    }
    
    k_opt <- as.numeric(names(sort(indices, decreasing = TRUE)[1]))
    
    return(list("k_opt" = k_opt, "indices" = indices))
  }
}




#' @title Validating k Prototypes Clustering: Silhouette index
#'
#' @description Calculating the Silhouette index for a k-Prototypes clustering with k clusters or 
#' computing the optimal number of clusters based on the Silhouette index for k-Prototype clustering.
#' 
#' @details \deqn{Silhouette = \frac{1}{n} \sum_{i=1}^n \frac{b(i)-a(i)}{max(a(i),b(i))}} \cr 
#' \eqn{a(i)} is the average dissimilarity of the i\emph{th} object to all other objects of the same/own cluster.
#' \eqn{b(i)=min(d(i,C))}, where \eqn{d(i,C)} is the average dissimilarity of the i\emph{th} object to all the other clusters except the own/same cluster.\cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#' 
#' @param object Object of class \code{kproto} resulting from a call with \code{kproto(..., keep.data=TRUE)}
#' @param data Original data; only required if \code{object == NULL}.
#' @param k Vector specifying the search range for optimum number of clusters; if \code{NULL} the range will set as \code{2:sqrt(n)}. Only required if \code{object == NULL}.
#' @param ... Further arguments passed to \code{\link[clustMixType]{kproto}}, like:
#'   \itemize{
#'     \item \code{nstart}: If > 1 repetetive computations of \code{kproto} with random initializations are computed.
#'     \item \code{lambda}: Factor to trade off between Euclidean distance of numeric variables and simple matching coefficient between categorical variables.
#'     \item \code{verbose}: Logical whether information about the cluster procedure should be given. Caution: If \code{verbose=FALSE}, the reduction of the number of clusters is not mentioned.
#'   }
#'
#' @return For computing the optimal number of clusters based on the Silhouette index for k-Prototype clustering the output contains:
#' @return \item{k_opt}{optimal number of clusters}
#' @return \item{indices}{calculated indices for \eqn{k=2,...,k_{max}}}
#' @return For computing the Silhouette index-value for a given k-Prototype clustering the output contains:
#' @return \item{index}{calculated index-value}
#'
#' @references \itemize{
#'     \item Charrad, M., Ghazzali, N., Boiteau, V., Niknafs, A. (2014): 
#'     NbClust: An R Package for Determining the Relevant Number of Clusters in a Data Set. 
#'     \href{http://www.jstatsoft.org/v61/i06/}{\emph{Journal of Statistical Software, Vol 61, Issue 6}}.
#'   }
#'
#' @seealso Other clustervalidation indices: \code{\link[clustMixType]{dunn_kproto}},
#' \code{\link[clustMixType]{dunn_kproto}}, \code{\link[clustMixType]{gamma_kproto}},
#' \code{\link[clustMixType]{gplus_kproto}}, \code{\link[clustMixType]{mcclain_kproto}},
#' \code{\link[clustMixType]{ptbiserial_kproto}}, \code{\link[clustMixType]{tau_kproto}}
#'
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 10
#' prb <- 0.99
#' muk <- 2.5
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' 
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' 
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' 
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' # apply k prototyps
#' kpres <- kproto(x, 4, keep.data=TRUE)
#' 
#' # calculate index-value
#' silhouette_value <- silhouette_kproto(object = kpres)
#' 
#' # calculate optimal number of cluster
#' k_opt <- silhouette_kproto(data = x, k = 3:5, nstart = 5, verbose = FALSE)
#' 
#' @author Rabea Aschenbruck
#' 
#' @rdname silhouette_kproto
#' 
#' @export
silhouette_kproto <- function(object = NULL, data = NULL, k = NULL, ...){
  
  if(is.null(data) && is.null(object)) stop("data or object muss be given!")
  
  if(!is.null(data) && !is.data.frame(data)) stop("data should be a data frame!")
  if(!is.null(data) && ncol(data) < 2) stop("For clustering data should contain at least two variables!")
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
    
    indices <- rep(NA,length(k))
    names(indices) <- k
    for(q in k){
      object <- kproto(x = data, k = q, ...)
      if(is.na(indices[length(object$size) == names(indices)])){
        indices[length(object$size) == names(indices)] <- silhouette_kproto(object = object)
      }
    }
    
    k_opt <- as.numeric(names(sort(indices, decreasing = TRUE)[1]))
    
    return(list("k_opt" = k_opt, "indices" = indices))
  }
}




#' @title Validating k Prototypes Clustering: Tau index
#'
#' @description Calculating the Tau index for a k-Prototypes clustering with k clusters or 
#' computing the optimal number of clusters based on the Tau index for k-Prototype clustering.
#' 
#' @details \deqn{Tau = \frac{s(+) - s(-)}{((\frac{N_t(N_t-1)}{2}-t)\frac{N_t(N_t-1)}{2})^{0.5}}} \cr 
#' Comparisons are made between all within-cluster dissimilarities and all between-cluster dissimilarities. 
#' \eqn{s(+)} is the number of concordant comparisons and \eqn{s(-)} is the number of discordant comparisons.
#' A comparison is named concordant (resp. discordant) if a within-cluster dissimilarity is strictly less 
#' (resp. strictly greater) than a between-cluster dissimilarity.\cr
#' \eqn{N_t} is the total number of distances \eqn{\frac{n(n-1)}{2}} and \eqn{t} is the number of comparisons 
#' of two pairs of objects where both pairs represent within-cluster comparisons or both pairs are between-cluster
#' comparisons. \cr
#' The maximum value of the index is used to indicate the optimal number of clusters.
#' 
#' @param object Object of class \code{kproto} resulting from a call with \code{kproto(..., keep.data=TRUE)}
#' @param data Original data; only required if \code{object == NULL}.
#' @param k Vector specifying the search range for optimum number of clusters; if \code{NULL} the range will set as \code{2:sqrt(n)}. Only required if \code{object == NULL}.
#' @param dists for internal purposes only
#' @param ... Further arguments passed to \code{\link[clustMixType]{kproto}}, like:
#'   \itemize{
#'     \item \code{nstart}: If > 1 repetetive computations of \code{kproto} with random initializations are computed.
#'     \item \code{lambda}: Factor to trade off between Euclidean distance of numeric variables and simple matching coefficient between categorical variables.
#'     \item \code{verbose}: Logical whether information about the cluster procedure should be given. Caution: If \code{verbose=FALSE}, the reduction of the number of clusters is not mentioned.
#'   }
#'
#' @return For computing the optimal number of clusters based on the Tau index for k-Prototype clustering the output contains:
#' @return \item{k_opt}{optimal number of clusters}
#' @return \item{indices}{calculated indices for \eqn{k=2,...,k_{max}}}
#' @return For computing the Tau index-value for a given k-Prototype clustering the output contains:
#' @return \item{index}{calculated index-value}
#'
#' @references \itemize{
#'     \item Charrad, M., Ghazzali, N., Boiteau, V., Niknafs, A. (2014): 
#'     NbClust: An R Package for Determining the Relevant Number of Clusters in a Data Set. 
#'     \href{http://www.jstatsoft.org/v61/i06/}{\emph{Journal of Statistical Software, Vol 61, Issue 6}}.
#'   }
#'
#' @seealso Other clustervalidation indices: \code{\link[clustMixType]{dunn_kproto}},
#' \code{\link[clustMixType]{dunn_kproto}}, \code{\link[clustMixType]{gamma_kproto}},
#' \code{\link[clustMixType]{gplus_kproto}}, \code{\link[clustMixType]{mcclain_kproto}},
#' \code{\link[clustMixType]{ptbiserial_kproto}}, \code{\link[clustMixType]{silhouette_kproto}}
#' 
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 10
#' prb <- 0.99
#' muk <- 2.5
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' 
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' 
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' 
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' # apply k prototyps
#' kpres <- kproto(x, 4, keep.data=TRUE)
#' 
#' # calculate index-value
#' tau_value <- tau_kproto(object = kpres)
#' 
#' # calculate optimal number of cluster
#' k_opt <- tau_kproto(data = x, k = 3:5, nstart = 5, verbose = FALSE)
#' 
#' @author Rabea Aschenbruck
#' 
#' @rdname tau_kproto
#' @importFrom stats na.omit
#' @importFrom utils combn
#' 
#' @export
tau_kproto <- function(object = NULL, data = NULL, k = NULL, dists = NULL, ...){
  
  if(is.null(data) && is.null(object)) stop("data or object muss be given!")
  
  if(!is.null(data) && !is.data.frame(data)) stop("data should be a data frame!")
  if(!is.null(data) && ncol(data) < 2) stop("For clustering data should contain at least two variables!")
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
    
    indices <- rep(NA,length(k))
    names(indices) <- k
    for(q in k){
      object <- kproto(x = data, k = q, ...)
      if(is.na(indices[length(object$size) == names(indices)])){
        indices[length(object$size) == names(indices)] <- tau_kproto(object = object, dists = dists)
      }
    }
    
    k_opt <- as.numeric(names(sort(indices, decreasing = TRUE)[1]))
    
    return(list("k_opt" = k_opt, "indices" = indices))
  }
}



















