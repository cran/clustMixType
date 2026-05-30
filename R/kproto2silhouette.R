###############################################
# According to a suggestion from Sam Sievertsen
###############################################

# add unit test

#' @title Compute individual Silhouette Widths
#' @description Computes an object of class \code{silhouette} from k-prototypes clustering for mixed-type data.
#' 
#' @details Information on the individual silhouette width for all observations is computed. 
#' The k-prototypes result of class \code{kproto} is turned into an object of class \code{silhouette} 
#' that can be further used with the functions \code{\link[cluster]{summary.silhouette}} or \code{\link[cluster]{plot.silhouette}}. 
#' 
#' @param object Object of class \code{kproto}.
#' 
#' @return Object of class \code{silhouette} to be further used with \code{\link[cluster]{summary.silhouette}} 
#' or \code{\link[cluster]{plot.silhouette}}. It consists in a matrix where each row corresponds to one observation 
#' and the columns contain the assigned cluster as well as the neighbor cluster as well as the individual silhouette 
#' width ob the observation.
#' 
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 50
#' prb <- 0.9
#' muk <- 1.5 
#' clusid <- rep(1:4, each = n)
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
#' # apply k-prototypes
#' kpres <- kproto(x, 4)
#' 
#' # compute silhouette widhts
#' # extract object of class sil
#' sils <- kproto2silhouette(kpres)
#' 
#' # apply functions from package cluster
#' library(cluster)
#' summary(sils)
#' plot(sils)
#' 
#' # compute average silhouette width  
#' mean(sils[,3])
#' 
#' # ...alternatively, compute average silhouette width directly using validation_kproto() 
#' cindex_value <- validation_kproto(method = "silhouette", object = kpres)
#' cindex_value
#' 
#' @author \email{gero.szepannek@@web.de}
#' 
#' 
#' @references \itemize{
#'     \item Szepannek, G. (2018): 
#'     clustMixType: User-Friendly Clustering of Mixed-Type Data in R, 
#'     {\emph{The R Journal 10/2}}, 200-208. 
#'     \doi{10.32614/RJ-2018-048}.
#'     
#'     \item Aschenbruck, R., Szepannek, G. (2020): 
#'     Cluster Validation for Mixed-Type Data. 
#'     \emph{Archives of Data Science, Series A, Vol 6, Issue 1}.
#'     \doi{10.5445/KSP/1000098011/02}.
#'     
#'     \item Kaufman, L., Rousseeuw, P. (1990): 
#'     \emph{Finding Groups in Data: An Introduction to Cluster Analysis.} 
#'      Wiley.
#'   }
#' @rdname kproto2silhouette
#' 
#' @importFrom cluster silhouette
#' 
#' @export
kproto2silhouette <- function(object){ 
  
  if(!(inherits(object, "kproto"))) stop("Object must be of class kproto!")
  type = object$type
  if(!("data" %in% names(object))) stop("Computation is not possible. The kproto object has to be vreated using the argument keep.data = TRUE!")
  data = object$data 
  k = length(object$withinss) 
  lambda = object$lambda 
  # kp_obj = "optimal", 
  # verbose = FALSE, ...
  # data_plus = NULL, 
  
  
  if(object$type == "gower"){  #...same as from validation_kproto.R ln. 1095--1166  
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
  
  
  #determine ai, bi and si as well as the neighbour cluster
  a <- numeric(nrow(object$data))
  b <- numeric(nrow(object$data))
  s <- numeric(nrow(object$data))
  neighbour <- numeric(nrow(object$data))
  
  for(i in 1:nrow(object$data)){
    if(is.na(object$cluster[i])){
      # special case: no cluster assignment cluster[i] (usually resulting of all variables NA in x[i])
      s[i] <- NA
    }else{
      a[i] <- cluster_dists[i, object$cluster[i]]
      b[i] <- min(cluster_dists[i, -object$cluster[i]])
      
      # computation of the neighbour is a bit complicated but hopefully also correct in case 
      # ...the cluster assignment has been sampled from two equally close clusters  
      neighbours.i <- order(cluster_dists[i,])
      neighbours.i <- neighbours.i[neighbours.i != object$cluster[i]]
      neighbour[i] <- neighbours.i[1]
      
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
  
  result <- cbind(object$cluster, neighbour, s)
  colnames(result) <- c("cluster", "neighbor", "sil_width")
  class(result) <- 'silhouette'
  attr(result, "Ordered") <- FALSE
  return(result)
}



# calc.dist <- function(type, lambda, data1, data2 = NULL, all.dists = FALSE, data_plus = NULL){
#   
#   if(type == "huang"){# calculation with standard distance of k-prototypes
#     
#     # if(is.null(lambda)){
#     #   lambda <- lambdaest(data, num.method = 1, fac.method = 1, outtype = "numeric", verbose = TRUE)
#     # }
#     
#     if(!is.null(data2)){#contains the input two datasets?
#       
#       #tbd: check if structure of data1 and data2 are equal
#       data <- rbind(data1, data2)
#       
#       if(nrow(data1) == 1 & nrow(data2) == 1){#calculation the distance between two observations (due to saving time)
#         
#         #numerical variables
#         d1 <- (data[1, data_plus$numvars, drop = FALSE]-data[2, data_plus$numvars, drop = FALSE])^2
#         d1[is.na(d1)] <- 0
#         if(length(lambda) == 1){
#           d1 <- sum(d1)
#         }else{ #length(lambda) > 1
#           d1 <- as.matrix(d1, nrow = 1) %*% matrix(lambda[data_plus$numvars], ncol = 1)
#         }
#         
#         #categorical variables
#         d2 <- sapply(which(data_plus$catvars), function(j) return(data[1,j] != data[2,j]))
#         d2[is.na(d2)] <- FALSE
#         if(length(lambda) == 1){
#           d2 <- lambda * sum(d2)
#         }else{#length(lambda) > 1
#           d2 <- matrix(d2, nrow = 1) %*% matrix(lambda[data_plus$catvars], ncol = 1)
#         }
#         
#         return(d1 + d2)
#       }else{#calculate distances between objects of two datasets
#         
#         #numerical variables
#         d1 <- (data1[rep(1:nrow(data1), each=nrow(data2)), data_plus$numvars, drop = FALSE] -
#                  data2[rep(1:nrow(data2), times=nrow(data1)), data_plus$numvars, drop = FALSE])^2
#         d1[is.na(d1)] <- 0
#         if(length(lambda) == 1){
#           d1 <- rowSums(d1)
#         }else{ #length(lambda) > 1
#           d1 <- as.matrix(d1) %*% lambda[data_plus$numvars]
#         }
#         
#         #categorical variables
#         d2 <- data1[rep(1:nrow(data1), each=nrow(data2)), data_plus$catvars, drop = FALSE] != 
#           data2[rep(1:nrow(data2), times = nrow(data1)), data_plus$catvars, drop = FALSE]
#         d2[is.na(d2)] <- FALSE
#         if(length(lambda) == 1){
#           d2 <- lambda * rowSums(d2)
#         }else{ #length(lambda) > 1
#           d2 <- as.matrix(d2) %*% lambda[data_plus$catvars]
#         }
#         
#         return(sum(d1 + d2))
#       }
#     }else{ #calculate distances between the objects of one dataset
#       
#       d1 <- (data1[rep(1:(nrow(data1)-1),times=(nrow(data1)-1):1, each=1), data_plus$numvars, drop = FALSE] -
#                data1[unlist(lapply(2:nrow(data1), seq, to=nrow(data1))), data_plus$numvars, drop = FALSE])^2
#       d1[is.na(d1)] <- 0
#       if(length(lambda) == 1) d1 <- rowSums(d1)
#       if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[data_plus$numvars]
#       
#       d2 <- data1[rep(1:(nrow(data1)-1),times=(nrow(data1)-1):1, each=1), data_plus$catvars, drop = FALSE] != 
#         data1[unlist(lapply(2:nrow(data1), seq, to=nrow(data1))), data_plus$catvars, drop = FALSE]
#       d2[is.na(d2)] <- FALSE
#       if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
#       if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[data_plus$catvars]
#       
#       if(all.dists == TRUE){
#         m <- matrix(rep(0,nrow(data1)^2), ncol = nrow(data1))
#         
#         m[lower.tri(m)] <- d1+d2
#         m <- t(m)
#         m[lower.tri(m)] <- d1+d2
#         return(m)
#       }else{
#         return(sum(d1 + d2))
#       }
#     }
#   }else{ #calculation with gower distance
#     
#     if(!is.null(data2)){#contains the input two datasets?
#       
#       #tbd: check if structure of data1 and data2 are equal
#       #data <- rbind(data1, data2)
#       
#       if(nrow(data1) == 1 & nrow(data2) == 1){#calculation the distance between two observations (due to saving time)
#         
#         d1 <- d2 <- d3 <- 0
#         
#         if(any(data_plus$numvars)){
#           d1 <- abs(data1[1, data_plus$numvars, drop = FALSE] - data2[1, data_plus$numvars, drop = FALSE])
#           for(jnum in 1:ncol(d1)) d1[,jnum] <- d1[,jnum] / data_plus$rgnums[jnum] 
#           d1[is.na(d1)] <- 0
#           if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[data_plus$numvars]
#           if(is.null(lambda)) d1 <- rowSums(d1)
#         }
#         
#         if(any(data_plus$catvars)){
#           d2 <- data1[1, data_plus$catvars, drop = FALSE] != data2[1, data_plus$catvars, drop = FALSE]
#           d2[is.na(d2)] <- FALSE
#           if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[data_plus$catvars]
#           if(is.null(lambda)) d2 <- rowSums(d2)
#         }
#         
#         if(any(data_plus$ordvars)){
#           d3 <- abs(data1[1, data_plus$ordvars, drop = FALSE] - data2[1, data_plus$ordvars, drop = FALSE])
#           for(jord in 1:ncol(d3)) d3[,jord] <- d3[,jord] / data_plus$rgords[jord] 
#           d3[is.na(d3)] <- 0
#           if(length(lambda) > 1) d3 <- as.matrix(d3) %*% lambda[data_plus$ordvars]
#           if(is.null(lambda)) d3 <- rowSums(d3)
#         }
#         
#         return(d1 + d2 + d3)
#         
#       }else{#calculate distances between objects of two datasets
#         
#         # in case of no numeric / factor / ordinal variables set:
#         d1 <- d2 <- d3 <- rep(0, nrow(data1)*nrow(data2))
#         
#         if(any(data_plus$numvars)){
#           d1 <- abs(data1[rep(1:nrow(data1), each=nrow(data2)), data_plus$numvars, drop = FALSE] - 
#                       data2[rep(1:nrow(data2), times=nrow(data1)), data_plus$numvars, drop = FALSE])
#           for(jnum in 1:ncol(d1)) d1[,jnum] <- d1[,jnum] / data_plus$rgnums[jnum] 
#           #time efficient alternative for codeline before:
#           # d1 <- d1 / rep(data_plus$rgnums, rep(nrow(d1), ncol(d1)))
#           d1[is.na(d1)] <- 0
#           if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[data_plus$numvars]
#           if(is.null(lambda)) d1 <- rowSums(d1)
#         }
#         
#         if(any(data_plus$catvars)){
#           d2 <- data1[rep(1:nrow(data1), each=nrow(data2)), data_plus$catvars, drop = FALSE] != 
#             data2[rep(1:nrow(data2), times = nrow(data1)), data_plus$catvars, drop = FALSE]
#           d2[is.na(d2)] <- FALSE
#           if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[data_plus$catvars]
#           if(is.null(lambda)) d2 <- rowSums(d2)
#         }
#         
#         if(any(data_plus$ordvars)){
#           d3 <- abs(data1[rep(1:nrow(data1), each=nrow(data2)), data_plus$ordvars, drop = FALSE] - 
#                       data2[rep(1:nrow(data2), times=nrow(data1)), data_plus$ordvars, drop = FALSE])
#           for(jord in 1:ncol(d3)) d3[,jord] <- d3[,jord] / data_plus$rgords[jord] 
#           #time efficient alternative for codeline before:
#           # d3 <- d3 / rep(data_plus$rgords, rep(nrow(d3), ncol(d3)))
#           d3[is.na(d3)] <- 0
#           if(length(lambda) > 1) d3 <- as.matrix(d3) %*% lambda[data_plus$ordvars]
#           if(is.null(lambda)) d3 <- rowSums(d3)
#         }
#         
#         return(sum(d1 + d2 + d3))
#       }
#     }else{ #calculate distances between the objects of one dataset
#       
#       d1 <- d2 <- d3 <- rep(0,length(rep(1:(nrow(data1)-1),times=(nrow(data1)-1):1, each=1)))
#       
#       if(any(data_plus$numvars)){
#         d1 <- abs(data1[rep(1:(nrow(data1)-1),times=(nrow(data1)-1):1, each=1), data_plus$numvars, drop = FALSE] - 
#                     data1[unlist(lapply(2:nrow(data1), seq, to=nrow(data1))), data_plus$numvars, drop = FALSE])
#         for(jnum in 1:ncol(d1)) d1[,jnum] <- d1[,jnum] / data_plus$rgnums[jnum] 
#         #time efficient alternative for codeline before:
#         # d1 <- d1 / rep(rgnums, rep(nrow(d1), ncol(d1)))
#         d1[is.na(d1)] <- 0
#         if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[data_plus$numvars]
#         if(is.null(lambda)) d1 <- rowSums(d1)
#       }
#       
#       if(any(data_plus$catvars)){
#         d2 <- data1[rep(1:(nrow(data1)-1),times=(nrow(data1)-1):1, each=1), data_plus$catvars, drop = FALSE] != 
#           data1[unlist(lapply(2:nrow(data1), seq, to=nrow(data1))), data_plus$catvars, drop = FALSE]
#         d2[is.na(d2)] <- FALSE
#         if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[data_plus$catvars]
#         if(is.null(lambda)) d2 <- rowSums(d2)
#       }
#       
#       if(any(data_plus$ordvars)){
#         d3 <- abs(data1[rep(1:(nrow(data1)-1),times=(nrow(data1)-1):1, each=1), data_plus$ordvars, drop = FALSE] - 
#                     data1[unlist(lapply(2:nrow(data1), seq, to=nrow(data1))), data_plus$ordvars, drop = FALSE])
#         for(jord in 1:ncol(d3)) d3[,jord] <- d3[,jord] / data_plus$rgords[jord] 
#         d3[is.na(d3)] <- 0
#         if(length(lambda) > 1) d3 <- as.matrix(d3) %*% lambda[data_plus$ordvars]
#         if(is.null(lambda)) d3 <- rowSums(d3)
#       }
#       
#       if(all.dists == TRUE){ #return distance matrix
#         m <- matrix(rep(0,nrow(data1)^2), ncol = nrow(data1))
#         
#         m[lower.tri(m)] <- d1+d2+d3
#         m <- t(m)
#         m[lower.tri(m)] <- d1+d2+d3
#         return(m)
#       }else{ #return the sum of all dists between objects
#         return(sum(d1 + d2 + d3))
#       }
#       
#     }
#   }
# }






# #####################################################################
# ### determination of distances, required for initialization functions
# dists_kproto <- function(x, y = NULL, lambda = NULL, verbose = FALSE){
#   
#   if(is.null(y)){
#     dat_part1 <- x[rep(c(1:nrow(x)), each = nrow(x)),]
#     dat_part2 <- x[rep(c(1:nrow(x)), times = nrow(x)),]
#   }else{
#     if(nrow(x) == 1 & nrow(y) == 1){
#       return(cbind(x, y, dist = dist_kproto(x,y,lambda = lambda)))
#     }else{
#       dat_part1 <- x[rep(c(1:nrow(x)), each = nrow(y)),]
#       dat_part2 <- y[rep(c(1:nrow(y)), times = nrow(x)),]
#     }
#   }
#   
#   # check for numeric and factor variables
#   numvars <- sapply(x, is.numeric)
#   anynum <- any(numvars)
#   catvars <- sapply(x, is.factor)
#   anyfact <- any(catvars)
#   
#   # determination of lambda
#   if(length(lambda) > 1) {if(length(lambda) != sum(c(numvars,catvars))) stop("If lambda is a vector, its length should be the sum of numeric and factor variables in the data frame!")}
#   if(is.null(lambda)){
#     if(anynum & anyfact){
#       vnum <- mean(sapply(x[,numvars, drop = FALSE], var, na.rm = TRUE))
#       vcat <- mean(sapply(x[,catvars, drop = FALSE], function(z) return(1-sum((table(z)/sum(!is.na(z)))^2))))
#       if (vnum == 0){
#         if(verbose) warning("All numerical variables have zero variance.")
#         anynum <- FALSE
#       } 
#       if (vcat == 0){
#         if(verbose) warning("All categorical variables have zero variance.")
#         anyfact <- FALSE
#       } 
#       if(anynum & anyfact){
#         lambda <- vnum/vcat
#         if(verbose) cat("Estimated lambda:", lambda, "\n\n")
#       }else{
#         lambda <- 1
#       }
#     }
#   }
#   
#   # compute distances 
#   nrows <- nrow(x)
#   d1 <- (dat_part1[,numvars, drop = FALSE] - dat_part2[,numvars, drop = FALSE])^2
#   if(length(lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
#   if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
#   d2 <- sapply(which(catvars), function(j) return(dat_part1[,j] != dat_part2[,j]))
#   d2[is.na(d2)] <- FALSE
#   if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
#   if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
#   
#   return(cbind(dat_part1, dat_part2, dist = as.vector(d1 + d2)))
# }
# 
# 
# 
# dist_kproto <- function(x, y, lambda, verbose = FALSE){
#   
#   x <- rbind(x,y)
#   
#   # check for numeric and factor variables
#   numvars <- sapply(x, is.numeric)
#   anynum <- any(numvars)
#   catvars <- sapply(x, is.factor)
#   anyfact <- any(catvars)
#   
#   # compute distances 
#   nrows <- nrow(x)
#   d1 <- (x[1, numvars, drop = FALSE] - x[2, numvars, drop = FALSE])^2
#   if(length(lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
#   if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
#   d2 <- sapply(which(catvars), function(j) return(x[1,j] != x[2,j]))
#   d2[is.na(d2)] <- FALSE
#   if(length(lambda) == 1) d2 <- lambda * rowSums(matrix(d2, nrow = 1))
#   if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
#   
#   return(as.numeric(d1 + d2))
# }
# 