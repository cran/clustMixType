# determination of the number of cluster based on stability values
# since the research showed no superior advantage over the internal validation indices,
#    there is no export of this function
# (cf. paper: Aschenbruck, Szepannek, Wilhelm: 
#             Stability of mixed-type cluster partitions for determination of the number of clusters,
#             proceedings of the IFCS 2022)
# 
# Input: - data    - data.frame
#        - k       - numeric vector, search range
#        - method  - character, one of these: luxburg, fowlkesmallows, rand, jaccard
#        - B       - numeric, number of bootstrap samples
#        - lambda  - numeric, weighting factor of k-Prototypes
#        - ...     - Further arguments passed to \code{\link[clustMixType]{kproto}}
# 
# Output: the stability based number of cluster and additionally the decison based on the adjusted values
#         - k_opt      - optimal number of clusters (sampled in case of ambiguity)
#         - index_opt  - stability value of the optimal clustering
#         - indices    - calculated stability for \eqn{k=2,...,k_{max}}
#         - kp_obj     - the kproto object of the optimal clustering

stability_det_k <- function(data, k, method, B = 100, lambda = NULL, ...){
  
  if(length(method) > 1)stop("input has to be exactly one of luxburg, fowlkesmallows, rand or jaccard")
  if(!(method %in% c("luxburg", "fowlkesmallows", "rand", "jaccard"))) stop("input has to be one of luxburg, fowlkesmallows, rand or jaccard")
  
  if(!is.null(k) && length(k) == 1){
    stop("k should be the search range for optimum number of clusters, e.g. c(2:sqrt(n))")
  }
  if(length(k) > 1){
    if(nrow(data) < max(k)) stop("Data frame has less observations than clusters!")
    if(any(k <= 1) | any(k >= nrow(data))) stop("Elements of k must be greater than 1 and strictly less than n!")
    if(all(!as.integer(k)==k)) stop("Elements of k must be type of integer")
  }
  
  n <- nrow(data)
  p <- ncol(data)
  
  if(is.null(lambda)){
    numvars <- sapply(data, is.numeric)
    anynum <- any(numvars)
    catvars <- sapply(data, is.factor)
    anyfact <- any(catvars)
    vnum <- mean(sapply(data[,numvars, drop = FALSE], var, na.rm = TRUE))
    vcat <- mean(sapply(data[,catvars, drop = FALSE], function(z) return(1-sum((table(z)/sum(!is.na(z)))^2))))
    lambda <- vnum/vcat
  }
  
  #calculate all kproto objects for k and determine the stability values based on given method
  object <- kproto(x = data, k = k[1], keep.data = TRUE, lambda = lambda, verbose = FALSE, ...)
  trace_kp <- list(list("index" = stability_kproto(method = method, object = object, B = B, verbose = FALSE,...)[[1]][[1]], 
                        "k" = length(object$size), "object" = object))
  for(q in k[-1]){
    object <- kproto(x = data, k = q, keep.data = TRUE, lambda = lambda, verbose = FALSE, ...)
    #save kproto object, if there isn't an object for this number of cluster
    if(!any(lapply(trace_kp, `[[`, 2) == length(object$size))){
      trace_kp <- c(trace_kp, list(list("index" = stability_kproto(method = method, object = object, B = B, verbose = FALSE,...)[[1]][[1]], 
                                        "k" = length(object$size), "object" = object)))
    }else{
      # save kproto object if there is a clusterpartition with same number of cluster but different validation index
      index_value <- stability_kproto(method = method, object = object, B = B, verbose = FALSE,...)[[1]][[1]]
      if(!(index_value %in% unlist(lapply(trace_kp, `[[`, 1))[which(unlist(lapply(trace_kp, `[[`, 2)) == length(object$size))])){
        trace_kp <- c(trace_kp, list(list("index" = index_value, "k" = length(object$size), "object" = object)))
      }
    }
  }
  
  # save all index values for comparison and as summary for output
  indices <- unlist(lapply(trace_kp, `[[`, 1))
  names(indices) <- lapply(trace_kp, `[[`, 2)
  
  # find the optimal k, if it is ambiguously: sample
  # distinguish between minimum and maximum stability values:
  if(method == "luxburg"){
    k_m <- which(indices == min(indices, na.rm = TRUE))
  }else{
    k_m <- which(indices == max(indices, na.rm = TRUE))
  }
  if(length(k_m)>1){k_m <- sample(k_m,1)}
  k_opt <- as.integer(names(indices[k_m]))
  index_opt <- indices[k_m]
  
  # apply adjustment of paper
  adjustment <- data.frame(k_eval = names(indices), imean = indices) %>%
    arrange(.data$k_eval) %>% 
    mutate(Nachbar_vor = ifelse(is.na(lag(.data$imean)), 1, lag(.data$imean)),
           Nachbar_nach = ifelse(is.na(lead(.data$imean)), .data$imean, lead(.data$imean)),
           ivalue = .data$imean * sqrt(.data$imean/.data$Nachbar_vor * .data$imean/.data$Nachbar_nach))
  indices_adj <- adjustment$ivalue
  names(indices_adj) <- adjustment$k_eval
  
  if(method == "luxburg"){
    k_m <- which(indices_adj == min(indices_adj, na.rm = TRUE))
  }else{
    k_m <- which(indices_adj == max(indices_adj, na.rm = TRUE))
  }
  if(length(k_m)>1){k_m <- sample(k_m,1)}
  k_opt_adj <- as.integer(names(indices_adj[k_m]))
  index_opt_adj <- indices_adj[k_m]
  
  
  output <- list("k_opt" = k_opt, "index_opt" = index_opt, "indices" = indices,
                 "k_opt_adj" = k_opt_adj, "index_opt_adj" = index_opt_adj, "indices_adj" = indices_adj,
                 "kp_obj" = trace_kp[[k_m]]$object)
  
  return(output)
}
# # not run:
# # generate toy data with factors and numerics
# n   <- 10
# prb <- 0.99
# muk <- 2.5
# 
# x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
# x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
# x1 <- as.factor(x1)
# x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
# x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
# x2 <- as.factor(x2)
# x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
# x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
# x <- data.frame(x1,x2,x3,x4)
# 
# stability_det_k(data = x, k = 2:6, B = 100, method = "fowlkesmallows")





#' @title Determination the stability of k Prototypes Clustering
#'
#' @description Calculating the stability for a k-Prototypes clustering with k clusters or computing the stability-based optimal number of clusters for k-Prototype clustering. Possible stability indices are: \code{Jaccard}, \code{Rand}, \code{Fowlkes \& Mallows} and \code{Luxburg}.
#' 
#' @param object Object of class \code{kproto} resulting from a call with \code{kproto(..., keep.data=TRUE)}
#' @param method character specifying the stability, either one or more of \code{luxburg}, \code{fowlkesmallows}, \code{rand} or/and \code{jaccard}.
#' @param B numeric, number of bootstrap samples
#' @param verbose Logical whether information about the bootstrap procedure should be given.
#' @param ... Further arguments passed to \code{\link[clustMixType]{kproto}}, like:
#'   \itemize{
#'     \item \code{nstart}: If > 1 repetitive computations of \code{kproto} with random initial prototypes are computed.
#'     \item \code{lambda}: Factor to trade off between Euclidean distance of numeric variables and simple matching coefficient between categorical variables.
#'   }
#' 
#' @return The output contains the stability for a given k-Prototype clustering in a list with two elements:
#' @return \item{kp_stab}{stability values for the given clustering}
#' @return \item{kp_bts_stab}{stability values for each bootstrap samples}
#'
#' @author Rabea Aschenbruck
#' 
#' @references \itemize{
#'     \item Aschenbruck, R., Szepannek, G., Wilhelm, A.F.X (2023): 
#'     Stability of mixed-type cluster partitions for determination of the number of clusters. 
#'     \emph{Submitted}.
#'     
#'     \item von Luxburg, U. (2010): 
#'      Clustering stability: an overview. 
#'     \emph{Foundations and Trends in Machine Learning, Vol 2, Issue 3}.
#'     \doi{10.1561/2200000008}.
#'     
#'     \item Ben-Hur, A., Elisseeff, A., Guyon, I. (2002): 
#'     A stability based method for discovering structure in clustered data. 
#'     \emph{Pacific Symposium on Biocomputing}.
#'     \doi{10/bhfxmf}.
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
#' #' # apply k-prototypes
#' kpres <- kproto(x, 4, keep.data = TRUE)
#' 
#' # calculate cluster stability
#' stab <- stability_kproto(method = c("luxburg","fowlkesmallows"), object = kpres)
#' 
#' }
#' 
#' @rdname stability_kproto
#' 
#' @importFrom stats na.omit
#' @importFrom combinat permn
#' @importFrom dplyr %>% 
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom dplyr lag
#' @importFrom dplyr lead
#' @importFrom rlang .data
#' 
#' 
#' @export


stability_kproto <- function(object, 
                             method = c("rand", "jaccard", "luxburg", "fowlkesmallows"), B = 100, 
                             verbose=FALSE, ...){
  
  if(!inherits(object, "kproto")) stop("object must be type of 'kproto'")
  if(!is.numeric(B)) stop("B has to be numeric, e. g. B = 100")
  
  if(is.null(object$data)) stop("object must contain clustered data (keep.data = TRUE)!")
  n <- nrow(object$data)
  
  if(is.character(method)){
    if(!any(!(method %in% c("luxburg", "fowlkesmallows", "rand", "jaccard")))){
      
      ### usage of bootstrap samples
      # object for the stability value for the whole cluster partition with respect to bts_res
      kp_bts_stab <- matrix(numeric(), ncol = length(method)+1, nrow = B)
      colnames(kp_bts_stab) <- c("B", "L", "FM", "R", "J")[c(TRUE, c("luxburg", "fowlkesmallows", "rand", "jaccard") %in% method)]
      kp_bts_stab[,"B"] <- 1:B
      
    }else{
      stop("method must only contain one or more of: luxburg, fowlkesmallows, rand or/and jaccard")
    }
  }else{
    stop("method must only contain character(s)!")
  }
  
  
  
  for(b in 1:B){ 
    # b <- 1
    if (verbose) cat(paste("boot ",b,"\n", sep = ""))
    
    # generate bootstrap samples:
    # AnmR.: which method? now: same approach as in fpc.default
    #        UPDATE: all samples data will be clustered and afterwards the subset is reduced for unique data
    bsamp <- sample(n, n, replace = TRUE)
    
    # cluster of bootstrap sample
    bts_res <- kproto(x = object$data[bsamp, ], k = length(object$size), verbose=FALSE, ...) 
    
    ### NEW: choose unique data
    # which objects are in the bootstrap sample and should be compared with object?
    ind_compare <- as.character(unique(bsamp))
    
    # cluster assignments for object and bts_res only with same subsample
    kp_cl <- object$cluster[unique(bsamp)]
    bts_cl <- bts_res$cluster[ind_compare]
    
    
    # Luxburg (2009) - "not normalized stability:
    if("luxburg" %in% method){
      # create all permutations of cluster group names and count minimum number of differences
      # short minimal example: after <- c(3,4,1,2); after[match(c(2,3,3,4,1), c(1,2,3,4))]
      kp_bts_stab[b,"L"] <- min(unlist(lapply(combinat::permn(unique(bts_cl)), 
                                                    FUN = function(x) sum(kp_cl != x[match(bts_cl, sort(unique(bts_cl)))]))))/length(ind_compare)
    }
    
    if(any(c("fowlkesmallows", "rand", "jaccard") %in% method)){
      # determination of number of pairs, which are in the same cluster
      L_kp <- sapply(kp_cl, FUN = function(x) as.numeric(x == kp_cl)) - diag(length(ind_compare))
      L_bts <- sapply(bts_cl, FUN = function(x) as.numeric(x == bts_cl)) - diag(length(ind_compare))
      
      if("fowlkesmallows" %in% method){
        # eqn. 3: similarity measure introduced by Fowlkes and Mallows:
        kp_bts_stab[b,"FM"] <- sum(L_kp * L_bts) / sqrt(sum(L_kp * L_kp) * sum(L_bts * L_bts) )
      }
      
      if("rand" %in% method){
        # eqn. 4: matching coefficient aka Rand:
        kp_bts_stab[b,"R"] <- 1 - 1/length(ind_compare)^2 * (sum(abs(L_kp - L_bts)^2)^(1/2))^2
      }

      if("jaccard" %in% method){
        # eqn. 5: jaccard coefficient
        kp_bts_stab[b,"J"] <- sum(L_kp * L_bts)/(sum(L_kp * L_kp) + sum(L_bts * L_bts) - sum(L_kp * L_bts))
      }
    }
    
  }
  
  # aggregation of stability values with respect to bootstrap samples
  #   to one value for object
  kp_stab <- colMeans(kp_bts_stab[,-1, drop=FALSE])
  
  return(list(kp_stab, kp_bts_stab))
}







