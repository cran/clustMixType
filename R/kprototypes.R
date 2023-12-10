###############################################
# For their feedback and suggestions thanks to:
# - Bridget Ssendagala 
# - Ben Feng 
# - Galen Terziysky
###############################################

#' @title k-Prototypes Clustering
#' @description Computes k-prototypes clustering for mixed-type data.
#' 
#' @details The algorithm like k-means iteratively recomputes cluster prototypes and reassigns clusters.
#' For \code{type = "standard"}
#' clusters are assigned using \eqn{d(x,y) =  d_{euclid}(x,y) + \lambda d_{simple\,matching}(x,y)}.
#' Cluster prototypes are computed as cluster means for numeric variables and modes for factors 
#' (cf. Huang, 1998).
#' Ordered factors variables are treated as categorical variables.
#' In case of \code{na.rm = FALSE}: for each observation variables with missings are ignored 
#' (i.e. only the remaining variables are considered for distance computation). 
#' In consequence for observations with missings this might result in a change of variable's weighting compared to the one specified
#' by \code{lambda}. For these observations distances to the prototypes will typically be smaller as they are based 
#' on fewer variables.
#' For \code{type = "gower"} cf. \code{\link{kproto_gower}}.
#' 
#' @keywords classif 
#' @keywords cluster
#' @keywords multivariate
#' 
#' 
#' @rdname kproto
#' @export 
kproto <- function (x, ...) 
  UseMethod("kproto")

#' @aliases kproto kproto.default 
#' 
#' @param x Data frame with both numerics and factors.
#' @param k Either the number of clusters, a vector specifying indices of initial prototypes, or a data frame of prototypes of the same columns as \code{x}.
#' @param lambda Parameter > 0 to trade off between Euclidean distance of numeric variables 
#' and simple matching coefficient between categorical variables. Also a vector of variable specific factors is possible where 
#' the order must correspond to the order of the variables in the data. In this case all variables' distances will be multiplied by 
#' their corresponding lambda value.
#' @param type Character, to specify the distance for clustering. Either \code{"standard"} (cf. details below) or \code{"gower"}. The latter calls \code{\link{kproto_gower}}.
#' @param iter.max Maximum number of iterations if no convergence before.
#' @param nstart If > 1 repetitive computations with random initializations are computed and the result with minimum tot.dist is returned.
#' @param na.rm Character; Either "yes" to strip NA values for complete case analysis, "no" to keep and ignore NA values, "imp.internal" to impute the NAs within the algorithm or "imp.onestep" to apply the algorithm ignoring the NAs and impute them after the partition is determined.
#' @param keep.data Logical whether original should be included in the returned object.
#' @param verbose Logical whether additional information about process should be printed. 
#' Caution: For \code{verbose=FALSE}, if the number of clusters is reduced during the iterations it will not mentioned.
#' @param init Character, to specify the initialization strategy. Either \code{"nbh.dens"}, \code{"sel.cen"} or \code{"nstart.m"}. Default is \code{"NULL"}, which results in nstart repetitive algorithm computations with random starting prototypes. Otherwise, \code{nstart} is not used. Argument k must be a number if a specific initialization strategy is choosen!
#' @param p_nstart.m Numeric, probability(=0.9 is default) for \code{init="nstart.m"}, where the strategy assures that with a probability of \code{p_nstart.m} at least one of the m sets of initial prototypes contains objects of every cluster group (cf. Aschenbruck et al. (2023): Random-based Initialization for clustering mixed-type data with the k-Prototypes algorithm. In: {\emph{Cladag 2023 Book of abstracts and short spapers}}, isbn: 9788891935632.).  
#' @param \dots Currently not used.
#
#' @return \code{\link{kmeans}} like object of class \code{kproto}:
#' @return \item{cluster}{Vector of cluster memberships.}
#' @return \item{centers}{Data frame of cluster prototypes.}
#' @return \item{lambda}{Distance parameter lambda.}
#' @return \item{size}{Vector of cluster sizes.}
#' @return \item{withinss}{Vector of within cluster distances for each cluster, i.e. summed distances of all observations belonging to a cluster to their respective prototype.}
#' @return \item{tot.withinss}{Target function: sum of all observations' distances to their corresponding cluster prototype.}
#' @return \item{dists}{Matrix with distances of observations to all cluster prototypes.}
#' @return \item{iter}{Prespecified maximum number of iterations.}
#' @return \item{trace}{List with two elements (vectors) tracing the iteration process: 
#' \code{tot.dists} and \code{moved} number of observations over all iterations.}
#' @return \item{inits}{Initial prototypes determined by specified initialization strategy, if init is either 'nbh.dens' or 'sel.cen'.}
#' @return \item{nstart.m}{only for 'init = nstart_m': determined number of randomly choosen sets.}
#' @return \item{data}{if 'keep.data = TRUE' than the original data will be added to the output list.}
#' @return \item{type}{Type argument of the function call.}
#' @return \item{stdization}{Only returned for \code{type = "gower"}: List of standardized ranks for ordinal variables 
#' and an additional element \code{num_ranges} with ranges of all numeric variables. Used by \code{\link{predict.kproto}}.}

#'   
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 100
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
#' clprofiles(kpres, x)
#' 
#' # in real world clusters are often not as clear cut
#' # by variation of lambda the emphasize is shifted towards factor / numeric variables    
#' kpres <- kproto(x, 2)
#' clprofiles(kpres, x)
#' 
#' kpres <- kproto(x, 2, lambda = 0.1)
#' clprofiles(kpres, x)
#' 
#' kpres <- kproto(x, 2, lambda = 25)
#' clprofiles(kpres, x)
#' 
#' @author \email{gero.szepannek@@web.de}
#' 
#' @references \itemize{
#'     \item Szepannek, G. (2018): clustMixType: User-Friendly Clustering of Mixed-Type Data in R, {\emph{The R Journal 10/2}}, 200-208, 
#'           \doi{10.32614/RJ-2018-048}.
#'     \item Aschenbruck, R., Szepannek, G., Wilhelm, A. (2022): Imputation Strategies for Clustering Mixed‑Type Data with Missing Values, 
#'     {\emph{Journal of Classification}}, \doi{10.1007/s00357-022-09422-y}. 
#'     \item Z.Huang (1998): 
#'           Extensions to the k-Means Algorithm for Clustering Large Data Sets with Categorical Variables, 
#'           Data Mining and Knowledge Discovery 2, 283-304.
#'   }
#' 
#' @rdname kproto
#' 
#' @importFrom stats complete.cases
#' @importFrom tibble is_tibble 
#' @importFrom stats qgeom
#' 
#' @method kproto default
#' @export 
kproto.default <- function(x, k, lambda = NULL, type = "standard", iter.max = 100, nstart = 1, na.rm = "yes", keep.data = TRUE, verbose = TRUE, init = NULL, p_nstart.m = 0.9, ...){
  
  # enable input of tibbles
  if(is_tibble(x) == TRUE){x <- as.data.frame(x)}
  
  if(nrow(x) == 1) stop("Only one observation clustering not meaningful.")
  k_input <- k # store input k for nstart > 1 as clusters can be merged 
  
  # treatment of missings
  NAcount <- apply(x, 2, function(z) sum(is.na(z)))
  if(verbose){
    cat("# NAs in variables:\n")
    print(NAcount)
  }
  if(any(NAcount == nrow(x))) stop(paste("Variable(s) have only NAs please remove them:",names(NAcount)[NAcount == nrow(x)],"!"))
  
  # Backward compatibility
  if(is.logical(na.rm)){
    if(na.rm){
      na.rm <- "yes"
    }else{
      na.rm <- "no"
    }
    message("Logical input for na.rm is deprecated. Please use either 'yes','no','imp.internal' or 'imp.onestep'.\n")
  }
  if(na.rm %in% c("yes","no","imp.internal", "imp.onestep")){
    if(na.rm == "yes") {
      miss <- apply(x, 1, function(z) any(is.na(z)))
      if(verbose){
        cat(sum(miss), "observation(s) with NAs.\n")
        if(sum(miss) > 0) message("Observations with NAs are removed.\n")
        cat("\n")
      } 
      x <- x[!miss,]
    } # remove missings
    
    if(na.rm != "yes"){
      allNAs <- apply(x,1,function(z) all(is.na(z)))
      if(sum(allNAs) > 0){
        if(verbose) cat(sum(allNAs), "observation(s) where all variables NA.\n")
        warning("No meaningful cluster assignment possible for observations where all variables NA.\n")
        if(verbose) cat("\n")
        
      }
    }
  }else{
    stop("Argument na.rm must be either 'yes','no','imp.internal' or 'imp.onestep'!")
  }
  
  if(!type %in% c("standard", "gower")) stop("Argument type must be either 'standard' or 'gower'!")
  if(type == "standard"){
    
    # save origin data, before starting the imputation process
    if(na.rm == "imp.internal"){origin <- x}
    
    # initial error checks
    if(!is.data.frame(x)) stop("x should be a data frame!")
    if(ncol(x) < 2) stop("For clustering x should contain at least two variables!")
    if(iter.max < 1 | nstart < 1) stop("iter.max and nstart must not be specified < 1!")
    if(!is.null(lambda)){
      if(any(lambda < 0)) stop("lambda must be specified >= 0!")
      if(!any(lambda > 0)) stop("lambda must be specified > 0 for at least two variable!")
    }
    # check for numeric and factor variables
    numvars <- sapply(x, is.numeric)
    anynum <- any(numvars)
    catvars <- sapply(x, is.factor)
    anyfact <- any(catvars)
    if(!anynum) stop("\n No numeric variables in x! Try using kmodes() from package klaR...\n\n")
    if(!anyfact) stop("\n No factor variables in x! Try using kmeans()...\n\n")
    
    # initialize prototypes
    if(!is.null(init)){ # determine prototypes with selected initialization strategy
      if(length(k) == 1){
        if(as.integer(k) != k){k <- as.integer(k); warning(paste("k has been set to", k,"!"))}
      }else{
        stop("Argument k must be a number if a specific initialization strategy is choosen!")
      }
      if(init == "nbh.dens"){
        k <- protos_init <- f_nbh_dens(dat = x, k = k)
      }else{
        if(init == "sel.cen"){
          k <- protos_init <- f_sel_cen(dat = x, k = k)
        }else{
          if(init == "nstart.m"){
            N <- sum(complete.cases(x))

            prob <- 1
            for(i in 0:(k-1)){prob <- prob * (N - i * N/k)/(N - i)}
            
            # pgeom describes x number of failures and x+1 number of failures until 1 success.. => therefore: qgeom + 1
            nstart <- init_nstart.m <- qgeom(p = p_nstart.m, prob = prob)+1
          }else{
            stop("Argument 'init' must be either 'nbh.dens','sel.cen','nstart.m' or NULL!")
          }
        }
      }
        
    }
    if(!is.data.frame(k)){
      if(length(k) == 1){
        if(as.integer(k) != k){k <- as.integer(k); warning(paste("k has been set to", k,"!"))}
        if(sum(complete.cases(x)) < k) stop("Data frame has less complete observations than clusters!")
        ids <- sample(row.names(x[complete.cases(x),]), k)
        protos <- x[ids,]
      }
      if(length(k) > 1){
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
      protos <- k
      if(any(!complete.cases(protos))) stop("Prototypes with missing values. Choose initial prototypes without missing values!")
      k <- nrow(protos)
    }
    if(k < 1) stop("Number of clusters k must not be smaller than 1!")
    
    # automatic calculation of lambda
    if(length(lambda) > 1){
      if(length(lambda) != sum(c(numvars,catvars))) stop("If lambda is a vector, its length should be the sum of numeric and factor variables in the data frame!")
      # warning for variable selection via lambda (which results in no numvars or no catvars)
      if(all(!as.logical(numvars*lambda))) warning("As a result of the choice of lambda: No numeric variables in x! Better try using kmodes() from package klaR...\n")
      if(all(!as.logical(catvars*lambda))) warning("As a result of the choice of lambda: No factor variables in x! Better try using kmeans()...\n")
    }else{
      if(length(lambda) == 1) {if(lambda == 0) stop("lambda has to be a value != 0. For automatic calculation use lambda = NULL (default setting)!")}
    }
    if(is.null(lambda)){
      lambda <- lambdaest(x = x, num.method = 1, fac.method = 1, outtype = "numeric", verbose = FALSE)
      if(verbose) cat("Estimated lambda:", lambda, "\n\n")
    }
    
    # deleting all incomplete observations
    if(na.rm == "yes"){
      miss <- apply(x, 1, function(z) any(is.na(z)))
      if(verbose){
        cat(sum(miss), "observation(s) with NAs.\n")
        if(sum(miss) > 0) message("Observations with NAs are removed.\n")
        cat("\n")
      }
      x <- x[!miss,]
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
          d1 <- sum((protos[l,numvars, drop = FALSE]-protos[m,numvars, drop = FALSE])^2) # euclidean for numerics
          d2 <- sum(protos[l,catvars, drop = FALSE] != protos[m,catvars, drop = FALSE]) # wtd simple matching for categorics 
          if((d1+d2) == 0) keep.protos[m] <- FALSE 
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
        #a0 <- proc.time()[3]      
        #d1 <- apply(x[,numvars],1, function(z) sum((z-protos[i,numvars])^2)) # euclidean for numerics
        d1 <- (x[,numvars, drop = FALSE] - matrix(rep(as.numeric(protos[i, numvars, drop = FALSE]), nrows), nrow=nrows, byrow=T))^2
        d1[is.na(d1)] <- 0
        if(length(lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
        if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
        #a1 <- proc.time()[3]      
        #d2 <- lambda * apply(x[,catvars],1, function(z) sum((z != protos[i,catvars]))) # wtd simple matching for categorics 
        d2 <- sapply(which(catvars), function(j) return(x[,j] != rep(protos[i,j], nrows)) )
        d2[is.na(d2)] <- FALSE
        if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
        if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
        #a2 <- proc.time()[3]      
        dists[,i] <- d1 + d2
        #cat(a1-a0, a2-a1, "\n")
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
          protos[which(remids == i), some_vals & numvars] <- sapply(x[clusters == i, some_vals & numvars, drop = FALSE], mean, na.rm = TRUE)
        }
        if(any(some_vals & catvars)){
          protos[which(remids == i), some_vals & catvars] <- sapply(x[clusters == i, some_vals & catvars, drop = FALSE], function(z) levels(z)[which.max(table(z))])
        }
      }
      
      # update missing values, if NAs should be imputed during the algorithm
      if(na.rm == "imp.internal"){
        x_old <- x
        x <- origin
        
        if(any(is.na(x[,numvars]))){
          x[,numvars][is.na(x[,numvars])] <- protos[clusters,numvars][is.na(x[,numvars])]
        }
        if(any(is.na(x[,catvars]))){
          x[,catvars][is.na(x[,catvars])] <- protos[clusters,catvars][is.na(x[,catvars])]
        }
        
        # observation with all values NA won't be imputed
        x[allNAs,] <- NA
      }
      
      if(k == 1){clusters <- rep(1, length(clusters)); size <- table(clusters); iter <- iter.max; break}
      
      # check for any equal prototypes and reduce cluster number in case of occurence
      if(iter == (iter.max-1)){ # REM: for last iteration equal prototypes are allowed. otherwise less prototypes than assigned clusters.
        keep.protos <- rep(TRUE,k)
        for(l in 1:(k-1)){
          for(m in (l+1):k){
            d1 <- sum((protos[l,numvars, drop = FALSE]-protos[m,numvars, drop = FALSE])^2) # euclidean for numerics
            d2 <- sum(protos[l,catvars, drop = FALSE] != protos[m,catvars, drop = FALSE]) # wtd simple matching for categorics 
            if((d1+d2) == 0) keep.protos[m] <- FALSE 
          }
        }
        if(!all(keep.protos)){
          protos <- protos[keep.protos,]
          k <- sum(keep.protos)
          if(verbose) cat("Equal prototypes merged. Cluster number reduced to:", k, "\n\n")      
        }
      }
      
      # add stopping rules
      if(moved[length(moved)] ==  0){
        if(na.rm != "imp.internal"){
          break
        }else{
          if(sum(!(x == x_old), na.rm = TRUE) == 0){
            break
          }
        }
      }
      
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
          protos[which(remids == i), some_vals & numvars] <- sapply(x[clusters == i, some_vals & numvars, drop = FALSE], mean, na.rm = TRUE)
        }
        if(any(some_vals & catvars)){
          protos[which(remids == i), some_vals & catvars] <- sapply(x[clusters == i, some_vals & catvars, drop = FALSE], function(z) levels(z)[which.max(table(z))])
        }
      }
      
      # compute distances 
      nrows <- nrow(x)
      dists <- matrix(NA, nrow=nrows, ncol = k)
      for(i in 1:k){
        d1 <- (x[,numvars, drop = FALSE] - matrix(rep(as.numeric(protos[i, numvars, drop = FALSE]), nrows), nrow=nrows, byrow=T))^2
        d1[is.na(d1)] <- 0
        if(length(lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
        if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
        d2 <- sapply(which(catvars), function(j) return(x[,j] != rep(protos[i,j], nrows)) )
        d2[is.na(d2)] <- FALSE
        if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
        if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
        dists[,i] <- d1 + d2
      }
      
      size          <- table(clusters)  
      min.dists     <- apply(cbind(clusters, dists), 1, function(z) z[z[1]+1])
      within        <- as.numeric(by(min.dists, clusters, sum))
      tot.within    <- sum(within)
    }
    
    # observations with all NA are not assigned to a cluster
    if(na.rm != "yes"){
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
                trace = list(tot.dists = tot.dists, moved = moved))
    
    # add the determined initial prototypes
    if(!is.null(init)){
      if(init %in% c("nbh.dens", "sel.cen")){res[["inits"]] <- protos_init}
    }
    
  }
  
  if(type == "gower"){
    if(na.rm %in% c("imp.internal", "imp.onestep")){
      stop("Argument na.rm must be either 'yes' or 'no', since imputation is not yet implemented for type = 'gower'!")
    }
    res <- kproto_gower(x=x, k=k_input, lambda = lambda, iter.max = iter.max, verbose=verbose, na.rm = na.rm)
  }
  
  # if keep.data == TRUE add data to list before nstart-loop otherwise the first determined x (inclusive imputed values)
  # is added and not necessarily the data of the x with the "best" imputed values (of clusterpartition with smallest tot.withinss)
  if(keep.data) res$data = x
  
  # loop: if nstart > 1:
  if(nstart > 1)
    for(j in 2:nstart){
      # start function kproto with the origin data in case of imputation
      if(na.rm == "imp.internal"){
        x <- origin
      }
      res.new <- kproto(x=x, k=k_input, lambda = lambda,  type = type, iter.max = iter.max, nstart=1, verbose=verbose, na.rm = na.rm, keep.data = keep.data)
      if(!is.null(init)){
        if(init == "nstart.m"){res.new[["nstart.m"]] <- init_nstart.m}
      } 
      if(res.new$tot.withinss < res$tot.withinss) res <- res.new
    }  
  
  if(na.rm == "imp.onestep"){
    x <- res$data
    if(any(is.na(x[,numvars]))){
      x[,numvars][is.na(x[,numvars])] <- protos[clusters,numvars][is.na(x[,numvars])]
    }
    if(any(is.na(x[,catvars]))){
      x[,catvars][is.na(x[,catvars])] <- protos[clusters,catvars][is.na(x[,catvars])]
    }
    res$data <- x
  }
  
  res$type <- type
  class(res) <- "kproto"
  return(res)
}



### determination of distances, required for initialization functions
dists_kproto <- function(x, y = NULL, lambda = NULL, verbose = FALSE){
  
  if(is.null(y)){
    dat_part1 <- x[rep(c(1:nrow(x)), each = nrow(x)),]
    dat_part2 <- x[rep(c(1:nrow(x)), times = nrow(x)),]
  }else{
    if(nrow(x) == 1 & nrow(y) == 1){
      return(cbind(x, y, dist = dist_kproto(x,y,lambda = lambda)))
    }else{
      dat_part1 <- x[rep(c(1:nrow(x)), each = nrow(y)),]
      dat_part2 <- y[rep(c(1:nrow(y)), times = nrow(x)),]
    }
  }
  
  # check for numeric and factor variables
  numvars <- sapply(x, is.numeric)
  anynum <- any(numvars)
  catvars <- sapply(x, is.factor)
  anyfact <- any(catvars)
  
  # determination of lambda
  if(length(lambda) > 1) {if(length(lambda) != sum(c(numvars,catvars))) stop("If lambda is a vector, its length should be the sum of numeric and factor variables in the data frame!")}
  if(is.null(lambda)){
    if(anynum & anyfact){
      vnum <- mean(sapply(x[,numvars, drop = FALSE], var, na.rm = TRUE))
      vcat <- mean(sapply(x[,catvars, drop = FALSE], function(z) return(1-sum((table(z)/sum(!is.na(z)))^2))))
      if (vnum == 0){
        if(verbose) warning("All numerical variables have zero variance.")
        anynum <- FALSE
      } 
      if (vcat == 0){
        if(verbose) warning("All categorical variables have zero variance.")
        anyfact <- FALSE
      } 
      if(anynum & anyfact){
        lambda <- vnum/vcat
        if(verbose) cat("Estimated lambda:", lambda, "\n\n")
      }else{
        lambda <- 1
      }
    }
  }
  
  # compute distances 
  nrows <- nrow(x)
  d1 <- (dat_part1[,numvars, drop = FALSE] - dat_part2[,numvars, drop = FALSE])^2
  if(length(lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
  if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
  d2 <- sapply(which(catvars), function(j) return(dat_part1[,j] != dat_part2[,j]))
  d2[is.na(d2)] <- FALSE
  if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
  if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
  
  return(cbind(dat_part1, dat_part2, dist = as.vector(d1 + d2)))
}



dist_kproto <- function(x, y, lambda, verbose = FALSE){
  
  x <- rbind(x,y)
  
  # check for numeric and factor variables
  numvars <- sapply(x, is.numeric)
  anynum <- any(numvars)
  catvars <- sapply(x, is.factor)
  anyfact <- any(catvars)
  
  # compute distances 
  nrows <- nrow(x)
  d1 <- (x[1, numvars, drop = FALSE] - x[2, numvars, drop = FALSE])^2
  if(length(lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
  if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
  d2 <- sapply(which(catvars), function(j) return(x[1,j] != x[2,j]))
  d2[is.na(d2)] <- FALSE
  if(length(lambda) == 1) d2 <- lambda * rowSums(matrix(d2, nrow = 1))
  if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
  
  return(as.numeric(d1 + d2))
}



### initial prototypes selection based on neighborhood density and distance
f_nbh_dens <- function(dat, k){
  
  # determine object to save initial prototypes
  protos_initial <- dat[0,]
  
  # determine distances between all pairs of objects:
  lambda <- lambdaest(x = dat, verbose = FALSE)
  all_dists <- dists_kproto(x = dat, lambda = lambda)
  
  # determine \bar{d}/2 (where \bar{d} is the average distance between all pairs of objects)
  nbh_radius <- (sum(all_dists$dist)/(nrow(all_dists)-nrow(dat)))/2
  
  # determine N_e, which is the sum of h() for every object 
  # Ne_nbh says the number of objects with distance smaller than nbh_radius for every object
  Ne_nbh <- numeric()
  for(i in 1:nrow(dat)){
    # noteR: -1 because the distance between X_i and X_i is included and always < nbh_radius...
    Ne_nbh[i] <- sum(all_dists[which(rowSums(mapply("==", dat[i,], all_dists[, 1:ncol(dat)])) == ncol(dat)), ]$dist <= nbh_radius) - 1
  }
  
  # determine first initial prototype
  arrange_index <- order(Ne_nbh, decreasing = TRUE)
  protos_initial[1,] <- dat[arrange_index[1],]
  
  # determination of initial prototypes, lowering L until we have enough initis found
  L <- nbh_radius*4 *2
  while(nrow(protos_initial) < k){
    # determine the other initial prototypes by comparison of objects to yet known prototypes
    L <- L/2
    
    for(i in arrange_index[-1]){
      dists <- dists_kproto(x = dat[i,, drop = FALSE], y = protos_initial, lambda = lambda)
      ## if(any(dists$dist > L)){ 
      # noteR: paper: !"...determine whether any M_i \in M has dist(Y_i,M_i) > L holds.
      #        If satisfied, then Y_i is put into the set M as next prototype..
      #        NOT MEANINGFUL, THE DISTANCE TO ALL PROTOS MUST BE > L...
      if(all(dists$dist > L)){ #the same but other code: if(min(dist%dist) > L)
        protos_initial <- rbind(protos_initial, dat[i,])
      }
    }
  }
  
  return(protos_initial[1:k,])
}



# Katsavounidis (1994), He: 
f_sel_cen <- function(dat, k){
  
  # determine object which has the highest 
  # sum of distances to all other objects as first prototypes
  lambda <- lambdaest(x = dat, verbose = FALSE)
  all_dists <- dists_kproto(x = dat, lambda = lambda)
  dists_matrix <- matrix(all_dists[,"dist"], ncol = nrow(dat), nrow = nrow(dat)) - diag(NA, nrow = nrow(dat))
  protos_initial <- dat[which.max(colSums(dists_matrix, na.rm = TRUE)),]
  
  # determine minimal distance to initial prototypes
  while(nrow(protos_initial) < k){
    Dmin <- numeric()
    for(j in 1:nrow(dat)){
      Dmin[j] <- min(dists_kproto(x = dat[j,, drop = FALSE], y = protos_initial, lambda = lambda)$dist)
    }
    
    # choose next initial prototype by highest distance to already choosen protos
    protos_initial <- rbind(protos_initial, dat[which.max(Dmin),])
    
  }
  
  return(protos_initial)
}


#' @title Assign k-Prototypes Clusters
#'
#' @description Predicts k-prototypes cluster memberships and distances for new data.
#'
#' @param object Object resulting from a call of \code{kproto}.
#' @param newdata New data frame (of same structure) where cluster memberships are to be predicted.
#' @param \dots Currently not used.
#
#' @return \code{\link{kmeans}} like object of class \code{kproto}:
#' @return \item{cluster}{Vector of cluster memberships.}
#' @return \item{dists}{Matrix with distances of observations to all cluster prototypes.}
# 
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 100
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
#' # apply k-prototyps
#' kpres <- kproto(x, 4)
#' predicted.clusters <- predict(kpres, x) 
#' 
#' 
#' @author \email{gero.szepannek@@web.de}
#' 
#' @rdname predict.kproto
#' @export
predict.kproto <-function(object, newdata, ...){
  lambda <- object$lambda
  k <- length(object$size)
  protos <- object$centers
  x <- newdata
  
  if(object$type == "standard"){
    # check for numeric and factor variables  
    numvars <- sapply(x, is.numeric)
    catvars <- sapply(x, is.factor)
    
    if(!all(names(object$centers) == names(x))) warning("Variable names of newdata do not match!")
    if(ncol(object$centers) != ncol(x)) stop("Column number of newdata does not match object!")
    if(!all(sapply(object$centers[,numvars, drop = FALSE], is.numeric))) stop("Numeric variables of object and newdata do not match!")
    if(!all(sapply(object$centers[,catvars, drop = FALSE], is.factor))) stop("Factor variables of object and newdata do not match!")
    
    # compute distances 
    dists <- matrix(NA, nrow=nrow(x), ncol = k)
    nrows <- nrow(x)
    for(i in 1:k){
      #a0 <- proc.time()[3]      
      #d1 <- apply(x[,numvars],1, function(z) sum((z-protos[i,numvars])^2)) # euclidean for numerics
      d1 <- (x[,numvars, drop = FALSE] - matrix(rep(as.numeric(protos[i,numvars, drop = FALSE]), nrows), nrow=nrows, byrow=T))^2
      if(length(lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
      if(length(lambda) > 1) d1<- as.matrix(d1) %*% lambda[numvars]
      #a1 <- proc.time()[3]      
      #d2 <- lambda * apply(x[,catvars],1, function(z) sum((z != protos[i,catvars]))) # wtd simple matching for categorics 
      d2 <- sapply(which(catvars), function(j) return(x[,j] != rep(protos[i,j], nrows)) )
      # fix for x data frame with only one observation according to G.Terziysky
      if (nrows == 1) d2 <- matrix(data = d2, nrow = 1, byrow = TRUE, dimnames = list(NULL, names(x)[catvars])) # <- FIX
      d2[is.na(d2)] <- FALSE
      if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
      if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
      #a2 <- proc.time()[3]      
      dists[,i] <- d1 + d2
      #cat(a1-a0, a2-a1, "\n")
    }
  }

  if(object$type == "gower"){
    # check for numeric, ordinal and factor variables
    numvars <- sapply(x, is.numeric)
    anynum <- any(numvars)
    ordvars <- sapply(x, is.ordered)
    anyord <- any(ordvars)
    catvars <- sapply(x, is.factor) & !ordvars
    anyfact <- any(catvars)

    # normalization numeric variables to unit range from lookup table 
    if(any(numvars)){
      for(jord in which(numvars)){
        n <- names(x)[jord]
        x[,jord] <- x[,jord] / object$stdization$num_ranges[[n]]
        protos[,jord] <- protos[,jord] / object$stdization$num_ranges[[n]]
      }
    } 
    # replace ordinal levels by standardized ranks from lookup table
    if(any(ordvars)){
      for(jord in which(ordvars)){
        val <- object$stdization[[names(x)[jord]]]$value
        names(val) <- object$stdization[[names(x)[jord]]]$level
        x[,jord] <- val[x[,jord]]
        protos[,jord] <- val[protos[,jord]]
      }  
    }

    # compute distances 
    nrows <- nrow(x)
    dists <- matrix(NA, nrow=nrows, ncol = k)
    for(i in 1:k){
      
      # in case of no numeric / factor / ordinale variables set:
      d1 <- d2 <- d3 <- rep(0, nrows)
      
      if(any(numvars)){
        d1 <- abs(x[, numvars, drop = FALSE] - matrix(rep(as.numeric(protos[i, numvars, drop = FALSE]), nrows), nrow=nrows, byrow=T))
        for(jnum in 1:ncol(d1)) d1[,jnum] <- d1[,jnum] / object$stdization$num_ranges[[names(d1)[jnum]]] 
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
        #for(jord in 1:ncol(d3)) d3[,jord] <- d3[,jord] / rgords[jord] # different to kproto_gower() not needed -- already standardized via exported lookup table in ln. 518ff
        d3[is.na(d3)] <- 0
        if(length(lambda) > 1) d3 <- as.matrix(d3) %*% lambda[ordvars]
        if(is.null(lambda)) d3 <- rowSums(d3)
      }
      
      dists[,i] <- d1 + d2 + d3
    }
    
  }
  
  # assign clusters
  clusters      <- apply(dists, 1, function(z) {a <- which.min(z); if (length(a)>1) a <- sample(a,1); return(a)}) # sample in case of multiple minima
  
  res <- list(cluster = clusters, dists = dists)
  return(res)
}


#' @title Profiling k-Prototypes Clustering
#'
#' @description Visualization of a k-prototypes clustering result for cluster interpretation.
#' 
#' @details For numerical variables boxplots and for factor variables barplots of each cluster are generated.
#' 
#' @param object Object resulting from a call of resulting \code{kproto}. Also other \code{kmeans} like objects with \code{object$cluster} and \code{object$size} are possible. 
#' @param x Original data.
#' @param vars Optional vector of either column indices or variable names.
#' @param col Palette of cluster colours to be used for the plots. As a default RColorBrewer's 
#' \code{brewer.pal(max(unique(object$cluster)), "Set3")} is used for k > 2 clusters and lightblue and orange else.
#'
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 100
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
#' # apply k-prototyps
#' kpres <- kproto(x, 4)
#' clprofiles(kpres, x)
#' 
#' # in real world clusters are often not as clear cut
#' # by variation of lambda the emphasize is shifted towards factor / numeric variables    
#' kpres <- kproto(x, 2)
#' clprofiles(kpres, x)

#' 
#' kpres <- kproto(x, 2, lambda = 0.1)
#' clprofiles(kpres, x)
#' 
#' kpres <- kproto(x, 2, lambda = 25)
#' clprofiles(kpres, x)

#' 
#' @author \email{gero.szepannek@@web.de}
#' 
#' @rdname clprofiles
#' @importFrom graphics barplot
#' @importFrom graphics boxplot
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
clprofiles <- function(object, x, vars = NULL, col = NULL){
  if(length(object$cluster) != nrow(x)) stop("Size of x does not match cluster result!")
  if(is.null(vars)) vars <- 1:ncol(x)
  if(!is.numeric(vars)) vars <- sapply(vars, function(z) return(which(colnames(x)==z)))
  if(length(vars) < 1) stop("Specified variable names do not match x!")
  if(is.null(col)){
    k <- max(unique(object$cluster))
    if(k > 2)  col <- brewer.pal(k, "Set3")
    if(k == 2) col <- c("lightblue","orange")
    if(k == 1) col <- "lightblue"
  }
  #clusids <- as.numeric(names(object$size)) # object size not named for kmeans
  clusids <- sort(unique(object$cluster)) 
  if(length(col) != max(clusids)) warning("Length of col should match number of clusters!")
  
  par(ask=TRUE)
  for(i in vars){
    if(is.numeric(x[,i])){
      boxplot(x[,i]~object$cluster, col = col, main = colnames(x)[i])
      legend("topright", legend=clusids, fill = col)
    } 
    if(is.factor(x[,i])){
      tab <- table(x[,i], object$cluster)
      for(j in 1:length(object$size)) tab[,j] <- tab[,j]/object$size[j]
      barplot(t(tab), beside = TRUE, main = colnames(x)[i], col = col)
    } 
  }
  par(ask=FALSE)
  invisible()
}


#' @title Compares Variability of Variables
#'
#' @description Investigation of the variables' variances/concentrations to support specification of lambda for k-prototypes clustering.
#' 
#' @details Variance (\code{num.method = 1}) or standard deviation (\code{num.method = 2}) of numeric variables 
#' and \eqn{1-\sum_i p_i^2} (\code{fac.method = 1}) or \eqn{1-\max_i p_i} (\code{fac.method = 2}) for factors is computed.
#' 
#' @param x Data.frame with both numerics and factors.
#' @param num.method Integer 1 or 2. Specifies the heuristic used for numeric variables.
#' @param fac.method Integer 1 or 2. Specifies the heuristic used for factor variables.
#' @param outtype Specifies the desired output: either 'numeric', 'vector' or 'variation'.
#' @param verbose Logical whether additional information about process should be printed. 
#' 
#' @return \item{lambda}{Ratio of averages over all numeric/factor variables is returned. 
#' In case of \code{outtype = "vector"} the separate lambda for all variables is returned as the inverse of the single variables' 
#' variation as specified by the \code{num.method} and \code{fac.method} argument. \code{outtype = "variation"} directly returns these quantities and is not meant to be 
#' passed directly to \code{kproto()}.}
#' 
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 100
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
#' lambdaest(x)
#' res <- kproto(x, 4, lambda = lambdaest(x))
#' 
#' @author \email{gero.szepannek@@web.de}
#' 
#' @rdname lambdaest
#' 
#' @importFrom stats var
#' @importFrom stats sd
#' @export
lambdaest <- function(x, num.method = 1, fac.method = 1, outtype = "numeric", verbose = TRUE){
  
  # enable input of tibbles
  if(is_tibble(x) == TRUE){x <- as.data.frame(x)}
  
  # initial error checks
  if(!is.data.frame(x)) stop("x should be a data frame!")
  if(nrow(x) == 1) stop("Determination for only one observation is not meaningful.")
  if(ncol(x) < 2) stop("x should contain at least two variables!")
  
  if(!num.method %in% 1:2) stop("Argument 'num.method' must be either 1 or 2!")
  if(!fac.method %in% 1:2) stop("Argument 'fac.method' must be either 1 or 2!")
  if(!outtype %in% c("numeric","vector","variation")) stop("Wrong specification of argument 'outtype'!")
    
  # check for numeric and factor variables
  numvars <- sapply(x, is.numeric)
  anynum <- any(numvars)
  catvars <- sapply(x, is.factor)
  anyfact <- any(catvars)
  if(!anynum) stop("\n No numeric variables in x! \n\n")
  if(!anyfact) stop("\n No factor variables in x! \n\n")
  
  if(anynum & num.method == 1) vnum <- sapply(x[,numvars, drop = FALSE], var, na.rm =TRUE)
  if(anynum & num.method == 2) vnum <- sapply(x[,numvars, drop = FALSE], sd, na.rm = TRUE)
  
  if(anyfact & fac.method == 1) vcat <- sapply(x[,catvars, drop = FALSE], function(z) return(1-sum((table(z)/sum(!is.na(z)))^2)))
  if(anyfact & fac.method == 2) vcat <- sapply(x[,catvars, drop = FALSE], function(z) return(1-max(table(z)/sum(!is.na(z)))))
  if (mean(vnum) == 0){
    if(verbose) warning("All numerical variables have zero variance.")
    anynum <- FALSE
  } 
  if (mean(vcat) == 0){
    if(verbose) warning("All categorical variables have zero variance.")
    anyfact <- FALSE
  } 
  
  if(verbose){
    if(num.method == 1) cat("Numeric variances:\n")
    if(num.method == 2) cat("Numeric standard deviations:\n")
    print(vnum)
    if(num.method == 1) cat("Average numeric variance:", mean(vnum), "\n\n")
    if(num.method == 2) cat("Average numeric standard deviation:", mean(vnum), "\n\n")
    
    cat(paste("Heuristic for categorical variables: (method = ",fac.method,") \n", sep = ""))
    print(vcat)
    cat("Average categorical variation:", mean(vcat), "\n\n")
  }
  
  # output depending on argument 'outtype'
  if(outtype == "numeric"){
    lambda <- mean(vnum)/mean(vcat)
    if(verbose) cat("Estimated lambda:", lambda, "\n\n")
  }
  if(outtype != "numeric"){
    lambda <- rep(0, ncol(x))
    names(lambda) <- names(x)
    lambda[numvars] <- vnum
    lambda[catvars] <- vcat
  }
  if(outtype == "vector"){
    lambda <- 1/lambda
  }
  
  return(lambda)
}



#' @export
print.kproto <- function(x, ...){
  cat("Distance type:", x$type, "\n\n")
  
  cat("Numeric predictors:", sum(sapply(x$centers, is.numeric)), "\n")
  numord <- sum(sapply(x$centers, is.factor))
  numcat <- sum(sapply(x$centers, is.factor))
  if(x$type == "gower"){
    cat("Ordinal predictors:", numord, "\n")
    cat("Categorical predictors:", numcat - numord, "\n")
  }
  if(x$type == "standard") cat("Categorical predictors:", numcat, "\n")
  if(!is.null(x$lambda)) cat("Lambda:", x$lambda, "\n\n")
  
  cat("Number of Clusters:", length(x$size), "\n")
  cat("Cluster sizes:", x$size, "\n")
  cat("Within cluster error:", x$withinss, "\n\n")
  
  cat("Cluster prototypes:\n")
  print(x$centers)
}

#' @title Summary Method for kproto Cluster Result
#'
#' @description Investigation of variances to specify lambda for k-prototypes clustering.
#' 
#' @details For numeric variables statistics are computed for each clusters using \code{summary()}. 
#' For categorical variables distribution percentages are computed.  
#' 
#' @param object Object of class \code{kproto}.
#' @param data Optional data set to be analyzed. If \code{!(is.null(data))} clusters for \code{data} are assigned by 
#' \code{predict(object, data)}. If not specified the clusters of the original data ara analyzed which is only possible if \code{kproto} 
#' has been called using \code{keep.data = TRUE}.
#' @param pct.dig Number of digits for rounding percentages of factor variables.
#' @param \dots Further arguments to be passed to internal call of \code{summary()} for numeric variables.
#'  
#' @return List where each element corresponds to one variable. Each row of any element corresponds to one cluster.
#' 
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 100
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
#' res <- kproto(x, 4)
#' summary(res)
#' 
#' @author \email{gero.szepannek@@web.de}
#' 
#' @rdname summary
#'
#' @importFrom stats predict 
#' @export
summary.kproto <- function(object, data = NULL, pct.dig = 3, ...){
  if(!inherits(object, "kproto")) stop("object must be of class kproto!")
  
  if(is.null(data)){
    data <- object$data
    cluster <- object$cluster
  }
  if(!is.null(data)) cluster <- predict(object, data)$cluster
  
  numvars <- sapply(data, is.numeric)
  #anynum <- any(numvars)
  catvars <- sapply(data, is.factor)
  #anyfact <- any(catvars)

  res <- NULL
  for(i in 1:ncol(data)){
    cat(names(data)[i],"\n")
    if(numvars[i]){
      resi <- by(data[,i], cluster, summary, ...)
      res[[i]] <- matrix(unlist(resi), nrow = length(unique(cluster)), byrow=TRUE)
      colnames(res[[i]]) <- names(resi[[1]])
      rownames(res[[i]]) <- sort(unique(cluster))
      }
    if(catvars[i])  res[[i]] <- round(prop.table(table(cluster, data[,i]),1), digits = pct.dig)
    print(res[[i]])
    cat("\n-----------------------------------------------------------------\n")
  }
  names(res) <- names(data)  

  #return(res)
  invisible(res)
}



# optilambda <- function(x, k, upper = 3*lambdaest(x), iter.max = 100, nstart=2, keep.data = FALSE){
#   
#   foo <- function(lamb, x = x, k = k, iter.max = iter.max, nstart=nstart, keep.data = keep.data){
#     res <- kproto(x = x, k = k, lambda = lamb, iter.max = iter.max, nstart=nstart, keep.data = keep.data)  
#     numvars <- sapply(x, is.numeric)
#     catvars <- sapply(x, is.factor)
#     nrows <- nrow(x)
#     protos <- res$centers
#     
#     dists <- ndists <- cdists <- matrix(NA, nrow=nrows, ncol = nrow(res$centers)) #dists <- rep(NA, nrows)
#     for(i in 1:nrow(res$centers)){
#       d1 <- (x[,numvars, drop = FALSE] - matrix(rep(as.numeric(protos[i,numvars, drop = FALSE]), nrows), nrow=nrows, byrow=T))^2
#       d2 <- sapply(which(catvars), function(j) return(x[,j] != rep(protos[i,j], nrows)) )
#       ndists[,i] <- rowSums(d1)
#       cdists[,i] <- lamb * rowSums(d2)
#       dists[,i]  <- ndists[,i] + cdists[,i]
#       }
#     
#     #min.dists     <- apply(cbind(res$cluster, dists), 1, function(z) z[z[1]+1])
#     min.ndists     <- apply(cbind(res$cluster, ndists), 1, function(z) z[z[1]+1])
#     min.cdists     <- apply(cbind(res$cluster, cdists), 1, function(z) z[z[1]+1])
#     crit           <- abs(sum(min.ndists) - sum(min.cdists))
#     return(crit)
#     }
#   
#     foo(lamb, x = x, k = k, iter.max = iter.max, nstart=nstart, keep.data = keep.data)
#     lambda <- optimize(foo, lower = 0, upper = upper, x = x, k = k, iter.max = iter.max, nstart=nstart, keep.data = keep.data)$minimum
#     return(lambda)  
#   }
# 
# optilambda(x=x, k = 4, nstart = 5)
# res <- kproto(x = x, k = k, lambda = 9.56, iter.max = iter.max, nstart=5, keep.data = keep.data) 


#' @title Assign k-Prototypes Clusters
#'
#' @description Plot distributions of the clusters across the variables.
#'
#' @details Wrapper around \code{\link{clprofiles}}. Only works for \code{kproto} object created with \code{keep.data = TRUE}. 
#' 
#' @param x Object resulting from a call of \code{kproto}.
#' @param \dots Additional arguments to be passet to \code{\link{clprofiles}} such as e.g. \code{vars}.
#'
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 100
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
#' # apply k-prototyps
#' kpres <- kproto(x, 4)
#' plot(kpres, vars = c("x1","x3")) 
#' 
#' 
#' @author \email{gero.szepannek@@web.de}
#' 
#' @rdname plot.kproto
#' @export
plot.kproto <- function(x, ...){
  if(!any("data" %in% names(x))) stop("plot method for kproto objects only works for keep.data = TRUE. Use clprofiles() instead.")
  clprofiles(x, x$data, ...)  
}
