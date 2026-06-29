#' @title Variable Importance
#' @description Computes two types of variable importance for an object of class \code{kproto}.
#' 
#' @details The contribution of the variables to \code{object$tot.withinss} is computed 
#' as well as those for a null model where all variables are in the same cluster. 
#' Variable importance is computed as variable-wise improvement of \code{object$tot.withinss}.
#' 
#' As proposed by Henning and Murphy (2023) the importance of a variable can also be assessed by 
#' comparing how much the resulting clusters will change if a variable is ignored: 
#' For \code{rand == TRUE} in addition variable importance is computed by reclustering the data 
#' where the variable of interest is replaced by a constant value.
#' Variable importance (VI) is computed by comparing the original partition 
#' with the one obtained after reclustering the data using adjusted Rand indices as \eqn{VI = 1-Rand_{adj}}. 
#' 
#' As reclustering for each variable may be time consuming, it may be descelected.
#' 
#' Note that for \code{na.rm = "imp.onestep"} and \code{"imp.internal"} the imputed data are used. 
#' For \code{na.rm = "yes"} observations with missing values are removed from the data and 
#' for \code{na.rm = "no"} \code{NA}s are ignored for distance computation just as in \code{kproto()}.
#'  
#' @param object Object of class \code{kproto}.
#' @param rand   Logical indicating whether additional variable importance based on rand indices should be computed.
#' 
#' @return List of three elements:
#' @return \item{vi_decomp}{Data frame where the rows are: the contribution of the variables to \code{object$tot.withinss}, 
#'                          the the contribution of the variables in a null model, the differences and releatiove differences}
#' @return \item{vi_rand}{Vector of variable importances based on the adjusted Rand index after reclustering.}
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
#' 
#' # compute different kinds of variable importance
#' importance_kproto(kpres, rand = TRUE)
#' 
#' # reducing the influence of the categorical variables (X1 and X2) on the clustering:
#' kpres <- kproto(x, 4, lambda = 0.01)
#' importance_kproto(kpres, rand = FALSE)
#' 
#' # reducing the importance of the numeric variables (X3 and X4) on the clustering:
#' kpres <- kproto(x, 4, lambda = 100)
#' importance_kproto(kpres, rand = FALSE)
#' 
#' # weighting all variables equally:
#' kpres <- kproto(x, 4, lambda = 1)
#' importance_kproto(kpres, rand = FALSE)
#' 
#' @author \email{gero.szepannek@@web.de}
#' 
#' @references \itemize{
#'     \item Szepannek, G. (2018):
#'     clustMixType: User-Friendly Clustering of Mixed-Type Data in R 
#'     {\emph{The R Journal 10/2}}, 200-208,  
#'     \doi{10.32614/RJ-2018-048}.
#'     \item Hennig, C., and K. Murphy (2023): 
#'     Quantyfying Variable Importance in Cluster Analysis, 
#'     {\emph{CLADAG} 2023 Book of Abstracts and Short Papers}, 515-518,  
#'     Pearson Education Resources.
#'     \item Aschenbruck, R., Szepannek, G., Wilhelm, A. (2022): 
#'     Imputation Strategies for Clustering Mixed‑Type Data with Missing Values, 
#'     {\emph{Journal of Classification}}, 
#'     \doi{10.1007/s00357-022-09422-y}. 
#'     }.
#'     
#' @rdname importance_kproto
#' 
#' @importFrom mclust adjustedRandIndex
#' 
#' @export


importance_kproto <- function(object, rand = TRUE){
  
  if(!(inherits(object, "kproto"))) stop("Object must be of class kproto!")
  if(!("data" %in% names(object)))  stop("Computation is not possible. The kproto object has to be created using the argument keep.data = TRUE!")
  
  # check for numeric and factor variables
  numvars <- sapply(object$data, is.numeric)
  anynum <- any(numvars)
  ordvars <- sapply(object$data, is.ordered)
  anyord <- any(ordvars)
  catvars <- sapply(object$data, is.factor) & !ordvars
  anyfact <- any(catvars)
  

  distance_decomposition <- function(object){
    
    type         <- object$type
    data         <- object$data 
    cl           <- object$cluster 
    tot.withinss <- object$tot.withinss
    k            <- length(object$withinss) 
    lambda       <- object$lambda
    
    # check for numeric and factor variables
    numvars <- sapply(data, is.numeric)
    anynum <- any(numvars)
    ordvars <- sapply(data, is.ordered)
    anyord <- any(ordvars)
    catvars <- sapply(data, is.factor) & !ordvars
    anyfact <- any(catvars)
    
    protos       <- object$centers
    
    # huang:
    ds <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
    if (type == "huang"){
      for(j in 1:ncol(data)){
        if (numvars[j])          ds[,j] <- (data[,j] - protos[cl,j])^2
        if (catvars[j])          ds[,j] <- as.numeric(data[,j] != protos[cl,j])
        if (length(lambda) > 1)  ds[,j] <- ds[,j] * lambda[j]
      }  
      if(length(lambda) == 1) ds[,catvars] <- ds[,catvars] * lambda 
    }
    
    # gower:
    if (type == "gower"){
      x <- data
      # vector of ranges for normalization  
      if(any(numvars)) rgnums <- sapply(x[, numvars, drop = FALSE], function(z) diff(range(z, na.rm = TRUE)))
      if(any(ordvars)){
        xord   <- x[, ordvars, drop = FALSE] # store original variables 
        # ...and replace ordered variables by their ranks
        for(jord in which(ordvars)) x[,jord] <- rank(x[,jord], na.last = "keep")
        rgords <- sapply(x[, ordvars, drop = FALSE], function(z) diff(range(z, na.rm = TRUE)))
      }
      
      for(j in 1:ncol(data)){
        if (numvars[j]) ds[,j] <- abs(data[,j] - protos[cl,j]) / rgnums[[colnames(data)[j]]]
        if (catvars[j]) ds[,j] <- (data[,j] != protos[cl,j])
        if (ordvars[j]) ds[,j] <- abs(as.numeric(data[,j]) - as.numeric(protos[cl,j])) / rgords[[colnames(data)[j]]]
        if (length(lambda) > 1)  ds[,j] <- ds[,j] * lambda[j]
      }
    }
    
    # NAs:
    ds[is.na(ds)] <- 0
    
    dj <- colSums(ds) # / tot.withinss
    return(dj)
  }
  
  
  compute_vi_rand <- function(object){
    arand <- NULL
    for(j in 1:ncol(object$data)){
      x <- object$data
      x[,j] <- x[1,j] # set variable j to constant value
      new.object <- kproto(x, k = object$k.init, lambda = object$lambda, type = object$type, iter.max = object$iter.max,
                           nstart = object$nstart, na.rm = object$na.rm, init = object$init, p_nstart.m = object$p_nstart.m, verbose = F) ### ...
      arand <- c(arand, 1 - mclust::adjustedRandIndex(object$cluster, new.object$cluster))
    }
    return(arand)  
  }


  # create null model (as a copy of the object) where all clusters have the same overall prototype. 
  object0   <- object
  protos0   <- object$centers
  x         <- object$data
  clusters  <- object$cluster
  
  some_vals <- sapply(x[, , drop = FALSE], function(z) !all(is.na(z))) # only update variables if not all values are NA
  if (object$type == "huang"){
    if(any(some_vals & numvars)){
      repl  <- sapply(x[, some_vals & numvars, drop = FALSE], mean, na.rm = TRUE)
      for(i in 1:nrow(protos0)) protos0[i, some_vals & numvars] <- repl
    }
    if(any(some_vals & catvars)){
      repl <- sapply(x[, some_vals & catvars, drop = FALSE], function(z) levels(z)[which.max(table(z))])
      for(i in 1:nrow(protos0)) protos0[i, some_vals & catvars] <- repl
    }
  }
  if (object$type == "gower"){
    if(any(some_vals & numvars)){
      repl  <- sapply(x[, some_vals & numvars, drop = FALSE], median, na.rm = TRUE)
      for(i in 1:nrow(protos0)) protos0[i, some_vals & numvars] <- repl
    }
    if(any(some_vals & catvars)){
      repl <- sapply(x[, some_vals & catvars, drop = FALSE], function(z) levels(z)[which.max(table(z))])
      for(i in 1:nrow(protos0)) protos0[i, some_vals & catvars] <- repl
    }
    if(any(some_vals & ordvars)){
      repl  <- sapply(x[, some_vals & ordvars, drop = FALSE], function(z) levels(z)[median(as.numeric(z), na.rm = TRUE)])
      for(i in 1:nrow(protos0)) protos0[i, some_vals & ordvars] <- repl
    }
  }
  object0$centers <- protos0 
  
  
  # compare contribution of all variables to tot.withinss for object and the null model 
  decomp      <- distance_decomposition(object)
  decomp0     <- distance_decomposition(object0)
  decomp_redu <- decomp0 - decomp
  vi_rel      <- decomp_redu / sum(decomp_redu)
  
  vi_decomp           <- rbind(decomp, decomp0, decomp_redu, vi_rel)
  vi_decomp           <- as.data.frame(vi_decomp)
  colnames(vi_decomp) <- colnames(object$data)
  
  vi_rand             <- NULL
  if(rand) vi_rand    <- compute_vi_rand(object)
  
  result <- list(vi_decomp = vi_decomp, vi_rand = vi_rand)
  #class(result) <- "kprotoimp"
  return(result)
}

