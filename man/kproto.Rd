% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kprototypes.R
\name{kproto}
\alias{kproto}
\alias{kproto.default}
\title{k-Prototypes Clustering}
\usage{
kproto(x, ...)

\method{kproto}{default}(
  x,
  k,
  lambda = NULL,
  type = "standard",
  iter.max = 100,
  nstart = 1,
  na.rm = "yes",
  keep.data = TRUE,
  verbose = TRUE,
  init = NULL,
  p_nstart.m = 0.9,
  ...
)
}
\arguments{
\item{x}{Data frame with both numerics and factors.}

\item{\dots}{Currently not used.}

\item{k}{Either the number of clusters, a vector specifying indices of initial prototypes, or a data frame of prototypes of the same columns as \code{x}.}

\item{lambda}{Parameter > 0 to trade off between Euclidean distance of numeric variables 
and simple matching coefficient between categorical variables. Also a vector of variable specific factors is possible where 
the order must correspond to the order of the variables in the data. In this case all variables' distances will be multiplied by 
their corresponding lambda value.}

\item{type}{Character, to specify the distance for clustering. Either \code{"standard"} (cf. details below) or \code{"gower"}. The latter calls \code{\link{kproto_gower}}.}

\item{iter.max}{Maximum number of iterations if no convergence before.}

\item{nstart}{If > 1 repetitive computations with random initializations are computed and the result with minimum tot.dist is returned.}

\item{na.rm}{Character; Either "yes" to strip NA values for complete case analysis, "no" to keep and ignore NA values, "imp.internal" to impute the NAs within the algorithm or "imp.onestep" to apply the algorithm ignoring the NAs and impute them after the partition is determined.}

\item{keep.data}{Logical whether original should be included in the returned object.}

\item{verbose}{Logical whether additional information about process should be printed. 
Caution: For \code{verbose=FALSE}, if the number of clusters is reduced during the iterations it will not mentioned.}

\item{init}{Character, to specify the initialization strategy. Either \code{"nbh.dens"}, \code{"sel.cen"} or \code{"nstart.m"}. Default is \code{"NULL"}, which results in nstart repetitive algorithm computations with random starting prototypes. Otherwise, \code{nstart} is not used. Argument k must be a number if a specific initialization strategy is choosen!}

\item{p_nstart.m}{Numeric, probability(=0.9 is default) for \code{init="nstart.m"}, where the strategy assures that with a probability of \code{p_nstart.m} at least one of the m sets of initial prototypes contains objects of every cluster group (cf. Aschenbruck et al. (2023): Random-based Initialization for clustering mixed-type data with the k-Prototypes algorithm. In: {\emph{Cladag 2023 Book of abstracts and short spapers}}, isbn: 9788891935632.).}
}
\value{
\code{\link{kmeans}} like object of class \code{kproto}:

\item{cluster}{Vector of cluster memberships.}

\item{centers}{Data frame of cluster prototypes.}

\item{lambda}{Distance parameter lambda.}

\item{size}{Vector of cluster sizes.}

\item{withinss}{Vector of within cluster distances for each cluster, i.e. summed distances of all observations belonging to a cluster to their respective prototype.}

\item{tot.withinss}{Target function: sum of all observations' distances to their corresponding cluster prototype.}

\item{dists}{Matrix with distances of observations to all cluster prototypes.}

\item{iter}{Prespecified maximum number of iterations.}

\item{trace}{List with two elements (vectors) tracing the iteration process: 
\code{tot.dists} and \code{moved} number of observations over all iterations.}

\item{inits}{Initial prototypes determined by specified initialization strategy, if init is either 'nbh.dens' or 'sel.cen'.}

\item{nstart.m}{only for 'init = nstart_m': determined number of randomly choosen sets.}

\item{data}{if 'keep.data = TRUE' than the original data will be added to the output list.}

\item{type}{Type argument of the function call.}

\item{stdization}{Only returned for \code{type = "gower"}: List of standardized ranks for ordinal variables 
and an additional element \code{num_ranges} with ranges of all numeric variables. Used by \code{\link{predict.kproto}}.}
}
\description{
Computes k-prototypes clustering for mixed-type data.
}
\details{
The algorithm like k-means iteratively recomputes cluster prototypes and reassigns clusters.
For \code{type = "standard"}
clusters are assigned using \eqn{d(x,y) =  d_{euclid}(x,y) + \lambda d_{simple\,matching}(x,y)}.
Cluster prototypes are computed as cluster means for numeric variables and modes for factors 
(cf. Huang, 1998).
Ordered factors variables are treated as categorical variables.
In case of \code{na.rm = FALSE}: for each observation variables with missings are ignored 
(i.e. only the remaining variables are considered for distance computation). 
In consequence for observations with missings this might result in a change of variable's weighting compared to the one specified
by \code{lambda}. For these observations distances to the prototypes will typically be smaller as they are based 
on fewer variables.
For \code{type = "gower"} cf. \code{\link{kproto_gower}}.
}
\examples{
# generate toy data with factors and numerics

n   <- 100
prb <- 0.9
muk <- 1.5 
clusid <- rep(1:4, each = n)

x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
x1 <- as.factor(x1)

x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
x2 <- as.factor(x2)

x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))

x <- data.frame(x1,x2,x3,x4)

# apply k-prototypes
kpres <- kproto(x, 4)
clprofiles(kpres, x)

# in real world clusters are often not as clear cut
# by variation of lambda the emphasize is shifted towards factor / numeric variables    
kpres <- kproto(x, 2)
clprofiles(kpres, x)

kpres <- kproto(x, 2, lambda = 0.1)
clprofiles(kpres, x)

kpres <- kproto(x, 2, lambda = 25)
clprofiles(kpres, x)

}
\references{
\itemize{
    \item Szepannek, G. (2018): clustMixType: User-Friendly Clustering of Mixed-Type Data in R, {\emph{The R Journal 10/2}}, 200-208, 
          \doi{10.32614/RJ-2018-048}.
    \item Aschenbruck, R., Szepannek, G., Wilhelm, A. (2022): Imputation Strategies for Clustering Mixed‑Type Data with Missing Values, 
    {\emph{Journal of Classification}}, \doi{10.1007/s00357-022-09422-y}. 
    \item Z.Huang (1998): 
          Extensions to the k-Means Algorithm for Clustering Large Data Sets with Categorical Variables, 
          Data Mining and Knowledge Discovery 2, 283-304.
  }
}
\author{
\email{gero.szepannek@web.de}
}
\keyword{classif}
\keyword{cluster}
\keyword{multivariate}
