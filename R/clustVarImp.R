# #############################
# ### ...under construction...
# #############################
# 
# 
# # library(clustMixType)
# # set.seed(42)
# # n   <- 100; prb <- 0.9; muk <- 1.5
# # clusid <- rep(1:4, each = n)
# # x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
# # x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
# # x1 <- as.factor(x1)
# # x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
# # x <- data.frame(x1,x3)
# # 
# # k <- 4
# # f <- "kproto"
# # a <- paste(f, "(x, k = ", k,")", sep = "")
# # 
# # obj <- eval(str2expression(a))
# # names(obj)
# # vars = NULL
# # data = obj$data
# # f = "kproto"
# # boot = 5
# # args = list()
# 
# clustVarImp <- function(obj, data = obj$data, f = "kproto", vars = NULL, boot = 100, args = list()){
#   cl.orig <- obj$cluster
#   k       <- length(unique(cl.orig))
#   n       <- nrow(data)
#   
#   if (is.null(vars)) vars <- 1:ncol(data)
#   if(is.numeric(vars)){
#     if(!any((1:ncol(data)) %in% vars)) stop("Indices in vars are not in data!")
#     if(!all(vars  %in% 1:ncol(data))) warning("Some indices in vars are not in data!")
#   }
#   if(is.character(vars)){
#     if(!any(colnames(data) %in% vars)) stop("No variable name in vars matches any colnames in data!")
#     if(!all(vars  %in% colnames(data))) warning("Some variable name(s) without match in data!")
#     vars <- which.max(colnames(data) %in% vars)
#   }
#   
#   varimps <- matrix(NA, nrow = boot, ncol = length(vars))
#   colnames(varimps) <- colnames(data)[vars]
#   
#   fcall <- paste(f, "(x, k = ", k,")", sep = "")
#   if (length(args) > 0) {
#     for(i in 1: length(args)){
#       fcall <- paste(fcall, ",", args[[i]])
#     }
#   }
# 
#     
#   for(j in 1:length(vars)){
#     for(b in 1:boot){
#       v      <- vars[j]
#       x      <- data 
#       x[,v]  <- sample(x[,v])
#       obj.vb <- eval(str2expression(fcall)) 
#       varimps[b,j] <- mclust::adjustedRandIndex(obj$cluster, obj.vb$cluster)     
#     }
#   }
# 
#   mean.VIs <- apply(varimps, 2, mean)
#   all.VIs  <- varimps
#   return(list(mean.VIs, all.VIs))
# }
