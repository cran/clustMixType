context("Basic tests.\n")

library(clustMixType)

# generate test data from example
set.seed(42)
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


kpres <- kproto(x, 4)
pred  <- predict(kpres, x)
test_that("resulting objects are as expected",{
  expect_is(kpres, "kproto")
  expect_is(pred$cluster, "integer")
  expect_equal(length(pred$cluster), nrow(x))
  expect_is(pred$dists, "matrix")
  expect_equal(nrow(pred$dists), nrow(x))
  expect_equal(ncol(pred$dists), length(kpres$size))
  expect_true(length(kpres$size) <= 4)
  expect_equal(sort(unique(kpres$cluster)), 1:length(kpres$size))
}
)

test_that("different initializations work properly",{
  expect_is(kproto(x, c(1,101,201,301)), "kproto")
  expect_error(kproto(x, c(-1,101,201,301)))            # initialization by vector of indizes
  expect_error(kproto(x, c(1,101,201,401)))
  expect_is(kproto(x, x[c(1,101,201,301),]), "kproto")  # initialization by df
  expect_message(kproto(x, x[c(1,1,2,3,3),]),"Equal prototypes merged. Cluster number reduced to:3")
  expect_is(kproto(x, x[1,]), "kproto")                 # one single observation
  expect_is(kproto(x, x[1,])$centers, "data.frame")
  expect_equal(nrow(kproto(x, x[1,])$centers), 1)
  expect_is(kproto(x, 1), "kproto")                 # k = 1
}
)

test_that("different ways of specifying lambda work properly",{
  expect_is(kproto(x, 4, lambda = 1), "kproto")
  expect_equal(kproto(x, 4, lambda = 0.834)$lambda, 0.834)
  expect_is(kproto(x, 4, lambda = 2:5), "kproto")           # vector
  expect_equal(kproto(x, 4, lambda = 2:5)$lambda, 2:5)
  expect_error(kproto(x, 4, lambda = 1:3))
  expect_error(kproto(x, 4, lambda = -1:2))
  expect_error(kproto(x, 4, lambda = c(0,0,0,0)))
  expect_error(kproto(x, 4, lambda = 0))
  expect_is(kproto(x, 4, lambda = c(0,0,0,1)), "kproto")
  expect_is(kproto(x, 4, lambda = c(1,0,0,0)), "kproto")
  expect_is(kproto(x, 4, lambda = c(0,0,1,1)), "kproto")
  expect_is(kproto(x, 4, lambda = c(1,1,0,0)), "kproto")
}
)


test_that("lambdaest does what it is supposed to do",{
  expect_is(lambdaest(x), "numeric")
  expect_is(lambdaest(x, outtype = "vector"), "numeric")
  expect_equal(length(lambdaest(x, outtype = "vector")), ncol(x))                         # correct dimension of lambda 
  expect_is(lambdaest(x, outtype = "variation"), "numeric")
  expect_equal(length(lambdaest(x, outtype = "variation")), ncol(x))
  expect_equal(lambdaest(x, outtype = "vector"), 1/lambdaest(x, outtype = "variation"))   # compare lambdas variations 
  expect_equal(kproto(x, 2, lambda = lambdaest(x))$lambda, lambdaest(x))                                          # estimated lambda used by kproto 
  expect_equal(kproto(x, 2, lambda = lambdaest(x, outtype = "vector"))$lambda, lambdaest(x, outtype = "vector"))  # estimated lambda used by kproto  
}
)

test_that("also works for 1D",{
    expect_is(kproto(x[,c(1,3)], 4), "kproto")
    expect_error(kproto(x[,1:2], 4))
    expect_error(kproto(x[,3:4], 4))
}
)

x2 <- data.frame(V1 = factor(c("A","A","B","B")), V2 = c(-3,-1,1,3)) 
kpres <- kproto(x2, x2[c(1,3),], lambda = 1)
test_that("distances are computed correctly",{
  expect_equal(as.numeric(kpres$dists[,1]), c(1,1,10,26))
}
)



# test missings and imputation

x3 <- x4 <- x
x3$x1[1] <- NA
x4$x1[1] <- x4$x2[1] <- x4$x3[1] <- x4$x4[1] <-  NA
test_that("message for NAs.",{
  expect_message(kproto(x3, 4),"Observations with NAs are removed.")
  expect_warning(kproto(x4, 4, na.rm = "no"),"No meaningful cluster assignment possible for observations where all variables NA.")
  expect_error(kproto(x3, 4, na.rm = "not"), "Argument na.rm must be either 'yes','no','imp.internal' or 'imp.onestep'!")
  expect_message(kproto(x3, 4, na.rm = TRUE),"Logical input for na.rm is deprecated. Please use either 'yes','no','imp.internal' or 'imp.onestep'.\n")
  expect_error(kproto(x3, 4, na.rm = "imp.internal", type = "gower"), "Argument na.rm must be either 'yes','no' or 'imp-onestep', since imp.internal is not yet implemented for type = 'gower'!")}
)

prototypes <- data.frame(V1 = factor(c("A","B")), V2 = c(-3,3)) 
x5 <- data.frame(V1 = factor(c(rep("A",10),rep("B",10))), V2 = c(rep(-3, 5), rep(-5, 5), rep(NA, 10))) 
kpres <- kproto(x = x5, k = prototypes, na.rm = "no")
test_that("handling all NAs in variable in a cluster.",{
  expect_equal(kpres$centers[2,2], prototypes[2,2])}
)

x6 <- data.frame(V1 = factor(c(rep("A",10),rep("B",10))), V2 = c(rep(-3, 5), rep(-5, 5), rep(5, 5), rep(NA, 5))) 
kpres_i <- kproto(x = x6, k = prototypes, na.rm = "imp.internal")
kpres_ios <- kproto(x = x6, k = prototypes, na.rm = "imp.onestep")
test_that("application of different imputation strategies",{
  expect_is(kpres_i, "kproto")
  expect_is(kpres_i$cluster, "integer")
  expect_is(kpres_i$dists, "matrix")
  expect_false(any(is.na(kpres_i$data)))
  expect_is(kpres_ios, "kproto")
  expect_is(kpres_ios$cluster, "integer")
  expect_is(kpres_ios$dists, "matrix")
  expect_false(any(is.na(kpres_ios$data)))}
)



# test initialization
kpres <- kproto(x, 4, init = "nbh.dens", verbose = FALSE)
test_that("testing kproto with selected initialization nbh.dens",{
  expect_is(kpres$nstart.m, "NULL")
  expect_is(kpres$inits, "data.frame")
  expect_equal(nrow(kpres$inits), 4)
  expect_equal(ncol(kpres$inits), ncol(x))
}
)
kpres <- kproto(x, 4, init = "nstart.m", verbose = FALSE)
test_that("testing kproto with selected initialization nstart.m",{
  expect_is(kpres$inits, "NULL")
  expect_is(kpres$nstart.m, "numeric")
}
)






# test gower extension
kpres <- kproto(x = x, k = 4, type = "gower")
test_that("Type = gower can be called instead of standard.",{
  expect_equal(kpres$type, "gower")}
)

clusid <- rep(1:4, each = n)
muk    <- 2
# numeric
mus <- c(rep(-muk, n), rep(-muk, n), rep(muk, n), rep(muk, n))
x1  <- rnorm(4*n) + mus
# ordered factor
mus <- c(rep(-muk, n),rep(muk, n),rep(-muk, n),rep(muk, n))
x2 <- rnorm(4*n) + mus
quants <- quantile(x2, seq(0, 1, length.out = (8+1)))
quants[1] <- -Inf
quants[length(quants)] <- Inf
x2 <- as.ordered(cut(x2, quants))
x <- data.frame(x1, x2)

kpres <- kproto(x = x, k = 4)
test_that("Standard still works if ordered factors.",{
  expect_equal(kpres$type, "huang")}
)

kpres <- kproto(x = x, k = 4, type = "standard")
test_that("Standard still works if ordered factors.",{
  expect_equal(kpres$type, "huang")}
)

kpres <- kproto(x = x, k = 4, type = "gower")
test_that("Type = gower will be called.",{
  expect_equal(kpres$type, "gower")}
)

xx <- x
xx[,3] <- factor(x[,2], ordered = F)
kpres <- kproto(x = xx, k = 4, type = "gower")
test_that("Prototypes for ordered factor variables using gower extension.",{
  expect_equal(class(kpres$centers[,2])[1], "ordered")}
)

test_that("Prototypes for (unordered) factor variables using gower extension.",{
  expect_equal(class(kpres$centers[,3]), "factor")}
)

data(iris)
model <- kproto(x = iris, k = 3)
pred  <- predict(model, iris[1, ])
test_that("Prediction for data frame with one single observation.",{
  expect_equal(class(pred$dists)[1], "matrix")}
)

iris2 <- iris
iris2$Species <- as.ordered(iris$Species)
model <- kproto(x = iris2, k = 3, type = "gower")
pred <- predict(model, iris2[1, ])
test_that("Prediction for data frame with one single observation for type = gower.",{
  expect_equal(class(pred$dists)[1], "matrix")}
)

kpres <- kproto(x = x, k = x[1:4,])
test_that("kproto_gower works for is.data.frame(k).",{
  expect_equal(class(kpres), "kproto")}
)


