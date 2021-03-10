context("clustervalidation - basic tests.\n")

library(clustMixType)

# generate test data from example
set.seed(42)
n   <- 10
prb <- 0.99
muk <- 2.5 
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

kpres <- kproto.default(x, 4, verbose = FALSE)
kpres_w <- kproto.default(x, 4, keep.data = FALSE, verbose = FALSE)
kpres_1 <- kproto.default(x, 1, keep.data = FALSE, verbose = FALSE)



test_that("checking input objects",{
  expect_error(validation_kproto(method = "other index", object = 2, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", object = 2, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", object = kpres, k = 4, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = NULL, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = data.frame("a" = 2,"b" = 2), verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = data.frame("a" = c(2, 2)), verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = x, k = -1:3, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = x, k = 5:43, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = NULL, object = kpres_w, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = NULL, object = kpres_1, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = x, k = 1, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = x, k = 1:3, verbose = FALSE))
}
)


context("clustervalidation - cindex.\n")
iv_cindex  <- validation_kproto(method = "cindex", object = kpres, verbose = FALSE)
ivk_cindex <- validation_kproto(method = "cindex", data = x, k = c(3:5), verbose = FALSE)
ivk_l_cindex <- validation_kproto(method = "cindex", object = NULL, data = x, k = c(3:5), 
                                  kp_obj = 'optimal', lambda = 2, verbose = FALSE)
test_that("resulting objects based on cindex are as expected",{
  
  expect_is(iv_cindex, "numeric")
  expect_true(length(iv_cindex) == 1)
  expect_true(iv_cindex > 0)
  expect_true(iv_cindex < 1)
  
  expect_is(ivk_cindex[[1]], "integer")
  expect_true(length(ivk_cindex[[1]]) == 1)
  expect_true(ivk_cindex[[1]] %in% c(3:5))
  expect_true(is.numeric(ivk_cindex[[2]]))
  expect_true(is.numeric(ivk_cindex[[3]]))
  expect_true(length(ivk_cindex[[3]]) <= 3)
  expect_equal(ivk_cindex[[1]], as.integer(names(which.min(ivk_cindex[[3]]))))
  
  expect_equal(ivk_l_cindex$kp_obj$lambda, 2)
}
)


context("clustervalidation - dunn.\n")
iv_dunn  <- validation_kproto(method = "dunn", object = kpres, verbose = FALSE)
ivk_dunn <- validation_kproto(method = "dunn", data = x, k = c(3:5), verbose = FALSE)
ivk_l_dunn <- validation_kproto(method = "dunn", object = NULL, data = x, k = c(3:5), 
                                  kp_obj = 'optimal', lambda = 2, verbose = FALSE)
test_that("resulting objects based on dunn index are as expected",{
  
  expect_is(iv_dunn, "numeric")
  expect_true(length(iv_dunn) == 1)
  expect_true(iv_dunn >= 0)
  
  expect_is(ivk_dunn[[1]], "integer")
  expect_true(length(ivk_dunn[[1]]) == 1)
  expect_true(ivk_dunn[[1]] %in% c(3:5))
  expect_true(is.numeric(ivk_dunn[[2]]))
  expect_true(is.numeric(ivk_dunn[[3]]))
  expect_true(length(ivk_dunn[[3]]) <= 3)
  expect_equal(ivk_dunn[[1]], as.integer(names(which.max(ivk_dunn[[3]]))))
  
  expect_equal(ivk_l_dunn$kp_obj$lambda, 2)
}
)


context("clustervalidation - gamma.\n")
iv_gamma  <- validation_kproto(method = "gamma", object = kpres, verbose = FALSE)
ivk_gamma <- validation_kproto(method = "gamma", data = x, k = c(3:5), verbose = FALSE)
ivk_l_gamma <- validation_kproto(method = "gamma", object = NULL, data = x, k = c(3:5), 
                                kp_obj = 'optimal', lambda = 2, verbose = FALSE)
test_that("resulting objects based on gamma index are as expected",{
  
  expect_is(iv_gamma, "numeric")
  expect_true(length(iv_gamma) == 1)
  expect_true(iv_gamma > -1)
  expect_true(iv_gamma <= 1)
  
  expect_is(ivk_gamma[[1]], "integer")
  expect_true(length(ivk_gamma[[1]]) == 1)
  expect_true(ivk_gamma[[1]] %in% c(3:5))
  expect_true(is.numeric(ivk_gamma[[2]]))
  expect_true(is.numeric(ivk_gamma[[3]]))
  expect_true(length(ivk_gamma[[3]]) <= 3)
  expect_equal(ivk_gamma[[1]], as.integer(names(which.max(ivk_gamma[[3]]))))
  
  expect_equal(ivk_l_gamma$kp_obj$lambda, 2)
}
)


context("clustervalidation - gplus.\n")
iv_gplus  <- validation_kproto(method = "gplus", object = kpres, verbose = FALSE)
ivk_gplus <- validation_kproto(method = "gplus", data = x, k = c(3:5), verbose = FALSE)
ivk_l_gplus <- validation_kproto(method = "gplus", object = NULL, data = x, k = c(3:5), 
                            kp_obj = 'optimal', lambda = 2, verbose = FALSE)
test_that("resulting objects based on gplus index are as expected",{
  
  expect_is(iv_gplus, "numeric")
  expect_true(length(iv_gplus) == 1)
  expect_true(iv_gplus >= 0)
  expect_true(iv_gplus < 1)
  
  expect_is(ivk_gplus[[1]], "integer")
  expect_true(length(ivk_gplus[[1]]) == 1)
  expect_true(ivk_gplus[[1]] %in% c(3:5))
  expect_true(is.numeric(ivk_gplus[[2]]))
  expect_true(is.numeric(ivk_gplus[[3]]))
  expect_true(length(ivk_gplus[[3]]) <= 3)
  expect_equal(ivk_gplus[[1]], as.integer(names(which.min(ivk_gplus[[3]]))))
  
  expect_equal(ivk_l_gplus$kp_obj$lambda, 2)
}
)


context("clustervalidation - mcclain.\n")
iv_mcclain  <- validation_kproto(method = "mcclain", object = kpres, verbose = FALSE)
ivk_mcclain <- validation_kproto(method = "mcclain", data = x, k = c(3:5), verbose = FALSE)
ivk_l_mcclain <- validation_kproto(method = "mcclain", object = NULL, data = x, k = c(3:5), 
                            kp_obj = 'optimal', lambda = 2, verbose = FALSE)
test_that("resulting objects based on mcclain index are as expected",{
  
  expect_is(iv_mcclain, "numeric")
  expect_true(length(iv_mcclain) == 1)
  expect_true(iv_mcclain >= 0)
  
  expect_is(ivk_mcclain[[1]], "integer")
  expect_true(length(ivk_mcclain[[1]]) == 1)
  expect_true(ivk_mcclain[[1]] %in% c(3:5))
  expect_true(is.numeric(ivk_mcclain[[2]]))
  expect_true(is.numeric(ivk_mcclain[[3]]))
  expect_true(length(ivk_mcclain[[3]]) <= 3)
  expect_equal(ivk_mcclain[[1]], as.integer(names(which.min(ivk_mcclain[[3]]))))
  
  expect_equal(ivk_l_mcclain$kp_obj$lambda, 2)
}
)


context("clustervalidation - ptbiserial.\n")
iv_ptbiserial  <- validation_kproto(method = "ptbiserial", object = kpres, verbose = FALSE)
ivk_ptbiserial <- validation_kproto(method = "ptbiserial", data = x, k = c(3:5), verbose = FALSE)
ivk_l_ptbiserial <- validation_kproto(method = "ptbiserial", object = NULL, data = x, k = c(3:5), 
                            kp_obj = 'optimal', lambda = 2, verbose = FALSE)
test_that("resulting objects based on ptbiserial index are as expected",{
  
  expect_is(iv_ptbiserial, "numeric")
  expect_true(length(iv_ptbiserial) == 1)
  
  expect_is(ivk_ptbiserial[[1]], "integer")
  expect_true(length(ivk_ptbiserial[[1]]) == 1)
  expect_true(ivk_ptbiserial[[1]] %in% c(3:5))
  expect_true(is.numeric(ivk_ptbiserial[[2]]))
  expect_true(is.numeric(ivk_ptbiserial[[3]]))
  expect_true(length(ivk_ptbiserial[[3]]) <= 3)
  expect_equal(ivk_ptbiserial[[1]], as.integer(names(which.max(ivk_ptbiserial[[3]]))))
  
  expect_equal(ivk_l_ptbiserial$kp_obj$lambda, 2)
}
)


context("clustervalidation - silhouette.\n")
iv_silhouette  <- validation_kproto(method = "silhouette", object = kpres, verbose = FALSE)
ivk_silhouette <- validation_kproto(method = "silhouette", data = x, k = c(3:5), verbose = FALSE)
ivk_l_silhouette <- validation_kproto(method = "silhouette", object = NULL, data = x, k = c(3:5), 
                            kp_obj = 'optimal', lambda = 2, verbose = FALSE)
test_that("resulting objects based on silhouette index are as expected",{
  
  expect_is(iv_silhouette, "numeric")
  expect_true(length(iv_silhouette) == 1)
  expect_true(iv_silhouette >= -1)
  expect_true(iv_silhouette <= 1)
  
  expect_is(ivk_silhouette[[1]], "integer")
  expect_true(length(ivk_silhouette[[1]]) == 1)
  expect_true(ivk_silhouette[[1]] %in% c(3:5))
  expect_true(is.numeric(ivk_silhouette[[2]]))
  expect_true(is.numeric(ivk_silhouette[[3]]))
  expect_true(length(ivk_silhouette[[3]]) <= 3)
  expect_equal(ivk_silhouette[[1]], as.integer(names(which.max(ivk_silhouette[[3]]))))
  
  expect_equal(ivk_l_silhouette$kp_obj$lambda, 2)
}
)


context("clustervalidation - tau.\n")
iv_tau  <- validation_kproto(method = "tau", object = kpres, verbose = FALSE)
ivk_tau <- validation_kproto(method = "tau", data = x, k = c(3:5), verbose = FALSE)
ivk_l_tau <- validation_kproto(method = "tau", object = NULL, data = x, k = c(3:5), 
                            kp_obj = 'optimal', lambda = 2, verbose = FALSE)
test_that("resulting objects based on tau index are as expected",{
  
  expect_is(iv_tau, "numeric")
  expect_true(length(iv_tau) == 1)
  expect_true(iv_tau >= -1)
  expect_true(iv_tau <= 1)
  
  expect_is(ivk_tau[[1]], "integer")
  expect_true(length(ivk_tau[[1]]) == 1)
  expect_true(ivk_tau[[1]] %in% c(3:5))
  expect_true(is.numeric(ivk_tau[[2]]))
  expect_true(is.numeric(ivk_tau[[3]]))
  expect_true(length(ivk_tau[[3]]) <= 3)
  expect_equal(ivk_tau[[1]], as.integer(names(which.max(ivk_tau[[3]]))))
  
  expect_equal(ivk_l_tau$kp_obj$lambda, 2)
}
)










