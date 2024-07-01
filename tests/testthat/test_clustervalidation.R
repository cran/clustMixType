context("clustervalidation - basic tests.\n")

#library(clustMixType)

# generate test data
set.seed(42)
n   <- 10
prb <- 0.98
muk <- 3
clusid <- rep(1:4, each = n)

x1 <- sample(c("A","B"), 4*n, replace = TRUE, prob = c(prb, 1-prb))
x1 <- as.factor(x1)
x2 <- sample(c("A","B"), n, replace = TRUE, prob = c(1-prb, prb))
x2 <- c(x2, sample(c("A","B"), 3*n, replace = TRUE, prob = c(prb, 1-prb)))
x2 <- as.factor(x2)

x3 <- rnorm(4*n, mean = -muk)
x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(2*n, mean = -muk))

x5 <- rnorm(4*n) - muk
x6 <- rnorm(4*n) + c(rep(-muk, n), rep(-muk, n), rep(muk, n), rep(-muk, n))
x5 <- as.ordered(cut(x5, c(-Inf, -2, 1, Inf)))
x6 <- as.ordered(cut(x6, c(-Inf, -2, 1, Inf)))


x_num_cat <- data.frame(x1,x2,x3,x4)
x_num_cat_ord <- data.frame(x1,x2,x3,x4,x5,x6)
remove(n, prb, muk, clusid, x1, x2, x3, x4, x5, x6)

kpres_nc_huang <- clustMixType::kproto(x_num_cat, 4, verbose = FALSE, type = "huang")
kpres_nc_huang_wo <- clustMixType::kproto(x_num_cat, 4, verbose = FALSE, type = "huang", keep.data = FALSE)
kpres_nc_huang_wo_1 <- clustMixType::kproto(x_num_cat, 1, verbose = FALSE, type = "huang", keep.data = FALSE)
kpres_nc_gower <- clustMixType::kproto(x_num_cat, 4, verbose = FALSE, type = "gower")
kpres_nco_gower <- clustMixType::kproto(x_num_cat_ord, 4, verbose = FALSE, type = "gower")
kpres_nco_gower_wo <- clustMixType::kproto(x_num_cat_ord, 4, verbose = FALSE, type = "gower", keep.data = FALSE)
kpres_nco_gower_wo_1 <- clustMixType::kproto(x_num_cat_ord, 1, verbose = FALSE, type = "gower", keep.data = FALSE)



test_that("checking input objects",{
  #input method
  expect_error(validation_kproto(method = "other index", object = 2, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", object = 2, verbose = FALSE))

  # input data
  expect_error(validation_kproto(method = "cindex", data = NULL, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = data.frame("a" = 2,"b" = 2), verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = data.frame("a" = c(2, 2)), verbose = FALSE))

  # input k
  expect_error(validation_kproto(method = "cindex", data = x_num_cat, k = 4, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = x_num_cat, k = -1:3, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = x_num_cat, k = 5:43, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = x_num_cat, k = 1, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = x_num_cat, k = 1:3, verbose = FALSE))

  #input object
  expect_error(validation_kproto(method = "cindex", data = NULL, object = kpres_nc_huang_wo, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = NULL, object = kpres_nco_gower_wo, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = NULL, object = kpres_nc_huang_wo_1, verbose = FALSE))
  expect_error(validation_kproto(method = "cindex", data = NULL, object = kpres_nco_gower_wo_1, verbose = FALSE))
}
)


context("clustervalidation - cindex.\n")
iv_nc_huang  <- validation_kproto(method = "cindex", object = kpres_nc_huang, verbose = FALSE)
iv_nc_gower  <- validation_kproto(method = "cindex", object = kpres_nc_gower, verbose = FALSE)
iv_nco_gower  <- validation_kproto(method = "cindex", object = kpres_nco_gower, verbose = FALSE)

ivk_nc_huang <- validation_kproto(method = "cindex", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "huang")
ivk_nc_gower <- validation_kproto(method = "cindex", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "gower")
ivk_nco_gower <- validation_kproto(method = "cindex", data = x_num_cat_ord, k = c(3:4), verbose = FALSE, type = "gower")

ivk_l_nc <- validation_kproto(method = "cindex", object = NULL, data = x_num_cat, k = c(3:4), kp_obj = 'optimal',
                                   lambda = 2, verbose = FALSE, type = "huang")
ivk_l_nco <- validation_kproto(method = "cindex", object = NULL, data = x_num_cat_ord, k = c(3:4), kp_obj = 'optimal',
                                       lambda = 2, verbose = FALSE, type = "gower")
ivk_huang_c <- validation_kproto(method = "cindex", object = NULL, data = x_num_cat, k = c(3:4),
                                        lambda = c(1,1,0,0), kp_obj = 'optimal', verbose = FALSE, type = "huang")
ivk_huang_n <- validation_kproto(method = "cindex", object = NULL, data = x_num_cat, k = c(3:4),
                                        lambda = c(0,0,1,1), kp_obj = 'optimal', verbose = FALSE, type = "huang")


test_that("resulting objects based on cindex are as expected",{
  ###validation index determination
  #nc_huang
  expect_is(iv_nc_huang, "numeric")
  expect_true(length(iv_nc_huang) == 1)
  expect_true(iv_nc_huang > 0)
  expect_true(iv_nc_huang < 1)
  #nc_gower
  expect_is(iv_nc_gower, "numeric")
  expect_true(length(iv_nc_gower) == 1)
  expect_true(iv_nc_gower > 0)
  expect_true(iv_nc_gower < 1)
  #nco_gower
  expect_is(iv_nco_gower, "numeric")
  expect_true(length(iv_nco_gower) == 1)
  expect_true(iv_nco_gower > 0)
  expect_true(iv_nco_gower < 1)

  ###determination index-optimal k
  #nc_huang
  expect_is(ivk_nc_huang[[1]], "integer")
  expect_true(length(ivk_nc_huang[[1]]) == 1)
  expect_true(ivk_nc_huang[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_huang[[2]]))
  expect_true(is.numeric(ivk_nc_huang[[3]]))
  expect_true(length(ivk_nc_huang[[3]]) <= 3)
  expect_equal(ivk_nc_huang[[1]], as.integer(names(which.min(ivk_nc_huang[[3]]))))
  #nc_gower
  expect_is(ivk_nc_gower[[1]], "integer")
  expect_true(length(ivk_nc_gower[[1]]) == 1)
  expect_true(ivk_nc_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_gower[[2]]))
  expect_true(is.numeric(ivk_nc_gower[[3]]))
  expect_true(length(ivk_nc_gower[[3]]) <= 3)
  expect_equal(ivk_nc_gower[[1]], as.integer(names(which.min(ivk_nc_gower[[3]]))))
  #nco_gower
  expect_is(ivk_nco_gower[[1]], "integer")
  expect_true(length(ivk_nco_gower[[1]]) == 1)
  expect_true(ivk_nco_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nco_gower[[2]]))
  expect_true(is.numeric(ivk_nco_gower[[3]]))
  expect_true(length(ivk_nco_gower[[3]]) <= 3)
  expect_equal(ivk_nco_gower[[1]], as.integer(names(which.min(ivk_nco_gower[[3]]))))

  #checking lambda
  expect_equal(ivk_l_nc$kp_obj$lambda, 2)
  expect_equal(ivk_l_nco$kp_obj$lambda, NULL)
  expect_true(length(ivk_huang_c) == 5)
  expect_true(length(ivk_huang_n) == 5)
}
)


context("clustervalidation - dunn.\n")
iv_nc_huang  <- validation_kproto(method = "dunn", object = kpres_nc_huang, verbose = FALSE)
iv_nc_gower  <- validation_kproto(method = "dunn", object = kpres_nc_gower, verbose = FALSE)
iv_nco_gower  <- validation_kproto(method = "dunn", object = kpres_nco_gower, verbose = FALSE)

ivk_nc_huang <- validation_kproto(method = "dunn", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "huang")
ivk_nc_gower <- validation_kproto(method = "dunn", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "gower")
ivk_nco_gower <- validation_kproto(method = "dunn", data = x_num_cat_ord, k = c(3:4), verbose = FALSE, type = "gower")

ivk_l_nc <- validation_kproto(method = "dunn", object = NULL, data = x_num_cat, k = c(3:4), kp_obj = 'optimal',
                              lambda = 2, verbose = FALSE, type = "huang")
ivk_l_nco <- validation_kproto(method = "dunn", object = NULL, data = x_num_cat_ord, k = c(3:4), kp_obj = 'optimal',
                               lambda = 2, verbose = FALSE, type = "gower")
ivk_huang_c <- validation_kproto(method = "dunn", object = NULL, data = x_num_cat, k = c(3:4),
                                 lambda = c(1,1,0,0), kp_obj = 'optimal', verbose = FALSE, type = "huang")
ivk_huang_n <- validation_kproto(method = "dunn", object = NULL, data = x_num_cat, k = c(3:4),
                                 lambda = c(0,0,1,1), kp_obj = 'optimal', verbose = FALSE, type = "huang")


test_that("resulting objects based on dunn index are as expected",{
  ###validation index determination
  #nc_huang
  expect_is(iv_nc_huang, "numeric")
  expect_true(length(iv_nc_huang) == 1)
  expect_true(iv_nc_huang >= 0)
  #nc_gower
  expect_is(iv_nc_gower, "numeric")
  expect_true(length(iv_nc_gower) == 1)
  expect_true(iv_nc_gower >= 0)
  #nco_gower
  expect_is(iv_nco_gower, "numeric")
  expect_true(length(iv_nco_gower) == 1)
  expect_true(iv_nco_gower >= 0)

  ###determination index-optimal k
  #nc_huang
  expect_is(ivk_nc_huang[[1]], "integer")
  expect_true(length(ivk_nc_huang[[1]]) == 1)
  expect_true(ivk_nc_huang[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_huang[[2]]))
  expect_true(is.numeric(ivk_nc_huang[[3]]))
  expect_true(length(ivk_nc_huang[[3]]) <= 3)
  expect_equal(ivk_nc_huang[[1]], as.integer(names(which.max(ivk_nc_huang[[3]]))))
  #nc_gower
  expect_is(ivk_nc_gower[[1]], "integer")
  expect_true(length(ivk_nc_gower[[1]]) == 1)
  expect_true(ivk_nc_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_gower[[2]]))
  expect_true(is.numeric(ivk_nc_gower[[3]]))
  expect_true(length(ivk_nc_gower[[3]]) <= 3)
  expect_equal(ivk_nc_gower[[1]], as.integer(names(which.max(ivk_nc_gower[[3]]))))
  #nco_gower
  expect_is(ivk_nco_gower[[1]], "integer")
  expect_true(length(ivk_nco_gower[[1]]) == 1)
  expect_true(ivk_nco_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nco_gower[[2]]))
  expect_true(is.numeric(ivk_nco_gower[[3]]))
  expect_true(length(ivk_nco_gower[[3]]) <= 3)
  expect_equal(ivk_nco_gower[[1]], as.integer(names(which.max(ivk_nco_gower[[3]]))))

  #checking lambda
  expect_equal(ivk_l_nc$kp_obj$lambda, 2)
  expect_equal(ivk_l_nco$kp_obj$lambda, NULL)
  expect_true(length(ivk_huang_c) == 5)
  expect_true(length(ivk_huang_n) == 5)
}
)


context("clustervalidation - gamma.\n")
iv_nc_huang  <- validation_kproto(method = "gamma", object = kpres_nc_huang, verbose = FALSE)
iv_nc_gower  <- validation_kproto(method = "gamma", object = kpres_nc_gower, verbose = FALSE)
iv_nco_gower  <- validation_kproto(method = "gamma", object = kpres_nco_gower, verbose = FALSE)

ivk_nc_huang <- validation_kproto(method = "gamma", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "huang")
ivk_nc_gower <- validation_kproto(method = "gamma", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "gower")
ivk_nco_gower <- validation_kproto(method = "gamma", data = x_num_cat_ord, k = c(3:4), verbose = FALSE, type = "gower")

ivk_l_nc <- validation_kproto(method = "gamma", object = NULL, data = x_num_cat, k = c(3:4), kp_obj = 'optimal',
                              lambda = 2, verbose = FALSE, type = "huang")
ivk_l_nco <- validation_kproto(method = "gamma", object = NULL, data = x_num_cat_ord, k = c(3:4), kp_obj = 'optimal',
                               lambda = 2, verbose = FALSE, type = "gower")
ivk_huang_c <- validation_kproto(method = "gamma", object = NULL, data = x_num_cat, k = c(3:4),
                                 lambda = c(1,1,0,0), kp_obj = 'optimal', verbose = FALSE, type = "huang")
ivk_huang_n <- validation_kproto(method = "gamma", object = NULL, data = x_num_cat, k = c(3:4),
                                 lambda = c(0,0,1,1), kp_obj = 'optimal', verbose = FALSE, type = "huang")

test_that("resulting objects based on gamma index are as expected",{
  ###validation index determination
  #nc_huang
  expect_is(iv_nc_huang, "numeric")
  expect_true(length(iv_nc_huang) == 1)
  expect_true(iv_nc_huang > -1)
  expect_true(iv_nc_huang <= 1)
  #nc_gower
  expect_is(iv_nc_gower, "numeric")
  expect_true(length(iv_nc_gower) == 1)
  expect_true(iv_nc_gower > -1)
  expect_true(iv_nc_gower <= 1)
  #nco_gower
  expect_is(iv_nco_gower, "numeric")
  expect_true(length(iv_nco_gower) == 1)
  expect_true(iv_nco_gower > -1)
  expect_true(iv_nco_gower <= 1)

  ###determination index-optimal k
  #nc_huang
  expect_is(ivk_nc_huang[[1]], "integer")
  expect_true(length(ivk_nc_huang[[1]]) == 1)
  expect_true(ivk_nc_huang[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_huang[[2]]))
  expect_true(is.numeric(ivk_nc_huang[[3]]))
  expect_true(length(ivk_nc_huang[[3]]) <= 3)
  expect_equal(ivk_nc_huang[[1]], as.integer(names(which.max(ivk_nc_huang[[3]]))))
  #nc_gower
  expect_is(ivk_nc_gower[[1]], "integer")
  expect_true(length(ivk_nc_gower[[1]]) == 1)
  expect_true(ivk_nc_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_gower[[2]]))
  expect_true(is.numeric(ivk_nc_gower[[3]]))
  expect_true(length(ivk_nc_gower[[3]]) <= 3)
  expect_equal(ivk_nc_gower[[1]], as.integer(names(which.max(ivk_nc_gower[[3]]))))
  #nco_gower
  expect_is(ivk_nco_gower[[1]], "integer")
  expect_true(length(ivk_nco_gower[[1]]) == 1)
  expect_true(ivk_nco_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nco_gower[[2]]))
  expect_true(is.numeric(ivk_nco_gower[[3]]))
  expect_true(length(ivk_nco_gower[[3]]) <= 3)
  expect_equal(ivk_nco_gower[[1]], as.integer(names(which.max(ivk_nco_gower[[3]]))))

  #checking lambda
  expect_equal(ivk_l_nc$kp_obj$lambda, 2)
  expect_equal(ivk_l_nco$kp_obj$lambda, NULL)
  expect_true(length(ivk_huang_c) == 5)
  expect_true(length(ivk_huang_n) == 5)
}
)


context("clustervalidation - gplus.\n")
iv_nc_huang  <- validation_kproto(method = "gplus", object = kpres_nc_huang, verbose = FALSE)
iv_nc_gower  <- validation_kproto(method = "gplus", object = kpres_nc_gower, verbose = FALSE)
iv_nco_gower  <- validation_kproto(method = "gplus", object = kpres_nco_gower, verbose = FALSE)

ivk_nc_huang <- validation_kproto(method = "gplus", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "huang")
ivk_nc_gower <- validation_kproto(method = "gplus", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "gower")
ivk_nco_gower <- validation_kproto(method = "gplus", data = x_num_cat_ord, k = c(3:4), verbose = FALSE, type = "gower")

ivk_l_nc <- validation_kproto(method = "gplus", object = NULL, data = x_num_cat, k = c(3:4), kp_obj = 'optimal',
                              lambda = 2, verbose = FALSE, type = "huang")
ivk_l_nco <- validation_kproto(method = "gplus", object = NULL, data = x_num_cat_ord, k = c(3:4), kp_obj = 'optimal',
                               lambda = 2, verbose = FALSE, type = "gower")
ivk_huang_c <- validation_kproto(method = "gplus", object = NULL, data = x_num_cat, k = c(3:4),
                                 lambda = c(1,1,0,0), kp_obj = 'optimal', verbose = FALSE, type = "huang")
ivk_huang_n <- validation_kproto(method = "gplus", object = NULL, data = x_num_cat, k = c(3:4),
                                 lambda = c(0,0,1,1), kp_obj = 'optimal', verbose = FALSE, type = "huang")

test_that("resulting objects based on gplus index are as expected",{
  ###validation index determination
  #nc_huang
  expect_is(iv_nc_huang, "numeric")
  expect_true(length(iv_nc_huang) == 1)
  expect_true(iv_nc_huang >= 0)
  expect_true(iv_nc_huang < 1)
  #nc_gower
  expect_is(iv_nc_gower, "numeric")
  expect_true(length(iv_nc_gower) == 1)
  expect_true(iv_nc_gower >= 0)
  expect_true(iv_nc_gower < 1)
  #nco_gower
  expect_is(iv_nco_gower, "numeric")
  expect_true(length(iv_nco_gower) == 1)
  expect_true(iv_nco_gower >= 0)
  expect_true(iv_nco_gower < 1)

  ###determination index-optimal k
  #nc_huang
  expect_is(ivk_nc_huang[[1]], "integer")
  expect_true(length(ivk_nc_huang[[1]]) == 1)
  expect_true(ivk_nc_huang[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_huang[[2]]))
  expect_true(is.numeric(ivk_nc_huang[[3]]))
  expect_true(length(ivk_nc_huang[[3]]) <= 3)
  expect_equal(ivk_nc_huang[[1]], as.integer(names(which.min(ivk_nc_huang[[3]]))))
  #nc_gower
  expect_is(ivk_nc_gower[[1]], "integer")
  expect_true(length(ivk_nc_gower[[1]]) == 1)
  expect_true(ivk_nc_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_gower[[2]]))
  expect_true(is.numeric(ivk_nc_gower[[3]]))
  expect_true(length(ivk_nc_gower[[3]]) <= 3)
  expect_equal(ivk_nc_gower[[1]], as.integer(names(which.min(ivk_nc_gower[[3]]))))
  #nco_gower
  expect_is(ivk_nco_gower[[1]], "integer")
  expect_true(length(ivk_nco_gower[[1]]) == 1)
  expect_true(ivk_nco_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nco_gower[[2]]))
  expect_true(is.numeric(ivk_nco_gower[[3]]))
  expect_true(length(ivk_nco_gower[[3]]) <= 3)
  expect_equal(ivk_nco_gower[[1]], as.integer(names(which.min(ivk_nco_gower[[3]]))))

  #checking lambda
  expect_equal(ivk_l_nc$kp_obj$lambda, 2)
  expect_equal(ivk_l_nco$kp_obj$lambda, NULL)
  expect_true(length(ivk_huang_c) == 5)
  expect_true(length(ivk_huang_n) == 5)
}
)


context("clustervalidation - mcclain.\n")
iv_nc_huang  <- validation_kproto(method = "mcclain", object = kpres_nc_huang, verbose = FALSE)
iv_nc_gower  <- validation_kproto(method = "mcclain", object = kpres_nc_gower, verbose = FALSE)
iv_nco_gower  <- validation_kproto(method = "mcclain", object = kpres_nco_gower, verbose = FALSE)

ivk_nc_huang <- validation_kproto(method = "mcclain", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "huang")
ivk_nc_gower <- validation_kproto(method = "mcclain", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "gower")
ivk_nco_gower <- validation_kproto(method = "mcclain", data = x_num_cat_ord, k = c(3:4), verbose = FALSE, type = "gower")

set.seed(42)
ivk_l_nc <- validation_kproto(method = "mcclain", object = NULL, data = x_num_cat, k = c(3:4), kp_obj = 'optimal', 
                              lambda = 2, verbose = FALSE, type = "huang")
ivk_l_nco <- validation_kproto(method = "mcclain", object = NULL, data = x_num_cat_ord, k = c(3:4), kp_obj = 'optimal', 
                               lambda = 2, verbose = FALSE, type = "gower")
ivk_huang_c <- validation_kproto(method = "mcclain", object = NULL, data = x_num_cat, k = c(3:4), 
                                 lambda = c(1,1,0,0), kp_obj = 'optimal', verbose = FALSE, type = "huang")
ivk_huang_n <- validation_kproto(method = "mcclain", object = NULL, data = x_num_cat, k = c(3:4), 
                                 lambda = c(0,0,1,1), kp_obj = 'optimal', verbose = FALSE, type = "huang")

test_that("resulting objects based on mcclain index are as expected",{
  ###validation index determination
  #nc_huang
  expect_is(iv_nc_huang, "numeric")
  expect_true(length(iv_nc_huang) == 1)
  expect_true(iv_nc_huang >= 0)
  #nc_gower
  expect_is(iv_nc_gower, "numeric")
  expect_true(length(iv_nc_gower) == 1)
  expect_true(iv_nc_gower >= 0)
  #nco_gower
  expect_is(iv_nco_gower, "numeric")
  expect_true(length(iv_nco_gower) == 1)
  expect_true(iv_nco_gower >= 0)
  
  ###determination index-optimal k
  #nc_huang
  expect_is(ivk_nc_huang[[1]], "integer")
  expect_true(length(ivk_nc_huang[[1]]) == 1)
  expect_true(ivk_nc_huang[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_huang[[2]]))
  expect_true(is.numeric(ivk_nc_huang[[3]]))
  expect_true(length(ivk_nc_huang[[3]]) <= 3)
  expect_equal(ivk_nc_huang[[1]], as.integer(names(which.min(ivk_nc_huang[[3]]))))
  #nc_gower
  expect_is(ivk_nc_gower[[1]], "integer")
  expect_true(length(ivk_nc_gower[[1]]) == 1)
  expect_true(ivk_nc_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_gower[[2]]))
  expect_true(is.numeric(ivk_nc_gower[[3]]))
  expect_true(length(ivk_nc_gower[[3]]) <= 3)
  expect_equal(ivk_nc_gower[[1]], as.integer(names(which.min(ivk_nc_gower[[3]]))))
  #nco_gower
  expect_is(ivk_nco_gower[[1]], "integer")
  expect_true(length(ivk_nco_gower[[1]]) == 1)
  expect_true(ivk_nco_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nco_gower[[2]]))
  expect_true(is.numeric(ivk_nco_gower[[3]]))
  expect_true(length(ivk_nco_gower[[3]]) <= 3)
  expect_equal(ivk_nco_gower[[1]], as.integer(names(which.min(ivk_nco_gower[[3]]))))
  
  #checking lambda
  expect_equal(ivk_l_nc$kp_obj$lambda, 2)
  expect_equal(ivk_l_nco$kp_obj$lambda, NULL)
  expect_true(length(ivk_huang_c) == 5)
  expect_true(length(ivk_huang_n) == 5)
}
)


context("clustervalidation - ptbiserial.\n")
iv_nc_huang  <- validation_kproto(method = "ptbiserial", object = kpres_nc_huang, verbose = FALSE)
iv_nc_gower  <- validation_kproto(method = "ptbiserial", object = kpres_nc_gower, verbose = FALSE)
iv_nco_gower  <- validation_kproto(method = "ptbiserial", object = kpres_nco_gower, verbose = FALSE)

ivk_nc_huang <- validation_kproto(method = "ptbiserial", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "huang")
ivk_nc_gower <- validation_kproto(method = "ptbiserial", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "gower")
ivk_nco_gower <- validation_kproto(method = "ptbiserial", data = x_num_cat_ord, k = c(3:4), verbose = FALSE, type = "gower")

ivk_l_nc <- validation_kproto(method = "ptbiserial", object = NULL, data = x_num_cat, k = c(3:4), kp_obj = 'optimal',
                              lambda = 2, verbose = FALSE, type = "huang")
ivk_l_nco <- validation_kproto(method = "ptbiserial", object = NULL, data = x_num_cat_ord, k = c(3:4), kp_obj = 'optimal',
                               lambda = 2, verbose = FALSE, type = "gower")
ivk_huang_c <- validation_kproto(method = "ptbiserial", object = NULL, data = x_num_cat, k = c(3:4),
                                 lambda = c(1,1,0,0), kp_obj = 'optimal', verbose = FALSE, type = "huang")
ivk_huang_n <- validation_kproto(method = "ptbiserial", object = NULL, data = x_num_cat, k = c(3:4),
                                 lambda = c(0,0,1,1), kp_obj = 'optimal', verbose = FALSE, type = "huang")

test_that("resulting objects based on ptbiserial index are as expected",{
  ###validation index determination
  #nc_huang
  expect_is(iv_nc_huang, "numeric")
  expect_true(length(iv_nc_huang) == 1)
  #nc_gower
  expect_is(iv_nc_gower, "numeric")
  expect_true(length(iv_nc_gower) == 1)
  #nco_gower
  expect_is(iv_nco_gower, "numeric")
  expect_true(length(iv_nco_gower) == 1)

  ###determination index-optimal k
  #nc_huang
  expect_is(ivk_nc_huang[[1]], "integer")
  expect_true(length(ivk_nc_huang[[1]]) == 1)
  expect_true(ivk_nc_huang[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_huang[[2]]))
  expect_true(is.numeric(ivk_nc_huang[[3]]))
  expect_true(length(ivk_nc_huang[[3]]) <= 3)
  expect_equal(ivk_nc_huang[[1]], as.integer(names(which.max(ivk_nc_huang[[3]]))))
  #nc_gower
  expect_is(ivk_nc_gower[[1]], "integer")
  expect_true(length(ivk_nc_gower[[1]]) == 1)
  expect_true(ivk_nc_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_gower[[2]]))
  expect_true(is.numeric(ivk_nc_gower[[3]]))
  expect_true(length(ivk_nc_gower[[3]]) <= 3)
  expect_equal(ivk_nc_gower[[1]], as.integer(names(which.max(ivk_nc_gower[[3]]))))
  #nco_gower
  expect_is(ivk_nco_gower[[1]], "integer")
  expect_true(length(ivk_nco_gower[[1]]) == 1)
  expect_true(ivk_nco_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nco_gower[[2]]))
  expect_true(is.numeric(ivk_nco_gower[[3]]))
  expect_true(length(ivk_nco_gower[[3]]) <= 3)
  expect_equal(ivk_nco_gower[[1]], as.integer(names(which.max(ivk_nco_gower[[3]]))))

  #checking lambda
  expect_equal(ivk_l_nc$kp_obj$lambda, 2)
  expect_equal(ivk_l_nco$kp_obj$lambda, NULL)
  expect_true(length(ivk_huang_c) == 5)
  expect_true(length(ivk_huang_n) == 5)
}
)


context("clustervalidation - silhouette.\n")
iv_nc_huang  <- validation_kproto(method = "silhouette", object = kpres_nc_huang, verbose = FALSE)
iv_nc_gower  <- validation_kproto(method = "silhouette", object = kpres_nc_gower, verbose = FALSE)
iv_nco_gower  <- validation_kproto(method = "silhouette", object = kpres_nco_gower, verbose = FALSE)

ivk_nc_huang <- validation_kproto(method = "silhouette", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "huang")
ivk_nc_gower <- validation_kproto(method = "silhouette", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "gower")
ivk_nco_gower <- validation_kproto(method = "silhouette", data = x_num_cat_ord, k = c(3:4), verbose = FALSE, type = "gower")

ivk_l_nc <- validation_kproto(method = "silhouette", object = NULL, data = x_num_cat, k = c(3:4), kp_obj = 'optimal',
                              lambda = 2, verbose = FALSE, type = "huang")
ivk_l_nco <- validation_kproto(method = "silhouette", object = NULL, data = x_num_cat_ord, k = c(3:4), kp_obj = 'optimal',
                               lambda = 2, verbose = FALSE, type = "gower")
ivk_huang_c <- validation_kproto(method = "silhouette", object = NULL, data = x_num_cat, k = c(3:4),
                                 lambda = c(1,1,0,0), kp_obj = 'optimal', verbose = FALSE, type = "huang")
ivk_huang_n <- validation_kproto(method = "silhouette", object = NULL, data = x_num_cat, k = c(3:4),
                                 lambda = c(0,0,1,1), kp_obj = 'optimal', verbose = FALSE, type = "huang")

test_that("resulting objects based on silhouette index are as expected",{
  ###validation index determination
  #nc_huang
  expect_is(iv_nc_huang, "numeric")
  expect_true(length(iv_nc_huang) == 1)
  expect_true(iv_nc_huang >= -1)
  expect_true(iv_nc_huang <= 1)
  #nc_gower
  expect_is(iv_nc_gower, "numeric")
  expect_true(length(iv_nc_gower) == 1)
  expect_true(iv_nc_gower >= -1)
  expect_true(iv_nc_gower <= 1)
  #nco_gower
  expect_is(iv_nco_gower, "numeric")
  expect_true(length(iv_nco_gower) == 1)
  expect_true(iv_nco_gower >= -1)
  expect_true(iv_nco_gower <= 1)

  ###determination index-optimal k
  #nc_huang
  expect_is(ivk_nc_huang[[1]], "integer")
  expect_true(length(ivk_nc_huang[[1]]) == 1)
  expect_true(ivk_nc_huang[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_huang[[2]]))
  expect_true(is.numeric(ivk_nc_huang[[3]]))
  expect_true(length(ivk_nc_huang[[3]]) <= 3)
  expect_equal(ivk_nc_huang[[1]], as.integer(names(which.max(ivk_nc_huang[[3]]))))
  #nc_gower
  expect_is(ivk_nc_gower[[1]], "integer")
  expect_true(length(ivk_nc_gower[[1]]) == 1)
  expect_true(ivk_nc_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_gower[[2]]))
  expect_true(is.numeric(ivk_nc_gower[[3]]))
  expect_true(length(ivk_nc_gower[[3]]) <= 3)
  expect_equal(ivk_nc_gower[[1]], as.integer(names(which.max(ivk_nc_gower[[3]]))))
  #nco_gower
  expect_is(ivk_nco_gower[[1]], "integer")
  expect_true(length(ivk_nco_gower[[1]]) == 1)
  expect_true(ivk_nco_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nco_gower[[2]]))
  expect_true(is.numeric(ivk_nco_gower[[3]]))
  expect_true(length(ivk_nco_gower[[3]]) <= 3)
  expect_equal(ivk_nco_gower[[1]], as.integer(names(which.max(ivk_nco_gower[[3]]))))

  #checking lambda
  expect_equal(ivk_l_nc$kp_obj$lambda, 2)
  expect_equal(ivk_l_nco$kp_obj$lambda, NULL)
  expect_true(length(ivk_huang_c) == 5)
  expect_true(length(ivk_huang_n) == 5)
}
)


context("clustervalidation - tau.\n")
iv_nc_huang  <- validation_kproto(method = "tau", object = kpres_nc_huang, verbose = FALSE)
iv_nc_gower  <- validation_kproto(method = "tau", object = kpres_nc_gower, verbose = FALSE)
iv_nco_gower  <- validation_kproto(method = "tau", object = kpres_nco_gower, verbose = FALSE)

ivk_nc_huang <- validation_kproto(method = "tau", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "huang")
ivk_nc_gower <- validation_kproto(method = "tau", data = x_num_cat, k = c(3:4), verbose = FALSE, type = "gower")
ivk_nco_gower <- validation_kproto(method = "tau", data = x_num_cat_ord, k = c(3:4), verbose = FALSE, type = "gower")

ivk_l_nc <- validation_kproto(method = "tau", object = NULL, data = x_num_cat, k = c(3:4), kp_obj = 'optimal',
                              lambda = 2, verbose = FALSE, type = "huang")
ivk_l_nco <- validation_kproto(method = "tau", object = NULL, data = x_num_cat_ord, k = c(3:4), kp_obj = 'optimal',
                               lambda = 2, verbose = FALSE, type = "gower")
ivk_huang_c <- validation_kproto(method = "tau", object = NULL, data = x_num_cat, k = c(3:4),
                                 lambda = c(1,1,0,0), kp_obj = 'optimal', verbose = FALSE, type = "huang")
ivk_huang_n <- validation_kproto(method = "tau", object = NULL, data = x_num_cat, k = c(3:4),
                                 lambda = c(0,0,1,1), kp_obj = 'optimal', verbose = FALSE, type = "huang")

test_that("resulting objects based on tau index are as expected",{
  ###validation index determination
  #nc_huang
  expect_is(iv_nc_huang, "numeric")
  expect_true(length(iv_nc_huang) == 1)
  expect_true(iv_nc_huang >= -1)
  expect_true(iv_nc_huang <= 1)
  #nc_gower
  expect_is(iv_nc_gower, "numeric")
  expect_true(length(iv_nc_gower) == 1)
  expect_true(iv_nc_gower >= -1)
  expect_true(iv_nc_gower <= 1)
  #nco_gower
  expect_is(iv_nco_gower, "numeric")
  expect_true(length(iv_nco_gower) == 1)
  expect_true(iv_nco_gower >= -1)
  expect_true(iv_nco_gower <= 1)

  ###determination index-optimal k
  #nc_huang
  expect_is(ivk_nc_huang[[1]], "integer")
  expect_true(length(ivk_nc_huang[[1]]) == 1)
  expect_true(ivk_nc_huang[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_huang[[2]]))
  expect_true(is.numeric(ivk_nc_huang[[3]]))
  expect_true(length(ivk_nc_huang[[3]]) <= 3)
  expect_equal(ivk_nc_huang[[1]], as.integer(names(which.max(ivk_nc_huang[[3]]))))
  #nc_gower
  expect_is(ivk_nc_gower[[1]], "integer")
  expect_true(length(ivk_nc_gower[[1]]) == 1)
  expect_true(ivk_nc_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nc_gower[[2]]))
  expect_true(is.numeric(ivk_nc_gower[[3]]))
  expect_true(length(ivk_nc_gower[[3]]) <= 3)
  expect_equal(ivk_nc_gower[[1]], as.integer(names(which.max(ivk_nc_gower[[3]]))))
  #nco_gower
  expect_is(ivk_nco_gower[[1]], "integer")
  expect_true(length(ivk_nco_gower[[1]]) == 1)
  expect_true(ivk_nco_gower[[1]] %in% c(3:4))
  expect_true(is.numeric(ivk_nco_gower[[2]]))
  expect_true(is.numeric(ivk_nco_gower[[3]]))
  expect_true(length(ivk_nco_gower[[3]]) <= 3)
  expect_equal(ivk_nco_gower[[1]], as.integer(names(which.max(ivk_nco_gower[[3]]))))

  #checking lambda
  expect_equal(ivk_l_nc$kp_obj$lambda, 2)
  expect_equal(ivk_l_nco$kp_obj$lambda, NULL)
  expect_true(length(ivk_huang_c) == 5)
  expect_true(length(ivk_huang_n) == 5)
}
)






context("checking stability_kproto.\n")

# stab_nc_huang <- stability_kproto(object = kpres_nc_huang, method = c("luxburg","fowlkesmallows"))
# stab_nc_gower <- stability_kproto(object = kpres_nc_gower, method = c("luxburg","fowlkesmallows"))
# stab_nco_gower <- stability_kproto(object = kpres_nco_gower, method = c("luxburg","fowlkesmallows"))

test_that("checking input objects of stability_kproto",{
  expect_error(stability_kproto(method = "stability", object = kpres_nc_huang, verbose = FALSE))
  expect_error(stability_kproto(method = c("luxburg", "jaccard"), object = 42, verbose = FALSE))
  expect_error(stability_kproto(method = "jaccard", object = kpres_nc_huang_wo, verbose = FALSE))
  expect_error(stability_kproto(method = "jaccard", object = kpres_nc_huang, verbose = FALSE, B = "forty-two"))

}
)

stab_1_nc_huang  <- stability_kproto(method = "jaccard", object = kpres_nc_huang, verbose = FALSE)
stab_1_nc_gower  <- stability_kproto(method = "jaccard", object = kpres_nc_gower, verbose = FALSE)
stab_1_nco_gower  <- stability_kproto(method = "jaccard", object = kpres_nco_gower, verbose = FALSE)
stab_2_nc_huang  <- stability_kproto(method = c("luxburg", "jaccard"), object = kpres_nc_huang, verbose = FALSE)
stab_2_nc_gower  <- stability_kproto(method = c("luxburg", "jaccard"), object = kpres_nc_gower, verbose = FALSE)
stab_2_nco_gower  <- stability_kproto(method = c("luxburg", "jaccard"), object = kpres_nco_gower, verbose = FALSE)
stab_4_nc_huang  <- stability_kproto(method = c("luxburg", "rand", "fowlkesmallows", "jaccard"), object = kpres_nc_huang, verbose = FALSE)
stab_4_nc_gower  <- stability_kproto(method = c("luxburg", "rand", "fowlkesmallows", "jaccard"), object = kpres_nc_gower, verbose = FALSE)
stab_4_nco_gower  <- stability_kproto(method = c("luxburg", "rand", "fowlkesmallows", "jaccard"), object = kpres_nco_gower, verbose = FALSE)

test_that("resulting stability objects are as expected",{

  #nc_huang
  expect_is(stab_1_nc_huang, "list")
  expect_true(length(stab_1_nc_huang) == 2)
  expect_true(length(stab_2_nc_huang) == 2)
  expect_true(length(stab_4_nc_huang) == 2)

  expect_is(stab_1_nc_huang[[1]], "numeric")
  expect_true(length(stab_1_nc_huang[[1]]) == 1)
  expect_true(length(stab_2_nc_huang[[1]]) == 2)
  expect_true(length(stab_4_nc_huang[[1]]) == 4)

  expect_true(dim(stab_1_nc_huang[[2]])[[2]] == 2)
  expect_true(dim(stab_2_nc_huang[[2]])[[2]] == 3)
  expect_true(dim(stab_4_nc_huang[[2]])[[2]] == 5)

  #nc_gower
  expect_is(stab_1_nc_gower, "list")
  expect_true(length(stab_1_nc_gower) == 2)
  expect_true(length(stab_2_nc_gower) == 2)
  expect_true(length(stab_4_nc_gower) == 2)

  expect_is(stab_1_nc_gower[[1]], "numeric")
  expect_true(length(stab_1_nc_gower[[1]]) == 1)
  expect_true(length(stab_2_nc_gower[[1]]) == 2)
  expect_true(length(stab_4_nc_gower[[1]]) == 4)

  expect_true(dim(stab_1_nc_gower[[2]])[[2]] == 2)
  expect_true(dim(stab_2_nc_gower[[2]])[[2]] == 3)
  expect_true(dim(stab_4_nc_gower[[2]])[[2]] == 5)

  #nco_gower
  expect_is(stab_1_nco_gower, "list")
  expect_true(length(stab_1_nco_gower) == 2)
  expect_true(length(stab_2_nco_gower) == 2)
  expect_true(length(stab_4_nco_gower) == 2)

  expect_is(stab_1_nco_gower[[1]], "numeric")
  expect_true(length(stab_1_nco_gower[[1]]) == 1)
  expect_true(length(stab_2_nco_gower[[1]]) == 2)
  expect_true(length(stab_4_nco_gower[[1]]) == 4)

  expect_true(dim(stab_1_nco_gower[[2]])[[2]] == 2)
  expect_true(dim(stab_2_nco_gower[[2]])[[2]] == 3)
  expect_true(dim(stab_4_nco_gower[[2]])[[2]] == 5)
}
)




















