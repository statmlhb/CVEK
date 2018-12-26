context("Model estimation")

test_that("type of output", {
  mode <- "loocv"
  strategy <- "erm"
  lambda <- exp(seq(-5, 5))
  sol <- estimation(fit$Y, fit$X1, fit$X2, fit$kern_list, 
                    mode, strategy, lambda)
  expect_is(sol$lambda, "numeric")
  expect_is(sol$beta, "matrix")
  expect_is(sol$alpha, "matrix")
  expect_is(sol$K, "matrix")
  expect_is(sol$u_hat, "numeric")
})

test_that("length of output", {
  expect_length(sol$lambda, 1)
  expect_equal(nrow(sol$beta), 1)
  expect_equal(nrow(sol$alpha), n)
  expect_equal(nrow(sol$K), n)
  expect_length(sol$u_hat, length(fit$kern_list))
})

test_that("warning message from tuning", {
  kern <- kern_list[[1]]
  K1_m <- kern(fit$X1, fit$X1)
  K2_m <- kern(fit$X2, fit$X2)
  K <- (K1_m + K2_m) / tr(K1_m + K2_m)
  expect_warning(tuning(fit$Y, K, mode = "loocv", lambda = rep(.5, 11)), 
                 "the smallest one")
})

