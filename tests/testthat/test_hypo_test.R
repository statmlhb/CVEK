context("Hypothesis testing")

test_that("range of output", {
  lambda <- exp(seq(-5, 5))
  formula_int <- Y ~ X1 * X2
  test <- "boot"
  B <- 100
  pvalue <- testing(formula_int, label_names, fit$Y, fit$X1, fit$X2, fit$kern_list, 
                    mode = "loocv", strategy, beta = 1, test, lambda, B)
  expect_lte(pvalue, 1)
  expect_gte(pvalue, 0)
})

test_that("warning message from tuning", {
  lambda = rep(.5, 11)
  expect_warning(testing(formula_int, label_names, fit$Y, fit$X1, fit$X2, fit$kern_list, 
                         mode = "loocv", strategy, beta = 1, test, lambda, B), 
                 "the smallest one")
})
