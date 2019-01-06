context("Hypothesis testing")

test_that("range of output", {
  formula <- Y ~ X + Z1 + Z2
  n <- 100
  label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4"))
  kern_par <- data.frame(method = c("rbf", "polynomial", "matern"), 
                         l = c(.5, 1, 1.5), d = 1:3)
  kern_par$method <- as.character(kern_par$method)
  data <- generate_data(n, fixed_num = 1, label_names, method = "rbf", 
                        int_effect = .3, l = 1, d = 2, eps = .01)
  fit <- define_model(formula, data, kern_par, fixed_num = 1, label_names)
  strategy <- "stack"
  lambda <- exp(seq(-10, 5))
  formula_int = Y ~ X + Z1 * Z2
  test <- "boot"
  B <- 100
  pvalue <- testing(formula_int, label_names, fit$Y, fit$X, fit$Z1, fit$Z2, fit$kern_list, 
                    mode = "loocv", strategy, beta_exp = 1, test, lambda, B)$pvalue

  expect_lte(pvalue, 1)
  expect_gte(pvalue, 0)
})

test_that("warning message from tuning", {
  formula <- Y ~ X + Z1 + Z2
  n <- 100
  label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4"))
  kern_par <- data.frame(method = c("rbf", "polynomial", "matern"), 
                         l = c(.5, 1, 1.5), d = 1:3)
  kern_par$method <- as.character(kern_par$method)
  data <- generate_data(n, fixed_num = 1, label_names, method = "rbf", 
                        int_effect = .3, l = 1, d = 2, eps = .01)
  fit <- define_model(formula, data, kern_par, fixed_num = 1, label_names)
  strategy <- "stack"
  formula_int = Y ~ X + Z1 * Z2
  test <- "boot"
  B <- 100
  lambda = rep(.5, 11)
  expect_warning(testing(formula_int, label_names, fit$Y, fit$X, fit$Z1, fit$Z2, fit$kern_list, 
                         mode = "loocv", strategy, beta_exp = 1, test, lambda, B), 
                 "the smallest one")
})
