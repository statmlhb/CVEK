context("Model estimation")

test_that("implementation of fixed effects", {
  names(longley) <- c("Y", paste0("z", 1:6))
  kern_par0 <- data.frame(method = "linear", l = 1, d = 2)
  kern_par0$method <- as.character(kern_par0$method)
  fit0 <- define_model(formula = Y ~ X + Z1 + Z2, 
                       longley, kern_par0, fixed_num = 0, 
                       label_names = list(Z1 = c("z1", "z2", "z3"), 
                                          Z2 = c("z4", "z5", "z6")))
  sol0 <- estimation(fit0$Y, fit0$X, fit0$Z1, fit0$Z2, 
                     fit0$kern_list, mode = "GCV", lambda = exp(seq(-10, 5)))
  y_fit0 <- sol0$K %*% sol0$alpha + sol0$beta[1, 1]
  mod0 <- lm.ridge(Y ~ ., longley, lambda = exp(seq(-10, 5)))
  whichIsBest <- which.min(mod0$GCV)
  y_mod0 <- as.matrix(cbind(const = 1, longley[, -1])) %*% coef(mod0)[whichIsBest,]
  diff <- sum((y_fit0 - y_mod0)^2) / length(y_fit0)
  expect_lte(diff, .01)
})


test_that("type of output", {
  formula <- Y ~ X + Z1 + Z2
  n <- 100
  label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4"))
  kern_par <- data.frame(method = c("rbf", "polynomial", "matern"), 
                         l = c(.5, 1, 1.5), d = 1:3)
  kern_par$method <- as.character(kern_par$method)
  data <- generate_data(n, fixed_num = 1, label_names, method = "rbf", 
                        int_effect = .3, l = 1, d = 2, eps = .01)
  fit <- define_model(formula, data, kern_par, fixed_num = 1, label_names)
  mode <- "loocv"
  strategy <- "stack"
  lambda <- exp(seq(-10, 5))
  sol <- estimation(fit$Y, fit$X, fit$Z1, fit$Z2, 
                    fit$kern_list, mode, strategy, lambda)
  expect_is(sol$lambda, "numeric")
  expect_is(sol$beta, "matrix")
  expect_is(sol$alpha, "matrix")
  expect_is(sol$K, "matrix")
  expect_is(sol$u_hat, "numeric")
})

test_that("length of output", {
  formula <- Y ~ X + Z1 + Z2
  n <- 100
  label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4"))
  kern_par <- data.frame(method = c("rbf", "polynomial", "matern"), 
                         l = c(.5, 1, 1.5), d = 1:3)
  kern_par$method <- as.character(kern_par$method)
  data <- generate_data(n, fixed_num = 1, label_names, method = "rbf", 
                        int_effect = .3, l = 1, d = 2, eps = .01)
  fit <- define_model(formula, data, kern_par, fixed_num = 1, label_names)
  mode <- "loocv"
  strategy <- "stack"
  lambda <- exp(seq(-10, 5))
  sol <- estimation(fit$Y, fit$X, fit$Z1, fit$Z2, 
                    fit$kern_list, mode, strategy, lambda)
  expect_length(sol$lambda, 1)
  expect_equal(nrow(sol$alpha), n)
  expect_equal(nrow(sol$K), n)
  expect_length(sol$u_hat, length(fit$kern_list))
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
  kern <- fit$kern_list[[1]]
  K1_m <- kern(fit$Z1, fit$Z1)
  K2_m <- kern(fit$Z2, fit$Z2)
  K <- (K1_m + K2_m) / tr(K1_m + K2_m)
  expect_warning(tuning(fit$Y, fit$X, K, mode = "loocv", lambda = rep(.5, 11)), 
                 "the smallest one")
})

