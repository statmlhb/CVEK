context("Model definition")

test_that("data and kernel parameters are dataframes", {
  n <- 50
  label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4"))
  data <- generate_data(n, fixed_num = 1, label_names, method = "rbf", 
                        int_effect = .3, l = 1, d = 2, eps = .01)
  expect_is(data, "data.frame")
  expect_is(kern_par, "data.frame")
})


test_that("type of output", {
  formula <- Y ~ X + Z1 + Z2
  n <- 50
  label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4"))
  kern_par <- data.frame(method = c("rbf", "polynomial", "matern"), 
                         l = c(.5, 1, 1.5), d = 1:3)
  kern_par$method <- as.character(kern_par$method)
  data <- generate_data(n, fixed_num = 1, label_names, method = "rbf", 
                        int_effect = .3, l = 1, d = 2, eps = .01)
  fit <- define_model(formula, data, kern_par, fixed_num = 1, label_names)
  expect_is(fit$Y, "numeric")
  expect_is(fit$Z1, "matrix")
  expect_is(fit$kern_list, "list")
})



