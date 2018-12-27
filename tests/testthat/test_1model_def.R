context("Model definition")

test_that("data and kernel parameters are dataframes", {
  n <- 100
  label_names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))
  data <- generate_data(n, label_names, method = "rbf", 
                        int_effect = .3, l = 1, eps = .01)
  expect_is(dora, "data.frame")
  expect_is(kern_par, "data.frame")
  expect_is(data, "data.frame")
})


test_that("type of output", {
  formula <- Y ~ X1 + X2
  fit <- define_model(formula, label_names, data = dora, kern_par)
  expect_is(fit$Y, "numeric")
  expect_is(fit$X1, "matrix")
  expect_is(fit$kern_list, "list")
})



