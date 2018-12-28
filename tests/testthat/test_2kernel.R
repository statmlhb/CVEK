context("Kernel Generation")

test_that("kernel implementation", {
  formula <- Y ~ X + Z1 + Z2
  n <- 100
  label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4"))
  kern_par <- data.frame(method = c("rbf", "polynomial", "matern"), 
                         l = c(.5, 1, 1.5), d = 1:3)
  kern_par$method <- as.character(kern_par$method)
  data <- generate_data(n, fixed_num = 1, label_names, method = "rbf", 
                        int_effect = .3, l = 1, d = 2, eps = .01)
  fit <- define_model(formula, data, kern_par, fixed_num = 1, label_names)
  K1 <- fit$kern_list[[1]]
  Kmat1 <- round(K1(Z1, Z1), 5)
  Ktest1 <- rbfdot(sigma = 2)
  Kmat_test1 <- round(kernelMatrix(Ktest1, Z1)@.Data, 5)
  expect_equal(Kmat1, Kmat_test1)
  
  K2 <- fit$kern_list[[2]]
  Kmat2 <- round(K2(Z1, Z1), 5)
  Ktest2 <- polydot(degree = 2, scale = 1, offset = 1)
  Kmat_test2 <- round(kernelMatrix(Ktest2, Z1)@.Data, 5)
  expect_equal(Kmat2, Kmat_test2)

})

