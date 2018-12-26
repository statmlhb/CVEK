context("Kernel Generation")

test_that("kernel implementation", {
  # kern_par <- data.frame(method = c("rbf", "polynomial", "matern"), 
  #                        Sigma = rep(0, 3), l = c(.5, 1, 1.5), p = 1:3)
  K1 <- fit$kern_list[[1]]
  Kmat1 <- round(K1(X1, X1), 5)
  Ktest1 <- rbfdot(sigma = 2)
  Kmat_test1 <- round(kernelMatrix(Ktest1, X1)@.Data, 5)
  expect_equal(Kmat1, Kmat_test1)
  
  K2 <- fit$kern_list[[2]]
  Kmat2 <- round(K2(X1, X1), 5)
  Ktest2 <- polydot(degree = 2, scale = 1, offset = 1)
  Kmat_test2 <- round(kernelMatrix(Ktest2, X1)@.Data, 5)
  expect_equal(Kmat2, Kmat_test2)
  
  # K3 <- fit$kern_list[[3]]
  # Kmat3 <- round(K3(X1, X1), 5)
  # Ktest3 <- rbfdot(sigma = 2)
  # Kmat_test3 <- round(kernelMatrix(Ktest3, X1), 5)
  # expect_equal(Kmat3, Kmat_test3)
})

