context("Kernel Generation")

test_that("intercept kernel", {
  point_wise <- function(xp, xq, l, d) 1
  intercept_kern <- function(X2, X1) 
    apply(X1, 1, function(xp){ 
      apply(X2, 1, function(xq){ 
        point_wise(xp, xq, l, d)
    })
  })
  Z1 <- rmvnorm(n = 50, mean = rep(0, 3), sigma = diag(3))
  Z2 <- rmvnorm(n = 52, mean = rep(0, 3), sigma = diag(3))
  Ktest <- round(intercept_kern(Z1, Z2), 5)
  intercept_fun <- generate_kernel("intercept")
  Kmat <- round(intercept_fun(Z1, Z2), 5)
})

test_that("linear kernel", {
  point_wise <- function(xp, xq, l, d) t(xp) %*% xq
  linear_kern <- function(X2, X1) 
    apply(X1, 1, function(xp){ 
      apply(X2, 1, function(xq){ 
        point_wise(xp, xq, l, d)
      })
    })
  Z1 <- rmvnorm(n = 50, mean = rep(0, 3), sigma = diag(3))
  Z2 <- rmvnorm(n = 52, mean = rep(0, 3), sigma = diag(3))
  Ktest <- round(linear_kern(Z1, Z2), 5)
  linear_fun <- generate_kernel("linear")
  Kmat <- round(linear_fun(Z1, Z2), 5)
})

test_that("polynomial kernel", {
  point_wise <- function(xp, xq, l, d) (t(xp) %*% xq + 1) ^ 3
  polynomial_kern <- function(X2, X1) 
    apply(X1, 1, function(xp){ 
      apply(X2, 1, function(xq){ 
        point_wise(xp, xq, l, d)
      })
    })
  Z1 <- rmvnorm(n = 50, mean = rep(0, 3), sigma = diag(3))
  Z2 <- rmvnorm(n = 52, mean = rep(0, 3), sigma = diag(3))
  Ktest <- round(polynomial_kern(Z1, Z2), 5)
  polynomial_fun <- generate_kernel("polynomial", d = 3)
  Kmat <- round(polynomial_fun(Z1, Z2), 5)
})


test_that("rbf kernel", {
  point_wise <- function(xp, xq, l, d) {
    exp(- sum((xp - xq) ^ 2) / (2 * 1 ^ 2))
    }
  rbf_kern <- function(X2, X1) 
    apply(X1, 1, function(xp){ 
      apply(X2, 1, function(xq){ 
        point_wise(xp, xq, l, d)
      })
    })
  Z1 <- rmvnorm(n = 50, mean = rep(0, 3), sigma = diag(3))
  Z2 <- rmvnorm(n = 52, mean = rep(0, 3), sigma = diag(3))
  Ktest <- round(rbf_kern(Z1, Z2), 5)
  rbf_fun <- generate_kernel("rbf", l = 1)
  Kmat <- round(rbf_fun(Z1, Z2), 5)
})

test_that("matern kernel", {
  point_wise <- function(xp, xq, l, d){
    l <- 1
    d <- 2
    r <- sqrt(sum((xp - xq) ^ 2))
    v <- d + 1 / 2
    s <- 0
    for (i in 0:d) {
      s <- s + factorial(d + i) / (factorial(i) * factorial(d - i)) *
        (sqrt(8 * v) * r / l) ^ (d - i)
    }
    k_v <- exp(-sqrt(2 * v) * r / l) * gamma(d + 1) / gamma(2 * d + 1) * s
    return(k_v)
  }
  matern_kern <- function(X2, X1) 
    apply(X1, 1, function(xp){ 
      apply(X2, 1, function(xq){ 
        point_wise(xp, xq, l, d)
      })
    })
  Z1 <- rmvnorm(n = 50, mean = rep(0, 3), sigma = diag(3))
  Z2 <- rmvnorm(n = 52, mean = rep(0, 3), sigma = diag(3))
  Ktest <- round(matern_kern(Z1, Z2), 5)
  matern_fun <- generate_kernel("matern", l = 1, d = 2)
  Kmat <- round(matern_fun(Z1, Z2), 5)
})


test_that("rational kernel", {
  point_wise <- function(xp, xq, l, d) {
    l <- 1
    d <- 2
    r <- sqrt(sum((xp - xq) ^ 2))
    k_a <- (1 + r ^ 2 / (2 * d * l ^ 2)) ^ (- d)
    return(k_a)
  }
  rational_kern <- function(X2, X1) 
    apply(X1, 1, function(xp){ 
      apply(X2, 1, function(xq){ 
        point_wise(xp, xq, l, d)
      })
    })
  Z1 <- rmvnorm(n = 50, mean = rep(0, 3), sigma = diag(3))
  Z2 <- rmvnorm(n = 52, mean = rep(0, 3), sigma = diag(3))
  Ktest <- round(rational_kern(Z1, Z2), 5)
  rational_fun <- generate_kernel("rational", l = 1, d = 2)
  Kmat <- round(rational_fun(Z1, Z2), 5)
})


test_that("nn kernel", {
  point_wise <- function(xp, xq, l, d) {
    xp <- c(1, xp)
    xq <- c(1, xq)
    s <- 2 * t(xp)  %*% xq / (sqrt((1 + 2 * t(xp) %*% xp) 
                                   * (1 + 2 * t(xq) %*% xq)))
    k_n <- asin(s)
    return(k_n)
  }
  nn_kern <- function(X2, X1) 
    apply(X1, 1, function(xp){ 
      apply(X2, 1, function(xq){ 
        point_wise(xp, xq, l, d)
      })
    })
  Z1 <- rmvnorm(n = 50, mean = rep(0, 3), sigma = diag(3))
  Z2 <- rmvnorm(n = 52, mean = rep(0, 3), sigma = diag(3))
  Ktest <- round(nn_kern(Z1, Z2), 5)
  nn_fun <- generate_kernel("nn")
  Kmat <- round(nn_fun(Z1, Z2), 5)
})


# test_that("kernel implementation", {
#   formula <- Y ~ X + Z1 + Z2
#   n <- 100
#   label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4"))
#   kern_par <- data.frame(method = c("rbf", "polynomial", "matern"), 
#                          l = c(.5, 1, 1.5), d = 1:3)
#   kern_par$method <- as.character(kern_par$method)
#   data <- generate_data(n, fixed_num = 1, label_names, method = "rbf", 
#                         int_effect = .3, l = 1, d = 2, eps = .01)
#   fit <- define_model(formula, data, kern_par, fixed_num = 1, label_names)
#   K1 <- fit$kern_list[[1]]
#   Kmat1 <- round(K1(Z1, Z1), 5)
#   Ktest1 <- rbfdot(sigma = 2)
#   Kmat_test1 <- round(kernelMatrix(Ktest1, Z1)@.Data, 5)
#   expect_equal(Kmat1, Kmat_test1)
#   
#   K2 <- fit$kern_list[[2]]
#   Kmat2 <- round(K2(Z1, Z1), 5)
#   Ktest2 <- polydot(degree = 2, scale = 1, offset = 1)
#   Kmat_test2 <- round(kernelMatrix(Ktest2, Z1)@.Data, 5)
#   expect_equal(Kmat2, Kmat_test2)
# })

