context("Hypothesis testing")

test_that("range of output", {
  kern_par <- data.frame(method = c("linear", "polynomial", "rbf"), 
                         l = c(.5, 1, 1.5), d = 1:3, stringsAsFactors = FALSE)
  # define kernel library
  kern_func_list <- list()
  for (j in 1:nrow(kern_par)) {
    kern_func_list[[j]] <- generate_kernel(kern_par[j,]$method, 
                                           kern_par[j,]$l, 
                                           kern_par[j,]$d)
  }
  n <- 100
  d <- 6
  formula <- y ~ x1 + x2 + k(x3, x4) + k(x5, x6)
  data <- as.data.frame(matrix(
    rnorm(n * d),
    ncol = d,
    dimnames = list(NULL, paste0("x", 1:d))
  ))
  lnr_kern_func <- generate_kernel(method = "linear")
  rbf_kern_func <- generate_kernel(method = "rbf", l = 1.25)
  kern_effect_lnr <- 
    parse_kernel_variable("k(x3, x4)", lnr_kern_func, data)
  kern_effect_rbf <- 
    parse_kernel_variable("k(x5, x6)", lnr_kern_func, data)
  beta_true <- c(1, .41, 2.37)
  alpha_lnr_true <- rnorm(n)
  alpha_rbf_true <- rnorm(n)
  
  kern_term_lnr <- kern_effect_lnr %*% alpha_lnr_true
  kern_term_rbf <- kern_effect_rbf %*% alpha_rbf_true
  
  data$y <- as.matrix(cbind(1, data[, c("x1", "x2")])) %*% beta_true + 
    kern_term_lnr + kern_term_rbf
  
  model_matrices <- parse_cvek_formula(formula, 
                                       kern_func_list = kern_func_list, 
                                       data = data, verbose = FALSE)
  K_int <- parse_kernel_terms(y ~ k(x1):k(x4, x6):x5, rbf_kern_func, data)
  pvalue <- testing(Y = model_matrices$y, X = model_matrices$X, 
                    K_list = model_matrices$K, K_int = K_int[[1]], 
                    mode = "loocv", strategy = "stack", 
                    beta_exp = 1, test = "boot", 
                    lambda = exp(seq(-10, 5)), B = 100)$pvalue

  expect_lte(pvalue, 1)
  expect_gte(pvalue, 0)
})

test_that("warning message from tuning", {
  kern_par <- data.frame(method = c("linear", "polynomial", "rbf"), 
                         l = c(.5, 1, 1.5), d = 1:3, stringsAsFactors = FALSE)
  # define kernel library
  kern_func_list <- list()
  for (j in 1:nrow(kern_par)) {
    kern_func_list[[j]] <- generate_kernel(kern_par[j,]$method, 
                                           kern_par[j,]$l, 
                                           kern_par[j,]$d)
  }
  n <- 100
  d <- 6
  formula <- y ~ x1 + x2 + k(x3, x4) + k(x5, x6)
  data <- as.data.frame(matrix(
    rnorm(n * d),
    ncol = d,
    dimnames = list(NULL, paste0("x", 1:d))
  ))
  lnr_kern_func <- generate_kernel(method = "linear")
  rbf_kern_func <- generate_kernel(method = "rbf", l = 1.25)
  kern_effect_lnr <- 
    parse_kernel_variable("k(x3, x4)", lnr_kern_func, data)
  kern_effect_rbf <- 
    parse_kernel_variable("k(x5, x6)", lnr_kern_func, data)
  beta_true <- c(1, .41, 2.37)
  alpha_lnr_true <- rnorm(n)
  alpha_rbf_true <- rnorm(n)
  
  kern_term_lnr <- kern_effect_lnr %*% alpha_lnr_true
  kern_term_rbf <- kern_effect_rbf %*% alpha_rbf_true
  
  data$y <- as.matrix(cbind(1, data[, c("x1", "x2")])) %*% beta_true + 
    kern_term_lnr + kern_term_rbf
  
  model_matrices <- parse_cvek_formula(formula, 
                                       kern_func_list = kern_func_list, 
                                       data = data, verbose = FALSE)
  K_int <- parse_kernel_terms(y ~ k(x1):k(x4, x6):x5, rbf_kern_func, data)
  lambda = rep(.5, 11)
  expect_warning(testing(Y = model_matrices$y, X = model_matrices$X, 
                         K_list = model_matrices$K, K_int = K_int[[1]], 
                         mode = "loocv", strategy = "stack", 
                         beta_exp = 1, test = "boot", 
                         lambda = lambda, B = 100),
                 "the smallest one")
})
