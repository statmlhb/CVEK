context("Model estimation")

test_that(desc = "single kernel term case",
          code = {
            # set up data
            kern_par <- data.frame(method = c("linear", "polynomial", "matern"), 
                                   l = c(.5, 1, 1.5), d = 1:3, stringsAsFactors = FALSE)
            # define kernel library
            kern_func_list <- list()
            for (j in 1:nrow(kern_par)) {
              kern_func_list[[j]] <- generate_kernel(kern_par[j,]$method, 
                                                     kern_par[j,]$l, 
                                                     kern_par[j,]$d)
            }
            n <- 100
            d <- 4
            formula <- y ~ x1 + x2 + k(x3, x4)
            data <- as.data.frame(matrix(
              rnorm(n * d),
              ncol = d,
              dimnames = list(NULL, paste0("x", 1:d))
            ))
            lnr_kern_func <- generate_kernel(method = "linear")
            kern_effect_mat <- 
              parse_kernel_variable("k(x3, x4)", lnr_kern_func, data)
            beta_true <- c(1, .41, 2.37)
            alpha_true <- rnorm(n)
            
            K <- kern_effect_mat
            data$y <- as.matrix(cbind(1, data[, c("x1", "x2")])) %*% beta_true + 
              K %*% alpha_true
            
            model_matrices <- parse_cvek_formula(formula, 
                                                 kern_func_list = kern_func_list, 
                                                 data = data, verbose = FALSE)
            
            result <- estimation(Y = model_matrices$y, 
                                 X = model_matrices$X, 
                                 K_list = model_matrices$K, 
                                 mode = "loocv", strategy = "stack", 
                                 beta_exp = 1, lambda = exp(seq(-10, 5)))
            
            # X shouldn't be standardized
            X <- model_matrices$X   
            y_new_ridge <- X %*% result$beta + result$K %*% result$alpha
            
            K <- kern_effect_mat
            K_scale <- tr(K)
            K <- K / K_scale
            lambda0 <- tuning(data$y, X, list(K), mode = "loocv", lambda = exp(seq(-10, 5)))
            V_inv <- ginv(K + lambda0 * diag(n))
            B_mat <- ginv(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv
            P_X <- X %*% B_mat
            beta_est <- B_mat %*% data$y
            alpha_est <- V_inv %*% (diag(n) - P_X) %*% data$y
            y_est <- X %*% beta_est + K %*% alpha_est
            diff1 <- euc_dist(data$y, y_est)^2 / length(y_est)
            diff2 <- euc_dist(y_new_ridge, y_est)^2 / length(y_est)
            expect_lte(diff1, .001)
            expect_lte(diff2, .001)
          })

test_that(desc = "two kernel terms case",
          code = {
            # set up data
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
            kern_effect_lnr1 <- 
              parse_kernel_variable("k(x3, x4)", lnr_kern_func, data)
            kern_effect_lnr2 <- 
              parse_kernel_variable("k(x5, x6)", lnr_kern_func, data)
            set.seed(1231)
            beta_true <- c(1, .41, 2.37)
            alpha_lnr1_true <- rnorm(n)
            alpha_lnr2_true <- rnorm(n)
            
            kern_term_lnr1 <- kern_effect_lnr1 %*% alpha_lnr1_true
            kern_term_lnr2 <- kern_effect_lnr2 %*% alpha_lnr2_true
            
            data$y <- as.matrix(cbind(1, data[, c("x1", "x2")])) %*% beta_true + 
              kern_term_lnr1 + kern_term_lnr2
            
            model_matrices <- parse_cvek_formula(formula, 
                                                 kern_func_list = kern_func_list, 
                                                 data = data, verbose = FALSE)
            
            result <- estimation(
              Y = model_matrices$y, 
              X = model_matrices$X, 
              K_list = model_matrices$K, 
              mode = "loocv", strategy = "stack", 
              beta_exp = 1, lambda = exp(seq(-10, 5)))
            
            X <- model_matrices$X   
            y_est <- X %*% result$beta + result$K %*% result$alpha
            y_est2 <- X %*% result$beta + result$kern_term_effect[, 1] + 
              result$kern_term_effect[, 2]
            
            diff1 <- euc_dist(data$y, y_est)^2 / length(y_est)
            diff2 <- euc_dist(y_est2, y_est)^2 / length(y_est)
            diff3 <- euc_dist(kern_term_lnr1, result$kern_term_effect[, 1])^2 / 
              length(kern_term_lnr1)
            diff4 <- euc_dist(kern_term_lnr2, result$kern_term_effect[, 2]) / 
              length(kern_effect_lnr2)
            expect_lte(diff1, .001)
            expect_lte(diff2, .001)
            expect_lte(diff3, .001)
            expect_lte(diff4, .001)
          })

test_that(desc = "individual term estimate consistent with projection matrix",
          code = {
            # set up data
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
            kern_effect_lnr1 <- 
              parse_kernel_variable("k(x3, x4)", lnr_kern_func, data)
            kern_effect_lnr2 <- 
              parse_kernel_variable("k(x5, x6)", lnr_kern_func, data)
            set.seed(1231)
            beta_true <- c(1, .41, 2.37)
            alpha_lnr1_true <- rnorm(n)
            alpha_lnr2_true <- rnorm(n)
            
            kern_term_lnr1 <- kern_effect_lnr1 %*% alpha_lnr1_true
            kern_term_lnr2 <- kern_effect_lnr2 %*% alpha_lnr2_true
            
            data$y <- as.matrix(cbind(1, data[, c("x1", "x2")])) %*% beta_true + 
              kern_term_lnr1 + kern_term_lnr2
            
            model_matrices <- parse_cvek_formula(formula, 
                                                 kern_func_list = kern_func_list, 
                                                 data = data, verbose = FALSE)
            
            result_byterm <- estimation(Y = model_matrices$y, 
                                        X = model_matrices$X, 
                                        K_list = model_matrices$K, 
                                        mode = "loocv", strategy = "stack", 
                                        beta_exp = 1, lambda = exp(seq(-10, 5)),
                                        compute_kernel_terms = TRUE)
            
            result_overall <- estimation(Y = model_matrices$y, 
                                         X = model_matrices$X, 
                                         K_list = model_matrices$K, 
                                         mode = "loocv", strategy = "stack", 
                                         beta_exp = 1, lambda = exp(seq(-10, 5)),
                                         compute_kernel_terms = FALSE)
            
            
            y_est_overall <- result_overall$kern_term_effect
            y_est_byterm <- result_byterm$kern_term_effect[, 1] + 
              result_byterm$kern_term_effect[, 2]
            
            diff <- euc_dist(y_est_overall, y_est_byterm)^2 / length(y_est)
            expect_lte(diff, .001)
          })

test_that("implementation of pure fixed effects", {
  names(longley) <- c("y", paste0("x", 1:6))
  formula <- y ~ x1 + x2 + x3 + x4 + x5 + x6
  model_matrices <- parse_cvek_formula(formula, 
                                       kern_func_list = NULL, 
                                       data = longley, verbose = FALSE)
  result <- estimation(Y = model_matrices$y, 
                       X = model_matrices$X, 
                       K_list = model_matrices$K, 
                       mode = "GCV", strategy = "stack", 
                       beta_exp = 1, lambda = exp(seq(-10, 5)))
  
  X <- model_matrices$X   
  y_est <- X %*% result$beta
  
  mod0 <- lm.ridge(y ~ ., longley, lambda = exp(seq(-10, 5)))
  whichIsBest <- which.min(mod0$GCV)
  y_mod0 <- as.matrix(cbind(const = 1, longley[, -1])) %*% coef(mod0)[whichIsBest,]
  diff <- euc_dist(y_mod0, y_est)^2 / length(y_est)
  expect_lte(diff, .5)
})


test_that("implementation of pure random effects", {
  names(longley) <- c("y", paste0("x", 1:6))
  kern_par <- data.frame(method = "linear", l = 1, d = 2,
                         stringsAsFactors = FALSE)
  kern_func_list <- list()
  for (j in 1:nrow(kern_par)) {
    kern_func_list[[j]] <- generate_kernel(kern_par[j,]$method,
                                           kern_par[j,]$l,
                                           kern_par[j,]$d)
  }
  
  formula <- y ~ k(x1, x2, x3, x4, x5, x6)
  model_matrices <- parse_cvek_formula(formula, 
                                       kern_func_list = kern_func_list, 
                                       data = longley, verbose = FALSE)
  result <- estimation(Y = model_matrices$y, 
                       X = model_matrices$X, 
                       K_list = model_matrices$K, 
                       mode = "GCV", strategy = "stack", 
                       beta_exp = 1, lambda = exp(seq(-10, 5)))
  
  X <- model_matrices$X   
  y_est <- X %*% result$beta + result$K %*% result$alpha
  
  mod0 <- lm.ridge(y ~ ., longley, lambda = exp(seq(-10, 5)))
  whichIsBest <- which.min(mod0$GCV)
  y_mod0 <- as.matrix(cbind(const = 1, longley[, -1])) %*% coef(mod0)[whichIsBest,]
  diff <- euc_dist(y_mod0, y_est)^2 / length(y_est)
  expect_lte(diff, .5)
})
