context("Model prediction")
n <- 50
d <- 4
formula <- y ~ x1 + x2 + k(x3, x4)
formula_test <- y ~ k(x1, x2) * k(x3, x4)
data <- as.data.frame(matrix(
  rnorm(n * d),
  ncol = d,
  dimnames = list(NULL, paste0("x", 1:d))
))
beta_true <- c(1, .41, 2.37)
lnr_kern_func <- generate_kernel(method = "linear")
kern_effect_lnr <- 
  parse_kernel_variable("k(x3, x4)", lnr_kern_func, data)
alpha_lnr_true <- rnorm(n)

data$y <- as.matrix(cbind(1, data[, c("x1", "x2")])) %*% beta_true + 
  kern_effect_lnr %*% alpha_lnr_true

test_that(desc = "single true kernel",
          code = {
            # set up data
            kern_par <- data.frame(method = c("linear"), 
                                   l = c(.5), d = 1, stringsAsFactors = FALSE)
            # define kernel library
            kern_func_list <- list()
            for (j in 1:nrow(kern_par)) {
              kern_func_list[[j]] <- generate_kernel(kern_par[j,]$method, 
                                                     kern_par[j,]$l, 
                                                     kern_par[j,]$d)
            }
            
            result <- cvek(formula,
                           kern_func_list,
                           data,
                           formula_test,
                           mode = "loocv",
                           strategy = "stack",
                           beta_exp = 1,
                           lambda = exp(seq(-10, 5)),
                           test = "boot",
                           alt_kernel_type = "linear",
                           B = 100,
                           verbose = FALSE)
            
            kern_effect_st <- kern_effect_lnr / tr(kern_effect_lnr)
            X <- result$model_matrices$X
            y_est <- X %*% result$beta + result$K %*% result$alpha
            y_est2 <- kern_effect_st %*% 
              ginv(kern_effect_st + result$lambda * diag(n)) %*% 
              (result$data$y - X %*% beta_true) + X %*% beta_true
            y_est3 <- predict(result, data)
            
            diff <- euc_dist(y_est3, y_est2)^2 / length(y_est2)
            expect_lte(diff, .1)
          })


test_that("multiple base linear kernels", 
          code = {
            # set up data
            kern_par <- data.frame(method = c("linear", "linear", "linear"), 
                                   l = c(.5, 1, 2), d = 1:3, stringsAsFactors = FALSE)
            # define kernel library
            kern_func_list <- list()
            for (j in 1:nrow(kern_par)) {
              kern_func_list[[j]] <- generate_kernel(kern_par[j,]$method, 
                                                     kern_par[j,]$l, 
                                                     kern_par[j,]$d)
            }
            
            result <- cvek(formula,
                           kern_func_list,
                           data,
                           formula_test,
                           mode = "GCV",
                           strategy = "stack",
                           beta_exp = 1,
                           lambda = exp(seq(-10, 5)),
                           test = "boot",
                           alt_kernel_type = "linear",
                           B = 100,
                           verbose = FALSE)
            
            X <- result$model_matrices$X
            y_est <- X %*% result$beta + result$K %*% result$alpha
            y_est3 <- predict(result, data)
            
            mod0 <- lm.ridge(y ~ ., data, lambda = exp(seq(-10, 5)))
            whichIsBest <- which.min(mod0$GCV)
            y_mod0 <- as.matrix(cbind(const = 1, data[, 1:4])) %*% coef(mod0)[whichIsBest,]
            diff <- euc_dist(y_mod0, y_est3)^2 / length(y_est3)
            expect_lte(diff, .1)
          })


test_that("multiple true kernels, train and test sets",
          code = {
            # set up dat
            kern_par <- data.frame(method = rep("rbf", 3),
                                   l = rep(3, 3), d = rep(2, 3), 
                                   stringsAsFactors = FALSE)
            # define kernel library
            kern_func_list <- list()
            for (j in 1:nrow(kern_par)) {
              kern_func_list[[j]] <- generate_kernel(kern_par[j,]$method,
                                                     kern_par[j,]$l,
                                                     kern_par[j,]$d)
            }

            n <- 100
            data <- as.data.frame(matrix(
              rnorm(n * d),
              ncol = d,
              dimnames = list(NULL, paste0("x", 1:d))
            ))
            beta_true <- c(1, .41, 2.37)
            lnr_kern_func <- generate_kernel(method = "rbf", l = 3)
            kern_effect_lnr <-
              parse_kernel_variable("k(x3, x4)", lnr_kern_func, data)
            alpha_lnr_true <- rnorm(n)

            data$y <- as.matrix(cbind(1, data[, c("x1", "x2")])) %*% beta_true +
              kern_effect_lnr %*% alpha_lnr_true

            data_train <- data[1:50, ]
            data_test <- data[51:100, ]

            result <- cvek(formula,
                           kern_func_list,
                           data_train,
                           formula_test,
                           mode = "loocv",
                           strategy = "stack",
                           beta_exp = 1,
                           lambda = exp(seq(-10, 5)),
                           test = "boot",
                           alt_kernel_type = "linear",
                           B = 100,
                           verbose = FALSE)

            y_est3 <- predict(result, data_test)
            SSE <- euc_dist(data_test[, "y"], y_est3)^2
            R_sq <- 1 - SSE / (var(data_test[, "y"]) * (length(y_est3) - 1))
            
            plot(data_test[, "y"], y_est3)
            y45 <- function(x) x
            curve(y45, -4, 11, add = TRUE)
            summary(lm(y_est3 ~ data_test[, "y"]))
          })

