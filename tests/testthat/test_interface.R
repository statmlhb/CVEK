context("Model estimation interface")

# set up data
n <- 100
d <- 7
data <- as.data.frame(matrix(
  rnorm(n * d),
  ncol = d,
  dimnames = list(NULL, paste0("x", 1:d))
))
data$y <- as.matrix(data) %*% rnorm(d)

data_new <- as.data.frame(matrix(
  rexp(n * d),
  ncol = d,
  dimnames = list(NULL, paste0("x", 1:d))
))

test_that(desc = "parse_kernel_variable computes kernel effects",
          code = {
            rbf_kern_func <- generate_kernel(method = "rbf", l = 1)
            lnr_kern_func <- generate_kernel(method = "linear")
            
            kern_effect_mat_singl <-
              parse_kernel_variable("k(x1)", rbf_kern_func, data)
            kern_effect_mat_singl_expected <-
              compute_expected_matrix("x1", rbf_kern_func, data)
            
            kern_effect_mat_multi <-
              parse_kernel_variable("k(x5,x7,x3)", rbf_kern_func, data)
            kern_effect_mat_multi_expected <-
              compute_expected_matrix(c("x5", "x7", "x3"), rbf_kern_func, data)
            
            fixd_effect_mat <-
              parse_kernel_variable("x4", rbf_kern_func, data)
            fixd_effect_mat_expected <-
              compute_expected_matrix("x4", lnr_kern_func, data)
            
            expect_identical(kern_effect_mat_singl, kern_effect_mat_singl_expected)
            expect_identical(kern_effect_mat_multi, kern_effect_mat_multi_expected)
            expect_identical(fixd_effect_mat, fixd_effect_mat_expected)
          })


test_that(desc = "parse_kernel_variable computes predictive kernel",
          code = {
            rbf_kern_func <- generate_kernel(method = "rbf", l = 1)
            lnr_kern_func <- generate_kernel(method = "linear")
            
            kern_effect_mat_singl <-
              parse_kernel_variable("k(x1)", rbf_kern_func, data, data_new)
            kern_effect_mat_singl_expected <-
              compute_expected_matrix("x1", rbf_kern_func, data, data_new)
            
            kern_effect_mat_multi <-
              parse_kernel_variable("k(x5,x7,x3)", rbf_kern_func, data, data_new)
            kern_effect_mat_multi_expected <-
              compute_expected_matrix(c("x5", "x7", "x3"), rbf_kern_func, data, data_new)
            
            fixd_effect_mat <-
              parse_kernel_variable("x4", rbf_kern_func, data, data_new)
            fixd_effect_mat_expected <-
              compute_expected_matrix("x4", lnr_kern_func, data, data_new)
            
            expect_identical(kern_effect_mat_singl, kern_effect_mat_singl_expected)
            expect_identical(kern_effect_mat_multi, kern_effect_mat_multi_expected)
            expect_identical(fixd_effect_mat, fixd_effect_mat_expected)
          })

test_that(desc = "parse_kernel_terms correctness of interaction term",
          code = {
            rbf_kern_func <- generate_kernel(method = "rbf", l = 1.25)
            lnr_kern_func <- generate_kernel(method = "linear")
            
            formula <- y ~ k(x1):k(x4,x6):x7
            
            kern_term_mat <- 
              parse_kernel_terms(formula, rbf_kern_func, data)
            kern_term_mat_expected <-
              compute_expected_matrix("x1", rbf_kern_func, data) * 
              compute_expected_matrix(c("x4", "x6"), rbf_kern_func, data) * 
              compute_expected_matrix("x7", lnr_kern_func, data)
            
            # should return exact one term
            expect_equal(length(kern_term_mat), 1)
            # evaluate correctness
            expect_identical(kern_term_mat[[1]], kern_term_mat_expected)
          })

test_that(desc = "parse_cvek_formula regular test",
          code = {
            formula <- y ~ x1 + k(x2) + k(x3, x4) + k(x5):x7
            kern_func_list <- setup_kernel_lib()
            
            model_mat_list <- parse_cvek_formula(formula, kern_func_list, data)
            
            # compute expected content 
            kern_check_id <- 3
            lib_kern_func <- kern_func_list[[kern_check_id]]
            lnr_kern_func <- generate_kernel(method = "linear")
            
            kern_term1_expected <- 
              compute_expected_matrix("x2", lib_kern_func, data)
            kern_term2_expected <- 
              compute_expected_matrix(c("x3", "x4"), lib_kern_func, data)
            kern_term3_expected <- 
              compute_expected_matrix("x5", lib_kern_func, data) * 
              compute_expected_matrix("x7", lnr_kern_func, data)
            
            # list structure
            expect_equal(length(model_mat_list), 3)
            expect_equal(names(model_mat_list), c("y", "X", "K"))
            # check list content: outcome and fixed-effect matrix
            expect_equal(model_mat_list$y, data$y)
            expect_equivalent(model_mat_list$X, cbind(1, data$x1))
            # check list content: kernel matrix
            expect_equal(length(model_mat_list$K), 3)
            expect_equivalent(model_mat_list$K[[kern_check_id]], 
                              list(kern_term1_expected,
                                   kern_term2_expected,
                                   kern_term3_expected))
            
          })


test_that(desc = "parse_cvek_formula with no intercept term",
          code = {
            formula <- y ~ -1 + x1 + k(x2) + k(x3, x4) + k(x5):x7
            kern_func_list <- setup_kernel_lib()
            
            model_mat_list <- parse_cvek_formula(formula, kern_func_list, data)
            expect_equal(colnames(model_mat_list$X), "x1")
          })

 test_that(desc = "parse_cvek_formula with only fixed effects",
           code = {
             formula <- y ~ x1 + x2 + x5*x7
             kern_func_list <- setup_kernel_lib()
             
             model_mat_list <- parse_cvek_formula(formula, kern_func_list, data)
             expect_equal(model_mat_list$K, vector("list", length=3))
             expect_equal(model_mat_list$X, 
                          model.matrix(formula, data = data))
             
           })

 test_that(desc = "parse_cvek_formula with only kernel effects",
           code = {
             formula <- y ~ k(x1) + k(x2) + k(x5, x7):x1
             kern_func_list <- setup_kernel_lib()
             
             model_mat_list <- parse_cvek_formula(formula, kern_func_list, data)
             expect_equal(sapply(model_mat_list$K, length), 
                          c(3, 3, 3))
             # should contain only intercept
             expect_equal(model_mat_list$X, 
                          model.matrix(~1, data = data))
             
           })
 
 test_that(desc = "parse_cvek_formula predict fixed effect",
           code = {
             formula <- y ~ x1 + x2 + x5*x7
             formula_new <- ~ x1 + x2 + x5*x7
             
             kern_func_list <- setup_kernel_lib()
             
             model_mat_list <- parse_cvek_formula(formula, kern_func_list, 
                                                  data, data_new)
             expect_equal(model_mat_list$X, 
                          model.matrix(formula_new, data = data_new))
             expect_equal(model_mat_list$y, NULL)
           })
 