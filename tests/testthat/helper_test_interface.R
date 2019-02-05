setup_kernel_lib <- function(){
  kern_par <- data.frame(
    method = c("rbf", "polynomial", "matern"),
    l = c(.5, 1, 1.5),
    p = 1:3,
    stringsAsFactors = FALSE
  )
  
  # set up kernel library
  kern_func_list <- list()
  for (d in 1:nrow(kern_par)) {
    kern_func_list[[d]] <- generate_kernel(kern_par[d, ]$method,
                                           kern_par[d, ]$l,
                                           kern_par[d, ]$p)
  }
  kern_func_list
}

# helper function for computing kernel matrix
compute_expected_matrix <- function(term_names, kern_func, data) {
  data_mat <- as.matrix(data[, term_names, drop = FALSE])
  kern_func(data_mat, data_mat)
}
