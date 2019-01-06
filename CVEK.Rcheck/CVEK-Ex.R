pkgname <- "CVEK"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "CVEK-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('CVEK')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("define_model")
### * define_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: define_model
### Title: Defining the Model
### Aliases: define_model

### ** Examples




kern_par <- data.frame(method = c("rbf", "polynomial", "matern"), 
l = c(.5, 1, 1.5), d = 1:3)
kern_par$method <- as.character(kern_par$method)
define_model(formula = Y ~ X + Z1 + Z2, data = mydata, kern_par, 
fixed_num = 1, label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4")))






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("define_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ensemble")
### * ensemble

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ensemble
### Title: Estimating Ensemble Kernel Matrices
### Aliases: ensemble

### ** Examples




ensemble(n = 100, kern_size = 3, strategy = "stack", beta_exp = 1, 
CVEK:::error_mat, CVEK:::A_hat)






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ensemble", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("estimate_base")
### * estimate_base

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: estimate_base
### Title: Estimating Projection Matrices
### Aliases: estimate_base
### Keywords: internal

### ** Examples




estimate_base(n = 100, kern_size = 3, CVEK:::Y, CVEK:::X, CVEK:::Z1, CVEK:::Z2, 
CVEK:::kern_list, mode = "loocv", lambda = exp(seq(-10, 5)))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("estimate_base", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("estimate_ridge")
### * estimate_ridge

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: estimate_ridge
### Title: Estimating a Single Model
### Aliases: estimate_ridge

### ** Examples




estimate_ridge(X = cbind(matrix(1, nrow = 100, ncol = 1), CVEK:::X), 
CVEK:::K_ens, CVEK:::Y, lambda = exp(seq(-10, 5)))






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("estimate_ridge", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("estimation")
### * estimation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: estimation
### Title: Conducting Gaussian Process Regression
### Aliases: estimation

### ** Examples




estimation(CVEK:::Y, CVEK:::X, CVEK:::Z1, CVEK:::Z2, CVEK:::kern_list, 
mode = "loocv", strategy = "stack", beta_exp = 1, lambda = exp(seq(-10, 5)))






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("estimation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generate_data")
### * generate_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generate_data
### Title: Generating Original Data
### Aliases: generate_data

### ** Examples




mydata <- generate_data(n = 100, fixed_num = 1, label_names =
list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4")),
method = "rbf", l = 1, d = 2, int_effect = 0, eps = .01)






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generate_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generate_formula")
### * generate_formula

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generate_formula
### Title: From Vectors to Single Variables
### Aliases: generate_formula

### ** Examples




generic_formula0 <- generate_formula(formula = Y ~ X + Z1 + Z2,
label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4")))






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generate_formula", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generate_kernel")
### * generate_kernel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generate_kernel
### Title: Generating A Single Kernel
### Aliases: generate_kernel

### ** Examples



kern_list <- list()
for (k in 1:nrow(kern_par)) {
  kern_list[[k]] <- generate_kernel(kern_par[k, ]$method,
                                    kern_par[k, ]$l,
                                    kern_par[k, ]$d)
}





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generate_kernel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("testing")
### * testing

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: testing
### Title: Conducting Score Tests for Interaction
### Aliases: testing

### ** Examples




testing(formula_int = Y ~ X + Z1 * Z2,
label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4")),
CVEK:::Y, CVEK:::X, CVEK:::Z1, CVEK:::Z2, CVEK:::kern_list, 
mode = "loocv", strategy = "stack",
beta_exp = 1, test = "boot", lambda = exp(seq(-10, 5)), B = 100)






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("testing", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tuning")
### * tuning

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tuning
### Title: Calculating Tuning Parameters
### Aliases: tuning

### ** Examples




lambda0 <- tuning(CVEK:::Y, CVEK:::X, K_mat = CVEK:::K_ens, 
mode = "loocv", lambda = exp(seq(-10, 5)))






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tuning", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
