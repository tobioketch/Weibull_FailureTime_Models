################################################################################
# This R script is executed as a function to solve
# values of interest from an RealDt (Real) Data set
################################################################################

# Loading Packages &Methods -----------------------------------------------

source(file.path("2--run_codes", "NOT_RUN_RCODES", "LoadMethods.R"))

# Parameters &Args - Data -------------------------------------------------

# Provide Data (Name of Data) Here As an Argument
# Args for RealDt (Application Data) Dt
ArgsApp <- this.path::progArgs()

# Subtitle: Name of the Original Data to Identify Precision Plot
CharArgs <- tail(ArgsApp, n=2)
SubTitle <- CharArgs[1]
Top_n <-  as.numeric(CharArgs[2])
n_nums <- NROW(ArgsApp)- NROW(CharArgs) # numeric Dt values
ArgsDt <- as.numeric(ArgsApp[1:n_nums])
DtNames <- gsub(" ", "", SubTitle)

LifeTime <- setNames(list(ArgsDt),  DtNames)

# MCMC Approach -----------------------------------------------------------

# Setting up options and multiple workers for the MCMC
rstan_options(auto_write = TRUE)
workers <- as.integer(parallel::detectCores(logical = FALSE)) - 2L
options(mc.cores = workers)

# ******************************************************************************************
# ModelMCMCRealDt <- parstan(x = LifeTime, m = 0, v = 25., fctor = 1., comparedefault = "exponential", tidyselectfun  = "starts_with")

ModelMCMCRealDt <- parstan(x = LifeTime, m = 0, v = 25., fctor = 1.)

# summarise parameters
MCMCStatsRealDt <- summarise_stan_pars(ModelMCMCRealDt) %>%
  # separates data name & methods into two columns
  separate_data_method_cols(drop = TRUE)

# %>%
#   efficient_model() # add pseudo pars and filter efficient models

# Classical Approach ------------------------------------------------------
# Fit Classical Models
classic.methods <- eval(formals(boot_weib_summaries)$method)
classic.methods ## 1] "Regression" "Moments"    "MLE"

# Fit a Classic Methods
ModelClassicRealDt <- purrr::map(classic.methods, \(m) {
  purrr::map(LifeTime, ~ boot_weib_summaries(.x, method = m, seed = 1333455234L) %>%
               magrittr::extract2(1) %>% dplyr::select(!contains(c("f.root", "iter", "estim.prec")) )) #names not in regression method
}) %>%
  stats::setNames(nm = classic.methods) %>%
  # Data in list format to facilitate computing for fisher Information and the associated Var-Covar matrix
  base::unlist(recursive = FALSE) %>%
  dplyr::bind_rows(.id = "method") %>%
  purrr::modify_if(is.numeric, round, digits = 4L) %>%
  tidyr::separate(col = "method", into = c("method", NA), sep = "[.]")
  # tidyr::separate_wider_delim(cols = "method", delim = ".", names = c("method", NA) )

# ModelClassicRealDt <- add_pseudo_pars(ModelClassicRealDt)

# COMBINE RESULTS ---------------------------------------------------------

RealDtPoolMetrics <- dplyr::bind_rows(MCMCStatsRealDt, ModelClassicRealDt) %>%
  method_to_priorcomb_cols()


save(RealDtPoolMetrics, file = "misc/app_misc_results/RealDtPoolMetrics" )


