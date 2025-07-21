# Loading Shape_Scale Prior Combinations
SHAPE_PRIORS <- c("LogNormal", "Exponential", "Gamma", "HalfNormal" )
SCALE_PRIORS <- c("Gamma", "HalfNormal", "LogNormal", "InverseGamma", "HalfCauchy", "Exponential")

# directory for the stan sampling models/ distributions for prior combinations
PATH_BASE_COMBS = file.path("2--run_codes/NOT_RUN_RCODES/StanSamplingModels")

.concat_priors <- function(shapeprior = SHAPE_PRIORS, scaleprior=SCALE_PRIORS) {
  # returns default prior combinations for We(param |a,b)
  df <- tidyr::expand_grid(shapeprior, scaleprior)
  # shape-scale priors combinations
  Out <- with(df, paste(df$shapeprior, df$scaleprior, sep = "_"))
  return(Out)
}


# base_prior - the name of the shape prior ( shape prior used as the main/ base name of the prior combinations)
# priorcombs - a list of prior_combinations 

.defaultShapeBasedSamplingModl <- function(baseprior, priorcombs = .concat_priors(), basepath = PATH_BASE_COMBS) {
  # A function that takes a name of a base prior, and reads the default stan scripts of the prior
  # combinations associated with the base prior

  # sampling models organized into folders
  # use base prior to identify the folder containing the prior combinations of the base prior
  SchmDir <- tidyselect::vars_select(dir(path = basepath), tidyselect::contains(baseprior))
  # character vector of priors combinations
  BasePrCombs <- tidyselect::vars_select(priorcombs, tidyselect::starts_with(baseprior))
  n <- length(BasePrCombs)
  combs <- purrr::map(
    seq_len(n),
    ~ vroom::vroom_lines(here::here(
      basepath, SchmDir, paste0(BasePrCombs[.x], ".stan")
    ))
  )
  names(combs) <- names(BasePrCombs)
  combs
}

.defaultSamplingModls <- function(shapeprior = SHAPE_PRIORS, prior_combs = .concat_priors(), base_path = PATH_BASE_COMBS) {
  # returns default sampling schemes
  mdl <- NULL
  for (prior in shapeprior) {
    mdl_prior <- .defaultShapeBasedSamplingModl(baseprior = prior, priorcombs = prior_combs, basepath = base_path)
    mdl <- c(mdl, mdl_prior)
  }
  mdl
}


# Load Rstan-Parameters
# Model specifications (list of specified sampler parameters)

.defaultTunePars <- function() {
  PrComb <- within(list(), {
    # Exp Shape (Base) Prior
    Exponential_Exponential <- list()
    Exponential_Gamma <- list()
    Exponential_HalfCauchy <- list()
    Exponential_HalfNormal <- list()
    Exponential_InverseGamma <- list()
    Exponential_LogNormal <- list()
    # LogNormal Base Priors
    LogNormal_Exponential <- list()
    LogNormal_Gamma <- list()
    LogNormal_HalfCauchy <- list()
    LogNormal_HalfNormal <- list()
    LogNormal_InverseGamma <- list()
    LogNormal_LogNormal <- list()
    # Gamma Base Priors
    Gamma_Exponential <- list()
    Gamma_Gamma <- list()
    Gamma_HalfCauchy <- list(iter = 4000)
    Gamma_HalfNormal <- list()
    Gamma_InverseGamma <- list()
    Gamma_LogNormal <- list()
    # HalfNormal Base Priors
    HalfNormal_Exponential <- list()
    HalfNormal_Gamma <- list()
    Halfnormal_HalfCauchy <- list()
    HalfNormal_HalfNormal <- list()
    Halfnormal_InverseGamma <- list()
    Halfnormal_Lognormal <- list()
  })
  return(PrComb)
}

# source(file = file.path(THIS_PATH__, "SamplerParameters.R")) 




