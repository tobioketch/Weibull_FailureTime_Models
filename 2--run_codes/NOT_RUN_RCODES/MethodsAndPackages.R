# Packages ----------------------------------------------------------------

if (!require(pacman)) {
  install.packages("pacman")
}
## Loading required package: pacman
suppressPackageStartupMessages(suppressWarnings(
  pacman::p_load(
    "rstan",
    "magrittr",
    "ggplot2",
    "pracma",
    "here",
    "kableExtra",
    # "NSM3",
    "purrr",
    "dplyr",
    "goft",
    "future",
    "furrr",
    "EnvStats",
    "rticles",
    "tidyr",
    "ggeasy",
    "knitr",
    "markdown",
    "mime",
    "rmarkdown",
    "stringi",
    "stringr",
    "tinytex",
    "xfun",
    "yaml",
    "boot",
    "cli",
    "progressr",
    "tidyselect",
    "renv",
    "this.path",
    "readr"
  )
))


# Methods -----------------------------------------------------------------

# Internal convenient function to facilitate sample size selection
# from the list of simulated data
is.nSmall <- function(lst) {
  # uses simulated sample names to
  # returns TRUE if sample size n less than 30
  nms <- names(lst)
  stringr::str_extract(nms, "n_\\d+") %>% readr::parse_number() < 30L
}

select_default_tuning_pars <- function(nm, tidyselectfun = NULL) {
  # Handy function to help select default tuning parameter
  # nm, a complete name of a prior distribution, such as exponential_gamma or
  # nm = exponential, in which case you specify tidyselectfun = "starts_with
  require("tidyselect")
  DefaultTuningPars <- .defaultTunePars()
  if (!missing(tidyselectfun)) {
    FUN <- match.fun(tidyselectfun)
    nm <- tidyselect::eval_select(FUN(nm), DefaultTuningPars)
    res <- DefaultTuningPars[nm]
  } else {
    res <- DefaultTuningPars[nm]
  }
  res
}

simulate_rweib <- function(n, shape, scale, seed = NULL) {
  # simulates a random sample from a Weibull distribution
  # n sample size, shape and scale parameters (parameters can be passed as vectors)
  # n = c(5, 10), shape = c(0.8, 3.0), scale = c(20, 1000)
  if (is.null(seed)) {
    set.seed(TRUE)
  } else {
    set.seed(seed = seed)
  }
  sim.pars <- base::expand.grid(n = n, shape = shape, scale = scale)
  Out <- dplyr::rowwise(sim.pars, n, shape, scale) %>%
    dplyr::summarise(sim_lifetime = base::list(stats::rweibull(n, shape, scale) %>% sort()), .groups = "keep") %>%
    dplyr::mutate(sim_lifetime = stats::setNames(sim_lifetime,
      nm = stringr::str_glue("n_{n}_shape_{shape}_scale_{scale}") %>% as.character()
    ))
  Out
}




lcasefold_name <- function(lst, upper = FALSE) {
  # Translates names of a list to a lower case by default
  # lst, a named list whose names to be translated
  # upper,  Boolean, TRUE, translate names to upper case, default = FALSE,
  nms <- names(lst)
  names(lst) <- casefold(nms, upper = upper)
  return(lst)
}

proper_case <- function(x) {
  # takes a string x and turns it to a proper case
  frst <- substr(x, start = 1, stop = 1)
  othr <- substr(x, start = 2, stop = 100000000L)
  resl <- paste0(toupper(frst), tolower(othr))
  resl
}

add_classic_to_priorCols <- function(x) {
  # SIMULATION DATA
  # adds rows of classic methods to the dashed ("-") rows of classic prior_comb
  x <- tidyr::separate(x,
    col = "prior_comb", into = c("shapePrior", "scalePrior"),
    fill = "right", sep = "_", remove = FALSE
  ) %>%
    dplyr::mutate(
      shapePrior = format_priors(shapePrior), scalePrior = format_priors(scalePrior),
      prior_comb = paste(shapePrior, scalePrior, sep = "-"), .keep = "unused"
    )
  x <- dplyr::mutate(x, method = format_methods(method))
  # In Simulation Method Specifies MCMC Vs Classic (equivalent of method_group in application)
  mask_classics <- x[["method"]] != "MCMC"
  x[mask_classics, "prior_comb"] <- x[mask_classics, "method"]
  x <- dplyr::select(x, -c("shapePrior", "scalePrior"))
  return(x)
}

# proper_classic_priors <- function(x, .drop = FALSE){
#   # APPLICATION DATA
#   # produces proper casing for methods(classical & priors), and prior_combination
#   # Assumes prior combination for the classical methods equals to dash ("-")
#   # It's a good idea to separate (column = methods) then process the prior combination
#   x <- tidyr::separate(x, col = "prior_comb", into = c("shapePrior", "scalePrior"), fill = "right", sep = "_") %>% 
#     dplyr::mutate(
#       shapePrior = format_priors(shapePrior), scalePrior = format_priors(scalePrior),
#       prior_comb = paste(shapePrior, scalePrior, sep = "-")
#     )
#   x <- dplyr::mutate(x, method = format_methods(method))
#   mask_classics <- x[["method"]] %in%  c("Regression", "Moments", "MLE")
#   x[mask_classics, "prior_comb"] <- x[mask_classics, "method"]
#   x[mask_classics, "scalePrior"] <- "-"
#   x["method"] <- x["prior_comb"]
#   if(.drop){
#     x <- dplyr::select(x, -"prior_comb" )
#   }
#   x
# }

methods_to_proper_case <- function(x, propcol, classicval = c("regression", "mle", "moments"), .drop = FALSE) {
  # APPLICATION DATA
  # produces proper casing for methods(classical & priors), and prior_combination
  # Also, creates shapePrior and scalePrior columns
  # propcol; column to process into proper case
  # classicval; classical values within the propcol
  # .drop = FALSE/ TRUE (drops prior_comb column if TRUE)
  x <- tidyr::separate(x, col = propcol, into = c("shapePrior", "scalePrior"), fill = "right", sep = "_", remove = FALSE) %>%
    dplyr::mutate(
      shapePrior = format_priors(shapePrior), scalePrior = format_priors(scalePrior),
      prior_comb = paste(shapePrior, scalePrior, sep = "-")
    )
  nm <- base::names(x[propcol])
  df <- stats::setNames(data.frame(format_methods(x[[propcol]])), nm = nm)
  x <- dplyr::bind_cols(df, dplyr::select(x, -tidyselect::all_of(propcol)))
  mask_classics <- x[[propcol]] %in% format_methods(classicval)
  x[mask_classics, "prior_comb"] <- x[mask_classics, propcol]
  x[mask_classics, "scalePrior"] <- x[mask_classics, "shapePrior"] <- "-"
  x[propcol] <- x["prior_comb"]
  if (.drop) {
    x <- dplyr::select(x, -"prior_comb")
  }
  x
}



#     working directory - is the the default directory
temp_plot_path <- here::here("manuscript/Figures")



# x   -----> character to be defused
# returns an expression of x without evaluating it.

defuse <- function(x) rlang::sym(x)



# Weibull Mean Residual Life
weib_mrl <- function(x, shape, scale, digits = 2L) {
  # calculates remaining lifetime of a system that has reached age x
  ll <- (x / scale)^shape
  Fx <- pgamma(ll, shape = 1 / shape, scale = 1, lower.tail = TRUE)
  mu.x <- (scale / shape) * exp(ll) * gamma(1 / shape) * (1 - Fx)
  round(mu.x, digits = digits)
}




# Utility Functions  - EXTRACT METHODS ------------------------------------

# EXTRACT FUNCTIONS - Internal Functions to help extract the a specific information
# extract shape param value
extract_shape <- function(x, shape_, scale_) {
  scl.pattern <- paste0("_", scale_, ".*")
  shp.pattern <- paste0("[", shape_, "_", "]")
  gsub(scl.pattern, "", x) %>%
    gsub(shp.pattern, "", .) %>%
    as.numeric()
}

# extract scale param value
extract_scale_mcmc <- function(x, scale_) {
  pattern <- paste0(".*", scale_, "_")
  gsub(pattern, "", x) %>%
    gsub(".[A-z].*", "", .) %>%
    as.numeric()
}

# extract prior information
extract_priors <- function(x, shape_, scale_) {
  pattern1 <- paste0(".*", scale_, "_\\d.")
  pattern2 <- paste0("\\.", shape_, "|", "\\.", scale_)
  gsub(pattern1, "", x) %>% gsub(pattern2, "", .)
}

# extract prior information- hazard data
extract_priors_hz <- function(x, scale_) {
  scl.pattern <- paste0(".*", scale_, "_\\d.")
  gsub(scl.pattern, "", x)
}


# extract names of parameters
extract_names_param <- function(x, shape_, scale_) {
  shape <- paste0("\\.", shape_)
  scale <- paste0("\\.", scale_)
  pattern <- paste0(shape, "|", scale)
  stringr::str_extract(x, pattern) %>% gsub("[.]", "", .)
}



# Plotting methods --------------------------------------------------------

# Customized theme for plots
customize_themePlot <- function() {
  theme(
    # strip.placement = "outside",
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    # panel.border = element_rect(linetype = "solid"),
    # strip.background = element_rect(colour = NA),
    plot.background = element_rect(linetype = "solid", colour = "black"), #linewidth = 0.5),
    axis.line.y = element_line(),
    axis.line.x = element_line(),
    strip.text = element_text(size = 20),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16.5),
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 22, face = "plain"),
    plot.subtitle = element_text(hjust = 0.5, size = 16)
  )
}

# paste("Shape", expression("\u03B2"), sep = ", ")

# plots analysis line graph
plot_lines <- function(data, x = Shape, y = Length, group = PriorCombOrMethod, facet.by = N, scales.facet = "free_y",
                       point.size = 3L, line.width = .7, title = NULL, lab.x = NULL, lab.y = NULL, nrow = 2) {
  if (missing(lab.y)) {
    lab.y <- eval_select(expr = enexpr(y), data) |> names()
  }
  if (missing(lab.x)) {
    lab.x <- eval_select(expr = enexpr(x), data) |> names()
  }
  quoted_group_nm <- rlang::enexpr(group) %>%
    tidyselect::eval_select(data)
  group_vals <- purrr::pluck(data, names(quoted_group_nm)) %>% unique()
  facet_vars <- rlang::enexpr(facet.by) %>%
    tidyselect::eval_select(data)

  ggplot(data = data, aes(x = {{ x }}, y = {{ y }}, group = {{ group }}, colour = {{ group }})) +
    geom_line(linewidth = line.width) +
    geom_point(aes(shape = {{ group }}), size = point.size) +
    facet_wrap(names(facet_vars), nrow = nrow, scales = scales.facet) +
    scale_shape_manual(values = seq_along(group_vals)) +
    customize_themePlot() +
    labs(title = title, x = lab.x, y = lab.y)
}

# saves a jpeg copy of a ggplot
# a wrapper function of ggsave
save_jpeg <- function(plot, file_name, path = "misc/Graphs/", height = 25, width = 35, ...) {
  ggplot2::ggsave(
    plot = plot,
    filename = paste0(file_name, ".jpeg"), device = "jpeg", path = path, units = "cm", height = height,
    width = width,
    limitsize = FALSE
  )
}



# Convenient function to draw tables of results
draw_table <- function(data, collapse_cols = 1:3, font_size = 8, ...) {
  # Draw tables
  kableExtra::kbl(data,
    longtable = FALSE, booktabs = TRUE, centering = TRUE, ...
  ) %>%
    kableExtra::collapse_rows(columns = collapse_cols, valign = "top", latex_hline = "major") %>%
    kableExtra::kable_styling(
      latex_options = c("hold_position", "repeat_header"),
      position = "center",
      font_size = font_size, full_width = FALSE
    )
}



# CLASSICAL METHODS - OUTPUT METHODS --------------------------------------


select_clean_cols <- function(x, param, ..., include = character(), exclude = character(), strict = TRUE) {
  pattn <- glue::glue(".{param}|{param}.")
  # contains used because the patterns could be jumbled (end of beginning of a column  name)
  nms <- tidyselect::vars_select(names(x), tidyselect::contains(param), ..., .include = include, .exclude = exclude, .strict = strict)
  # select columns with the patterns
  dplyr::select(x, tidyselect::all_of(nms)) %>%
    # remove pattern on the columns names (so as to remain with same column names for merging)
    dplyr::rename_with(function(x) gsub(pattn, "", x)) %>%
    # indicate the pattern column
    dplyr::mutate(Param = param)
}

# ..., additional parameter (.include) passed on to var_select function
# .include = "the common variable to be included in both splits of shape and scale portions of dt"

# Works similar to gather long except handling messy column names
# only works with two separation (patterns)
tidy_long_messy_cols <- function(dt, pattn1 = "shape", pattn2 = "scale", ...) {
  # takes dt of classic ( ML, Moments, and Rank Regression) results
  # split the dt by parameters and bind the rows together
  set1 <- select_clean_cols(dt, param = pattn1, ...)
  set2 <- select_clean_cols(dt, param = pattn2, ...)
  resl <- dplyr::bind_rows(set1, set2)
  resl
}




# Classical Methods -Parameter Estimation ---------------------------------

# functions for classical approaches
# to estimate parameters of Weibull distributions

# RANK REGRESSION METHOD
#
# Estimates parameters of the Weibull distribution using linear regression function
# x, time to failure data
# benard, approximation of the cdf (default False)

weib_regression <- function(x, index, benard = FALSE, ...) {
  # process failure data
  if (missing(index)) {
    index <- seq_along(x)
  }
  D <- x[index]
  D <- sort(D, decreasing = FALSE)
  n <- length(D)
  dln <- log(D)
  inds <- seq_along(D)
  if (!benard) Fx <- (inds - 0.5) / n
  if (benard) Fx <- (inds - 0.3) / (n + 0.4)
  y <- log(-log(1 - Fx))
  # fit reg model
  fit <- lm(y ~ dln)
  slope <- fit$coefficients[[2]]
  intercept <- fit$coefficients[[1]]
  out <- c(shape = slope, scale = exp(-(intercept / slope))) # , r.squared = broom::glance(fit)$r.squared )

  return(out)
}

# Method of Moments -------------------------------------------------------
weib_moments <- function(x, index, interval = c(0.1, 5), extendInt = "yes", ...) {
  # method of moments approach
  # uses stats::uniroot function to solve the zeros of the high dimensional beta (shape model)
  # x, sample data
  # ..., additional parameters to be passed to the stats::uniroot function
  if (missing(index)) {
    index <- seq_along(x)
  }
  D <- x[index] # bootstrap style data
  # Encloses an environment for the shape model
  ConstructShapeMOM <- function(shape0 = FALSE, x) {
    b <- shape0 # initial beta
    mu <- mean(x)
    sigma <- sd(x)
    function(b) {
      (sigma^2 + mu^2) / mu^2 - gamma(1 + 2 / b) / gamma(1 + 1 / b)^2
      # log(gamma(1 + 2 / b)) - 2 * log(gamma(1 + 1 / b)) - log(sigma^2 + mu^2) + 2 * log(mu)
    }
  }
  # estimate scale param
  WeibullScaleMOM <- function(shape, x) {
    mu <- mean(x)
    mu / gamma(1 + 1 / shape)
  }
  fun <- ConstructShapeMOM(x = D)
  out <- stats::uniroot(f = fun, interval = interval, extendInt = extendInt, ...)
  bhat <- out$root
  ahat <- WeibullScaleMOM(bhat, D)
  return(c(
    shape = bhat, scale = ahat, f.root = out[["f.root"]],
    iter = out[["iter"]], estim.prec = out[["estim.prec"]]
  ))
}

weib_summary <- function(b, a) {
  # Estimates theoretical moments of Weibull distribution
  # b , estimated shape parameter
  # a, estimated scale parameter
  Mean <- a * gamma(1 + 1 / b)
  Var <- a^2 * gamma(1 + 2 / b) - Mean^2
  return(list(mean = Mean, var = Var))
}

# METHOD, MLE -------------------------------------------------------------

weib_mle <- function(x, index, interval = c(0.01, 10), extendInt = "yes", ...) {
  # mle method to estimate weibull parameters
  # uses stats::uniroot function to solve the zeros of the high dimensional beta (shape model)
  # x, sample data
  # ..., additional parameters to be passed to the stats::uniroot function
  if (missing(index)) {
    index <- seq_along(x)
  }
  D <- x[index]
  ConstructShapeMLE <- function(shape0 = FALSE, x) {
    # shape beta model enclosure
    # Takes data x
    # used in a stats::uniroot function to return beta hat
    n <- length(x)
    shp <- shape0
    function(shp) {
      # estimates shape
      xshp <- x^shp
      lnx <- log(x)
      A <- base::sum(xshp * lnx) / base::sum(xshp)
      B <- base::sum(lnx) / n

      return(A - (1 / shp + B))
    }
  }
  ScaleMLE <- function(b, x) {
    # estimate scale
    n <- length(x)
    xb <- x^b
    ahat <- sum(xb) / n
    ahat <- ahat^(1 / b)
    return(ahat)
  }
  bbfun <- ConstructShapeMLE(x = D)
  out <- stats::uniroot(f = bbfun, interval = interval, extendInt = extendInt, ...)
  # estimated shape
  bhat <- out$root
  ahat <- ScaleMLE(bhat, x = D)
  return(c(
    shape = bhat, scale = ahat, f.root = out[["f.root"]],
    iter = out[["iter"]], estim.prec = out[["estim.prec"]]
  ))
}

# Boot Models -------------------------------------------------------------

boot_summary <- function(bootOut, alpha = 0.05, P, ...) {
  # P number of parameters in the model
  t0 <- bootOut$t0[1:P]
  nm <- names(t0[[1]])
  t1 <- matrix(bootOut$t[, 1:P], nrow = bootOut$R, ncol = P, dimnames = list(NULL, nm))
  t <- apply(t1, 2L, mean)
  stdErr <- apply(t1, 2L, sd, na.rm = TRUE)
  # Basic Bootstrap Confidence Interval
  aa <- alpha / 2
  bCI <- list()
  for (inds in seq_len(P)) {
    qq <- unname(quantile(sort(bootOut$t[, inds]), probs = c(1 - aa, aa)))
    bCI[[inds]] <- 2 * bootOut$t0[inds] - qq
  }
  bCI <- matrix(do.call("rbind", bCI), nrow = P, ncol = 2, dimnames = list(nm, c("LL", "UL")))
  return(cbind(t0 = t0, t = t, stdErr = stdErr, bCI))
}


boot_weibOr <- function(x, method = c("mle", "moments", "regression"), parallel = c("snow", "multicore", "no"),
                        R = 1000, P = 2, ncpus = getOption("boot.ncpus", 2L), alpha = 0.05, ...) {
  # x, sample data
  # ..., additional parameters to be passed to boot::boot function
  # P, number of parameters to be estimated, default = 2
  method <- match.arg(method)
  fun <- get(paste("weib", method, sep = "_"))
  out <- boot(data = x, statistic = fun, R = R, parallel = parallel[1], ncpus = ncpus, ...) |>
    boot_summary(alpha = alpha, P = P, ...)
  out
}


# MCMC Approach/ Method -----------------------------------------------

# returns stan model data
# ..., additional parameters passed on to boot & uniroot functions
# B, number of bootstraps used to summarize parameter estimates
# weib_stanx_weib
weib_stanx <- function(x, m, v, fctor, B = 1000, seed = TRUE, interval = c(0.1, 10), extendInt = "yes") {
  # x, lifetime
  # hyper-parameters m -mean, s - sigma
  Fun <- match.fun(weib_mle)
  set.seed(seed)
  bootOut <- boot::boot(x, statistic = Fun, R = B, interval = interval, extendInt = extendInt)
  prior_mean <- apply(bootOut$t, 2L, mean)
  prior_stdErr <- apply(bootOut$t, 2L, sd)
  # model data
  Dt <- within(list(), {
    n <- length(x)
    y_n <- x
    mu_shape <- prior_mean[1]
    var_shape <- prior_stdErr[1]^2
    mu_scale <- prior_mean[2]
    var_scale <- prior_stdErr[2]^2
    m <- as.numeric(m)
    v <- as.numeric(v)
    fctor <- as.numeric(fctor)
  })
  Dt
}

.checkTuningParamSamplingModel <- function(tunepars, comparedefault, usermodel, addmodel, tidyselectfun) {
  # checks conditions to define tuning parameters and sampling model(s)
  if (all(tunepars != "default" & (is.null(comparedefault) || !is.null(comparedefault)) & (is.null(usermodel) || !is.null(usermodel)))) {
    tunepars <- lcasefold_name(tunepars) # define tuning parameters
  } else {
    tunepars <- lcasefold_name(.defaultTunePars())
  }
  if (!is.null(usermodel)) {
    SmplngMdl <- lcasefold_name(usermodel, upper = FALSE) # define sampling models
  } else {
    SmplngMdl <- .defaultSamplingModls() # named list of default models
    SmplngMdl <- lcasefold_name(SmplngMdl, upper = FALSE) # format names to lower cases
  }
  if (all(tunepars == "default" & is.null(usermodel) & !addmodel)) {
    # From condition 1) # addmodel =FALSE (default)
    if (is.null(comparedefault)) {
      tunepars <- tunepars
      SmplngMdl <- SmplngMdl
      cat("\n runs all default models on all default parameterizations\n")
    } else {
      # # 3) shorter lists for tuning parameters and Sampling Models
      if (!missing(tidyselectfun)) {
        FUN <- get(tidyselectfun)
        nm.compredefalt <- tidyselect::eval_select(FUN(comparedefault), SmplngMdl)
        nmTunepars <- tidyselect::eval_select(FUN(comparedefault), tunepars)
      } else {
        nm.compredefalt <- tidyselect::eval_select(comparedefault, SmplngMdl) # makes sense (it's not in if condition)
        nmTunepars <- tidyselect::eval_select(comparedefault, tunepars)
      }
      SmplngMdl <- SmplngMdl[nm.compredefalt]
      tunepars <- tunepars[nmTunepars]
      checkTrue <- all(nm.compredefalt %in% names(tunepars))
      if (!checkTrue) { # for user with no complete list of default models to be specified under "comparedefault"
        stop(paste("Here's the list of all default models\n", names(tunepars), sep = "--"))
      }
    }
  }

  if (all(tunepars != "default" & is.null(comparedefault))) {
    if (is.null(usermodel)) {
      ## 5 ) define new parameters for a few selected models
      if (addmodel) {
        stop("Supply and add a user model")
      } else {
        # check if all tunepars are in default
        nmTunepars <- names(tunepars)
        nmSmplngMdl <- names(SmplngMdl)
        checkNames <- all(nmTunepars %in% names(SmplngMdl))
        if (!checkNames) {
          stop(paste("Here's the list of all default models\n", names(SmplngMdl), sep = "--"))
        }
        SmplngMdl <- SmplngMdl[nmTunepars]
      }
    } else {
      # 6 )user model provided and newly specified parameters (default model)
      if (addmodel) {
        SmplngMdl <- c(usermodel, SmplngMdl)
        # check if tunepars defined for all SmplngMdl
        nmTunepars <- names(tunepars)
        if (!all(nmTunepars) %in% names(SmplngMdl)) {
          stop(paste("some parameters not for default models. Here are the default models", names(SmplngMdl)))
        }
      } else {
        # when addmodel is false # just the provided user model
        SmplngMdl <- lcasefold_name(usermodel)
        if (!all(names(tunepars) %in% names(SmplngMdl))) {
          stop("specify parameters for each model")
        }
      }
    }
  }
  if (all(tunepars != "default" & !is.null(comparedefault))) {
    if (is.null(usermodel)) {
      if (addmodel) {
        stop("User model is missing")
      } else {
        # User model missing & addmodel is FALSE
        comparedefault <- casefold(comparedefault) # takes care of the case for comparisons
        if (!missing(tidyselectfun)) {
          FUN <- get(tidyselectfun)
          nm.compredefalt <- tidyselect::eval_select(FUN(comparedefault), SmplngMdl)
          nmTunepars <- tidyselect::eval_select(FUN(comparedefault), tunepars)
        } else {
          nm.compredefalt <- tidyselect::eval_select(comparedefault, SmplngMdl) # makes sense (it's not in if condition)
          nmTunepars <- tidyselect::eval_select(comparedefault, tunepars)
        }
        SmplngMdl <- SmplngMdl[nm.compredefalt]
        tunepars <- tunepars[nmTunepars]
        if (!all(names(tunepars) %in% names(SmplngMdl))) {
          stop("specify parameters for each model")
        }
      }
    } else {
      # # 8 ) user model provided
      if (addmodel) {
        SmplngMdl <- c(comparedefault, SmplngMdl) # add user model to compare model
        tunepars <- tunepars
        if (!all(names(tunepars) %in% names(SmplngMdl))) {
          stop("specify parameters for each model")
        }
      } else {
        # not adding user model to compare model when both are specifies
        stop("Change Addmodel= FALSE (default) to admodel = TRUE")
      }
    }
  }

  Out <- list(tunepars = tunepars, SmplngMdl = SmplngMdl)
}


# implements parallel stan mcmc, a two parameter weibull likelihood by default
# update, define a user model
# x, a list of observed failure time
# comparedefault = NULL (default) to compare combinations from two or more base priors (priors for the shape)
# tunepars = "default" , tuning parameters - otherwise use within(list(),{"specify prior_comb = list(pars)"})
# to specify parameters when using user models or when adding user models to comparedefault ( compare with a few selected models)
# m = 0, v = 1, fctor = 1 passed on to weib_stanx (stan data model)

parstan <- function(x, tunepars = "default", chains = 4L, xfun = weib_stanx, usermodel = NULL, addmodel = FALSE, trtitle = NULL,
                    comparedefault = NULL, tidyselectfun = NULL, summarise = TRUE, probs = c(0.025, 0.975), m, v, fctor = 1,
                    parspattn = "beta|alpha", vrcovpattn = "varCov", infopattn = "fisher") {
  isNumeric <- all(sapply(x, \(x)is.numeric(x)))
  if (length(x) == 0 || !isNumeric) {
    stop("provide a list of numeric vector for x")
  }
  stopifnot(is.list(x))
  if (all(tunepars != "default" & !is.list(tunepars))) {
    stop("tunepars should be defined as tunepars = within(list(),{\"specify prior_comb = list(pars)\"})")
  }
  WngUserModl <- "define user models and the parameters to add or specify addmodel = FALSE (default)"
  Wrng.UserModlPar <- "Specify parameters of the user model"

  # 1) missing user model to be added
  if (all(tunepars == "default" & is.null(comparedefault) & is.null(usermodel) & addmodel)) {
    stop(WngUserModl)
  }
  if (all((tunepars == "default" & !is.null(usermodel)) & (!is.null(comparedefault) || is.null(comparedefault)))) {
    stop(Wrng.UserModlPar) # 4 & 2 # parameters required for the user specified model
  }
  ParsModel <- .checkTuningParamSamplingModel(
    tunepars = tunepars, comparedefault = comparedefault,
    usermodel = usermodel, addmodel = addmodel, tidyselectfun = tidyselectfun
  )
  tunepars <- ParsModel$tunepars
  SmplngMdl <- ParsModel$SmplngMdl
  # return(tunepars)
  future::plan(future::multisession, workers = getOption("mc.cores", default = 2L)) # Setting up multiple workers
  stantryError <- try(
    Modl <- purrr::map(x, \(D) {
      tunepars %>%
        furrr::future_map2(., names(.), \(param, nmPrior) {
          # Parameters for the data function should be more open, to cover for any possible user defined model
          Dt <- do.call(xfun, list(x = D, m = m, v = v, fctor = fctor))
          do.call(as.function(rstan::stan), c(list(data = Dt, model_code = SmplngMdl[[nmPrior]], chains = chains), as.list(param)))
        }, .options = furrr::furrr_options(seed = TRUE))
    }) %>% unlist(recursive = FALSE),
    silent = TRUE
  )
  if (inherits(stantryError, "try-error")) {
    matrx <- matrix(ncol = 2, nrow = 2)
    warning("The model did not run smoothly\n")
    print(stantryError)
    Out <- list(pars = rep(NA, 2), fisherInfo = matrx, VarCov = matrx, x = x)
    return(Out)
  }
  # Closing down multiple workers
  future::plan(future::sequential)
  if (summarise) {
    # Diagnostic plots
    trc <- purrr::map(Modl, \(m) trace_plt(m, title = trtitle, chains = chains))
    div <- purrr::map(Modl, \(m) sampler_divergence(m))
    # Model Summary
    ModlSumry <- purrr::map(Modl, \(m) rstan::summary(m, probs = probs)$summary)
    pars <- purrr::map(ModlSumry, \(s) extract_parameter(s, pattern = parspattn) %>%
      tibble::rownames_to_column(var = "param"))
    fisher <- map(ModlSumry, \(s)extract_parameter(s, pattern = infopattn) %>%
      to_matrix_measure(pattern = infopattn))
    VrCv <- purrr::map(ModlSumry, \(s) extract_parameter(s, pattern = vrcovpattn) %>%
      to_matrix_measure(pattern = vrcovpattn))
    Out <- list(estimatedPars = pars, fisherInfo = fisher, VarCov = VrCv, trace = trc, divergence = div)
  } else {
    Out <- list(model = Modl, trace = trc, divergence = div)
  }
  return(Out)
}

# parstanobj, parstan object
# shape = beta (default), aliase for the shape parameter in the stanpar object
# scale = "alpha (default), aliase for the shape parmeter (stanpar object)
summarise_stan_pars <- function(parstanobj, shape = "beta", scale = "alpha") {
  # Summarizes the parameters against measures of model fit
  # Total Sampling Variance and Parameter Coverage
  ParStats <- parstanobj$estimatedPars %>%
    purrr::map(\(x) dplyr::select(x, -c("se_mean", "n_eff", "Rhat")) %>%
      pivot_wider(names_from = "param", values_from = -"param")) %>%
    dplyr::bind_rows(.id = "method") %>%
    dplyr::rename_with(\(x) gsub(shape, "shape", x) %>% gsub(scale, "scale", .)) %>%
    dplyr::mutate(total_sampling_var = sd_shape^2 + sd_scale^2)
  # Total Fisher Information
  TotFisher <- parstanobj$fisherInfo %>%
    purrr::map(\(f)c(total_fisher_info = sum(diag(f)) - 2 * f[row(f) == col(f) + 1])) %>%
    dplyr::bind_rows(.id = "method")
  # Total Variance
  TotVar <- parstanobj$VarCov %>%
    purrr::map(\(v)c(total_var = sum(diag(v)) - 2 * v[row(v) == col(v) + 1])) %>%
    dplyr::bind_rows(.id = "method")
  Out <- right_join(right_join(ParStats, TotFisher, by = "method"), TotVar, by = "method") %>%
    dplyr::rename_with(\(x) gsub("alpha", "scale", x) %>% gsub("beta", "shape", .)) %>%
    dplyr::relocate(contains("total"), .after = "mean_scale") %>%
    dplyr::relocate(contains("%_shape"), .before = contains("%_scale"))
  return(Out)
}

weighted_rel_eff <- function(ssobj) {
  # ssobj, summarized stan pars object
  nms <- names(ssobj)
  Out <- ssobj %>%
    dplyr::mutate(
      wgt_tot_var = sum(total_var),
      wgt_tot_samplng = sum(total_sampling_var),
      prop_tot_var = total_var / wgt_tot_var,
      prop_tot_samplng = total_sampling_var / wgt_tot_samplng,
      wgt_rel_eff = prop_tot_samplng / prop_tot_var
    ) %>%
    dplyr::select(tidyselect::all_of(nms), "wgt_rel_eff") 
  return(Out)
}

method_to_priorcomb_cols <- function(sstnobj) {
  # Creates two more columns (prior_comb, method_group) from the method column
  Out <- dplyr::mutate(sstnobj,
    prior_comb = method,
    prior_comb = dplyr::case_when(prior_comb %in% c("regression", "mle", "moments") ~ "-", TRUE ~ prior_comb),
    # method_group = dplyr::case_when(method %in% c("regression", "mle", "moments") ~ "classic", .default = "mcmc"), .after = "method"
    method_group = dplyr::case_when(method %in% c("regression", "mle", "moments") ~ "classic", TRUE ~ "mcmc"), .after = "method"
  ) %>% dplyr::relocate(prior_comb, .after = "method")
  return(Out)
}


# sstnobj, summarized_stan_par object
# returns, a method column with mcmc entries,
# simulated data entries: n (sample size), values for the shape and scale parameters
separate_method_pars_cols <- function(sstnobj) {
  # process method column of summarized simulated stan par object
  resl <- dplyr::mutate(
    sstnobj,
    # creates three columns: a prior combination,
    # sim.params (separated further into n, shape & scale),
    # method = 'mcmc' columns
    prior_comb = gsub(".*_scale_\\d+.", "", method),
    sim.params = stringr::str_extract(method, ".*_scale_\\d+"), method = "mcmc", .keep = "unused"
  ) %>%
    tidyr::separate(col = "sim.params", into = c(NA, "n", NA, "shape", NA, "scale"), convert = TRUE, remove = TRUE, sep = "\\_") %>%
    dplyr::relocate(c("method", "prior_comb", "n", "shape", "scale"), .before = "mean_shape")
  return(resl)
}


separate_data_method_cols <- function(sstanobj, drop = FALSE) {
  # # separates summarized stan object's method column into
  # # data and method columns
  # Out <- tidyr::separate_wider_delim(sstanobj,
  #   cols = "method",
  #   delim = ".", names = c("Data", "method")
  # )
  # if (drop) {
  #   Out <- dplyr::select(Out, everything(), -Data)
  # }
  Out <- tidyr::separate(sstanobj, col = "method", into = c(NA, "method"), sep = "[.]", remove = drop)
  Out
}


# stan_obj, stan_fit (model object)
# n_chains, number of chains used to generate the stan_obj
# title, desired plot title
# ..., other arguments passed to stan_trace, and
#       optional plot arguments (linetype, size, alpha) passed to geom_path

trace_plt <- function(stan_obj, chains, title, add_burn_in = TRUE, ...) {
  suppressMessages(
    rstan::stan_trace(stan_obj, inc_warmup = add_burn_in, ...) +
      ggtitle({{ title }}) +
      scale_color_manual(values = rep("black", chains)) +
      theme_classic() +
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
  )
}

# Sampler parameters and Divergence
sampler_pars <- function(fitObj) {
  # lists by number chains used
  # combine chains together into a tbl
  samplerPar <- rstan::get_sampler_params(fitObj, inc_warmup = FALSE) %>% do.call("rbind", .)
  logPost <- rstan::get_logposterior(fitObj, inc_warmup = FALSE) %>% unlist(recursive = TRUE)
  samplerPar <- dplyr::bind_cols(samplerPar,
    logPosterior = logPost,
    Divergent = dplyr::if_else(rstan::get_divergent_iterations(fitObj), "Divergent", "Not Divergent")
  )
  samplerPar
}

# plots Mean Metropolis acceptance rate Vs divergence (Divergent / Not Divergent)
# samplerObj, generated by sampler_pars
#
sampler_divergence <- function(stanObj) {
  sampParam <- sampler_pars(fitObj = stanObj)
  ggplot2::ggplot(sampParam, aes(x = Divergent, y = accept_stat__)) +
    geom_violin() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"), axis.title = element_text(face = "bold")) +
    labs(title = "Divergence", y = ("Mean Metropolis Accepatance Rate"), x = NULL)
}

# ..., other parameters passed on to underlying print method of the monitor
# extract: inc_warmup = FALSE (default) if TRUE adjust/ indicate the number of warmups,
#           include = TRUE (default), Logical, indicates whether parameter specified in pars
#                     should be included in the summary or not
#           pars - can be specified - otherwise all parameters are summarized

# iter:  the number of iterations after warmup ( when warmup = 0) /Or
#       the total number of iterations including warmup (when warmup not equal to zero)

# monitor : print = FALSE (default)
stan_monitor <- function(x, iter, chains, probs = c(0.025, 0.975), inc_warmup = FALSE, warmup = 0, print = FALSE, ...) {
  sims <- rstan::extract(x, permuted = FALSE, inc_warmup = inc_warmup)
  resl <- rstan::monitor(sims, probs = probs, warmup = warmup, print = print, ...)

  resl <- tibble::rownames_to_column(as.data.frame(resl), var = "param")

  tibble::as_tibble(resl) %>%
    dplyr::mutate(
      n_eff_valid = dplyr::case_when((n_eff / iter) > 0.01 ~ 1, TRUE ~ 0),
      B_ESS_valid = dplyr::case_when((Bulk_ESS / chains) > 100 ~ 1, TRUE ~ 0),
      T_ESS_valid = dplyr::case_when((Bulk_ESS / chains) > 100 ~ 1, TRUE ~ 0),
      ESS_valid = dplyr::case_when(B_ESS_valid == 1 & T_ESS_valid == 1 ~ 1, B_ESS_valid == 0 | T_ESS_valid == 0 ~ 0)
    )
}

###
# sumryObj, a stan summary stanObj
# pattern, regex type expression
# ..., additional arguments passed on to str_detect function (negate = FALSE)
extract_parameter <- function(stanObj, pattern, ...) {
  # to extracts parameters by fixed string matching
  df <- as.data.frame(stanObj)
  df[stringr::str_detect(rownames(df), pattern, ...), ]
}
# Converts MCMCM object for 2 parameter Weibull distribution
# to matrix format of the Fisher Information & Var-Cov parameters
to_matrix_measure <- function(stanObj, pattern) {
  df <- as.data.frame(stanObj)
  Dt <- df[stringr::str_detect(rownames(df), pattern), ]
  Dt <- Dt$mean
  mtrx <- matrix(nrow = 2, ncol = 2)
  mtrx[1, 1] <- Dt[1]
  mtrx[1, 2] <- mtrx[2, 1] <- Dt[2]
  mtrx[2, 2] <- Dt[4]
  return(mtrx)
}


# applies mean summary to a data frame of MCMC objects
## ... additional parameters passed on to as.data.frame
summarise_stan_chains <- function(stanObj, FUN, ...) {
  # permutes and mixes chains for each parameter together
  fn <- match.fun(FUN)
  df <- as.data.frame(stanObj, ...)
  Out <- sapply(df, fn)
  Out
}

# uses bootstrap approach on the methods of regression, moments and mle
# to estimates mean and standard error of parameters ( beta and alpha) of the weibull model
# x, vector of data sample
# iter, number of bootstrap iterations
# method, for estimating Weibull parameters
# probs, for quantile intervals of the estimated parameters
# ..., additional parameters for the function of methods

boot_weib_summaries <- function(x, iter = 1000, method = c("mle", "regression", "moments"), probs = c(0.025, 0.975), seed = NULL, ...) {
  case_name <- match.arg(method)
  FUN <- get(paste("weib", method, sep = "_"))
  set.seed(seed)
  # bootstrap data
  BtDt <- purrr::map(base::seq_len(iter), ~ base::sample(x, replace = TRUE) %>% base::sort())
  # bootstrap iterations of parameters
  BtPars <- purrr::map(BtDt, FUN, ...)
  TotFisher <- purrr::map(BtPars, \(f){
    tfInfo <- weib2par_fisherInfo(n = length(x), shape = f[["shape"]], scale = f[["scale"]])
    Out <- c(total_fisher_info = sum(diag(tfInfo)) - 2 * tfInfo[row(tfInfo) == col(tfInfo) + 1])
    Out
  })
  TotVar <- purrr::map(BtPars, \(v){
    totVar <- weib2par_varCov(x = x, shape = v[["shape"]], scale = v[["scale"]])
    Out <- c(total_var = sum(diag(totVar)) - 2 * totVar[row(totVar) == col(totVar) + 1])
    Out
  })
  BtPars <- BtPars %>%
    dplyr::bind_rows() %>%
    purrr::modify_if(is.numeric, round, digits = 4L)
  # Quantile intervals
  Quant <- apply(BtPars, 2, quantile, probs = probs)
  Quant <- tibble::rownames_to_column(as.data.frame(Quant), var = "limit") %>%
    pivot_wider(names_from = "limit", values_from = -"limit")
  # re-aligning names for back compatibility
  ShapDt <- dplyr::select(Quant, tidyselect::contains("shape")) # Shape
  nms1 <- paste(names(ShapDt), "shape", sep = "_") %>% gsub("shape_", "", .)
  ShapDt <- setNames(ShapDt, nms1)
  SclDt <- dplyr::select(Quant, tidyselect::contains("scale")) # Scale
  nms2 <- paste(names(SclDt), "scale", sep = "_") %>% gsub("scale_", "", .)
  SclDt <- setNames(SclDt, nms2)
  Quant <- dplyr::bind_cols(ShapDt, SclDt)
  # total fisher and variances (also a tibble with rows as long as boot iterations)
  Info <- cbind(total_fisher_info = TotFisher, total_var = TotVar) %>%
    as_tibble() %>%
    unnest(cols = c("total_fisher_info", "total_var"))
  # bootstrap iterations of parameters and measures of information
  Bts <- bind_cols(as.data.frame(BtPars), Info) %>%
    dplyr::relocate(c("total_fisher_info", "total_var"), .after = "scale") %>%
    dplyr::rename(mean_shape = "shape", mean_scale = "scale")
  StatsDt <- Bts[, c("mean_shape", "mean_scale", "total_fisher_info", "total_var")]
  Avg <- apply(StatsDt, 2, mean) %>%
    data.matrix() %>%
    t()
  sd_pars <- apply(StatsDt[, c("mean_shape", "mean_scale")], 2, sd) %>% t()
  sd_pars <- data.frame(sd_pars) %>%
    rename_with(\(x)gsub("mean", "sd", x))
  total_sampling_var <- sum(sd_pars^2)
  stats <- data.frame(Avg, total_sampling_var, sd_pars)
  stats <- dplyr::bind_cols(stats, Quant) %>%
    dplyr::relocate(tidyselect::contains("sd"), .before = tidyselect::contains("total"))
  Out <- list(stats = stats, boots = Bts)
  return(Out)
}


asympCI <- function(stat, std_error, alpha) {
  # computes asymptotic CI
  q <- qnorm(alpha / 2, lower.tail = FALSE)
  data.frame(LL = stat - q * std_error, UL = stat + q * std_error)
}

# Exact Fisher Information and VAR-COV ------------------------------------

weib2par_varCov <- function(x, shape, scale) {
  Out <- matrix(nrow = 2, ncol = 2)
  n <- length(x)
  Out[2, 2] <- 1.1087 * (scale / shape)^2
  Out[1, 2] <- Out[2, 1] <- 0.2570 * scale
  Out[1, 1] <- 0.6079 * shape^2
  Out <- Out / n
  return(Out)
}

weib2par_fisherInfo <- function(x = NULL, n = NULL, shape, scale) {
  if (missing(x) && missing(n)) {
    stop("x or n is required inorder to proceed", call. = FALSE)
  }
  if (!missing(x)) {
    n <- length(x)
  } else {
    n <- n
  }
  Out <- matrix(ncol = 2L, nrow = 2L)
  Out[1, 1] <- 1.823680 * (n / shape^2)
  Out[2, 2] <- n * (shape / scale)^2
  Out[1, 2] <- Out[2, 1] <- -0.422784 * n / scale
  return(Out)
}


# Convenient Table Processing Functions -----------------------------------

# obj, summary_stan_par object
# shape(shape value),      scale (scale value)
# llshape (lower limit shape), ulshape (upper limit shape)
# llscale (lower limit scale), upscale (upper limit scale) parameter(s)
weib2par_coverage <- function(obj, shape, scale, llshape, ulshape, llscale, ulscale) {
  # Deals with coverage of estimated parameters
  resl <- mutate(obj,
    coverage = (.data[[shape]] > .data[[llshape]] & .data[[shape]] < .data[[ulshape]]) &
      (.data[[scale]] > .data[[llscale]] & .data[[scale]] < .data[[ulscale]])
  )
  return(resl)
}

# Coverage for both shape and the scale parameter
# head(Args.dt[, -c(2:5,8:9)] )%>%
# Convenient function to format priors
format_priors <- function(x){
  x = proper_case(x)
  x = case_when(
    x == "Halfcauchy" ~ "HalfCauchy",
    x == "Inversegamma" ~ "InvGamma",
    x == "Lognormal" ~ "LogNormal",
    x == "Halfnormal" ~ "HalfNormal", TRUE ~ x
    # x == "Exponential" ~ "Exp", TRUE ~ x
  )
}
format_methods <- function(x) {
  x <- proper_case(x)
  resl <- ifelse(x == "Mcmc", "MCMC", ifelse(x == "Mle", "MLE", x))
  resl
}

# **********************************************************
# Prepare Analysis data for plotting and tabulation
# **********************************************************


# ceiling for efficiency
ceil_eff <- function(x, eff){
  Out <- dplyr::mutate(x, eff_state = .data[[eff]] <= 1.01) %>%
    dplyr::filter(eff_state == TRUE)  %>%
    dplyr::select(-eff_state)
  Out
}

rank_methods <- function(.data, rank.by, n, descend = TRUE, ...) {
  # .data, dataframe/or tibble
  # rank.by, variable to use for ranking data
  # n, top n rows of data to return
  # descend, to order ranks either ascending (descend = FALSE) or descending before returning to n rows
  df <- base::within(.data, Ranks <- base::rank(.data[rank.by], ...))
  if (descend) {
    df <- dplyr::arrange(df, desc(Ranks))
  } else {
    df <- dplyr::arrange(df, Ranks)
  }
  df <- utils::head(df, n = n)
  Out <- subset(df, select = -Ranks)
  return(Out)
}




# Simulation Data ---------------------------------------------------------

.Mean_Eff <- function(tib){
  Out <- tib %>%
    dplyr::group_by(hazard, n, prior_comb) %>%
    dplyr::mutate(
      mean_rel_wgt_eff = mean(rel_wgt_eff)
    ) %>%
    dplyr::arrange(hazard, n, prior_comb) %>%
    dplyr::ungroup() 
  return(Out)
}

.Grouped_Mean_Eff <- function(tib, ...){
  Out <- tib %>%
    # hazard, n, priors_methods
    dplyr::group_by(...) %>%
    dplyr::mutate(
      # mean_fisher_info = mean(total_fisher_info),
      mean_rel_wgt_eff = mean(rel_wgt_eff)
    ) %>%
    # hazard, n, priors_methods
    dplyr::arrange(...) %>%
    dplyr::ungroup() 
  return(Out)
}


# mutating weighted efficiency
.Mutate_WeightedEff <- function(tib){
  nms <- names(tib)
  Out <- tib %>%
    # all methods at each simulated data
    dplyr::group_by(n, shape) %>%
    dplyr::mutate(
      wgt_tot_var = sum(total_var),
      wgt_tot_samplng = sum(total_sampling_var),
      prop_tot_var = total_var / wgt_tot_var,
      prop_tot_samplng = total_sampling_var / wgt_tot_samplng,
      rel_wgt_eff = prop_tot_samplng / prop_tot_var
    ) %>%
    dplyr::arrange(n, shape) %>%
    dplyr::ungroup() %>% 
    dplyr::select(tidyselect::all_of(nms), "rel_wgt_eff")
  return(Out)
}


.Rank_PlotDt <- function(.data, rank.by, top.n, descend = TRUE) {
  # top.n, top ranked methods
  # group, grouping variable for ranking by groups
  # rank.by, ranking variable
  # descend
  # methods = mcmc, regression, mle, and moments
  # Basically, ranking all mcmc at each n & shape
  Out <- split(.data, list(.data[["method"]], .data[["n"]], .data[["shape"]]) ) %>%
    purrr::map(~ rank_methods(.x, rank.by = rank.by, n = top.n, descend = descend)) %>%
    purrr::reduce(rbind)
  Out
}


.Rank_TableDt <- function(.data, rank.by, top.n, descend = TRUE) {
  # top.n, top ranked methods
  # group, grouping variable for ranking by groups
  # rank.by, ranking variable
  # descend
  # methods = mcmc, regression, mle, and moments
  # Basically, ranking all mcmc at each n & shape
  Out <- split(.data, list(.data[["hazard"]], .data[["n"]], .data[["method"]]) ) %>%
    purrr::map(~ rank_methods(.x, rank.by = rank.by, n = top.n, descend = descend)) %>%
    purrr::reduce(rbind)
  Out
}




as_analysis_sim_df <- function(ssobj, top.n= 5, table = FALSE, descend = TRUE){
  # WORKED FOR TABLE = FALSE
  # ssobj, summarized stan object
  # top.n = 5 (top methods)
  # table = FALSE (produces a df for line plots) otherwise concise table
  # descend = TRUE (methods ranked in descending order otherwise ascending)
  ssobj <- add_classic_to_priorCols(ssobj)
  OutWgt <- .Mutate_WeightedEff(ssobj)
  MCMC <- dplyr::filter(OutWgt, method=="MCMC") %>%
    ceil_eff(eff = "rel_wgt_eff")
  Classic <- dplyr::filter(OutWgt, method!="MCMC")
  OutDt <- bind_rows(MCMC, Classic)
  OutDt <- .Rank_PlotDt(OutDt, rank.by = "rel_wgt_eff", top.n = top.n, descend = descend)
  if(table){
    TabDf <- dplyr::mutate(OutDt, hazard = dplyr::case_when(shape<1 ~ "DHR", TRUE ~ "IHR") )
    
    OutDt <- .Mean_Eff(TabDf) %>%
      .Rank_TableDt(rank.by = "mean_rel_wgt_eff", top.n = top.n, descend = descend)
    
    OutDt <- dplyr::select(OutDt, n, hazard, prior_comb, mean_rel_wgt_eff) %>%
      base::unique() %>%
      dplyr::arrange(n, hazard, desc(mean_rel_wgt_eff)) %>%
      purrr::modify_if(is.numeric, base::round, digits = 3L)
  }
  OutDt
}






percent_change <- function(from, to) {
  out <- ((to - from) / from) * 100
  out
}
