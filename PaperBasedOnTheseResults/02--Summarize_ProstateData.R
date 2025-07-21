
base::rm( list = ls( all.names = TRUE )  )

# Loading Functions, Packages and Prostate Cancer Analysis Data -----------

base::source(file = "PaperBasedOnTheseResults/MethodsAndPackages.R")


ProstCancer <- c(
  rep(0, 3), 2, 3, 4, 6, rep(7, 2), 8, rep(9, 2), rep(11, 3), rep(12, 3), rep(15, 2),
  rep(16, 3), rep(17, 2), 18, rep(19, 2), 20, 21, rep(22, 2), 23, 24, rep(25, 2), rep(26, 3),
  rep(27, 2), rep(28, 2), rep(29, 2), 30, 31, rep(32, 3), rep(33, 2), 34, 35, 36, rep(37, 2),
  38, 40, rep(41, 2), rep(42, 2), 43, rep(45, 3), 46, rep(47, 2), rep(48, 2), 51, rep(53, 2),
  rep(54, 2), 57, 60, 61, rep(62, 2), 67, 69, 87, rep(97, 2), 100, 145, 158
)

# t = 0, Taken as a Diagnosis Date
ProstCancer <- ProstCancer[ProstCancer > 0]

LifeTime <- list(ProstCancer = ProstCancer)


base::load(file = "PaperBasedOnTheseResults/misc_PaperBasedOnTheseResults/app_misc_results/RealDtPoolMetrics_PaperBasedOnThisData")


# MCMC workflow Prostate Data
# 1. compute weighted_rel_eff (combined methods - MCMC & Classic Methods)
# 2. filter ceil_eff methods (MCMC)
# 3. Rank top methods (MCMC)
# 4. Combine  (MCMC and Classic Output Data)

MetricsProstDt <- weighted_rel_eff(RealDtPoolMetrics) #(combined methods - MCMC & Classic Methods)
MetricsProstDt

# MCMC Methods
MetricsProstDt_mcmc <- MetricsProstDt %>% 
  dplyr::filter(method_group=="mcmc") %>% 
  # filter out MCMC-Methods with WRE greater than 1 (Non Efficient Methods)
  ceil_eff(eff = "wgt_rel_eff") 
  
# Classic Methods
MetricsProstDt_class <- dplyr::filter(MetricsProstDt, method_group!="mcmc") %>% 
  dplyr::arrange(wgt_rel_eff)

MetricsProstDt <- dplyr::bind_rows(MetricsProstDt_mcmc, MetricsProstDt_class) %>%
  # proper_classic_priors()
  methods_to_proper_case(propcol = "method")

Top3Methods <- base::split(MetricsProstDt, ~ method_group) %>% 
  purrr::map(\(x)rank_methods(x, rank.by = "wgt_rel_eff", n = 3,descend = TRUE)) %>% 
  purrr::reduce(rbind) 

# Top3Methods

ErrorBarPlot <- ggplot(Top3Methods, aes(x = reorder(method, wgt_rel_eff), y = mean_shape)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_shape - sd_shape, ymax = mean_shape + sd_shape), width = 0.2, linewidth = 0.5) +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  customize_themePlot() +
  labs(
    title = "Precision, Shape Parameter",
    subtitle = 'Prostate Cancer Data', # From Args- SubTitle
    x = "Prior Combinations/ Methods", y = paste0("Shape (", expression("\u03B2"), ") Parameter")
  )

ErrorBarPlot

save_jpeg(ErrorBarPlot,
          file_name = "Fig3",
          path = "PaperBasedOnTheseResults/demoPlots_PaperBasedonThesePlots/"
)


# Mean Residual Life (MRL) Plot  ------------------------------------------

# Parameters for MEAN RESIDUAL LIFE, MRL(x)
ParamsMRL <- Top3Methods %>%
  select(method, contains(c("mean"))) %>%
  tibble::add_row(
    method = "Integrated", mean_shape = mean(RealDtPoolMetrics$mean_shape),
    mean_scale = mean(RealDtPoolMetrics$mean_scale)
  )
# A sample of lifetime data
set.seed(TRUE)
SampleLifeTime <- LifeTime %>%
  purrr::map(\(x)sample(x, size = 10L, replace = FALSE) %>% sort())

MRLifeDt <- ParamsMRL %>%
  dplyr::rowwise(mean_shape, mean_scale) %>%
  dplyr::mutate(time = SampleLifeTime, meanResidLife = list(weib_mrl(time, shape = mean_shape, scale = mean_scale))) %>%
  tidyr::unnest_longer(col = c("time", "meanResidLife")) %>%
  dplyr::rename_with(\(x)gsub("mean_", "", x)) %>%
  dplyr::select(method, time, meanResidLife)

# MRL Plot
MRLPlot <- with(MRLifeDt, {
  ggplot(MRLifeDt, aes(x = time, y = meanResidLife, color = method, group = method)) +
    geom_line(linewidth = 0.5) +
    geom_point(aes(shape = method)) +
    scale_shape_manual(values = seq_along(unique(method))) +
    customize_themePlot() +
    labs(x = "Time (Months) ", y = "Average Remaining Lifetime", title = "Mean Residual Life", subtitle = "Prostate Cancer Data")
})

print(MRLPlot)

save_jpeg(
  plot = MRLPlot,
  # DtNames -- From Args
  file_name = "Fig4",
  path = "PaperBasedOnTheseResults/demoPlots_PaperBasedonThesePlots/"
)


# Sampled Mean Residual Life-Table ------------------------------------------

SampleMRLifeTable <- MRLifeDt %>%
  tidyr::pivot_wider(names_from = method, values_from = meanResidLife, values_fn = list) %>%
  tidyr::unnest(cols = -time) %>%
  dplyr::select(time, Integrated, everything()) %>%
  dplyr::relocate(MLE, .before = "Regression") %>%
  dplyr::relocate(tidyr::contains(c("Gamma", "Exponential") ), .before = "MLE")

nms <- setdiff(names(SampleMRLifeTable), c("time", "Integrated"))

SampleMRLifeTable <- mutate(SampleMRLifeTable, across(tidyselect::all_of(nms), ~ ((.x / Integrated - 1) * 100) %>%
                                                        scales::number(accuracy = 0.01, suffix = "%", style_positive = "plus", big.mark = ",")))


# Table
CaptionLifeTab <- stringr::str_glue(
  "Mean residual lifetime of top three MCMC-based and the classic Weibull model of a randomly sampled prostate cancer survival times.
  The Integrated model uses the actual parameter values {expression(\u03B2)} = {round(mean(RealDtPoolMetrics$mean_shape),2L)} and {expression(\u03B1)} ={round(mean(RealDtPoolMetrics$mean_scale),2L)}, 
  pooled from the results of all the models under the study.
  MCMC and Classic results represent percent deviations from the true mean residual lifetime of the integrated method",
  notation = "number"
)

# Remaining Life of the 10 randomly selected lifetimes
lftab <- SampleMRLifeTable %>%
  kbl(longtable = FALSE, format = "html", caption = CaptionLifeTab, align = "c") %>%
  collapse_rows(valign = "top", latex_hline = "none") %>%
  kable_styling(
    latex_options = c("hold_position", "repeat_header"),
    position = "center", font_size = 8, full_width = FALSE
  ) %>%
  add_header_above(c(" " = 2, "MCMC" = 3, "Classics" = 3))

print(lftab)


# Table -Top Ten Methods --------------------------------------------------

Top10Methods <- base::split(MetricsProstDt, ~ method_group) %>% 
  # ranks methods by groups (MCMC and Classics separately)
  purrr::map(\(x)rank_methods(x, rank.by = "wgt_rel_eff", n = 10,descend = TRUE)) %>% 
  purrr::reduce(rbind)

ResultingTable <- Top10Methods %>%
  mutate(
    sample_size = length(LifeTime$ProstCancer), 
    method = case_when(!(method %in% c("Moments", "MLE", "Regression")) ~ "MCMC", TRUE ~ method)
  ) %>% 
  dplyr::select(sample_size, method, shapePrior, scalePrior, mean_shape, mean_scale, wgt_rel_eff) %>% 
  purrr::modify_at(c("mean_shape", "mean_scale", "wgt_rel_eff"), round, digits = 3L)

# Top 10 MCMC
CAPTION <- stringr::str_glue(
  "The best combination of prior distributions versus classic methods for accurately estimating the Weibull model’s parameters for the Prostrate cancer data. 
  The actual parameter values are {expression(\u03B2)} = {round(mean(RealDtPoolMetrics$mean_shape),2L)} and {expression(\u03B1)} ={round(mean(RealDtPoolMetrics$mean_scale),2L)}, 
  pooled from the results of the 27 fitted models."
)

ResultingTable <- ResultingTable %>%
  draw_table(collapse_cols = 1:2, format = "html", caption = CAPTION) %>%
  add_header_above(c(" " = 2, "Priors" = 2, "Estimated Parameters" = 2, "Model Efficiency" = 1))
# add_footnote(str_glue("The true parameter values are {expression(\u03B2)} = {unique(ApplicationPoolMetrics$TruePar_B)}  and {expression(\u03B1)} ={unique(ApplicationPoolMetrics$TruePar_A)},
# # the average of all parameter estimates from the {length(unique(ApplicationPoolMetrics$PriorComb) )} fitted models."))


print(ResultingTable)





























