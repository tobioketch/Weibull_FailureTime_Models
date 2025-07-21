# Data --------------------------------------------------------------------
base::rm( list = ls( all.names = TRUE )  )

source(file = "PaperBasedOnTheseResults/MethodsAndPackages.R")

# Small Sample Size
base::load(file = "PaperBasedOnTheseResults/misc_PaperBasedOnTheseResults/sim_misc_results/ResultsSimMethods_SmallSmpls")

# Large Sample Size
base::load(file = "PaperBasedOnTheseResults/misc_PaperBasedOnTheseResults/sim_misc_results/ResultsSimMethods_LargeSmpls")



# Add Plots for Large Sample Sizes & Format Plots and Data for DHR to start from 0.4

ResultsSimMethods <- dplyr::bind_rows(ResultsSimMethods_SmallSmpls, ResultsSimMethods_LargeSmpls) %>%
  dplyr::filter(shape > 0.3) %>% 
  dplyr::arrange(n) 

# Efficiency Trend Plot ---------------------------------------------------

xlab <- paste0("Shape (", expression("\u03B2"),") Parameter")
ylab <- paste0("Weighted Relative Efficiency, WRE(", expression("\u03B8"),")")

# Plot Data, n = 15
PltDt <- ResultsSimMethods %>%
  as_analysis_sim_df(table = FALSE) %>%
  dplyr::mutate(n = factor(paste("n", n, sep = ", "), levels = paste("n", sort(unique(n)), sep = ", ")) )
      
# Plot
Plt <- PltDt %>%
  plot_lines(
    x = shape, y = rel_wgt_eff, group = prior_comb, facet.by = n,
    lab.x = xlab, lab.y = ylab, title = "Model Efficiency"
  )

Plt

save_jpeg(
  plot = Plt,
  file_name = "Fig2",
  path = "PaperBasedOnTheseResults/demoPlots_PaperBasedonThesePlots/"
)

# Excerpt Tables for Efficiency Trends ------------------------------------

PltDtSectionTab <- PltDt %>%
  select(n, shape, method, prior_comb, rel_wgt_eff) %>%
  mutate(prior_comb = case_when(prior_comb %in% c("Regression", "Moments", "MLE") ~ "-", TRUE ~ prior_comb),
         n = as.character(n), n = readr::parse_number(n), rel_wgt_eff = round(rel_wgt_eff, digits = 3L)) %>%
  separate(col = "prior_comb", into = c("shape_prior", "scale_prior"), fill = "right", sep = "-") %>%
  mutate(across(c("shape_prior", "scale_prior"), \(x) na_if(x, "") %>% replace_na("-")))  %>%
  group_by(shape) %>%
  arrange(desc(rel_wgt_eff), .by_group = TRUE) %>%
  ungroup()

# Table excerpt (n, 15) for different shape values (0.4, 0.5, 0.6)

PltDtSectionTab %>% 
filter(n==15, shape >= 0.4, shape <= 0.6) %>% 
  draw_table(collapse_cols = 1:3, format = "html", 
             caption = "A section of results showing weighted relative efficiency for a small sample size of 15 units.") %>% 
  add_header_above(c(" " = 3, "Prior Distribution" = 2, "Model Efficiency" = 1))

# Table excerpt (n, 100) for different shape values (1.4, 1.5, 1.6)

PltDtSectionTab %>% 
filter(n==100, shape >1.3, shape < 1.8) %>% 
  draw_table(collapse_cols = 1:3, format = "html", 
             caption = "A section of results showing weighted relative efficiency for a large sample size of 100 units.") %>% 
  add_header_above(c(" " = 3, "Prior Distribution" = 2, "Model Efficiency" = 1))


# Average weighted efficiency tables --------------------------------------

AvgWeightedREff <- ResultsSimMethods %>%
  as_analysis_sim_df(table = TRUE, top.n = 25L) %>%
  tidyr::separate(col = "prior_comb", into = c("shape_prior", "scale_prior"), fill = "right", sep = "-") %>%
  dplyr::mutate(
    method = dplyr::case_when(!(shape_prior %in% c("Moments", "MLE", "Regression")) ~ "MCMC", TRUE ~ shape_prior),
    shape_prior = dplyr::case_when((shape_prior %in% c("Moments", "MLE", "Regression")) ~ "-", TRUE ~ shape_prior),
    scale_prior = case_when(is.na(scale_prior) ~ "-", TRUE ~ scale_prior)
  ) %>%
  dplyr::relocate(method, .after = "hazard") %>% 
  dplyr::rename(AWRE = "mean_rel_wgt_eff")

# Sample Size=15, 25
Avg_Eff_SmallSmpl <-   AvgWeightedREff %>% 
  dplyr::filter(n < 55) 

# Small Samples Avg Weighted Relative Efficiency Table
Avg_Eff_SmallSmpl %>% 
   draw_table(collapse_cols = 1:3, font_size = 12, format = "html", caption = "Average weighted relative efficiency small samples") %>% 
  add_header_above(c(" " = 3, "Prior Distribution" = 2, "Model Efficiency" = 1))

# Sample Size=55, 100
Avg_Eff_LargeSmpl <-   AvgWeightedREff %>% 
  dplyr::filter(n >= 55) 

# Large Samples Avg Weighted Relative Efficiency Table
Avg_Eff_LargeSmpl %>% 
  draw_table(collapse_cols = 1:3, font_size = 12, format = "html", caption = "Average weighted relative efficiency large samples") %>% 
  add_header_above(c(" " = 3, "Prior Distribution" = 2, "Model Efficiency" = 1))

# Conclusion and Discussion
######################################################################################

# Small Samples Recommended Prior
RecommendedSmallSamples <- Avg_Eff_SmallSmpl %>% 
  dplyr::filter(method=="MCMC") %>% 
  dplyr::group_by(hazard, method, shape_prior, scale_prior) %>% 
  dplyr::mutate(AWRE= mean(AWRE), .keep = "unused") %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-c(n, method) ) %>% 
  base::unique() %>% 
  dplyr::arrange(hazard, desc(AWRE)) %>% 
  dplyr::mutate(n = "n = (15, 25)", .before = "hazard")
# Large Samples Recommended Prior
RecommendedLargeSamples <- Avg_Eff_LargeSmpl %>% 
  dplyr::filter(method=="MCMC") %>% 
  dplyr::group_by(hazard, method, shape_prior, scale_prior) %>% 
  dplyr::mutate(AWRE= mean(AWRE), .keep = "unused") %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-c(n, method)) %>% 
  base::unique() %>% 
  dplyr::arrange(hazard, desc(AWRE)) %>% 
  dplyr::mutate(n = "n = (55, 100)", .before = "hazard")

RecommendedMethods <- dplyr::bind_rows(RecommendedSmallSamples, RecommendedLargeSamples)

RecommendedMethods %>% 
  dplyr::filter(hazard=="DHR") %>% 
  draw_table(collapse_cols = c(1:2, 4), format = "html", 
  caption = "A list of recommended combinations of prior distributions for small and large
Weibull-distributed data sets with decreasing hazard rate properties. The priors are ranked based
on the AWRE.")

RecommendedMethods %>% 
  dplyr::filter(hazard=="IHR") %>% 
  draw_table(collapse_cols = c(1:2, 4), format = "html", 
  caption = "A list of recommended combinations of prior distributions for small and large
Weibull-distributed data sets with increasing hazard rate properties. The priors are ranked based on
the AWRE.")


# Constructed from prior methods with higher Avg_WREff over I/DHR for small and large datasets
RecomendPriorsDHR <- data.frame(
  shape.nsmall = c(
    "HalfNormal",
    "LogNormal",
    "Gamma",
    "HalfNormal",
    "Gamma",
    "LogNormal"
  ),
  scale.nsmall = c(rep("HalfCauchy", 3), "LogNormal", "Gamma", "Gamma"),
  shape.nlarge = c(
    "Exponential",
    "HalfNormal",
    "Gamma",
    "LogNormal",
    "LogNormal",
    "HalfNormal"
  ),
  scale.nlarge = c(rep("HalfCauchy", 4), rep("Gamma", 2))
)

RecomendPriorsDHR %>%
  draw_table(
    collapse_cols = c(2, 4), format = "html",
    caption = "A list of recommended combinations of prior distributions\n
    for small and large Weibull-distributed data sets with decreasing hazard rate properties."
  ) %>%
  add_header_above(c("n = c(15, 25)" = 2, "n = c(55, 100)" = 2))


RecomendPriorsIHR <- data.frame(
  shape.nsmall = c("Gamma", "HalfNormal", "LogNormal", "-"),
  scale.nsmall = c(rep("HalfCauchy", 3), "-"),
  shape.nlarge = c("Exponential", "Lognormal", "Gamma", "HalfNormal"),
  scale.nlarge = rep("HalfCauchy", 4)
)


RecomendPriorsIHR %>%
  draw_table(
    collapse_cols = c(2, 4), format= "html",
    caption = "A list of recommended combinations of prior \n
    distributions for small and large Weibull-distributed data sets with increasing hazard rate properties."
  ) %>%
  add_header_above(c("n = c(15, 25)" = 2, "n = c(55, 100)" = 2))







