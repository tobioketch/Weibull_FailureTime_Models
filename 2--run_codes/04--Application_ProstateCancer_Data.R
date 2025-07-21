
# if(!require("NSM3"))
#   install.packages("NSM3")

# Prostate Cancer Patients
ProstCancer <- c(
  rep(0, 3), 2, 3, 4, 6, rep(7, 2), 8, rep(9, 2), rep(11, 3), rep(12, 3), rep(15, 2),
  rep(16, 3), rep(17, 2), 18, rep(19, 2), 20, 21, rep(22, 2), 23, 24, rep(25, 2), rep(26, 3),
  rep(27, 2), rep(28, 2), rep(29, 2), 30, 31, rep(32, 3), rep(33, 2), 34, 35, 36, rep(37, 2),
  38, 40, rep(41, 2), rep(42, 2), 43, rep(45, 3), 46, rep(47, 2), rep(48, 2), 51, rep(53, 2),
  rep(54, 2), 57, 60, 61, rep(62, 2), 67, 69, 87, rep(97, 2), 100, 145, 158
)

# t = 0, Taken as a Diagnosis Date
ProstCancer <- ProstCancer[ProstCancer > 0]

# MLE Results Suggests Increasing Failure Rate Prostrate Cancer Data
# Large Sample Approximation Test; n > 9

# # H0: Constant Failure Rate
# # H1: Increasing Failure Rate
# EpsteinTest <- NSM3::epstein(x = ProstCancer, alternative = "ifr")
# base::rbind(EpsteinTest)


this.path::withArgs(
  # Calls "SolverProstateCancerDt.R" script to fit MCMC and Classic Models. 
  # Returns
  #     1. Pooled/ combined results of the fitted models that's 
  #        saved in the path file = "misc/app_misc_results/" as RealDtPoolMetrics.RData
  #     2. The saved data set is later used to plot graphs and tables for visual results
  
  source("2--run_codes/NOT_RUN_RCODES/SolverProstateCancerDt.R", verbose = FALSE),
  # ###################### Args3 (currently not implemented)
  # Arg1: Data,     
  # Arg2 : a name used for Subtitle & PlotFileName
  # # Arg3 : Top n (10) MCMC models(Tabulated Results)
  #        : We plot top 3 of the MCMC models vs. Classical methods (MLE, OLS, & MOM)
  ProstCancer, "Prostate Cancer Data", 10L
)

 







