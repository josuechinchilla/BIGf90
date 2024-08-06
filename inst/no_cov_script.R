#Set seed, working directory and load packages needed
set.seed(101919)
library("BIGf90")
setwd("~/Desktop/weight_length_2024/weight_2022phenos_2022ped/no_cov")
#set path to executables
path_2_execs = system.file("exec_windows",package = "BIGf90")

####    1) Process raw parameter file with renumf90    ####
run_renum(path_2_execs, 
          output_files_dir = "results", 
          raw_par_files = "data/weight_2022_no_cov.par")

####    2) Estimate  variance components with gibbsf90+ and postgibbsf90    ####

run_gibbs(path_2_execs, 
          input_files_dir = "results",
          gibbs_iter = 250000,
          gibbs_burn = 50000,
          gibbs_keep = 1)

run_postgibbs(path_2_execs, input_files_dir = "results",
              postgibbs_burn = 1,
              postgibbs_keep = 100)

####    3) Run a cross-validation analysis (CVA)  to estimate the accuracy of your EBV predictions  ####

bf90_cv(
  path_2_execs,
  missing_value_code = -999,
  random_effect_col = 3,
  h2 = 0.95,
  num_runs = 10,
  num_folds = 5,
  renf90_ped_name = "results/renadd03.ped",
  input_files_dir = "results",
  output_files_dir = "bf90_cv_results/",
  output_table_name = "CVA_results_weight_2022phenos_2022ped_no_cov.txt"
  )

# Check in your version if the files generated have the same results below (with some tolerance)
yhat <- read.table("bf90_cv_results/yhat_residual")
dim(yhat)
sum(yhat$V2) # 89.84

cv_results <- read.table("bf90_cv_results/CVA_results_weight_2022phenos_2022ped_no_cov.txt", skip = 1, nrows = 30)
sum(cv_results$V4) # 25.746
