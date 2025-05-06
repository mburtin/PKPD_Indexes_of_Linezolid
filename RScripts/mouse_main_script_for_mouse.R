# Record execution time
library(tictoc)
tic()

# Check depencies and load them if necessary
require(mrgsolve)       # Simulations
require(dplyr)          # Data manipulation
require(tidyr)          # Data manipulation
require(purrr)          # Data manipulation
require(tibble)         # Data manipulation
require(future)         # Multithreading
require(future.apply)   # Multithreading

#############################################################################
###########                   MrgSolve Simulations                ###########                       
#############################################################################

#######
##  Generic functions
######
clean_sim_table <- function(df) { 
  df |>
    group_by(ID) |>
    filter(!(time == 0 & duplicated(time))) |> # MrgSolve doesn't allow to filter different EVID output, so we get by default two T0 (one for the dose administration and another for the first observation)
    select(STRN, MIC, DoseGroup, ID, AMT, time, B0, Log10CFU, 
           C_CENTRAL, Cmax_CENTRAL, AUC_CENTRAL, TOVER_MIC_CENTRAL) |>
    ungroup()
}

#######
##  Fractioned simulations
######
simulate_fractioned_dose <- function(i, fractioned_daily_AMT) {
  strain_name <- names(model_list[i])
  
  q12 <- model_list[[i]] |>
    mrgsolve::idata_set(i_data) |>
    mrgsolve::ev(amt = fractioned_daily_AMT / 2, ii = 12, tinf = 0.5, addl = 1) |>
    mrgsolve::mrgsim(delta = 0.1, end = 24) |>
    dplyr::mutate(STRN = strain_name, AMT = fractioned_daily_AMT, DoseGroup = "fractioned_12h") |>
    clean_sim_table()
  
  q6 <- model_list[[i]] |>
    mrgsolve::idata_set(i_data) |>
    mrgsolve::ev(amt = fractioned_daily_AMT / 4, ii = 6, tinf = 0.5, addl = 3) |>
    mrgsolve::mrgsim(delta = 0.1, end = 24) |>
    dplyr::mutate(STRN = strain_name, AMT = fractioned_daily_AMT, DoseGroup = "fractioned_6h") |>
    clean_sim_table()
  
  q3 <- model_list[[i]] |>
    mrgsolve::idata_set(i_data) |>
    mrgsolve::ev(amt = fractioned_daily_AMT / 8, ii = 3, tinf = 0.5, addl = 7) |>
    mrgsolve::mrgsim(delta = 0.1, end = 24) |>
    dplyr::mutate(STRN = strain_name, AMT = fractioned_daily_AMT, DoseGroup = "fractioned_3h") |>
    clean_sim_table()
  
  q1half <- model_list[[i]] |>
    mrgsolve::idata_set(i_data) |>
    mrgsolve::ev(amt = fractioned_daily_AMT / 16, ii = 1.5, tinf = 0.5, addl = 15) |>
    mrgsolve::mrgsim(delta = 0.1, end = 24) |>
    dplyr::mutate(STRN = strain_name, AMT = fractioned_daily_AMT, DoseGroup = "fractioned_1-5h") |>
    clean_sim_table()
  
  return(rbind(q12, q6, q3, q1half) |> tibble::as_tibble())
}

# Print a progress message in the console
message("Simulations start...")

# Start multithreading session
plan(multisession, workers = thread_number)

# Multithreading loop: one thread per strain
fractioned_results <- future_lapply(1:length(model_list), function(i) {
  result_list <- lapply(fractioned_daily_AMT, function(fractioned_daily_AMT_value) {
    simulate_fractioned_dose(i, fractioned_daily_AMT_value)
  })
  bind_rows(result_list)
}, future.seed = TRUE)

# Give the proper strain name to each list element
names(fractioned_results) <- names(model_list)

fractioned_24h_full <- bind_rows(fractioned_results) |> filter(time == 24.0)

#############################################################################
###########                   Data pre-processing                 ###########                       
#############################################################################
# Print a progress message in the console
message("Pre-process data...")

# Merge the data for fractioned doses and constant perfusion
# Transform and pivot also indexes
# Drop unused colunms
sim_formated_data <- bind_rows(fractioned_24h_full) |>
  mutate(ID = row_number()) |>
  group_by(STRN, DoseGroup, AMT, ID) |>
  mutate(
    CENTRAL_Cmax = Cmax_CENTRAL/MIC,
    CENTRAL_AUC_MIC = AUC_CENTRAL/MIC,
    CENTRAL_ToverMIC = TOVER_MIC_CENTRAL/max(time)*100
  ) |>
  pivot_longer(
    cols = c(CENTRAL_Cmax, CENTRAL_AUC_MIC, CENTRAL_ToverMIC),
    names_to = "PKPD_Index",
    values_to = "value"
  ) |>
  select(-c(Cmax_CENTRAL, AUC_CENTRAL, TOVER_MIC_CENTRAL)) |>
  ungroup()

# Split sim_formated_data for each target, and reorder columns
sim_data <- sim_formated_data |>
  group_by(STRN) |>
  mutate(MaxCFU = max(Log10CFU))

#############################################################################
###########          Estimation of correlation parameters         ###########                       
#############################################################################

# From simulated data, we estimate an Emax model for each PKPD indexes. 
# This model will be used to generate the correlation curve. We will also 
# determine R-square and Target attainment for each indexes by DeltaMethods.

source("RScripts/mouse_pkpd_fitting_script.R")

generate_PKPD_fit_data <- function(data) {
  
  correlation_data <- data |>
    unnest(c(data)) |> # Unwrap the main dataset
    select(STRN, MIC, DoseGroup, ID, AMT, time, B0, MaxCFU, PKPD_Index, value, deltaLog10CFU, Rsq) |>
    # Facet_grid order by name, so this fix AUC/MIC displayed before Cmax
    mutate(PKPD_Index = case_when(
      PKPD_Index %in% c("CENTRAL_Cmax", "CSF_Cmax") ~ "I4_Cmax",
      PKPD_Index %in% c("CENTRAL_AUC_MIC", "CSF_AUC_MIC") ~ "I2_AUC",
      PKPD_Index %in% c("CENTRAL_ToverMIC", "CSF_ToverMIC") ~ "I3_ToverMIC",
      TRUE ~ PKPD_Index
    )) |>
    ungroup()
  
  return(correlation_data)
}

# Print progress messages in the console
message("Estimate correlation parameters...")
corrCurve_data <- index_PKPD_fit_curve(sim_data)
sim_data <- generate_PKPD_fit_data(corrCurve_data)


#############################################################################
###########                Generate predicted data                ###########                       
#############################################################################

# Print a progress message in the console
message("Generate predicted data & Obs mean data...")

# These data are used to generate the correlation curve from the estimate parameters precendently.
# Generic function to generate predicted data from the correlation curve
generate_pred_data <- function(data) {
  
  pred_data <- data |>
    select(-c("data", "nls")) |> 
    group_by(STRN, PKPD_Index, I0, Imax, IC50, H)  |>
    expand(value = seq(0.1, 1000, 0.1)) |>
    mutate(pred = Emax_model(value, I0, Imax, IC50, H)) |>
    mutate(PKPD_Index = case_when(
      PKPD_Index %in% c("CENTRAL_Cmax", "CSF_Cmax") ~ "I4_Cmax",
      PKPD_Index %in% c("CENTRAL_AUC_MIC", "CSF_AUC_MIC") ~ "I2_AUC",
      PKPD_Index %in% c("CENTRAL_ToverMIC", "CSF_ToverMIC") ~ "I3_ToverMIC",
      TRUE ~ PKPD_Index
    )) |> filter(!((PKPD_Index %in% c("I3_ToverMIC")) & value > 100)) |>
    ungroup()
  
  return(pred_data)
}

pred_data <- generate_pred_data(corrCurve_data)

#############################################################################
########                    Environment cleaning                     ########                       
#############################################################################
# Reorganize all data in lists
sim_results_fractioned <- list(fractioned_results = fractioned_results,
                               fractioned_24h_full = fractioned_24h_full)
data <- list(
  sim_data = sim_data,
  correlation_data = corrCurve_data,
  pred_data = pred_data
)

mrgsolve_model <- list(model_file = model_file, model_list = model_list)

# Stop recording time execution
elapsed <- toc(quiet = TRUE)
message("Total execution time: ", round(elapsed$toc - elapsed$tic, 0), " seconds")

# Clean the environment
rm(fractioned_results, fractioned_24h_full, 
   i_data, 
   sim_data,
   pred_data,
   corrCurve_data,
   model_file, model_list, elapsed)

gc()