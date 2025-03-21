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
    select(STRN, MIC, DoseGroup, ID, AMT, time, B0, Log10CFU_CSF, 
           C_CENTRAL, Cmax_CENTRAL, AUCCENTRAL, TOVER_MIC_CENTRAL,
           C_CSF, Cmax_CSF, AUCCSF, TOVER_MIC_CSF) |>
    ungroup()
}

#######
##  Fractioned simulations
######
simulate_fractioned_dose <- function(i, fractioned_daily_AMT) {
  strain_name <- names(model_list[i])
  
  q24 <- model_list[[i]] |>
    mrgsolve::idata_set(i_data) |>
    mrgsolve::ev(amt = fractioned_daily_AMT, tinf = 0.5, ii = 24, addl = 0) |>
    mrgsolve::mrgsim(delta = 1, end = 24) |>
    dplyr::mutate(STRN = strain_name, AMT = fractioned_daily_AMT, DoseGroup = "fractioned_24h") |>
    clean_sim_table()
  
  q12 <- model_list[[i]] |>
    mrgsolve::idata_set(i_data) |>
    mrgsolve::ev(amt = fractioned_daily_AMT / 2, tinf = 0.5, ii = 12, addl = 1) |>
    mrgsolve::mrgsim(delta = 1, end = 24) |>
    dplyr::mutate(STRN = strain_name, AMT = fractioned_daily_AMT, DoseGroup = "fractioned_12h") |>
    clean_sim_table()
  
  q8 <- model_list[[i]] |>
    mrgsolve::idata_set(i_data) |>
    mrgsolve::ev(amt = fractioned_daily_AMT / 3, tinf = 0.5, ii = 8, addl = 2) |>
    mrgsolve::mrgsim(delta = 1, end = 24) |>
    dplyr::mutate(STRN = strain_name, AMT = fractioned_daily_AMT, DoseGroup = "fractioned_8h") |>
    clean_sim_table()
  
  q4 <- model_list[[i]] |>
    mrgsolve::idata_set(i_data) |>
    mrgsolve::ev(amt = fractioned_daily_AMT / 6, tinf = 0.5, ii = 4, addl = 5) |>
    mrgsolve::mrgsim(delta = 1, end = 24) |>
    dplyr::mutate(STRN = strain_name, AMT = fractioned_daily_AMT, DoseGroup = "fractioned_4h") |>
    clean_sim_table()
  
  return(rbind(q24, q12, q8, q4) |> tibble::as_tibble())
}

# Print a progress message in the console
message("Simulations start:")
message(" - Fractioned doses... (1/3)")

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
#######
##  Continuous infusion simulations
######
# Print a progress message in the console
message(" - Continuous infusions... (2/3)")

# Initialize an empty dataframe for continuous infusion
continuous_inf_results <- tibble()

# Recover the average concentrations of fractioned regimens to calculate the infusion dose and loading dose
continuous_inf_amt <- fractioned_24h_full |>
  group_by(AMT) |>
  summarize(AMT_Inf = mean(AUCCENTRAL)*(CL/FU),
            AMT_Loading = mean(AUCCENTRAL)/24*V1/FU,
            Conc_mean = mean(AUCCENTRAL)/24) |>
  select(-c(AMT))

# Create event table for both doses
create_dosing_events <- function(id, loading_amt, inf_amt) {
  loading_dose <- ev(time = 0, 
                     amt = loading_amt,
                     tinf = 0,
                     cmt = 1,
                     ID = id)
  
  infusion_dose <- ev(time = 0,
                      amt = inf_amt,
                      tinf = 24,
                      cmt = 1,
                      ID = id)
  
  return(c(loading_dose, infusion_dose))
}

# Loop over all the models to simulate the continuous infusions
continuous_inf_results <- future_lapply(1:length(model_list), function(i) {
  
  # Initialize an empty list to store results
  model_results <- list()
  
  # Loop over all the continuous infusion doses
  for (j in 1:nrow(continuous_inf_amt)) {
    
    # Create dosing events
    inf_dose <- mapply(function(id) {
      create_dosing_events(
        id = id,
        loading_amt = continuous_inf_amt$AMT_Loading[j],
        inf_amt = continuous_inf_amt$AMT_Inf[j]
      )
    }, 1:nrow(i_data), SIMPLIFY = FALSE)
    
    inf_dose <- do.call(c, inf_dose)
    
    # MrgSolve simulations
    continuous_inf <- model_list[[i]] |>
      mrgsolve::data_set(inf_dose) |>
      mrgsolve::mrgsim(delta = 1, end = 24) |>
      mutate(STRN = names(model_list)[i], 
             AMT = round(continuous_inf_amt$AMT_Inf[j], 2), 
             DoseGroup = "continuous_infusion") |>
      clean_sim_table()
    
    # Add results to the list
    model_results[[j]] <- continuous_inf
  }
  
  # Return result in one dataframe
  return(bind_rows(model_results))
  
}, future.seed = TRUE)

# Combine all continuous infusion dataframes in one
continuous_inf_results <- bind_rows(continuous_inf_results)

# Keep only values at 24h
continuous_inf_24h <- continuous_inf_results |> filter(time == 24.0)

#######
##  Control simulations
######
# Print a progress message in the console
message(" - Control... (3/3)")

# Loop over all the models to simulate
control_results_list <- future_lapply(1:length(model_list), function(i) {
  
  # MrgSolve simulations
  sim_result <- model_list[[i]] |>
    mrgsolve::idata_set(i_data) |>
    mrgsolve::ev(amt = 0, tinf = 0.5, ii = 24, addl = 0) |>
    mrgsolve::mrgsim(delta = 1, end = 24) |>
    dplyr::mutate(STRN = names(model_list[i]), AMT = 0, DoseGroup = "control") |>
    clean_sim_table()
  
  # Return result in one dataframe
  return(sim_result)
  
}, future.seed = TRUE)

# Combine all continuous infusion dataframes in one
control_results <- bind_rows(control_results_list) |> tibble::as_tibble()

# Keep only values at 24h
control_24h <- control_results |> filter(time == 24.0)

# Stop multithreading session
plan(sequential)

#############################################################################
###########                   Data pre-processing                 ###########                       
#############################################################################
# Print a progress message in the console
message("Pre-process data...")

# Merge the data for fractioned doses and constant perfusion
# Transform and pivot also indexes
# Drop unused colunms
sim_formated_data <- bind_rows(fractioned_24h_full, continuous_inf_24h, control_24h) |>
  mutate(ID = row_number()) |>
  group_by(STRN, DoseGroup, AMT, ID) |>
  mutate(
    CENTRAL_Cmax = Cmax_CENTRAL/MIC,
    CENTRAL_AUC_MIC = AUCCENTRAL/MIC,
    CENTRAL_ToverMIC = TOVER_MIC_CENTRAL/max(time)*100,
    CSF_Cmax = Cmax_CSF/MIC,
    CSF_AUC_MIC = AUCCSF/MIC,
    CSF_ToverMIC = TOVER_MIC_CSF/max(time)*100
  ) |>
  pivot_longer(
    cols = c(CENTRAL_Cmax, CENTRAL_AUC_MIC, CENTRAL_ToverMIC,
             CSF_Cmax, CSF_AUC_MIC, CSF_ToverMIC),
    names_to = "PKPD_Index",
    values_to = "value"
  ) |>
  select(-c(Cmax_CENTRAL, AUCCENTRAL, TOVER_MIC_CENTRAL,
            Cmax_CSF, AUCCSF, TOVER_MIC_CSF)) |>
  ungroup()

# Split sim_formated_data for each target, and reorder columns
csf_sim_data <- sim_formated_data |> rename("Log10CFU" = Log10CFU_CSF) |> filter(PKPD_Index %in% c("CSF_Cmax", "CSF_AUC_MIC", "CSF_ToverMIC"))
mix_sim_data <- sim_formated_data |> rename("Log10CFU" = Log10CFU_CSF) |> filter(PKPD_Index %in% c("CENTRAL_Cmax", "CENTRAL_AUC_MIC", "CENTRAL_ToverMIC"))

#############################################################################
###########          Estimation of correlation parameters         ###########                       
#############################################################################

# From simulated data, we estimate an Emax model for each PKPD indexes. 
# This model will be used to generate the correlation curve. We will also 
# determine R-square and Target attainment for each indexes by DeltaMethods.

source("RScripts/pkpd_fitting_script.R")

generate_PKPD_fit_data <- function(data) {
  
  correlation_data <- data |>
    unnest(c(data)) |> # Unwrap the main dataset
    select(STRN, MIC, DoseGroup, ID, AMT, time, B0, PKPD_Index, value, Log10CFU, deltaLog10CFU, Rsq) |>
    filter(ID > 0) |>
    # Facet_grid order by name, so this fix AUC/MIC displayed before Cmax
    mutate(PKPD_Index = case_when(
      PKPD_Index %in% c("CENTRAL_Cmax", "CSF_Cmax") ~ "I1_Cmax",
      PKPD_Index %in% c("CENTRAL_AUC_MIC", "CSF_AUC_MIC") ~ "I2_AUC",
      PKPD_Index %in% c("CENTRAL_ToverMIC", "CSF_ToverMIC") ~ "I3_ToverMIC",
      TRUE ~ PKPD_Index
    )) |>
    ungroup()
  
  return(correlation_data)
}

# Print progress messages in the console
message("Estimate correlation parameters:")

message(" - CSF... (1/2)")
csf_corrCurve_data <- index_PKPD_fit_curve(csf_sim_data)
csf_sim_data <- generate_PKPD_fit_data(csf_corrCurve_data)

message(" - CSF with Plasma indexes... (2/2)")
mix_corrCurve_data <- index_PKPD_fit_curve(mix_sim_data)
mix_sim_data <- generate_PKPD_fit_data(mix_corrCurve_data)

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
    group_by(PKPD_Index, I0, Imax, IC50, H)  |>
    expand(value = seq(0.1, 1000, 0.1)) |>
    mutate(pred = Emax_model(value, I0, Imax, IC50, H)) |>
    mutate(PKPD_Index = case_when(
      PKPD_Index %in% c("CENTRAL_Cmax", "CSF_Cmax") ~ "I1_Cmax",
      PKPD_Index %in% c("CENTRAL_AUC_MIC", "CSF_AUC_MIC") ~ "I2_AUC",
      PKPD_Index %in% c("CENTRAL_ToverMIC", "CSF_ToverMIC") ~ "I3_ToverMIC",
      TRUE ~ PKPD_Index
    )) |> filter(!((PKPD_Index %in% c("I3_ToverMIC")) & value > 100)) |>
    ungroup()
  
  return(pred_data)
}

csf_pred_data <- generate_pred_data(csf_corrCurve_data)
mix_pred_data <- generate_pred_data(mix_corrCurve_data)

#############################################################################
########          Generate a dataset with observation mean           ########                       
#############################################################################

csf_obs_mean <- csf_sim_data |>
  group_by(STRN, DoseGroup, AMT, PKPD_Index) |>
  summarize( 
    value = mean(value),
    deltaLog10CFU = mean(deltaLog10CFU),
    Rsq = unique(Rsq),
  )

mix_obs_mean <- mix_sim_data |>
  group_by(STRN, DoseGroup, AMT, PKPD_Index) |>
  summarize( 
    value = mean(value),
    deltaLog10CFU = mean(deltaLog10CFU),
    Rsq = unique(Rsq),
  )

#############################################################################
########                    Environment cleaning                     ########                       
#############################################################################
# Reorganize all data in lists
sim_results_fractioned <- list(fractioned_results = fractioned_results,
                               fractioned_24h_full = fractioned_24h_full)

sim_results_continuous_inf <- list(continuous_inf_results = continuous_inf_results, 
                                   continuous_inf_24h = continuous_inf_24h, 
                                   continuous_inf_amt = continuous_inf_amt)

sim_results_control <- list(control_results = control_results, 
                            control_24h = control_24h)

csf_data <- list(
  sim_data = csf_sim_data,
  correlation_data = csf_corrCurve_data,
  pred_data = csf_pred_data,
  obs_mean = csf_obs_mean
)

csf_with_plasma_index_data <- list(
  sim_data = mix_sim_data,
  correlation_data = mix_corrCurve_data,
  pred_data = mix_pred_data,
  obs_mean = mix_obs_mean
)

mrgsolve_model <- list(model_file = model_file, model_list = model_list)

# Stop recording time execution
elapsed <- toc(quiet = TRUE)
message("Total execution time: ", round(elapsed$toc - elapsed$tic, 0), " seconds")

# Clean the environment
rm(fractioned_results, fractioned_24h_full, 
   continuous_inf_results, continuous_inf_24h, continuous_inf_amt, 
   control_results, control_24h, i_data, 
   csf_sim_data, mix_sim_data,
   csf_pred_data, mix_pred_data,
   csf_obs_mean, mix_obs_mean,
   csf_corrCurve_data, mix_corrCurve_data,
   model_file, model_list, elapsed)

gc()