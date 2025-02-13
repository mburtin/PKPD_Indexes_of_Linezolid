# Record execution time using tictoc for more precise timing
library(tictoc)
tic()

# Check depencies
require(mrgsolve)       # Simulation package
require(dplyr)          # Data manipulation
require(tidyr)          # Data manipulation
require(purrr)          # Data manipulation
require(tibble)         # Data manipulation
require(future)         # Multithreading
require(future.apply)   # Multithreading
require(ggplot2)        # Plotting
require(ggh4x)          # Plotting (allow to set custom scales for each facet)

#############################################################################
###########                   MrgSolve Simulations                ###########                       
#############################################################################

#######
##  Generic functions
######
clean_sim_table <- function(df) { 
  df |>
    group_by(ID) |>
    filter(!(time == 0 & duplicated(time))) |>
    select(STRN, MIC, DoseGroup, ID, AMT, time, B0, Log10CFU_CENTRAL, Log10CFU_CSF, 
           Cmax_CENTRAL, AUCCENTRAL, TOVER_MIC_CENTRAL, TOVER_4MIC_CENTRAL, TOVER_10MIC_CENTRAL, 
           Cmax_CSF, AUCCSF, TOVER_MIC_CSF, TOVER_4MIC_CSF, TOVER_10MIC_CSF, C_CENTRAL, C_CSF) |>
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

# Multithreading loop, one threading do all simulations for one strain
fractioned_results <- future_lapply(1:length(model_list), function(i) {
  result_list <- lapply(fractioned_daily_AMT, function(fractioned_daily_AMT_value) {
    simulate_fractioned_dose(i, fractioned_daily_AMT_value)
  })
  # Combine all results
  do.call(rbind, result_list)
}, future.seed = TRUE)

# Given the strain name to each list elements
names(fractioned_results) <- names(model_list)

# Filter the results to keep only 24h values
fractioned_24h <- future_lapply(fractioned_results, function(df) {
  df[df$time == 24.0, ]
})

# Stop multithreading session
plan(sequential)

# Combine all fractioned results in one dataframe
fractioned_24h_full <- do.call(rbind, fractioned_24h)

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
  mutate(DoseGroup = 0) |>
  select(AMT_Inf, AMT_Loading, Conc_mean, DoseGroup)

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


plan(multisession, workers = thread_number)

continuous_inf_results <- future_lapply(1:length(model_list), function(i) {
  
  # Initialisation de la liste pour stocker les résultats de chaque modèle
  model_results <- list()
  
  # Boucle interne sur les doses dans continuous_inf_amt
  for (j in 1:nrow(continuous_inf_amt)) {
    
    # Créer les événements de dosage en une seule étape sans do.call + lapply
    inf_dose <- mapply(function(id) {
      create_dosing_events(
        id = id,
        loading_amt = continuous_inf_amt$AMT_Loading[j],
        inf_amt = continuous_inf_amt$AMT_Inf[j]
      )
    }, 1:nrow(i_data), SIMPLIFY = FALSE)
    
    inf_dose <- do.call(c, inf_dose)
    
    # Simuler avec mrgsolve
    continuous_inf <- model_list[[i]] |>
      mrgsolve::data_set(inf_dose) |>
      mrgsolve::mrgsim(delta = 1, end = 24) |>
      mutate(STRN = names(model_list)[i], 
             AMT = round(continuous_inf_amt$AMT_Inf[j], 2), 
             DoseGroup = "continuous_infusion") |>
      clean_sim_table()
    
    # Ajouter les résultats à la liste
    model_results[[j]] <- continuous_inf
  }
  
  # Retourner les résultats du modèle en une seule data frame
  return(bind_rows(model_results))
}, future.seed = TRUE)

plan(sequential)

# Combiner tous les résultats des différents modèles en un seul data frame
continuous_inf_results <- bind_rows(continuous_inf_results)
  
continuous_inf_24h <- continuous_inf_results |> filter(time == 24.0)

#######
##  Control simulations
######
# Print a progress message in the console
message(" - Control... (3/3)")
## Simulate a control with AMT = 0
control_results <- do.call(rbind, lapply(1:length(model_list), function(i) {
  model_list[[i]] |>
    mrgsolve::idata_set(i_data) |>
    mrgsolve::ev(amt = 0, tinf = 0.5, ii = 24, addl = 0) |>
    mrgsolve::mrgsim(delta = 1, end = 24) |>
    dplyr::mutate(STRN = names(model_list[i]), AMT = 0, DoseGroup = "control") |>
    clean_sim_table()
})) |>
  tibble::as_tibble()

control_24h <- control_results |> filter(time == 24.0)

#############################################################################
###########                   Data pre-processing                 ###########                       
#############################################################################
# Print a progress message in the console
message("Pre-process data...")

# Merge the data for fractioned doses and constant perfusion
# Transform and pivot also indexes
# Drop unused colunms
sim_formated_data <- rbind(fractioned_24h_full, continuous_inf_24h, control_24h) |>
  select(-c(ID)) |>
  mutate(ID = row_number()) |>
  group_by(STRN, DoseGroup, AMT, ID) |>
  mutate(
    CENTRAL_Cmax = Cmax_CENTRAL/MIC,
    CENTRAL_AUC_MIC = AUCCENTRAL/MIC,
    CENTRAL_ToverMIC = TOVER_MIC_CENTRAL/max(time)*100,
    CENTRAL_ToverMIC_4X = TOVER_4MIC_CENTRAL/max(time)*100,
    CENTRAL_ToverMIC_10X = TOVER_10MIC_CENTRAL/max(time)*100,
    CSF_Cmax = Cmax_CSF/MIC,
    CSF_AUC_MIC = AUCCSF/MIC,
    CSF_ToverMIC = TOVER_MIC_CSF/max(time)*100,
    CSF_ToverMIC_4X = TOVER_4MIC_CSF/max(time)*100,
    CSF_ToverMIC_10X = TOVER_10MIC_CSF/max(time)*100
  ) |>
  pivot_longer(
    cols = c(CENTRAL_Cmax, CENTRAL_AUC_MIC, CENTRAL_ToverMIC, CENTRAL_ToverMIC_4X, CENTRAL_ToverMIC_10X,
             CSF_Cmax, CSF_AUC_MIC, CSF_ToverMIC, CSF_ToverMIC_4X, CSF_ToverMIC_10X),
    names_to = "PKPD_Index",
    values_to = "value"
  ) |>
  select(-c(Cmax_CENTRAL, AUCCENTRAL, TOVER_MIC_CENTRAL, TOVER_4MIC_CENTRAL, TOVER_10MIC_CENTRAL,
            Cmax_CSF, AUCCSF, TOVER_MIC_CSF, TOVER_4MIC_CSF, TOVER_10MIC_CSF)) |>
  ungroup()

# Split sim_formated_data for each target, and reorder columns
csf_sim_data <- sim_formated_data |> rename("Log10CFU" = Log10CFU_CSF) |> filter(PKPD_Index %in% c("CSF_Cmax", "CSF_AUC_MIC", "CSF_ToverMIC", "CSF_ToverMIC_4X", "CSF_ToverMIC_10X"))

plasma_sim_data <- sim_formated_data |> rename("Log10CFU" = Log10CFU_CENTRAL) |> filter(PKPD_Index %in% c("CENTRAL_Cmax", "CENTRAL_AUC_MIC", "CENTRAL_ToverMIC", "CENTRAL_ToverMIC_4X", "CENTRAL_ToverMIC_10X"))

mix_sim_data <- sim_formated_data |> rename("Log10CFU" = Log10CFU_CSF) |> filter(PKPD_Index %in% c("CENTRAL_Cmax", "CENTRAL_AUC_MIC", "CENTRAL_ToverMIC", "CENTRAL_ToverMIC_4X", "CENTRAL_ToverMIC_10X"))

#############################################################################
###########          Estimation of correlation parameters         ###########                       
#############################################################################

# From simulated data, we estimate an Emax model for each PKPD indexes. 
# This model will be used to generate the correlation curve. We will also 
# determine R-square and Target attainment for each indexes by DeltaMethods.

source("RScripts/pkpd_fitting_script.R")

generate_PKPD_fit_data <- function(corrData) {
  
  correlation_data <- corrData |>
    unnest(c(data)) |> # Unwrap the main dataset
    select(STRN, MIC, DoseGroup, ID, AMT, time, B0, PKPD_Index, value, Log10CFU, deltaLog10CFU, Rsq, Adj.Rsq) |>
    filter(ID > 0) |>
    # Facet_grid order by name, so this fix AUC/MIC displayed before Cmax
    mutate(PKPD_Index = case_when(
      PKPD_Index %in% c("CENTRAL_Cmax", "CSF_Cmax") ~ "I1_Cmax",
      PKPD_Index %in% c("CENTRAL_AUC_MIC", "CSF_AUC_MIC") ~ "I2_AUC",
      PKPD_Index %in% c("CENTRAL_ToverMIC", "CSF_ToverMIC") ~ "I3_ToverMIC",
      PKPD_Index %in% c("CENTRAL_ToverMIC_4X", "CSF_ToverMIC_4X") ~ "I4_ToverMIC_4X",
      PKPD_Index %in% c("CENTRAL_ToverMIC_10X", "CSF_ToverMIC_10X") ~ "I5_ToverMIC_10X",
      TRUE ~ PKPD_Index
    ))
  
  return(correlation_data)
}

# Print progress messages in the console
message("Estimate correlation parameters:")

message(" - Plasma... (1/3)")
plasma_corrCurve_data <- index_PKPD_fit_curve(plasma_sim_data)
plasma_sim_data <- generate_PKPD_fit_data(plasma_corrCurve_data)

message(" - CSF... (2/3)")
csf_corrCurve_data <- index_PKPD_fit_curve(csf_sim_data)
csf_sim_data <- generate_PKPD_fit_data(csf_corrCurve_data)

message(" - CSF with Plasma indexes... (3/3)")
mix_corrCurve_data <- index_PKPD_fit_curve(mix_sim_data)
mix_sim_data <- generate_PKPD_fit_data(mix_corrCurve_data)

#############################################################################
###########                Generate predicted data                ###########                       
#############################################################################

# Print a progress message in the console
message("Generate predicted data & Obs mean data...")

# These data are used to generate the correlation curve from the estimate parameters precendently.
# Generic function to generate predicted data from the correlation curve
generate_pred_data <- function(data, target) {
  
  # Filter by STRN, recover B0 for each strain, and expand the dataset to have a value from 0.1 to 1000 
  # and calculate the pred values
  pred_data <- data |>
    select(-c("data", "nls")) |> 
    group_by(PKPD_Index, I0, Imax, IC50, H)  |>
    expand(value = seq(0.1, 1000, 0.1)) |>
    mutate(pred = Imax_model(value, I0, Imax, IC50, H))
  
  # Depending of the target, we need to filter the data outside of x scaling in plot 
  # and align indexes names with sim_data
  if(target == "CSF") {
    pred_data <- pred_data |>
      filter(!((PKPD_Index %in% c("CSF_ToverMIC", "CSF_ToverMIC_4X", "CSF_ToverMIC_10X")) & value > 100)) |>
      # Rename indexes colunms to be align with csf_pkpd_data
      mutate(PKPD_Index = case_when(
        PKPD_Index == "CSF_Cmax" ~ "I1_Cmax",
        PKPD_Index == "CSF_AUC_MIC" ~ "I2_AUC",
        PKPD_Index == "CSF_ToverMIC" ~ "I3_ToverMIC",
        PKPD_Index == "CSF_ToverMIC_4X" ~ "I4_ToverMIC_4X",
        PKPD_Index == "CSF_ToverMIC_10X" ~ "I5_ToverMIC_10X",
        TRUE ~ PKPD_Index
      ))
    return(pred_data)
    
  } else if (target == "Plasma") {
    pred_data <- pred_data |>
      filter(!((PKPD_Index %in% c("CENTRAL_ToverMIC", "CENTRAL_ToverMIC_4X", "CENTRAL_ToverMIC_10X")) & value > 100)) |>
      # Rename indexes colunms to be align with csf_pkpd_data
      mutate(PKPD_Index = case_when(
        PKPD_Index == "CENTRAL_Cmax" ~ "I1_Cmax",
        PKPD_Index == "CENTRAL_AUC_MIC" ~ "I2_AUC",
        PKPD_Index == "CENTRAL_ToverMIC" ~ "I3_ToverMIC",
        PKPD_Index == "CENTRAL_ToverMIC_4X" ~ "I4_ToverMIC_4X",
        PKPD_Index == "CENTRAL_ToverMIC_10X" ~ "I5_ToverMIC_10X",
        TRUE ~ PKPD_Index
      ))
    return(pred_data)
    
  } else if (target == "Mix") {
    pred_data <- pred_data |>
      filter(!((PKPD_Index %in% c("CENTRAL_ToverMIC", "CENTRAL_ToverMIC_4X", "CENTRAL_ToverMIC_10X")) & value > 100)) |>
      # Rename indexes colunms to be align with csf_pkpd_data
      mutate(PKPD_Index = case_when(
        PKPD_Index == "CENTRAL_Cmax" ~ "I1_Cmax",
        PKPD_Index == "CENTRAL_AUC_MIC" ~ "I2_AUC",
        PKPD_Index == "CENTRAL_ToverMIC" ~ "I3_ToverMIC",
        PKPD_Index == "CENTRAL_ToverMIC_4X" ~ "I4_ToverMIC_4X",
        PKPD_Index == "CENTRAL_ToverMIC_10X" ~ "I5_ToverMIC_10X",
        TRUE ~ PKPD_Index
      ))
    return(pred_data)
    
  } else
    warning("Target doesn't exist !")
}

csf_pred_data <- generate_pred_data(csf_corrCurve_data, "CSF")
plasma_pred_data <- generate_pred_data(plasma_corrCurve_data, "Plasma")
mix_pred_data <- generate_pred_data(mix_corrCurve_data, "Mix")

#############################################################################
########          Generate a dataset with observation mean           ########                       
#############################################################################

csf_obs_mean <- csf_sim_data |>
  group_by(STRN, MIC, DoseGroup, AMT, time, B0, PKPD_Index) |>
  summarize( 
    value = mean(value),
    deltaLog10CFU = mean(deltaLog10CFU),
    Rsq = mean(Rsq),
    Adj.Rsq = mean(Adj.Rsq)
  )

plasma_obs_mean <- plasma_sim_data |>
  group_by(STRN, MIC, DoseGroup, AMT, time, B0, PKPD_Index) |>
  summarize( 
    value = mean(value),
    deltaLog10CFU = mean(deltaLog10CFU),
    Rsq = mean(Rsq),
    Adj.Rsq = mean(Adj.Rsq)
  )

mix_obs_mean <- mix_sim_data |>
  group_by(STRN, MIC, DoseGroup, AMT, time, B0, PKPD_Index) |>
  summarize( 
    value = mean(value),
    deltaLog10CFU = mean(deltaLog10CFU),
    Rsq = mean(Rsq),
    Adj.Rsq = mean(Adj.Rsq)
  )

#############################################################################
########                    Environement cleaning                    ########                       
#############################################################################

sim_results_fractioned <- list(fractioned_results = fractioned_results, 
                               fractioned_24h = fractioned_24h, 
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

plasma_data <- list(
  sim_data = plasma_sim_data,
  correlation_data = plasma_corrCurve_data,
  pred_data = plasma_pred_data,
  obs_mean = plasma_obs_mean
)

csf_with_plasma_index_data <- list(
  sim_data = mix_sim_data,
  correlation_data = mix_corrCurve_data,
  pred_data = mix_pred_data,
  obs_mean = mix_obs_mean
)

mrgsolve_model <- list(model_file = model_file, model_list = model_list)

elapsed <- toc(quiet = TRUE)
message("Total execution time: ", round(elapsed$toc - elapsed$tic, 2), " seconds")

rm(fractioned_results, fractioned_24h, fractioned_24h_full, 
   continuous_inf_results, continuous_inf_24h, continuous_inf_amt, 
   control_results, control_24h, i_data, 
   csf_sim_data, plasma_sim_data, mix_sim_data,
   csf_pred_data, plasma_pred_data, mix_pred_data,
   csf_obs_mean, plasma_obs_mean, mix_obs_mean,
   csf_corrCurve_data, plasma_corrCurve_data, mix_corrCurve_data,
   model_file, model_list, elapsed)

gc()