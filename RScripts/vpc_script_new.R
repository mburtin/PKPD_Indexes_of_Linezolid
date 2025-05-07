# Record execution time
library(tictoc)
tic()

#######
##  Generic functions
######
generate_obs_data <- function(data, pkpd_index, n_bins = 20) {
  x <- data |>
    ungroup() |>
    filter(DoseGroup != "control") |>
    select(STRN, DoseGroup, AMT, ID, PKPD_Index, value, deltaLog10CFU) |>
    filter(PKPD_Index == pkpd_index)
  
  # Use cut for ToverMIC, ntile for others
  if (pkpd_index == "I3_ToverMIC") {
    x <- x |>
      mutate(bin = cut(value, breaks = seq(0, 100, by = 20), right = FALSE))
  } else {
    x <- x |>
      mutate(bin = ntile(value, n_bins))
  }
  
  x <- x |>
    group_by(bin) |>
    summarize(
      x_mean = mean(value),
      p5 = quantile(deltaLog10CFU, 0.05),
      p50 = quantile(deltaLog10CFU, 0.5),
      p95 = quantile(deltaLog10CFU, 0.95)
    )
  
  if (any(is.na(x))) {
    stop("NA values in obs_data")
  }
  
  message(paste0(" - Empirical data completed for ", pkpd_index))
  
  return(x)
}

generate_pred_data <- function(sim_data, obs_data, pkpd_index, n_bins = 20) {
  # Prepare observed data
  obs_data <- obs_data |>
    ungroup() |>
    filter(DoseGroup != "control") |>
    select(STRN, DoseGroup, AMT, ID, PKPD_Index, value, deltaLog10CFU) |>
    filter(PKPD_Index == pkpd_index)
  
  # Prepare simulated data
  x <- sim_data |> 
    mutate(PKPD_Index = case_when(
      PKPD_Index %in% c("CENTRAL_Cmax", "CSF_Cmax") ~ "I1_Cmax",
      PKPD_Index %in% c("CENTRAL_AUC_MIC", "CSF_AUC_MIC") ~ "I2_AUC",
      PKPD_Index %in% c("CENTRAL_ToverMIC", "CSF_ToverMIC") ~ "I3_ToverMIC",
      TRUE ~ PKPD_Index
    )) |>
    select(-c(data, nls, Rsq, Adj.Rsq, Target_stasis, Target_1LogKill, Target_2LogKill)) |>
    filter(PKPD_Index == pkpd_index) |>
    ungroup()
  
  x <- x |> 
    left_join(obs_data, by = c("PKPD_Index")) |>
    mutate(
      PRED = Emax_model(value, I0, Imax, IC50, H),
      IPRED = Emax_model(value, I0, Imax, IC50, H) + rnorm(n(), mean = 0, sd = Sd_Residuals)
    )
  
  # Use cut for ToverMIC, ntile for others
  if (pkpd_index == "I3_ToverMIC") {
    x <- x |>
      mutate(bin = cut(value, breaks = seq(0, 100, by = 20), right = FALSE))
  } else {
    x <- x |>
      mutate(bin = ntile(value, n_bins))
  }
  
  x <- x |>
    group_by(bin) |>
    summarize(
      x_mean = mean(value),
      p5 = quantile(IPRED, 0.05),
      p50 = quantile(IPRED, 0.5),
      p95 = quantile(IPRED, 0.95)
    )
   
  if (any(is.na(x))) {
    stop("NA values in sim_data")
  }
   
  return(x)
}

generate_pred_data_replicates <- function(sim_data, obs_data, pkpd_index, n_replicates = 500, n_bins = 20) {
  # Generate all replicates
  replicates <- map(1:n_replicates, ~{
    generate_pred_data(sim_data, obs_data, pkpd_index, n_bins) |>
      mutate(REP = .x)
  }) |>
    bind_rows()
  
  # Calculate percentiles between each percentiles of all bins/replicates
  result <- replicates |>
    pivot_longer(cols = c(p5, p50, p95), names_to = "percentile", values_to = "value") |>
    group_by(percentile, bin, x_mean) |>
    summarize(
      ic_05 = quantile(value, 0.05),
      ic_50 = quantile(value, 0.5),
      ic_95 = quantile(value, 0.95)
    )
  
  message(paste0(" - Predictions data completed for ", pkpd_index))
  
  return(result)
}

#######
##  Dataset generation
######
message("Starting VPC data generation...")

csf_vpc_obs_data <- list(
    cmax = generate_obs_data(csf_data[["sim_data"]], pkpd_index = "I1_Cmax"),
    auc = generate_obs_data(csf_data[["sim_data"]], pkpd_index = "I2_AUC"),
    tmic = generate_obs_data(csf_data[["sim_data"]], pkpd_index = "I3_ToverMIC")
  )

csf_vpc_pred_data <- list(
  cmax_sim = generate_pred_data_replicates(sim_data = csf_data[["correlation_data"]], 
                                          obs_data = csf_data[["sim_data"]],
                                          pkpd_index = "I1_Cmax",
                                          n_replicates = 100,
                                          n_bins = 20),
  auc_sim = generate_pred_data_replicates(sim_data = csf_data[["correlation_data"]], 
                                         obs_data = csf_data[["sim_data"]],
                                         pkpd_index = "I2_AUC",
                                         n_replicates = 100,
                                         n_bins = 20),
  tmic_sim = generate_pred_data_replicates(sim_data = csf_data[["correlation_data"]], 
                                          obs_data = csf_data[["sim_data"]],
                                          pkpd_index = "I3_ToverMIC",
                                          n_replicates = 100,
                                          n_bins = 20)
)

# Stop recording time execution
elapsed <- toc(quiet = TRUE)
message("VPC data generated in: ", round(elapsed$toc - elapsed$tic, 0), " seconds")
