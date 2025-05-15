#######
##  Generic functions
######
generate_obs_data <- function(data) {
  x <- data |>
    ungroup() |>
    filter(DoseGroup != "control") |>
    select(STRN, DoseGroup, AMT, ID, PKPD_Index, value, deltaLog10CFU)

  if (any(is.na(x))) {
    stop("NA values in obs_data")
  }

  return(x)
}

generate_pred_data <- function(sim_data, obs_data, pkpd_index, n_sim) {
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
    left_join(obs_data |> filter(PKPD_Index == pkpd_index), by = c("PKPD_Index")) |>
    slice(rep(1:n(), each = n_sim)) |>
    mutate(
      REP = rep(1:n_sim, length.out = n()),
      PRED = Emax_model(value, I0, Imax, IC50, H),
      IPRED = Emax_model(value, I0, Imax, IC50, H) + rnorm(n(), mean = 0, sd = Sd_Residuals)
    ) |>
    select(STRN, DoseGroup, AMT, ID, PKPD_Index, REP, value, deltaLog10CFU, PRED, IPRED)

  if (any(is.na(x))) {
    stop("NA values in sim_data")
  }

  return(x)
}

ntile_obs_bins <- function(data, n_bins) {
  data |>
    mutate(bin = ntile(value, n_bins)) |>
    group_by(bin) |>
    summarize(
      x_mean = mean(value),
      p5 = quantile(deltaLog10CFU, 0.05),
      p50 = quantile(deltaLog10CFU, 0.5),
      p95 = quantile(deltaLog10CFU, 0.95)
    )
}

ntile_sim_bins <- function(data, n_bins) {
  data |>
    group_by(REP) |>
    mutate(bin = ntile(value, n_bins)) |>
    group_by(REP, bin) |>
    summarize(
      x_mean = mean(value),
      p5 = quantile(IPRED, 0.05),
      p50 = quantile(IPRED, 0.5),
      p95 = quantile(IPRED, 0.95)
    ) |>
    pivot_longer(cols = c(p5, p50, p95), names_to = "percentile", values_to = "value") |>
    group_by(percentile, bin, x_mean) |>
    summarize(
      ic_05 = quantile(value, 0.05),
      ic_50 = quantile(value, 0.5),
      ic_95 = quantile(value, 0.95)
    )
}

cut_obs_bins <- function(data, seq_start, seq_end, seq_steps) {
  data |>
    mutate(bin = cut(value, breaks = seq(seq_start, seq_end, by = seq_steps), right = FALSE)) |>
    group_by(bin) |>
    summarize(
      x_mean = mean(value),
      p5 = quantile(deltaLog10CFU, 0.05),
      p50 = quantile(deltaLog10CFU, 0.5),
      p95 = quantile(deltaLog10CFU, 0.95)
    )
}

cut_sim_bins <- function(data, seq_start, seq_end, seq_steps) {
  data |>
    group_by(REP) |>
    mutate(bin = cut(value, breaks = seq(seq_start, seq_end, by = seq_steps), right = FALSE)) |>
    group_by(REP, bin) |>
    summarize(
      x_mean = mean(value),
      p5 = quantile(IPRED, 0.05),
      p50 = quantile(IPRED, 0.5),
      p95 = quantile(IPRED, 0.95)
    ) |>
    pivot_longer(cols = c(p5, p50, p95), names_to = "percentile", values_to = "value") |>
    group_by(percentile, bin, x_mean) |>
    summarize(
      ic_05 = quantile(value, 0.05),
      ic_50 = quantile(value, 0.5),
      ic_95 = quantile(value, 0.95)
    )
}

#######
##  Dataset generation
######
csf_vpc_data <- list(
  obs_data = generate_obs_data(csf_data[["sim_data"]]),
  cmax_sim = generate_pred_data(
    sim_data = csf_data[["correlation_data"]],
    obs_data = generate_obs_data(csf_data[["sim_data"]]),
    pkpd_index = "I1_Cmax",
    n_sim = 100
  ),
  auc_sim = generate_pred_data(
    sim_data = csf_data[["correlation_data"]],
    obs_data = generate_obs_data(csf_data[["sim_data"]]),
    pkpd_index = "I2_AUC",
    n_sim = 100
  ),
  tmic_sim = generate_pred_data(
    sim_data = csf_data[["correlation_data"]],
    obs_data = generate_obs_data(csf_data[["sim_data"]]),
    pkpd_index = "I3_ToverMIC",
    n_sim = 100
  )
)

csf_vpc_with_plasma_index_data <- list(
  obs_data = generate_obs_data(csf_with_plasma_index_data[["sim_data"]]),
  cmax_sim = generate_pred_data(
    sim_data = csf_with_plasma_index_data[["correlation_data"]],
    obs_data = generate_obs_data(csf_with_plasma_index_data[["sim_data"]]),
    pkpd_index = "I1_Cmax",
    n_sim = 100
  ),
  auc_sim = generate_pred_data(
    sim_data = csf_with_plasma_index_data[["correlation_data"]],
    obs_data = generate_obs_data(csf_with_plasma_index_data[["sim_data"]]),
    pkpd_index = "I2_AUC",
    n_sim = 100
  ),
  tmic_sim = generate_pred_data(
    sim_data = csf_with_plasma_index_data[["correlation_data"]],
    obs_data = generate_obs_data(csf_with_plasma_index_data[["sim_data"]]),
    pkpd_index = "I3_ToverMIC",
    n_sim = 100
  )
)

gc()
