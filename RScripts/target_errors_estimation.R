temp <- csf_data[["correlation_data"]] |>
  select(c(
    PKPD_Index,
    Target_stasis,
    Target_1LogKill,
    Target_2LogKill,
  )) |>
  pivot_longer(-c(PKPD_Index),
    names_to = "EfficacyCriteria",
    values_to = "TargetIndexValue"
  ) |>
  mutate(
    TargetIndexValue =
      map_dbl(
        TargetIndexValue,
        ~ if (anyNA(.x)) {
          return(NA)
        } else {
          round(.x[1, 1], 2)
        }
      )
  ) |>
  pivot_wider(names_from = EfficacyCriteria, values_from = TargetIndexValue)

target_error_csf <- csf_data[["correlation_data"]] |>
  select(-c(data, nls, I0_RSE, Imax_RSE, IC50_RSE, H_RSE, Rsq, Adj.Rsq, Target_stasis, Target_1LogKill, Target_2LogKill))

# bind result and target_csf_error by PKPD_Index
target_error_csf <- target_error_csf |>
  left_join(temp, by = "PKPD_Index") |>
  slice(rep(1:n(), each = 1000)) |>
  mutate(
    REP = rep(1:1000, length.out = n()),
    Stasis_pred = Imax_model(Target_stasis, I0, Imax, IC50, H) + rnorm(n = n(), mean = 0, sd = Sd_Residuals),
    T1Log_pred = Imax_model(Target_1LogKill, I0, Imax, IC50, H) + rnorm(n = n(), mean = 0, sd = Sd_Residuals),
    T2Log_pred = Imax_model(Target_2LogKill, I0, Imax, IC50, H) + rnorm(n = n(), mean = 0, sd = Sd_Residuals)
  ) |>
  group_by(PKPD_Index) |>
  summarise(
    Stasis_target = unique(Target_stasis),
    Stasis_error = sum(Stasis_pred > 0) / 1000 * 100,
    T1Log_target = unique(Target_1LogKill),
    T1Log_error = sum(T1Log_pred > -1) / 1000 * 100,
    T1Log_stase_error = sum(T1Log_pred > 0) / 1000 * 100,
    T2Log_target = unique(Target_2LogKill),
    T2Log_error = sum(T2Log_pred > -2) / 1000 * 100,
    T2Log_stase_error = sum(T2Log_pred > 0) / 1000 * 100
  ) |>
  filter(PKPD_Index != "CSF_ToverMIC_4X" & PKPD_Index != "CSF_ToverMIC_10X")

### CSF With plasma index data
temp <- csf_with_plasma_index_data[["correlation_data"]] |>
  select(c(
    PKPD_Index,
    Target_stasis,
    Target_1LogKill,
    Target_2LogKill,
  )) |>
  pivot_longer(-c(PKPD_Index),
    names_to = "EfficacyCriteria",
    values_to = "TargetIndexValue"
  ) |>
  mutate(
    TargetIndexValue =
      map_dbl(
        TargetIndexValue,
        ~ if (anyNA(.x)) {
          return(NA)
        } else {
          round(.x[1, 1], 2)
        }
      )
  ) |>
  pivot_wider(names_from = EfficacyCriteria, values_from = TargetIndexValue)

target_error_csf_with_plasma_index <- csf_with_plasma_index_data[["correlation_data"]] |>
  select(-c(data, nls, I0_RSE, Imax_RSE, IC50_RSE, H_RSE, Rsq, Adj.Rsq, Target_stasis, Target_1LogKill, Target_2LogKill))

target_error_csf_with_plasma_index <- target_error_csf_with_plasma_index |>
  left_join(temp, by = "PKPD_Index") |>
  slice(rep(1:n(), each = 1000)) |>
  mutate(
    REP = rep(1:1000, length.out = n()),
    Stasis_pred = Imax_model(Target_stasis, I0, Imax, IC50, H) + rnorm(n = n(), mean = 0, sd = Sd_Residuals),
    T1Log_pred = Imax_model(Target_1LogKill, I0, Imax, IC50, H) + rnorm(n = n(), mean = 0, sd = Sd_Residuals),
    T2Log_pred = Imax_model(Target_2LogKill, I0, Imax, IC50, H) + rnorm(n = n(), mean = 0, sd = Sd_Residuals)
  ) |>
  group_by(PKPD_Index) |>
  summarise(
    Stasis_target = unique(Target_stasis),
    Stasis_error = sum(Stasis_pred > 0) / 1000 * 100,
    T1Log_target = unique(Target_1LogKill),
    T1Log_error = sum(T1Log_pred > -1) / 1000 * 100,
    T1Log_stase_error = sum(T1Log_pred > 0) / 1000 * 100,
    T2Log_target = unique(Target_2LogKill),
    T2Log_error = sum(T2Log_pred > -2) / 1000 * 100,
    T2Log_stase_error = sum(T2Log_pred > 0) / 1000 * 100
  ) |>
  filter(PKPD_Index != "CSF_ToverMIC_4X" & PKPD_Index != "CSF_ToverMIC_10X")

rm(temp)

## Generic function to calculate the target error
calculate_target_error <- function(index_value, I0, Imax, IC50, H, Sd_Residuals, n_sim) {
  results <- tibble()

  for (i in 1:n_sim) {
    results[i, 1] <- Imax_model(index_value, I0, Imax, IC50, H) + rnorm(n = 1, mean = 0, sd = Sd_Residuals)
  }

  error_percentage <- sum(results > 0) / n_sim * 100

  return(error_percentage)
}

calculate_target_error(
  index_value = 80,
  I0 = csf_data[["correlation_data"]][2, "I0"],
  Imax = csf_data[["correlation_data"]][2, "Imax"],
  IC50 = csf_data[["correlation_data"]][2, "IC50"],
  H = csf_data[["correlation_data"]][2, "H"],
  Sd_Residuals = csf_data[["correlation_data"]][2, "Sd_Residuals"] |> pull(),
  n_sim = 10000
)

calculate_target_error(
  index_value = 80,
  I0 = csf_with_plasma_index_data[["correlation_data"]][2, "I0"],
  Imax = csf_with_plasma_index_data[["correlation_data"]][2, "Imax"],
  IC50 = csf_with_plasma_index_data[["correlation_data"]][2, "IC50"],
  H = csf_with_plasma_index_data[["correlation_data"]][2, "H"],
  Sd_Residuals = csf_with_plasma_index_data[["correlation_data"]][2, "Sd_Residuals"] |> pull(),
  n_sim = 10000
)
