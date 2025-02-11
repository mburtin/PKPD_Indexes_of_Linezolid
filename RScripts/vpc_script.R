generate_obs_data <- function(data) {
  
  x <- data |>
    ungroup() |>
    filter(DoseGroup != "control") |>
    select(-c(Log10CFU, Rsq, Adj.Rsq, B0, time)) |>
    rename(REP = ID,
           DV = deltaLog10CFU,
           Index_value = value) |>
    group_by(DoseGroup, AMT) |>
    mutate(ID = cur_group_id()) |>
    ungroup() |>
    select(STRN, MIC, ID, REP, DoseGroup, AMT, PKPD_Index, Index_value, DV)
  
  return(x)
}

generate_pred_data <- function(data) {
  x <- data |> 
    mutate(PKPD_Index = case_when(
      PKPD_Index %in% c("CENTRAL_Cmax", "CSF_Cmax") ~ "I1_Cmax",
      PKPD_Index %in% c("CENTRAL_AUC_MIC", "CSF_AUC_MIC") ~ "I2_AUC",
      PKPD_Index %in% c("CENTRAL_ToverMIC", "CSF_ToverMIC") ~ "I3_ToverMIC",
      PKPD_Index %in% c("CENTRAL_ToverMIC_4X", "CSF_ToverMIC_4X") ~ "I4_ToverMIC_4X",
      PKPD_Index %in% c("CENTRAL_ToverMIC_10X", "CSF_ToverMIC_10X") ~ "I5_ToverMIC_10X",
      TRUE ~ PKPD_Index
    )) |>
    select(-c(data, nls, I0_RSE, Imax_RSE, IC50_RSE, H_RSE, Rsq, Adj.Rsq, Target_stasis, Target_1LogKill, Target_2LogKill)) |>
    ungroup()
  
  return(x)
}

# Obsbs generate df
obs_plasma_data <- generate_obs_data(plasma_data[["sim_data"]])
obs_csf_data <- generate_obs_data(csf_data[["sim_data"]])

if(any(is.na(obs_plasma_data))) {
  stop("NA values in obs_plasma_data")
} else if (any(is.na(obs_csf_data))) {
  stop("NA values in obs_csf_data")
} else {
  message("Obs data generated successfully")
}

# Plasma sim generation df
isim_plasma_data <- generate_pred_data(plasma_data[["correlation_data"]])

isim_plasma_data <- left_join(obs_plasma_data, isim_plasma_data, by = "PKPD_Index") |>
  group_by(STRN, ID, REP, DoseGroup, PKPD_Index) |>
  mutate(PRED  = Imax_model(Index_value, I0, Imax, IC50, H),
         IPRED = Imax_model(Index_value, I0, Imax, IC50, H) + sample(Residuals[[1]] |> pull(), size = 1, replace = TRUE)) |>
  select(-c(Residuals))

if(any(is.na(isim_plasma_data))) {
  stop("NA values in isim_plasma_data")
} else {
  message("Plasma sim data generated successfully")
}

# CSF sim generation df
isim_csf_data <- generate_pred_data(csf_data[["correlation_data"]])

isim_csf_data <- isim_csf_data |>
  left_join(obs_csf_data, isim_csf_data, by = "PKPD_Index") |>
  group_by(STRN, ID, REP, DoseGroup, PKPD_Index) |>
  mutate(PRED  = Imax_model(Index_value, I0, Imax, IC50, H),
         IPRED = Imax_model(Index_value, I0, Imax, IC50, H) + sample(Residuals[[1]] |> pull(), size = 1, replace = TRUE)) |>
  select(-c(Residuals))

if(any(is.na(isim_csf_data))) {
  stop("NA values in isim_csf_data")
} else {
  message("CSF sim data generated successfully")
}