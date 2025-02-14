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

generate_pred_data <- function(data, obs_data, n_sim) {
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
  
   x <- x |> left_join(obs_data, by = c("PKPD_Index")) |>
     slice(rep(1:n(), each = n_sim)) |>
     mutate(REP = rep(1:n_sim, length.out = n()),
            PRED  = Imax_model(value, I0, Imax, IC50, H),
            IPRED = Imax_model(value, I0, Imax, IC50, H) + rnorm(n(), mean = 0, sd = Sd_Residuals)) |>
     select(STRN, DoseGroup, AMT, ID, PKPD_Index, REP, value, deltaLog10CFU, PRED, IPRED)
   
   if (any(is.na(x))) {
     stop("NA values in sim_data")
   }
   
  return(x)
}