require(patchwork)

PKPD_target_plots_600mg <- function(obs_data, pred_data, effect, auc_target, tmic_target, title) {
  # Check if library are loaded
  require(ggplot2)        # Plotting
  require(ggh4x)          # Plotting (allow to set custom scales for each facet)
  
  auc_atteignment = (nrow(obs_data |> filter(PKPD_Index == "I2_AUC") |> filter(value >= auc_target)) / nrow(obs_data |> filter(PKPD_Index == "I2_AUC")))*100
  tmic_atteignment = (nrow(obs_data |> filter(PKPD_Index == "I3_ToverMIC") |> filter(value >= tmic_target)) / nrow(obs_data |> filter(PKPD_Index == "I3_ToverMIC")))*100
  
  # Generate the plot
  ggplot(obs_data |> filter(PKPD_Index %in% c("I2_AUC", "I3_ToverMIC")), aes(x=value, y=deltaLog10CFU)) +
    geom_point(size = 0.25) +
    geom_point(data = obs_data |> filter(PKPD_Index == "I2_AUC") |> filter(value >= auc_target), size = 0.1, color = "blue") +
    geom_point(data = obs_data |> filter(PKPD_Index == "I3_ToverMIC") |> filter(value >= tmic_target), size = 0.1, color = "blue") +
    facet_wrap(~ PKPD_Index, scales = "free_x", strip.position = "bottom", nrow = 1, 
               labeller = labeller(
                 PKPD_Index = c(
                   "I2_AUC" = "fAUC/MIC",
                   "I3_ToverMIC" = "fT>MIC (%)"))) +
    facetted_pos_scales(
      x = list(
        PKPD_Index == "I2_AUC" ~ scale_x_continuous(limits = c(1, 1000), 
                                                    breaks = c(1, 10, 100, 1000),
                                                    trans="log10",
                                                    name = expression(italic(f)*AUC/MIC)),
        PKPD_Index == "I3_ToverMIC" ~ scale_x_continuous(limits = c(0, 100),
                                                         breaks = c(0, 25, 50, 75, 100),
                                                         name = expression(T['>MIC']~'(%)'))
      ),
    ) +
    scale_y_continuous(
      breaks = c(-4, -2, 0, 2, 4),
      limits = c(-5, 5),
      name = expression(Delta~log[10]~'CFU/ml'),
    ) +
    geom_hline(yintercept = effect, linetype = "dotted") +
    geom_vline(data = obs_data |> filter(PKPD_Index == "I2_AUC"), aes(xintercept = auc_target), linetype = "dotted") +
    geom_vline(data = obs_data |> filter(PKPD_Index == "I3_ToverMIC"), aes(xintercept = tmic_target), linetype = "dotted") +
    geom_text(data = obs_data |> filter(PKPD_Index == "I3_ToverMIC"), 
              aes(x = 85, y = 4, label =  paste(round(tmic_atteignment, 0), "%")), size = 4, hjust = 0) +
    geom_text(data = obs_data |> filter(PKPD_Index == "I2_AUC"), 
              aes(x = 400, y = 4, label = paste(round(auc_atteignment, 0), "%")), size = 4, hjust = 0) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour="black"),
          strip.background = element_rect(fill = "white", colour = "transparent"),
          strip.placement = "outside",
          axis.title.y = element_blank(),
          strip.text.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()
    ) +
    geom_line(data = pred_data |> filter(PKPD_Index %in% c("I2_AUC", "I3_ToverMIC")), mapping = aes(x=value, y=pred), linewidth = 0.8) +
    labs(x = title)
}


PKPD_effect_plots_600mg <- function(obs_data, pred_data, effect, auc_target, tmic_target, title) {
  # Check if library are loaded
  require(ggplot2)        # Plotting
  require(ggh4x)          # Plotting (allow to set custom scales for each facet)
  
  auc_effect = (nrow(obs_data |> filter(PKPD_Index == "I2_AUC") |> filter(deltaLog10CFU <= effect)) / nrow(obs_data |> filter(PKPD_Index == "I2_AUC")))*100
  tmic_effect = (nrow(obs_data |> filter(PKPD_Index == "I3_ToverMIC") |> filter(deltaLog10CFU <= effect)) / nrow(obs_data |> filter(PKPD_Index == "I3_ToverMIC")))*100
  
  # Generate the plot
  ggplot(obs_data |> filter(PKPD_Index %in% c("I2_AUC", "I3_ToverMIC")), aes(x=value, y=deltaLog10CFU)) +
    geom_point(size = 0.1) +
    geom_point(data = obs_data |> filter(PKPD_Index == "I2_AUC") |> filter(deltaLog10CFU <= effect), size = 0.1, color = "orange") +
    geom_point(data = obs_data |> filter(PKPD_Index == "I3_ToverMIC") |> filter(deltaLog10CFU <= effect), size = 0.1, color = "orange") +
    facet_wrap(~ PKPD_Index, scales = "free_x", strip.position = "bottom", nrow = 1, 
               labeller = labeller(
                 PKPD_Index = c(
                   "I2_AUC" = "fAUC/MIC",
                   "I3_ToverMIC" = "fT>MIC (%)"))) +
    facetted_pos_scales(
      x = list(
        PKPD_Index == "I2_AUC" ~ scale_x_continuous(limits = c(1, 1000), 
                                                    breaks = c(1, 10, 100, 1000),
                                                    trans="log10",
                                                    name = expression(italic(f)*AUC/MIC)),
        PKPD_Index == "I3_ToverMIC" ~ scale_x_continuous(limits = c(0, 100),
                                                         breaks = c(0, 25, 50, 75, 100),
                                                         name = expression(T['>MIC']~'(%)'))
      ),
    ) +
    scale_y_continuous(
      breaks = c(-4, -2, 0, 2, 4),
      limits = c(-5, 5),
      name = expression(Delta~log[10]~'CFU/ml')
    ) +
    geom_hline(yintercept = effect, linetype = "dotted") +
    geom_vline(data = obs_data |> filter(PKPD_Index == "I2_AUC"), aes(xintercept = auc_target), linetype = "dotted") +
    geom_vline(data = obs_data |> filter(PKPD_Index == "I3_ToverMIC"), aes(xintercept = tmic_target), linetype = "dotted") +
    geom_text(data = obs_data |> filter(PKPD_Index == "I3_ToverMIC"), 
              aes(x = 85, y = 4, label =  paste(round(tmic_effect, 0), "%")), size = 4, hjust = 0) +
    geom_text(data = obs_data |> filter(PKPD_Index == "I2_AUC"), 
              aes(x = 400, y = 4, label = paste(round(auc_effect, 0), "%")), size = 4, hjust = 0) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour="black"),
          strip.background = element_rect(fill = "white", colour = "transparent"),
          strip.placement = "outside",
          plot.margin = margin(0, 0, 0, 0),
          axis.title.y = element_blank(),
          strip.text.x = element_text(size = 12),
          axis.title.x = element_blank(),
    ) +
    geom_line(data = pred_data |> filter(PKPD_Index %in% c("I2_AUC", "I3_ToverMIC")), mapping = aes(x=value, y=pred), linewidth = 0.8)
}