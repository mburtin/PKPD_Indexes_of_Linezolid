generate_PKPD_plots <- function(obs_data, pred_data, title, filename) {
  # Check if library are loaded
  require(ggplot2)        # Plotting
  require(ggh4x)          # Plotting (allow to set custom scales for each facet)
  
  # Generate the plot
  ggplot(obs_data |> filter(PKPD_Index %in% c("I2_AUC", "I3_ToverMIC")), aes(x=value, y=deltaLog10CFU)) +
    geom_point(size=0.25) +
    geom_point(data = obs_data |> filter(PKPD_Index == "I2_AUC") |> filter(value >= 73), size=0.25, color = "blue") +
    geom_point(data = obs_data |> filter(PKPD_Index == "I3_ToverMIC") |> filter(value >= 91), size=0.25, color = "blue") +
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
    geom_hline(yintercept = -2, linetype = "dotted") +
    geom_vline(data = obs_data |> filter(PKPD_Index == "I2_AUC"), aes(xintercept = 73), linetype = "dotted") +
    geom_vline(data = obs_data |> filter(PKPD_Index == "I3_ToverMIC"), aes(xintercept = 91), linetype = "dotted") +
    geom_text(data = obs_data |> filter(PKPD_Index == "I2_AUC"), aes(x = 500, y = 3.5, label = "16%"), size = 4) +
    geom_text(data = obs_data |> filter(PKPD_Index == "I3_ToverMIC"), aes(x = 78, y = 3.5, label = "47%"), size = 4) +
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


generate_PKPD_plots2 <- function(obs_data, pred_data, title, filename) {
  # Check if library are loaded
  require(ggplot2)        # Plotting
  require(ggh4x)          # Plotting (allow to set custom scales for each facet)
  
  # Generate the plot
  ggplot(obs_data |> filter(PKPD_Index %in% c("I2_AUC", "I3_ToverMIC")), aes(x=value, y=deltaLog10CFU)) +
    geom_point(size=0.25) +
    geom_point(data = obs_data |> filter(PKPD_Index == "I2_AUC") |> filter(deltaLog10CFU <= -2), size=0.25, color = "orange") +
    geom_point(data = obs_data |> filter(PKPD_Index == "I3_ToverMIC") |> filter(deltaLog10CFU <= -2), size=0.25, color = "orange") +
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
    geom_hline(yintercept = -2, linetype = "dotted") +
    geom_vline(data = obs_data |> filter(PKPD_Index == "I2_AUC"), aes(xintercept = 73), linetype = "dotted") +
    geom_vline(data = obs_data |> filter(PKPD_Index == "I3_ToverMIC"), aes(xintercept = 91), linetype = "dotted") +
    geom_text(data = obs_data |> filter(PKPD_Index == "I2_AUC"), aes(x = 500, y = 3.5, label = "26%"), size = 4) +
    geom_text(data = obs_data |> filter(PKPD_Index == "I3_ToverMIC"), aes(x = 78, y = 3.5, label = "26%"), size = 4) +
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