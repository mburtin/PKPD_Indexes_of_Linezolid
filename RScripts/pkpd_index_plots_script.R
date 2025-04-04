PKPD_index_plots <- function(obs_data, pred_data, title, file_name) {
  # Check if library are loaded
  require(ggplot2)        # Plotting
  require(ggh4x)          # Plotting (allow to set custom scales for each facet)
  
  # Generate the plot
  ggplot(obs_data, aes(x=value, y=deltaLog10CFU, color = DoseGroup)) +  # Ajout de color = DoseGroup
    geom_point(size=0.25) +
    facet_wrap(~ PKPD_Index, scales = "free_x", strip.position = "bottom", nrow = 1, 
               labeller = labeller(
                 PKPD_Index = c(
                   "I1_Cmax" = "fCmax/MIC",
                   "I2_AUC" = "fAUC/MIC",
                   "I3_ToverMIC" = "fT>MIC (%)"))) +
    facetted_pos_scales(
      x = list(
        PKPD_Index == "I1_Cmax" ~ scale_x_continuous(limits = c(1, 100),
                                                     breaks = c(1, 10, 100),
                                                     trans="log10",
                                                     name = expression(italic(f)*Cmax/MIC)),
        PKPD_Index == "I2_AUC" ~ scale_x_continuous(limits = c(10, 1000), 
                                                    breaks = c(10, 100, 1000),
                                                    trans="log10",
                                                    name = expression(italic(f)*AUC/MIC)),
        PKPD_Index == "I3_ToverMIC" ~ scale_x_continuous(limits = c(0, 100),
                                                         breaks = c(0, 20, 40, 60, 80, 100),
                                                         name = expression(T['>MIC']~'(%)'))
      ),
    ) +
    scale_y_continuous(
      breaks = c(0, 1, 2, 3, 4, 5, 6, 7),
      limits = c(-1, 7),
      name = expression(Delta~log[10]~'CFU/ml')
    ) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour="black"),
          strip.background = element_rect(fill = "white", colour = "transparent"),
          strip.placement = "outside",
          axis.title.y = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          axis.title.x = element_blank()
    ) +
    ggtitle(title)
  
  # To keep organize the workspace, generate a result folder if it doesn't exist
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  # Save the plot
  ggsave(file.path(getwd(), "results", file_name), width = 12, height = 3, units = "in", dpi = 1200)
}