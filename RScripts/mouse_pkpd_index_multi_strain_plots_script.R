PKPD_index_plots <- function(obs_data, pred_data, title, file_name) {
  # Check if library are loaded
  require(ggplot2)        # Plotting
  require(ggh4x)          # Plotting (allow to set custom scales for each facet)
  
  # Generate the plot
  ggplot(obs_data, aes(x=value, y=deltaLog10CFU)) +  # Ajout de color = DoseGroup
    geom_point(size=1, aes(colour = as.factor(DoseGroup))) +
    facet_wrap(STRN ~ PKPD_Index, scales = "free_x", strip.position = "bottom", nrow = 4, 
               labeller = labeller(
                 PKPD_Index = c(
                   "I4_Cmax" = "Cmax/MIC",
                   "I2_AUC" = "AUC/MIC",
                   "I3_ToverMIC" = "T>MIC (%)"))) +
    facetted_pos_scales(
      x = list(
        PKPD_Index == "I4_Cmax" ~ scale_x_continuous(limits = c(0.75, 1300),
                                                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9,
                                                                10, 20, 30, 40, 50, 60, 70, 80, 90, 
                                                                100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
                                                     trans="log10",
                                                     labels = function(x) ifelse(x %in% c(1, 10, 100, 1000), as.character(x), ""),
                                                     name = expression(Cmax/MIC)),
        PKPD_Index == "I2_AUC" ~ scale_x_continuous(limits = c(8, 1300), 
                                                    breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 
                                                               100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
                                                    labels = function(x) ifelse(x %in% c(10, 100, 1000), as.character(x), ""),
                                                    trans="log10",
                                                    name = expression(AUC/MIC)),
        PKPD_Index == "I3_ToverMIC" ~ scale_x_continuous(limits = c(0, 100),
                                                         breaks = c(0, 20, 40, 60, 80, 100),
                                                         name = expression(T['>MIC']~'(%)'))
      ),
    ) +
    geom_hline(aes(yintercept = 3.27), linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_line(data = pred_data, aes(x=value, y=pred), linewidth=0.5) +
    geom_text(aes(x = Inf, y = Inf, label = paste("R² = ", round(Rsq, 2))), 
              hjust = 1.2, vjust = 2, size = 4, fontface = "italic", color = "black") +
    scale_y_continuous(
      breaks = c(0, 1, 2, 3, 4, 5, 6, 7),
      limits = c(-0.5, 7),
      labels = c("0", "-1", "-2", "-3", "-4", "-5", "-6", "-7"),
      name = "Change in log10 CFU/ml"
    ) +
    labs(color = "Dosing regimens") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour="black"),
          strip.background = element_rect(fill = "white", colour = "transparent"),
          strip.placement = "outside",
          axis.title.y = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.position = "bottom",
          panel.spacing = unit(1.5, "lines")
    ) +
    ggtitle(title)
  
  # To keep organize the workspace, generate a result folder if it doesn't exist
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  # Save the plot
  ggsave(file.path(getwd(), "results/mouse", file_name), width = 7.5, height = 16, units = "in", dpi = "retina")
}