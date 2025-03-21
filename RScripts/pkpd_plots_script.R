generate_PKPD_plots <- function(obs_data, pred_data, title, filename) {
  # Check if library are loaded
  require(ggplot2)        # Plotting
  require(ggh4x)          # Plotting (allow to set custom scales for each facet)
  
  # Generate the plot
  ggplot(obs_data, aes(x=value, y=deltaLog10CFU)) +
    geom_point(size=0.25) +
    facet_wrap(~ PKPD_Index, scales = "free_x", strip.position = "bottom", nrow = 1, 
               labeller = labeller(
                 PKPD_Index = c(
                   "I1_Cmax" = "fCmax/MIC",
                   "I2_AUC" = "fAUC/MIC",
                   "I3_ToverMIC" = "fT>MIC (%)"))) +
    facetted_pos_scales(
      x = list(
        PKPD_Index == "I1_Cmax" ~ scale_x_continuous(limits = c(0.1, 100),
                                                     breaks = c(0.1, 1, 10, 100),
                                                     trans="log10",
                                                     name = expression(italic(f)*Cmax/MIC)),
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
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour="black"),
          strip.background = element_rect(fill = "white", colour = "transparent"),
          strip.placement = "outside",
          axis.title.y = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          axis.title.x = element_blank()
    ) +
    geom_line(data = pred_data, mapping = aes(x=value, y=pred), linewidth = 0.8) +
    geom_text(aes(x = Inf, y = Inf, label = paste("R² = ", round(Rsq, 2))), 
              hjust = 1.2, vjust = 2, size = 4, fontface = "italic", color = "black") +
    ggtitle(title)
  
  # To keep organize the workspace, generate a result folder if it doesn't exist
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  # Save the plot
  ggsave(file.path(getwd(), "results", filename), width = 12, height = 3, units = "in", dpi = 1200)
}