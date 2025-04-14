PKPD_index_plots <- function(obs_data, pred_data, title, file_name) {
  # Check if library are loaded
  require(ggplot2)        # Plotting
  require(ggh4x)          # Plotting (allow to set custom scales for each facet)
  
  # Generate the plot
  ggplot(obs_data, aes(x=value, y=deltaLog10CFU)) +  # Ajout de color = DoseGroup
    geom_point(size=1, aes(color = as.factor(DoseGroup))) +
    facet_wrap(~ PKPD_Index, scales = "free_x", strip.position = "bottom", nrow = 1, 
               labeller = labeller(
                 PKPD_Index = c(
                   "I4_Cmax" = "fCmax/MIC",
                   "I2_AUC" = "fAUC/MIC",
                   "I3_ToverMIC" = "fT>MIC (%)"))) +
    facetted_pos_scales(
      x = list(
        PKPD_Index == "I4_Cmax" ~ scale_x_continuous(limits = c(1, 25),
                                                     breaks = c(1, 10),
                                                     trans = "log10",
                                                     name = expression(italic(f)*Cmax/MIC)),
        PKPD_Index == "I2_AUC" ~ scale_x_continuous(limits = c(1, 125), 
                                                    breaks = c(1, 10, 100),
                                                    trans = "log10",
                                                    name = expression(italic(f)*AUC/MIC)),
        PKPD_Index == "I3_ToverMIC" ~ scale_x_continuous(limits = c(1, 125),
                                                         breaks = c(1, 10, 100),
                                                         trans = "log10",
                                                         name = expression(T['>MIC']~'(%)'))
      ),
    ) +
    geom_line(data = pred_data |> filter(PKPD_Index != "I4_Cmax"), aes(x=value, y=pred), linewidth=0.5) +
    geom_text(aes(x = Inf, y = Inf, label = paste("R² = ", round(Rsq, 2))), 
              hjust = 1.2, vjust = 2, size = 4, fontface = "italic", color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_y_continuous(
      breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
      limits = c(-5, 5),
      name = expression(Delta~log[10]~'CFU/ml')
    ) +
    # Rename legends groups
    scale_color_manual(
      values = c("#377EB8", "purple", "#4DAF4A", "#FF7F00", "red"),  # Assigner des couleurs spécifiques
      labels = c("q4h", "q6h", "q8h", "q12h", "q24h"),  # Nom des groupes
      name = "Dosing regimens"  # Titre de la légende
    ) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour="black"),
          strip.background = element_rect(fill = "white", colour = "transparent"),
          strip.placement = "outside",
          axis.title.y = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.position = "bottom",
    ) +
    ggtitle(title)
  
  # To keep organize the workspace, generate a result folder if it doesn't exist
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  # Save the plot
  ggsave(file.path(getwd(), "results", file_name), width = 9, height = 3, units = "in", dpi = 1200)
}