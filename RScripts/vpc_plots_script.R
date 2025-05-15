VPC_Plots <- function(data, file_name, Rsq_data) {
  require(ggplot2)
  require(patchwork)

  # Ajouter une légende pour les couleurs dans les géométries
  p1 <- ggplot() +
    geom_smooth(
      data = data$cmax_obs, aes(x = x_mean, y = p5, linetype = "p5_obs"),
      method = "loess", color = "black", linewidth = 0.5, se = FALSE, span = 1
    ) +
    geom_smooth(
      data = data$cmax_obs, aes(x = x_mean, y = p50, linetype = "p50_obs"),
      method = "loess", color = "black", linewidth = 0.75, se = FALSE, span = 1
    ) +
    geom_smooth(
      data = data$cmax_obs, aes(x = x_mean, y = p95, linetype = "p95_obs"),
      method = "loess", color = "black", linewidth = 0.5, se = FALSE, span = 1
    ) +
    geom_smooth(
      data = data$cmax_sim |> filter(percentile == "p5"), aes(x = x_mean, y = ic_50, color = "p5_sim"),
      method = "loess", linetype = "solid", linewidth = 0.6, se = TRUE, span = 1
    ) +
    geom_smooth(
      data = data$cmax_sim |> filter(percentile == "p50"), aes(x = x_mean, y = ic_50, color = "p50_sim"),
      method = "loess", linetype = "solid", linewidth = 0.6, se = TRUE, span = 1
    ) +
    geom_smooth(
      data = data$cmax_sim |> filter(percentile == "p95"), aes(x = x_mean, y = ic_50, color = "p95_sim"),
      method = "loess", linetype = "solid", linewidth = 0.6, se = TRUE, span = 1
    ) +
    scale_x_continuous(
      limits = c(0.1, 100),
      breaks = c(0.1, 1, 10, 100),
      trans = "log10",
      name = expression(italic(f) * Cmax / MIC)
    ) +
    scale_y_continuous(
      limits = c(-5.5, 6),
      breaks = c(4, 2, 0, -2, -4),
    ) +
    scale_linetype_manual(
      name = "Empirical percentiles",
      values = c("p5_obs" = "dashed", "p50_obs" = "solid", "p95_obs" = "dashed"),
      labels = c("p5_obs" = "5%", "p50_obs" = "50%", "p95_obs" = "95%")
    ) +
    scale_color_manual(
      name = "Predicted percentiles",
      values = c("p5_sim" = "blue", "p50_sim" = "red", "p95_sim" = "blue"),
      labels = c("p5_sim" = "5%", "p50_sim" = "50%", "p95_sim" = "95%")
    ) +
    guides(
      color = guide_legend(order = 2, ncol = 3),
      linetype = guide_legend(order = 1, ncol = 3)
    ) +
    labs(x = "fCmax/MIC", y = expression(Delta ~ log[10] ~ "CFU/ml")) +
    geom_text(
      aes(
        x = Inf, y = Inf,
        label = paste("R² = ", round(Rsq_data$correlation_data |> filter(PKPD_Index %in% c("CENTRAL_Cmax", "CSF_Cmax")) |> pull(Rsq), 2))
      ),
      hjust = 1.1, vjust = 3, size = 4, fontface = "italic", color = "black"
    ) +
    theme_minimal() +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.3)
    )

  # Ajouter une légende pour les couleurs dans les géométries
  p2 <- ggplot() +
    geom_smooth(
      data = data$auc_obs, aes(x = x_mean, y = p5, linetype = "p5_obs"),
      method = "loess", color = "black", linewidth = 0.5, se = FALSE, span = 1
    ) +
    geom_smooth(
      data = data$auc_obs, aes(x = x_mean, y = p50, linetype = "p50_obs"),
      method = "loess", color = "black", linewidth = 0.75, se = FALSE, span = 1
    ) +
    geom_smooth(
      data = data$auc_obs, aes(x = x_mean, y = p95, linetype = "p95_obs"),
      method = "loess", color = "black", linewidth = 0.5, se = FALSE, span = 1
    ) +
    geom_smooth(
      data = data$auc_sim |> filter(percentile == "p5"), aes(x = x_mean, y = ic_50, color = "p5_sim"),
      method = "loess", linetype = "solid", linewidth = 0.6, se = TRUE, span = 1
    ) +
    geom_smooth(
      data = data$auc_sim |> filter(percentile == "p50"), aes(x = x_mean, y = ic_50, color = "p50_sim"),
      method = "loess", linetype = "solid", linewidth = 0.6, se = TRUE, span = 1
    ) +
    geom_smooth(
      data = data$auc_sim |> filter(percentile == "p95"), aes(x = x_mean, y = ic_50, color = "p95_sim"),
      method = "loess", linetype = "solid", linewidth = 0.6, se = TRUE, span = 1
    ) +
    scale_x_continuous(
      limits = c(1, 1000),
      breaks = c(1, 10, 100, 1000),
      trans = "log10",
      name = expression(italic(f) * AUC / MIC)
    ) +
    scale_y_continuous(
      limits = c(-5.5, 6),
      breaks = c(4, 2, 0, -2, -4),
    ) +
    scale_linetype_manual(
      name = "Empirical percentiles",
      values = c("p5_obs" = "dashed", "p50_obs" = "solid", "p95_obs" = "dashed"),
      labels = c("p5_obs" = "5%", "p50_obs" = "50%", "p95_obs" = "95%")
    ) +
    scale_color_manual(
      name = "Predicted percentiles",
      values = c("p5_sim" = "blue", "p50_sim" = "red", "p95_sim" = "blue"),
      labels = c("p5_sim" = "5%", "p50_sim" = "50%", "p95_sim" = "95%")
    ) +
    guides(
      color = guide_legend(order = 2, ncol = 3),
      linetype = guide_legend(order = 1, ncol = 3)
    ) +
    labs(x = "fAUC/MIC", y = NULL) +
    geom_text(
      aes(
        x = Inf, y = Inf,
        label = paste("R² = ", round(Rsq_data$correlation_data |> filter(PKPD_Index %in% c("CENTRAL_AUC_MIC", "CSF_AUC_MIC")) |> pull(Rsq), 2))
      ),
      hjust = 1.1, vjust = 3, size = 4, fontface = "italic", color = "black"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.3)
    )

  # Ajouter une légende pour les couleurs dans les géométries
  p3 <- ggplot() +
    geom_smooth(
      data = data$tmic_obs, aes(x = x_mean, y = p5, linetype = "p5_obs"),
      method = "loess", color = "black", linewidth = 0.5, se = FALSE, span = 1
    ) +
    geom_smooth(
      data = data$tmic_obs, aes(x = x_mean, y = p50, linetype = "p50_obs"),
      method = "loess", color = "black", linewidth = 0.75, se = FALSE, span = 1
    ) +
    geom_smooth(
      data = data$tmic_obs, aes(x = x_mean, y = p95, linetype = "p95_obs"),
      method = "loess", color = "black", linewidth = 0.5, se = FALSE, span = 1
    ) +
    geom_smooth(
      data = data$tmic_sim |> filter(percentile == "p5"), aes(x = x_mean, y = ic_50, color = "p5_sim"),
      method = "loess", linetype = "solid", linewidth = 0.6, se = TRUE, span = 1
    ) +
    geom_smooth(
      data = data$tmic_sim |> filter(percentile == "p50"), aes(x = x_mean, y = ic_50, color = "p50_sim"),
      method = "loess", linetype = "solid", linewidth = 0.6, se = TRUE, span = 1
    ) +
    geom_smooth(
      data = data$tmic_sim |> filter(percentile == "p95"), aes(x = x_mean, y = ic_50, color = "p95_sim"),
      method = "loess", linetype = "solid", linewidth = 0.6, se = TRUE, span = 1
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      name = expression(italic(f) * T[">MIC"] ~ "(%)")
    ) +
    scale_y_continuous(
      limits = c(-5.5, 6),
      breaks = c(4, 2, 0, -2, -4),
    ) +
    scale_linetype_manual(
      name = "Empirical percentiles",
      values = c("p5_obs" = "dashed", "p50_obs" = "solid", "p95_obs" = "dashed"),
      labels = c("p5_obs" = "5%", "p50_obs" = "50%", "p95_obs" = "95%")
    ) +
    scale_color_manual(
      name = "Predicted percentiles",
      values = c("p5_sim" = "blue", "p50_sim" = "red", "p95_sim" = "blue"),
      labels = c("p5_sim" = "5%", "p50_sim" = "50%", "p95_sim" = "95%")
    ) +
    guides(
      color = guide_legend(order = 2, ncol = 3),
      linetype = guide_legend(order = 1, ncol = 3)
    ) +
    labs(x = "fT>MIC (%)", y = NULL) +
    geom_text(
      aes(
        x = Inf, y = Inf,
        label = paste("R² = ", round(Rsq_data$correlation_data |> filter(PKPD_Index %in% c("CENTRAL_ToverMIC", "CSF_ToverMIC")) |> pull(Rsq), 2))
      ),
      hjust = 1.1, vjust = 3, size = 4, fontface = "italic", color = "black"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.3)
    )

  # Merge the 3 VPCs and collect the legend at the top
  p1 + p2 + p3 +
    plot_layout(ncol = 3, guides = "collect") +
    plot_annotation(theme = theme(legend.position = "top"))

  ggsave(file.path(getwd(), "results", file_name), width = 9, height = 3, dpi = "retina")
}
