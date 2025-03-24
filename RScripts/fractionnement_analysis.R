q24h <- csf_data$sim_data |>
  filter(DoseGroup == "fractioned_24h") |>
  select(-c(Rsq))
q24h_corrCruve <- index_PKPD_fit_curve(q24h)
q24h_pred_data <- generate_pred_data(q24h_corrCruve)

q12h <- csf_data$sim_data |>
  filter(DoseGroup == "fractioned_12h") |>
  select(-c(Rsq))
q12h_corrCruve <- index_PKPD_fit_curve(q12h)
q12h_pred_data <- generate_pred_data(q12h_corrCruve)

q8h <- csf_data$sim_data |>
  filter(DoseGroup == "fractioned_8h") |>
  select(-c(Rsq))
q8h_corrCruve <- index_PKPD_fit_curve(q8h)
q8h_pred_data <- generate_pred_data(q8h_corrCruve)

q4h <- csf_data$sim_data |>
  filter(DoseGroup == "fractioned_4h") |>
  select(-c(Rsq))
q4h_corrCruve <- index_PKPD_fit_curve(q4h)
q4h_pred_data <- generate_pred_data(q4h_corrCruve)

# Define your data for plotting
global_pred_dat = csf_data$pred_data |> filter(PKPD_Index == "I2_AUC")

pred_q24_data = q24h_pred_data |> filter(PKPD_Index == "I2_AUC")
corr_q24_data = q24h_corrCruve |> filter(PKPD_Index == "I2_AUC")
  
pred_q12_data = q12h_pred_data |> filter(PKPD_Index == "I2_AUC")
corr_q12_data = q12h_corrCruve |> filter(PKPD_Index == "I2_AUC")
  
pred_q8_data = q8h_pred_data |> filter(PKPD_Index == "I2_AUC")
corr_q8_data = q8h_corrCruve |> filter(PKPD_Index == "I2_AUC")
  
pred_q4_data = q4h_pred_data |> filter(PKPD_Index == "I2_AUC")
corr_q4_data = q4h_corrCruve |> filter(PKPD_Index == "I2_AUC")

## AUC
p1 <- ggplot() +
  scale_x_continuous(limits = c(1, 1000), 
                     breaks = c(1, 10, 100, 1000),
                     trans="log10",
                     name = expression(italic(f)*AUC/MIC)) +
  scale_y_continuous(
    breaks = c(-4, -2, 0, 2, 4),
    limits = c(-5, 5),
    name = expression(Delta~log[10]~'CFU/ml')
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(data = global_pred_dat, mapping = aes(x=value, y=pred), linewidth = 0.8) +
  geom_line(data = pred_q24_data, 
            mapping = aes(x = value, y = pred), 
            linewidth = 0.8, 
            color = "red") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        strip.background = element_rect(fill = "white", colour = "transparent"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0)
      ) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R² = ", round(corr_q24_data |> pull(Rsq), 2))), 
            hjust = 1.2, vjust = 2, size = 4, fontface = "italic", color = "black") +
  geom_text(aes(x = 2, y = -4, label = "q24h"), hjust = 0, vjust = 0, size = 4, fontface = "italic", color = "black")

## q12h
p2 <- ggplot() +
  scale_x_continuous(limits = c(1, 1000), 
                     breaks = c(1, 10, 100, 1000),
                     trans="log10",
                     name = expression(italic(f)*AUC/MIC)) +
  scale_y_continuous(
    breaks = c(-4, -2, 0, 2, 4),
    limits = c(-5, 5),
    name = expression(Delta~log[10]~'CFU/ml')
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(data = global_pred_dat, mapping = aes(x=value, y=pred), linewidth = 0.8) +
  geom_line(data = pred_q12_data, 
            mapping = aes(x = value, y = pred), 
            linewidth = 0.8, 
            color = "red") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        strip.background = element_rect(fill = "white", colour = "transparent"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0)
  ) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R² = ", round(corr_q12_data |> pull(Rsq), 2))), 
            hjust = 1.2, vjust = 2, size = 4, fontface = "italic", color = "black") +
  geom_text(aes(x = 2, y = -4, label = "q12h"), hjust = 0, vjust = 0, size = 4, fontface = "italic", color = "black")

## q8h
p3 <- ggplot() +
  scale_x_continuous(limits = c(1, 1000), 
                     breaks = c(1, 10, 100, 1000),
                     trans="log10",
                     name = expression(italic(f)*AUC/MIC)) +
  scale_y_continuous(
    breaks = c(-4, -2, 0, 2, 4),
    limits = c(-5, 5),
    name = expression(Delta~log[10]~'CFU/ml')
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(data = global_pred_dat, mapping = aes(x=value, y=pred), linewidth = 0.8) +
  geom_line(data = pred_q8_data, 
            mapping = aes(x = value, y = pred), 
            linewidth = 0.8, 
            color = "red") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        strip.background = element_rect(fill = "white", colour = "transparent"),
        plot.margin = margin(0, 0, 0, 0),
        strip.placement = "outside",
        strip.text.x = element_text(size = 12)
  ) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R² = ", round(corr_q8_data |> pull(Rsq), 2))), 
            hjust = 1.2, vjust = 2, size = 4, fontface = "italic", color = "black") +
  geom_text(aes(x = 2, y = -4, label = "q8h"), hjust = 0, vjust = 0, size = 4, fontface = "italic", color = "black")


## qh
p4 <- ggplot() +
  scale_x_continuous(limits = c(1, 1000), 
                     breaks = c(1, 10, 100, 1000),
                     trans="log10",
                     name = expression(italic(f)*AUC/MIC)) +
  scale_y_continuous(
    breaks = c(-4, -2, 0, 2, 4),
    limits = c(-5, 5),
    name = expression(Delta~log[10]~'CFU/ml')
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(data = global_pred_dat, mapping = aes(x=value, y=pred), linewidth = 0.8) +
  geom_line(data = pred_q4_data, 
            mapping = aes(x = value, y = pred), 
            linewidth = 0.8, 
            color = "red") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        strip.background = element_rect(fill = "white", colour = "transparent"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        strip.placement = "outside",
        strip.text.x = element_text(size = 12)
  ) +
  geom_text(aes(x = Inf, y = Inf, label = paste("R² = ", round(corr_q4_data |> pull(Rsq), 2))), 
            hjust = 1.2, vjust = 2, size = 4, fontface = "italic", color = "black") +
  geom_text(aes(x = 2, y = -4, label = "q4h"), hjust = 0, vjust = 0, size = 4, fontface = "italic", color = "black")


library(patchwork)
p1 + p2 + p3 + p4 +
  plot_layout(ncol = 2, guides = "collect")

ggsave(file.path(getwd(), "results", "Impact du fractionnement - AUC.png"), width = 8, height = 6, units = "in", dpi = 1200)
