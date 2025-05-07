## Note:
#
#
# 
##

x = sim_formated_data |> filter(DoseGroup != "control")

CSF_AUC_AIO <- x |> 
  filter(PKPD_Index == "CSF_AUC_MIC") |>
  ungroup() |>
  group_by(DoseGroup, AMT) |> 
  summarize(
    DoseGroup = unique(DoseGroup),
    AMT = unique(AMT),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 30) & ((Log10CFU_CSF - B0) > 0)) / n() * 100, 2),
    Stase_FN = round(sum((value < 30) & ((Log10CFU_CSF - B0) <= 0)) / n() * 100, 2),
    Stase_TP = round(sum((value >= 30) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 30) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 73) & ((Log10CFU_CSF - B0) > -2)) / n() * 100, 2),
    Bactericidal_FN = round(sum((value < 73) & ((Log10CFU_CSF - B0) <= -2)) / n() * 100, 2),
    Bactericidal_TP = round(sum((value >= 73) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 73) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CSF_AUC_DoseGroup <- x |> 
  filter(PKPD_Index == "CSF_AUC_MIC") |>
  ungroup() |>
  group_by(DoseGroup) |> 
  summarize(
    DoseGroup = unique(DoseGroup),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 30) & ((Log10CFU_CSF - B0) > 0)) / n() * 100, 2),
    Stase_FN = round(sum((value < 30) & ((Log10CFU_CSF - B0) <= 0)) / n() * 100, 2),
    Stase_TP = round(sum((value >= 30) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 30) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 73) & ((Log10CFU_CSF - B0) > -2)) / n() * 100, 2),
    Bactericidal_FN = round(sum((value < 73) & ((Log10CFU_CSF - B0) <= -2)) / n() * 100, 2),
    Bactericidal_TP = round(sum((value >= 73) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 73) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CSF_AUC_AMT <- x |> 
  filter(PKPD_Index == "CSF_AUC_MIC") |>
  ungroup() |>
  group_by(AMT) |> 
  summarize(
    AMT = unique(AMT),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 30) & ((Log10CFU_CSF - B0) > 0)) / n() * 100, 2),
    Stase_FN = round(sum((value < 30) & ((Log10CFU_CSF - B0) <= 0)) / n() * 100, 2),
    Stase_TP = round(sum((value >= 30) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 30) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 73) & ((Log10CFU_CSF - B0) > -2)) / n() * 100, 2),
    Bactericidal_FN = round(sum((value < 73) & ((Log10CFU_CSF - B0) <= -2)) / n() * 100, 2),
    Bactericidal_TP = round(sum((value >= 73) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 73) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CENTRAL_AUC_AIO <- x |> 
  filter(PKPD_Index == "CENTRAL_AUC_MIC") |>
  ungroup() |>
  group_by(DoseGroup, AMT) |> 
  summarize(
    DoseGroup = unique(DoseGroup),
    AMT = unique(AMT),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 38) & ((Log10CFU_CSF - B0) > 0)) / n() * 100, 2),
    Stase_FN = round(sum((value < 38) & ((Log10CFU_CSF - B0) <= 0)) / n() * 100, 2),
    Stase_TP = round(sum((value >= 38) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 38) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) > -2)) / n() * 100, 2),
    Bactericidal_FN = round(sum((value < 91) & ((Log10CFU_CSF - B0) <= -2)) / n() * 100, 2),
    Bactericidal_TP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 91) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CENTRAL_AUC_DoseGroup <- x |> 
  filter(PKPD_Index == "CENTRAL_AUC_MIC") |>
  ungroup() |>
  group_by(DoseGroup) |> 
  summarize(
    DoseGroup = unique(DoseGroup),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 38) & ((Log10CFU_CSF - B0) > 0)) / n() * 100, 2),
    Stase_FN = round(sum((value < 38) & ((Log10CFU_CSF - B0) <= 0)) / n() * 100, 2),
    Stase_TP = round(sum((value >= 38) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 38) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) > -2)) / n() * 100, 2),
    Bactericidal_FN = round(sum((value < 91) & ((Log10CFU_CSF - B0) <= -2)) / n() * 100, 2),
    Bactericidal_TP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 91) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CENTRAL_AUC_AMT <- x |> 
  filter(PKPD_Index == "CENTRAL_AUC_MIC") |>
  ungroup() |>
  group_by(AMT) |> 
  summarize(
    AMT = unique(AMT),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 38) & ((Log10CFU_CSF - B0) > 0)) / n() * 100, 2),
    Stase_FN = round(sum((value < 38) & ((Log10CFU_CSF - B0) <= 0)) / n() * 100, 2),
    Stase_TP = round(sum((value >= 38) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 38) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) > -2)) / n() * 100, 2),
    Bactericidal_FN = round(sum((value < 91) & ((Log10CFU_CSF - B0) <= -2)) / n() * 100, 2),
    Bactericidal_TP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 91) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CSF_ToverMIC_AIO <- x |> 
  filter(PKPD_Index == "CSF_ToverMIC") |>
  ungroup() |>
  group_by(DoseGroup, AMT) |> 
  summarize(
    DoseGroup = unique(DoseGroup),
    AMT = unique(AMT),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 63) & ((Log10CFU_CSF - B0) > 0)) / n() * 100, 2),
    Stase_FN = round(sum((value < 63) & ((Log10CFU_CSF - B0) <= 0)) / n() * 100, 2),
    Stase_TP = round(sum((value >= 63) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 63) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) > -2)) / n() * 100, 2),
    Bactericidal_FN = round(sum((value < 91) & ((Log10CFU_CSF - B0) <= -2)) / n() * 100, 2),
    Bactericidal_TP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 91) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CSF_ToverMIC_DoseGroup <- x |> 
  filter(PKPD_Index == "CSF_ToverMIC") |>
  ungroup() |>
  group_by(DoseGroup) |> 
  summarize(
    DoseGroup = unique(DoseGroup),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 63) & ((Log10CFU_CSF - B0) > 0)) / n() * 100, 2),
    Stase_FN = round(sum((value < 63) & ((Log10CFU_CSF - B0) <= 0)) / n() * 100, 2),
    Stase_TP = round(sum((value >= 63) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 63) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) > -2)) / n() * 100, 2),
    Bactericidal_FN = round(sum((value < 91) & ((Log10CFU_CSF - B0) <= -2)) / n() * 100, 2),
    Bactericidal_TP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 91) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CSF_ToverMIC_AMT <- x |> 
  filter(PKPD_Index == "CSF_ToverMIC") |>
  ungroup() |>
  group_by(AMT) |> 
  summarize(
    AMT = unique(AMT),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 63) & ((Log10CFU_CSF - B0) > 0)) / n() * 100, 2),
    Stase_FN = round(sum((value < 63) & ((Log10CFU_CSF - B0) <= 0)) / n() * 100, 2),
    Stase_TP = round(sum((value >= 63) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 63) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) > -2)) / n() * 100, 2),
    Bactericidal_FN = round(sum((value < 91) & ((Log10CFU_CSF - B0) <= -2)) / n() * 100, 2),
    Bactericidal_TP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 91) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

#### Effect
ggplot(
  data = CSF_AUC_AIO,
  aes(
    x = AMT,
    y = Stase_effect,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  labs(
    title = "Stase effect through AMT/DoseGroup for fAUC/MIC",
    x = "AMT",
    y = "Stase effect (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CSF_AUC_AIO,
  aes(
    x = AMT,
    y = Bactericidal_effect,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  labs(
    title = "Bactericidal effect through AMT/DoseGroup for fAUC/MIC",
    x = "AMT",
    y = "Bactericide effect (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

#### False positive - Stase effect
ggplot(
  data = CSF_AUC_AIO,
  aes(
    x = AMT,
    y = Stase_FP,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "False positive for CSF fAUC/MIC with stasis effect",
    x = "AMT",
    y = "Percent of false positive (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CENTRAL_AUC_AIO,
  aes(
    x = AMT,
    y = Stase_FP,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "False positive for plasma fAUC/MIC with stasis effect",
    x = "AMT",
    y = "Percent of false positive (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CSF_ToverMIC_AIO,
  aes(
    x = AMT,
    y = Stase_FP,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "False positive for CSF fT>MIC with stasis effect",
    x = "AMT",
    y = "Percent of false positive (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")


#### False negative - Stase effect
ggplot(
  data = CSF_AUC_AIO,
  aes(
    x = AMT,
    y = Stase_FN,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "False negative for CSF fAUC/MIC with stasis effect",
    x = "AMT",
    y = "Percent of false negative (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CENTRAL_AUC_AIO,
  aes(
    x = AMT,
    y = Stase_FN,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "False negative for plasma fAUC/MIC with stasis effect",
    x = "AMT",
    y = "Percent of false negative (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CSF_ToverMIC_AIO,
  aes(
    x = AMT,
    y = Stase_FN,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "False negative for CSF fT>MIC with stasis effect",
    x = "AMT",
    y = "Percent of false negative (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

#### True positive - Stase effect
ggplot(
  data = CSF_AUC_AIO,
  aes(
    x = AMT,
    y = Stase_TP,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "True positive for CSF fAUC/MIC with stasis effect",
    x = "AMT",
    y = "Percent of true positive (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CENTRAL_AUC_AIO,
  aes(
    x = AMT,
    y = Stase_TP,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "True positive for plasma fAUC/MIC with stasis effect",
    x = "AMT",
    y = "Percent of true positive (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CSF_ToverMIC_AIO,
  aes(
    x = AMT,
    y = Stase_TP,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "True positive for CSF fT>MIC with stasis effect",
    x = "AMT",
    y = "Percent of true positive (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

#### False positive - Bactericidal effect
ggplot(
  data = CSF_AUC_AIO,
  aes(
    x = AMT,
    y = Bactericidal_FP,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "False positive for CSF fAUC/MIC with bactericidal effect",
    x = "AMT",
    y = "Percent of false positive (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CENTRAL_AUC_AIO,
  aes(
    x = AMT,
    y = Bactericidal_FP,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "False positive for plasma fAUC/MIC with bactericidal effect",
    x = "AMT",
    y = "Percent of false positive (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CSF_ToverMIC_AIO,
  aes(
    x = AMT,
    y = Bactericidal_FP,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "False positive for CSF fT>MIC with bactericidal effect",
    x = "AMT",
    y = "Percent of false positive (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")


#### False negative - Bactericidal effect
ggplot(
  data = CSF_AUC_AIO,
  aes(
    x = AMT,
    y = Bactericidal_FN,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "False negative for CSF fAUC/MIC with bactericidal effect",
    x = "AMT",
    y = "Percent of false negative (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CENTRAL_AUC_AIO,
  aes(
    x = AMT,
    y = Bactericidal_FN,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "False negative for plasma fAUC/MIC with bactericidal effect",
    x = "AMT",
    y = "Percent of false negative (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CSF_ToverMIC_AIO,
  aes(
    x = AMT,
    y = Bactericidal_FN,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "False negative for CSF fT>MIC with bactericidal effect",
    x = "AMT",
    y = "Percent of false negative (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

#### True positive - Bactericidal effect
ggplot(
  data = CSF_AUC_AIO,
  aes(
    x = AMT,
    y = Bactericidal_TP,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "True positive for CSF fAUC/MIC with bactericidal effect",
    x = "AMT",
    y = "Percent of true positive (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CENTRAL_AUC_AIO,
  aes(
    x = AMT,
    y = Bactericidal_TP,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "True positive for plasma fAUC/MIC with bactericidal effect",
    x = "AMT",
    y = "Percent of true positive (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")

ggplot(
  data = CSF_ToverMIC_AIO,
  aes(
    x = AMT,
    y = Bactericidal_TP,
    group = DoseGroup,
    color = DoseGroup
  )
) +
  geom_line() +
  geom_point() +
  scale_x_log10(limits = c(100, 10000)) +
  labs(
    title = "True positive for CSF fT>MIC with bactericidal effect",
    x = "AMT",
    y = "Percent of true positive (%)"
  ) +
  theme_minimal() +
  annotation_logticks(base = 10, sides = "b", short = unit(0.1, "cm"), mid = unit(0.1, "cm"),long = unit(0.1, "cm")) + 
  theme(legend.position = "bottom")