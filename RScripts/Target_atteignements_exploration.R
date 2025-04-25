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
    Stase_FP = round(sum((value >= 30) & ((Log10CFU_CSF - B0) > 0)) / sum(value >= 30) * 100, 2),
    Stase_FN = round(sum((value < 30) & ((Log10CFU_CSF - B0) <= 0)) / sum((Log10CFU_CSF - B0) <= 0) * 100, 2),
    Stase_TP = round(sum((value >= 30) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 30) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 73) & ((Log10CFU_CSF - B0) > -2)) / sum(value >= 73) * 100, 2),
    Bactericidal_FN = round(sum((value < 73) & ((Log10CFU_CSF - B0) <= -2)) / sum((Log10CFU_CSF - B0) <= -2) * 100, 2),
    Bactericidal_TP = round(sum((value >= 73) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 73) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CSF_AUC_DoseGroup <- x |> 
  filter(PKPD_Index == "CSF_AUC_MIC") |>
  ungroup() |>
  group_by(DoseGroup) |> 
  summarize(
    DoseGroup = unique(DoseGroup),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 30) & ((Log10CFU_CSF - B0) > 0)) / sum(value >= 30) * 100, 2),
    Stase_FN = round(sum((value < 30) & ((Log10CFU_CSF - B0) <= 0)) / sum((Log10CFU_CSF - B0) <= 0) * 100, 2),
    Stase_TP = round(sum((value >= 30) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 30) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 73) & ((Log10CFU_CSF - B0) > -2)) / sum(value >= 73) * 100, 2),
    Bactericidal_FN = round(sum((value < 73) & ((Log10CFU_CSF - B0) <= -2)) / sum((Log10CFU_CSF - B0) <= -2) * 100, 2),
    Bactericidal_TP = round(sum((value >= 73) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 73) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CSF_AUC_AMT <- x |> 
  filter(PKPD_Index == "CSF_AUC_MIC") |>
  ungroup() |>
  group_by(AMT) |> 
  summarize(
    AMT = unique(AMT),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 30) & ((Log10CFU_CSF - B0) > 0)) / sum(value >= 30) * 100, 2),
    Stase_FN = round(sum((value < 30) & ((Log10CFU_CSF - B0) <= 0)) / sum((Log10CFU_CSF - B0) <= 0) * 100, 2),
    Stase_TP = round(sum((value >= 30) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 30) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 73) & ((Log10CFU_CSF - B0) > -2)) / sum(value >= 73) * 100, 2),
    Bactericidal_FN = round(sum((value < 73) & ((Log10CFU_CSF - B0) <= -2)) / sum((Log10CFU_CSF - B0) <= -2) * 100, 2),
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
    Stase_FP = round(sum((value >= 38) & ((Log10CFU_CSF - B0) > 0)) / sum(value >= 38) * 100, 2),
    Stase_FN = round(sum((value < 38) & ((Log10CFU_CSF - B0) <= 0)) / sum((Log10CFU_CSF - B0) <= 0) * 100, 2),
    Stase_TP = round(sum((value >= 38) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 38) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) > -2)) / sum(value >= 91) * 100, 2),
    Bactericidal_FN = round(sum((value < 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((Log10CFU_CSF - B0) <= -2) * 100, 2),
    Bactericidal_TP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 91) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CENTRAL_AUC_DoseGroup <- x |> 
  filter(PKPD_Index == "CENTRAL_AUC_MIC") |>
  ungroup() |>
  group_by(DoseGroup) |> 
  summarize(
    DoseGroup = unique(DoseGroup),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 38) & ((Log10CFU_CSF - B0) > 0)) / sum(value >= 38) * 100, 2),
    Stase_FN = round(sum((value < 38) & ((Log10CFU_CSF - B0) <= 0)) / sum((Log10CFU_CSF - B0) <= 0) * 100, 2),
    Stase_TP = round(sum((value >= 38) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 38) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) > -2)) / sum(value >= 91) * 100, 2),
    Bactericidal_FN = round(sum((value < 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((Log10CFU_CSF - B0) <= -2) * 100, 2),
    Bactericidal_TP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 91) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CENTRAL_AUC_AMT <- x |> 
  filter(PKPD_Index == "CENTRAL_AUC_MIC") |>
  ungroup() |>
  group_by(AMT) |> 
  summarize(
    AMT = unique(AMT),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 38) & ((Log10CFU_CSF - B0) > 0)) / sum(value >= 38) * 100, 2),
    Stase_FN = round(sum((value < 38) & ((Log10CFU_CSF - B0) <= 0)) / sum((Log10CFU_CSF - B0) <= 0) * 100, 2),
    Stase_TP = round(sum((value >= 38) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 38) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) > -2)) / sum(value >= 91) * 100, 2),
    Bactericidal_FN = round(sum((value < 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((Log10CFU_CSF - B0) <= -2) * 100, 2),
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
    Stase_FP = round(sum((value >= 63) & ((Log10CFU_CSF - B0) > 0)) / sum(value >= 63) * 100, 2),
    Stase_FN = round(sum((value < 63) & ((Log10CFU_CSF - B0) <= 0)) / sum((Log10CFU_CSF - B0) <= 0) * 100, 2),
    Stase_TP = round(sum((value >= 63) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 63) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) > -2)) / sum(value >= 91) * 100, 2),
    Bactericidal_FN = round(sum((value < 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((Log10CFU_CSF - B0) <= -2) * 100, 2),
    Bactericidal_TP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 91) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CSF_ToverMIC_DoseGroup <- x |> 
  filter(PKPD_Index == "CSF_ToverMIC") |>
  ungroup() |>
  group_by(DoseGroup) |> 
  summarize(
    DoseGroup = unique(DoseGroup),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 63) & ((Log10CFU_CSF - B0) > 0)) / sum(value >= 63) * 100, 2),
    Stase_FN = round(sum((value < 63) & ((Log10CFU_CSF - B0) <= 0)) / sum((Log10CFU_CSF - B0) <= 0) * 100, 2),
    Stase_TP = round(sum((value >= 63) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 63) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) > -2)) / sum(value >= 91) * 100, 2),
    Bactericidal_FN = round(sum((value < 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((Log10CFU_CSF - B0) <= -2) * 100, 2),
    Bactericidal_TP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 91) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

CSF_ToverMIC_AMT <- x |> 
  filter(PKPD_Index == "CSF_ToverMIC") |>
  ungroup() |>
  group_by(AMT) |> 
  summarize(
    AMT = unique(AMT),
    Stase_effect = round(sum((Log10CFU_CSF - B0) <= 0) / n() * 100, 2),
    Stase_FP = round(sum((value >= 63) & ((Log10CFU_CSF - B0) > 0)) / sum(value >= 63) * 100, 2),
    Stase_FN = round(sum((value < 63) & ((Log10CFU_CSF - B0) <= 0)) / sum((Log10CFU_CSF - B0) <= 0) * 100, 2),
    Stase_TP = round(sum((value >= 63) & ((Log10CFU_CSF - B0) <= 0)) / sum((value >= 63) | ((Log10CFU_CSF - B0) <= 0)) * 100, 2),
    Bactericidal_effect = round(sum((Log10CFU_CSF - B0) <= -2) / n() * 100, 2),
    Bactericidal_FP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) > -2)) / sum(value >= 91) * 100, 2),
    Bactericidal_FN = round(sum((value < 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((Log10CFU_CSF - B0) <= -2) * 100, 2),
    Bactericidal_TP = round(sum((value >= 91) & ((Log10CFU_CSF - B0) <= -2)) / sum((value >= 91) | ((Log10CFU_CSF - B0) <= -2)) * 100, 2),
  )

