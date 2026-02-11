library(caret)
library(MLeval)
library(tidyverse)
library(xgboost)
library(pROC)
source('main_ml_func.R')

keys <- c(
  "ID",
  "class",
  "DPF3_rs4903064",
  "ZEB",
  "condition",
  "PVT1_rs35252396",
  "EPAS1_rs7579899",
  "ITPR2_rs10771279",
  "ITPR2_rs1049380"
)


dominant_ml <- purrr::reduce(dominant_long_filt, full_join, by = keys) %>%
  na.omit()
codominant_ml <- purrr::reduce(codominant_long_filt, full_join, by = keys) %>%
  na.omit()
recessive_ml <- purrr::reduce(recessive_long_filt, full_join, by = keys) %>%
  na.omit()


ml_datalist <- list(
  dominant = dominant_ml,
  codominant = codominant_ml,
  recessive = recessive_ml
)

ml_datalist <- lapply(ml_datalist, function(df) {
  df %>%
    mutate(across(
      where(is.character) & !any_of(c("ID", "condition")),
      as.factor
    ))
})


ml_results_dominant <- ml_models(
  data = ml_datalist$dominant,
  y = 'condition',
  number = 100,
  ctrl = 'boot632',
  metric = 'ROC'
)

ml_results_recessive <- ml_models(
  data = ml_datalist$recessive,
  y = 'condition',
  number = 100,
  ctrl = 'boot632',
  metric = 'ROC'
)


ml_results_codominant <- ml_models(
  data = ml_datalist$codominant,
  y = 'condition',
  number = 100,
  ctrl = 'boot632',
  metric = 'ROC'
)


ml_metrics_dominant <- ml_metrics(
  modelist = ml_results_dominant$models,
  y = 'condition',
  train_data = ml_results_dominant$train_data,
  test_data = ml_results_dominant$test_data
)

ml_metrics_recessive <- ml_metrics(
  modelist = ml_results_recessive$models,
  y = 'condition',
  train_data = ml_results_recessive$train_data,
  test_data = ml_results_recessive$test_data
)

ml_metrics_codominant <- ml_metrics(
  modelist = ml_results_codominant$models,
  y = 'condition',
  train_data = ml_results_codominant$train_data,
  test_data = ml_results_codominant$test_data
)

ml_metrics_list <- list(
  dominant = ml_metrics_dominant,
  recessive = ml_metrics_recessive,
  codominant = ml_metrics_codominant
)


library(openxlsx)

cm_to_df <- function(cm) {
  as.data.frame(cm$table)
}

lapply(names(ml_metrics_list), function(model_type) {
  metrics <- ml_metrics_list[[model_type]]
  file_name <- paste0(model_type, "_ml_metrics.xlsx")

  wb <- createWorkbook()

  for (metric_name in names(metrics)) {
    metric <- metrics[[metric_name]]

    if (inherits(metric, "confusionMatrix")) {
      metric <- cm_to_df(metric)
    }

    if (is.data.frame(metric)) {
      addWorksheet(wb, sheetName = metric_name)
      writeData(wb, sheet = metric_name, x = metric)
    }
  }

  saveWorkbook(wb, file = file_name, overwrite = TRUE)
})


mtx_surv <- read.delim(
  'mtx_surv.txt',
  header = T,
  sep = '\t',
  dec = "."
) %>%
  select(., -DIAG, -SEG) %>%
  mutate(across(everything(), ~ ifelse(. == "si", "y", .))) %>%
  mutate(across(everything(), ~ ifelse(. == "no", "n", .))) %>%
  mutate(across(everything(), ~ ifelse(. == "1\n", "1", .))) %>%
  mutate(across(
    c(
      'ITPR2_rs1049380',
      'ITPR2_rs10771279',
      'DPF3_rs4903064',
      'PVT1_rs35252396',
      'EPAS1_rs7579899',
      'CLASS',
      'TRT',
      'MTX_DIAG',
      'MTX_POST'
    ),
    as.factor
  ))


mtx_surv_dominant <- mtx_surv %>%
  mutate(
    ITPR2_rs1049380 = case_when(
      ITPR2_rs1049380 == "T/T" ~ "T/T",
      TRUE ~ "G/T + G/G"
    ),
    ITPR2_rs10771279 = case_when(
      ITPR2_rs10771279 == "T/T" ~ "T/T",
      TRUE ~ "C/T + C/C"
    ),
    DPF3_rs4903064 = case_when(
      DPF3_rs4903064 == "T/T" ~ "T/T",
      TRUE ~ "C/T + C/C"
    ),
    EPAS1_rs7579899 = case_when(
      EPAS1_rs7579899 == "A/A" ~ "A/A",
      TRUE ~ "A/G + G/G"
    ),
    PVT1_rs35252396 = case_when(
      PVT1_rs35252396 == "CG/CG" ~ "CG/CG",
      TRUE ~ "AC/CG + AC/AC"
    )
  )


mtx_surv_recessive <- mtx_surv %>%
  mutate(
    ITPR2_rs10771279 = case_when(
      ITPR2_rs10771279 == "C/C" ~ "C/C",
      TRUE ~ "C/T + T/T"
    ),
    ITPR2_rs1049380 = case_when(
      ITPR2_rs1049380 == "G/G" ~ "G/G",
      TRUE ~ "G/T + T/T"
    ),
    DPF3_rs4903064 = case_when(
      DPF3_rs4903064 == "C/C" ~ "C/C",
      TRUE ~ "C/T + T/T"
    ),
    EPAS1_rs7579899 = case_when(
      EPAS1_rs7579899 == "G/G" ~ "G/G",
      TRUE ~ "A/G + A/A"
    ),
    PVT1_rs35252396 = case_when(
      PVT1_rs35252396 == "AC/AC" ~ "AC/AC",
      TRUE ~ "AC/CG + CG/CG"
    )
  )


mtx_surv_list <- list(
  codominant = mtx_surv,
  dominant = mtx_surv_dominant,
  recessive = mtx_surv_recessive
)
mtx_surv_list <- lapply(mtx_surv_list, function(x) {
  x <- x %>%
    select(matches("^(PVT1|EPAS1|MYC|ITPR2|DPF3|ZEB)"), MTX_DIAG) %>%
    select(-ITPR2_rs1049380, MTX_DIAG) %>%
    na.omit()
})


ml_mtx_dominant <- ml_models(
  data = mtx_surv_list$dominant,
  y = 'MTX_DIAG',
  number = 100,
  ctrl = 'boot632',
  metric = 'ROC'
)

ml_mtx_recessive <- ml_models(
  data = mtx_surv_list$recessive,
  y = 'MTX_DIAG',
  number = 100,
  ctrl = 'boot632',
  metric = 'ROC'
)


ml_mtx_codominant <- ml_models(
  data = mtx_surv_list$codominant,
  y = 'MTX_DIAG',
  number = 100,
  ctrl = 'boot632',
  metric = 'ROC'
)


ml_mtx_metrics_dominant <- ml_metrics(
  modelist = ml_mtx_dominant$models,
  y = 'MTX_DIAG',
  train_data = ml_mtx_dominant$train_data,
  test_data = ml_mtx_dominant$test_data
)

ml_mtx_metrics_recessive <- ml_metrics(
  modelist = ml_mtx_recessive$models,
  y = 'MTX_DIAG',
  train_data = ml_mtx_recessive$train_data,
  test_data = ml_mtx_recessive$test_data
)

ml_mtx_metrics_codominant <- ml_metrics(
  modelist = ml_mtx_codominant$models,
  y = 'MTX_DIAG',
  train_data = ml_mtx_codominant$train_data,
  test_data = ml_mtx_codominant$test_data
)


save.image(file = "ml_expr.RData")
