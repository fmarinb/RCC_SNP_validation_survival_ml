library(survival)
library(survminer)
library(tidyverse)
library(ggeffects)

surv_genot <- read.delim(
  'survival.txt',
  header = T,
  sep = '\t',
  dec = "."
) %>%
  mutate(TRT = str_remove(TRT, "\\n$")) %>%
  mutate(across(
    c(
      'ITPR2_rs1049380',
      'ITPR2_rs10771279',
      'DPF3_rs4903064',
      'PVT1_rs35252396',
      'EPAS1_rs7579899',
      'TRT',
      'Sex',
      'Age.group',
      'Metastasis.at.diagnosis',
      'MTX_POST'
    ),
    as.factor
  ))


surv_genot_dominant <- surv_genot %>%
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
  ) %>%
  mutate(
    across(
      c(
        ITPR2_rs1049380,
        ITPR2_rs10771279,
        DPF3_rs4903064,
        EPAS1_rs7579899,
        PVT1_rs35252396
      ),
      as.factor
    )
  )


surv_genot_recessive <- surv_genot %>%
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
  ) %>%
  mutate(
    across(
      c(
        ITPR2_rs1049380,
        ITPR2_rs10771279,
        DPF3_rs4903064,
        EPAS1_rs7579899,
        PVT1_rs35252396
      ),
      as.factor
    )
  )


surv_list <- list(
  codominant = surv_genot,
  dominant = surv_genot_dominant,
  recessive = surv_genot_recessive
)

surv_list <- lapply(surv_list, function(x) {
  x <- x %>%
    select(., -ID, -MTX_POST, -DATE1, -DATE2)
})

surv_objects <- lapply(surv_list, function(x) {
  "Surv(OS, CLASS)"
})


cox_s <- function(data, s) {
  genot_vars <- c(
    "ITPR2_rs1049380",
    "ITPR2_rs10771279",
    "DPF3_rs4903064",
    "EPAS1_rs7579899",
    "PVT1_rs35252396"
  )

  covars <- c(
    "Metastasis.at.diagnosis",
    "Age.group",
    'Sex'
  )

  met = "Metastasis.at.diagnosis"

  rhs_int <- paste0(
    "(",
    paste(genot_vars, collapse = "+"),
    ")*",
    met,
    " + ",
    paste(covars, collapse = " + ")
  )
  fml_int <- as.formula(paste(s, "~", rhs_int))
  cox_ph_int <- coxph(
    fml_int,
    data = data,
    na.action = na.omit,
    ties = "efron",
    control = coxph.control(iter.max = 2000, eps = 1e-09)
  )

  rhs_strata <- paste0(
    paste(genot_vars, collapse = "+"),
    " + ",
    paste(covars[-1], collapse = " + "),
    " + ",
    sprintf("strata(%s)", met)
  )

  fml_strata <- as.formula(paste(s, "~", rhs_strata))
  cox_ph_strata <- coxph(
    fml_strata,
    data = data,
    na.action = na.omit,
    ties = "efron",
    control = coxph.control(iter.max = 2000, eps = 1e-09)
  )

  rhs_basic <- paste0(
    paste(covars, collapse = " + "),
    " + ",
    paste(genot_vars, collapse = " + ")
  )
  fml_basic <- as.formula(paste(s, "~", rhs_basic))
  cox_ph <- coxph(
    fml_basic,
    data = data,
    na.action = na.omit,
    ties = "efron",
    control = coxph.control(iter.max = 2000, eps = 1e-09)
  )

  return(list(
    cox_basic = cox_ph,
    cox_int = cox_ph_int,
    cox_strata = cox_ph_strata
  ))
}

cox_genot_adjs <- mapply(
  cox_s,
  data = surv_list,
  s = surv_objects,
  SIMPLIFY = FALSE
)


coxzph_dominant <- lapply(cox_genot_adjs$dominant, function(x) {
  xzph <- cox.zph(x)
  plots <- ggcoxzph(xzph)
  return(list(results = xzph, plots = plots))
})

coxzph_codominant <- lapply(cox_genot_adjs$codominant, function(x) {
  xzph <- cox.zph(x)
  plots <- ggcoxzph(xzph)
  return(list(results = xzph, plots = plots))
})

coxzph_recessive <- lapply(cox_genot_adjs$recessive, function(x) {
  xzph <- cox.zph(x)
  plots <- ggcoxzph(xzph)
  return(list(results = xzph, plots = plots))
})


surv_list$dominant <- surv_list$dominant %>%
  mutate(
    ITPR2_rs1049380_num = case_when(
      ITPR2_rs1049380 == "T/T" ~ 1,
      TRUE ~ 0
    )
  )


genot_vars2 <- c(
  "ITPR2_rs10771279",
  "DPF3_rs4903064",
  "EPAS1_rs7579899",
  "PVT1_rs35252396"
)

covars <- c(
  "Metastasis.at.diagnosis",
  "Age.group",
  'Sex'
)

met = "Metastasis.at.diagnosis"

rhs_strata2 <- paste0(
  paste(genot_vars2, collapse = "+"),
  " + ",
  paste(covars[-1], collapse = " + "),
  " + ",
  sprintf("strata(%s)", met),
  " + ",
  "tt(ITPR2_rs1049380_num)"
)

fml_strata2 <- as.formula(paste("Surv(OS, CLASS)", "~", rhs_strata2))

cox_time_dep <- coxph(
  fml_strata2,
  data = surv_list$dominant,
  tt = function(x, t, ...) x * sqrt(t)
)


rhs_strata3 <- paste0(
  paste(genot_vars2, collapse = "+"),
  " + ",
  paste(covars[-1], collapse = " + "),
  " + ",
  sprintf("strata(%s)", met),
  " + ",
  "tt(ITPR2_rs1049380_num)"
)
fml_strata3 <- as.formula(paste("Surv(OS, CLASS)", "~", rhs_strata3))
cox_time_pspline <- coxph(
  fml_strata3,
  data = surv_list$dominant,
  tt = function(x, t, ...) pspline(x + t / 365.25)
)


dominant_strata_models <- list(
  cox_genot_adjs$dominant$cox_strata,
  cox_time_pspline,
  cox_time_dep
)

lapply(dominant_strata_models, summary)

cox_curve_dominant <- lapply(cox_genot_adjs$dominant, function(x) {
  ggadjustedcurves(
    x,
    data = surv_list$dominant,
    variable = 'ITPR2_rs1049380',
    palette = "npg",
    curv.size = 3
  )
})


termplot(
  cox_genot_adjs$dominant$cox_strata,
  se = TRUE,
  terms = "ITPR2_rs1049380"
)

print(forestmodel::forest_model(
  cox_genot_adjs$dominant$cox_strata
))


dominant_strata_fit <- survfit(
  Surv(OS, CLASS) ~ ITPR2_rs1049380 + strata(Metastasis.at.diagnosis),
  data = surv_list$dominant
)

ggsurvplot(
  dominant_strata_fit,
  data = surv_list$dominant,
  pval = TRUE,
  palette = "npg",
  risk.table = TRUE,
  risk.table.height = 0.3,
  risk.table.y.text.col = TRUE,
  legend.labs = c(
    "G/T + G/G non-met",
    "G/T + G/G met",
    "T/T  non-met",
    "T/T met"
  )
)

ggsurvplot(
  dominant_strata_fit,
  data = surv_list$dominant,
  fun = "event",
  palette = "npg",
  risk.table = TRUE,
  risk.table.height = 0.3,
  risk.table.y.text.col = TRUE,
  legend.labs = c(
    "G/T + G/G non-met",
    "G/T + G/G met",
    "T/T  non-met",
    "T/T met"
  )
)

ggsurvplot(
  dominant_strata_fit,
  data = surv_list$dominant,
  fun = "cumhaz",
  palette = "npg",
  risk.table = TRUE,
  risk.table.height = 0.3,
  risk.table.y.text.col = TRUE,
  legend.labs = c(
    "G/T + G/G non-met",
    "G/T + G/G met",
    "T/T  non-met",
    "T/T met"
  )
)


surv_list$codominant <- surv_list$codominant %>%
  mutate(
    ITPR2_rs1049380_num = case_when(
      ITPR2_rs1049380 == "T/T" ~ 1,
      ITPR2_rs1049380 == "C/T" ~ 2,
      TRUE ~ 0
    )
  )

genot_vars2 <- c(
  "ITPR2_rs10771279",
  "DPF3_rs4903064",
  "EPAS1_rs7579899",
  "PVT1_rs35252396"
)

covars <- c(
  "Metastasis.at.diagnosis",
  "Age.group",
  'Sex'
)

met = "Metastasis.at.diagnosis"


rhs_strata3_codominant <- paste0(
  paste(genot_vars2, collapse = "+"),
  " + ",
  paste(covars[-1], collapse = " + "),
  " + ",
  sprintf("strata(%s)", met),
  " + ",
  "tt(ITPR2_rs1049380_num)"
)
fml_strata3_codominant <- as.formula(paste(
  "Surv(OS, CLASS)",
  "~",
  rhs_strata3_codominant
))
cox_time_pspline_codominant <- coxph(
  fml_strata3_codominant,
  data = surv_list$codominant,
  tt = function(x, t, ...) pspline(x + t / 365.25)
)


print(forestmodel::forest_model(
  cox_genot_adjs$codominant$cox_strata
))


surv_5y_list <- lapply(surv_list, function(x) {
  x <- x %>%
    mutate(
      OS_5y = pmin(OS, 60),
      CLASS_5y = case_when(
        OS <= 60 & CLASS == 1 ~ 1,
        OS <= 60 & CLASS == 0 & DATE2 == "01/01/2025" ~ 0,
        OS <= 60 & CLASS == 0 ~ 0,
        OS > 60 ~ 0
      )
    )
})


surv_5y_objects <- lapply(surv_list, function(x) {
  "Surv(OS_5y, CLASS_5y)"
})

cox_genot_adjs_5y <- mapply(
  cox_s,
  data = surv_5y_list,
  s = surv_5y_objects,
  SIMPLIFY = FALSE
)

coxzph_dominant_5y <- lapply(cox_genot_adjs_5y$dominant, function(x) {
  xzph <- cox.zph(x)
  plots <- ggcoxzph(xzph)
  return(list(results = xzph, plots = plots))
})

coxzph_codominant_5y <- lapply(cox_genot_adjs_5y$codominant, function(x) {
  xzph <- cox.zph(x)
  plots <- ggcoxzph(xzph)
  return(list(results = xzph, plots = plots))
})

coxzph_recessive_5y <- lapply(cox_genot_adjs_5y$recessive, function(x) {
  xzph <- cox.zph(x)
  plots <- ggcoxzph(xzph)
  return(list(results = xzph, plots = plots))
})


ggadjustedcurves(
  cox_genot_adjs_5y$dominant$cox_strata,
  data = surv_5y_list$dominant,
  variable = 'ITPR2_rs1049380',
  palette = "npg",
  curv.size = 3
)

ggadjustedcurves(
  cox_genot_adjs_5y$dominant$cox_strata,
  data = surv_5y_list$dominant,
  variable = 'PVT1_rs35252396',
  palette = "npg",
  curv.size = 3
)

termplot(
  cox_genot_adjs_5y$dominant$cox_strata,
  se = TRUE,
  terms = "ITPR2_rs1049380"
)

termplot(
  cox_genot_adjs_5y$dominant$cox_strata,
  se = TRUE,
  terms = "PVT1_rs35252396",
  col.term =
)


print(forestmodel::forest_model(
  cox_genot_adjs_5y$dominant$cox_strata
))

print(forestmodel::forest_model(
  cox_genot_adjs_5y$codominant$cox_strata,
  covariates = c(
    "ITPR2_rs1049380",
    "ITPR2_rs10771279",
    "EPAS1_rs7579899",
    "PVT1_rs35252396",
    "Metastasis.at.diagnosis",
    "Age.group"
  )
))

dominant_5y_strata_fit_104 <- survfit(
  Surv(OS_5y, CLASS_5y) ~ ITPR2_rs1049380 + strata(Metastasis.at.diagnosis),
  data = surv_5y_list$dominant
)

dominant_5y_strata_fit_352 <- survfit(
  Surv(OS_5y, CLASS_5y) ~ PVT1_rs35252396 + strata(Metastasis.at.diagnosis),
  data = surv_5y_list$dominant
)


ggsurvplot(
  dominant_5y_strata_fit_104,
  data = surv_5y_list$dominant,
  pval = TRUE,
  palette = "npg",
  risk.table = TRUE,
  risk.table.height = 0.3,
  risk.table.y.text.col = TRUE,
  legend.labs = c(
    "G/T + G/G non-met",
    "G/T + G/G met",
    "T/T  non-met",
    "T/T met"
  )
)

ggsurvplot(
  dominant_5y_strata_fit_104,
  data = surv_5y_list$dominant,
  fun = "event",
  palette = "npg",
  risk.table.height = 0.3,
  risk.table.y.text.col = TRUE,
  legend.labs = c(
    "G/T + G/G non-met",
    "G/T + G/G met",
    "T/T  non-met",
    "T/T met"
  )
)

ggsurvplot(
  dominant_5y_strata_fit_104,
  data = surv_5y_list$dominant,
  fun = "cumhaz",
  palette = "npg",
  risk.table.height = 0.3,
  risk.table.y.text.col = TRUE,
  legend.labs = c(
    "G/T + G/G non-met",
    "G/T + G/G met",
    "T/T  non-met",
    "T/T met"
  )
)


ggsurvplot(
  dominant_5y_strata_fit_352,
  data = surv_5y_list$dominant,
  pval = TRUE,
  palette = "npg",
  risk.table = TRUE,
  risk.table.height = 0.3,
  risk.table.y.text.col = TRUE,
  legend.labs = c(
    "AC/CG + AC/AC non-met",
    "AC/CG + AC/AC met",
    "CG/CG non-met",
    "CG/CG met"
  )
)

ggsurvplot(
  dominant_5y_strata_fit_352,
  data = surv_5y_list$dominant,
  fun = "event",
  palette = "npg",
  risk.table.height = 0.3,
  risk.table.y.text.col = TRUE,
  legend.labs = c(
    "AC/CG + AC/AC non-met",
    "AC/CG + AC/AC met",
    "CG/CG non-met",
    "CG/CG met"
  )
)

ggsurvplot(
  dominant_5y_strata_fit_352,
  data = surv_5y_list$dominant,
  fun = "cumhaz",
  palette = "npg",
  risk.table.height = 0.3,
  risk.table.y.text.col = TRUE,
  legend.labs = c(
    "AC/CG + AC/AC non-met",
    "AC/CG + AC/AC met",
    "CG/CG non-met",
    "CG/CG met"
  )
)


ggadjustedcurves(
  cox_genot_adjs_5y$codominant$cox_strata,
  data = surv_5y_list$codominant,
  variable = 'ITPR2_rs1049380',
  palette = "npg",
  curv.size = 3
)

termplot(
  cox_genot_adjs_5y$codominant$cox_strata,
  se = TRUE,
  terms = "ITPR2_rs1049380"
)

print(forestmodel::forest_model(
  cox_genot_adjs_5y$codominant$cox_strata,
  covariates = c(
    "ITPR2_rs1049380",
    "ITPR2_rs10771279",
    "EPAS1_rs7579899",
    "PVT1_rs35252396",
    "Metastasis.at.diagnosis",
    "Age.group"
  )
))


codominant_5y_strata_fit_104 <- survfit(
  Surv(OS_5y, CLASS_5y) ~ ITPR2_rs1049380 + strata(Metastasis.at.diagnosis),
  data = surv_5y_list$codominant
)

ggsurvplot(
  codominant_5y_strata_fit_104,
  data = surv_5y_list$codominant,
  pval = TRUE,
  palette = "npg",
  risk.table = TRUE,
  risk.table.height = 0.3,
  risk.table.y.text.col = TRUE,
  legend.labs = c(
    "G/G met",
    "G/G non-met",
    "G/T non-met",
    "G/T met",
    "T/T non-met",
    "T/T met"
  )
)

ggsurvplot(
  codominant_5y_strata_fit_104,
  data = surv_5y_list$codominant,
  fun = "event",
  palette = "npg",
  risk.table.height = 0.3,
  risk.table.y.text.col = TRUE,
  legend.labs = c(
    "G/G met",
    "G/G non-met",
    "G/T non-met",
    "G/T met",
    "T/T non-met",
    "T/T met"
  )
)

ggsurvplot(
  codominant_5y_strata_fit_104,
  data = surv_5y_list$codominant,
  fun = "cumhaz",
  palette = "npg",
  risk.table.height = 0.3,
  risk.table.y.text.col = TRUE,
  legend.labs = c(
    "G/G met",
    "G/G non-met",
    "G/T non-met",
    "G/T met",
    "T/T non-met",
    "T/T met"
  )
)


write.xlsx(
  coxzph_dominant_5y$dominant$results$table,
  file = "5yr_coxzph_dominant_5y.xlsx",
  col.names = TRUE,
  row.names = TRUE,
  append = FALSE
)

write.xlsx(
  coxzph_codominant_5y$codominant$results$table,
  file = "5yr_coxzph_codominant_5y.xlsx",
  col.names = TRUE,
  row.names = TRUE,
  append = FALSE
)

write.xlsx(
  coxzph_recessive_5y$recessive$results$table,
  file = "5yr_coxzph_recessive_5y.xlsx",
  col.names = TRUE,
  row.names = TRUE,
  append = FALSE
)

write.xlsx(
  lapply(cox_strata_list, broom::tidy),
  file = "5yr_cox_strata_list.xlsx",
  col.names = TRUE,
  row.names = TRUE,
  append = FALSE
)


write.xlsx(
  coxzph_dominant$dominant$results$table,
  file = "global_coxzph_dominant.xlsx",
  col.names = TRUE,
  row.names = TRUE,
  append = FALSE
)

write.xlsx(
  coxzph_codominant$codominant$results$table,
  file = "global_coxzph_codominant.xlsx",
  col.names = TRUE,
  row.names = TRUE,
  append = FALSE
)

write.xlsx(
  coxzph_recessive$recessive$results$table,
  file = "global_coxzph_recessive.xlsx",
  col.names = TRUE,
  row.names = TRUE,
  append = FALSE
)

write.xlsx(
  lapply(cox_strata_list, broom::tidy),
  file = "5yr_cox_strata_list.xlsx",
  col.names = TRUE,
  row.names = TRUE,
  append = FALSE
)
