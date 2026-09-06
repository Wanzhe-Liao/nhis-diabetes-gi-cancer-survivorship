# rev1_p1_main.R —— 修回 P1 主分析（确定性部分）
# 依据 revision/analysis_spec_v1.md（v1.1：cholesterol 移出嵌套主序列；
# v4 中 58.5% 为构念性缺失）
# R 4.6.1 / 22-year primary output: outputs/revision_round1_v4_r461_sens2007/
options(survey.lonely.psu = "fail")
source(file.path("scripts", "rev1_runtime.R"))
source(file.path("scripts", "rev1_domain_design.R"))

suppressPackageStartupMessages({
  library(survey); library(survival)
})

set.seed(20260829)
SUPPORT_OUT <- rev1_support_dir()
OUT <- rev1_output_dir()
COHORT_PATH <- rev1_cohort_path()
log_con <- file(file.path(OUT, "run_log.txt"), open = "wt")
sink(log_con, split = TRUE)
cat("rev1_p1_main.R |", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("cohort:", COHORT_PATH, "\n")
cat("support:", SUPPORT_OUT, "\n")
cat("output:", OUT, "\n")
print(rev1_runtime_info(), row.names = FALSE)

## ---------- 0. 数据准备 ----------
cohort <- read_rev1_parquet(COHORT_PATH)
POVERTY_PATH <- file.path(SUPPORT_OUT, "income_repair", "poverty_repaired.csv")
pov <- read.csv(POVERTY_PATH,
                stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN"))
cohort$year <- as.integer(cohort$year)
cohort$hhx_k <- as.character(as.integer(cohort$hhx))
cohort$fmx_k <- as.character(as.integer(cohort$fmx))
pov$hhx_k <- as.character(as.integer(pov$hhx))
pov$fmx_k <- as.character(as.integer(pov$fmx))
pov$year <- as.integer(pov$year)
cohort <- rev1_left_join(
  cohort,
  pov[, c("year", "hhx_k", "fmx_k", "poverty_3cat_rep", "poverty_4cat_rep"), drop = FALSE],
  by = c("year", "hhx_k", "fmx_k")
)
cat("收入修复合并缺失率:", mean(is.na(cohort$poverty_3cat_rep)), "\n")

# 生存变量（沿用 round4 规则）
cohort$time_months <- pmin(cohort$followup_years * 12, 120)
cohort$event <- as.numeric(cohort$mortstat == 1 & cohort$followup_years <= 10)
# Domain-correct revision: 行过滤推迟到派生完成之后（保留 zero-time 行供 G 敏感性）。
cohort <- rev1_period_identifiers(cohort)

# 新暴露层级
cohort$trajectory_5cat_rev <- as.character(cohort$trajectory_6cat)
cohort$trajectory_5cat_rev[cohort$trajectory_5cat_rev %in% c("dm_to_gi_2_10y", "dm_to_gi_gt10y")] <- "established_pre_cancer_dm"
cohort$trajectory_5cat_rev[cohort$trajectory_5cat_rev == "gi_to_dm"] <- "post_cancer_dm"
cohort$trajectory_5cat_rev[cohort$trajectory_5cat_rev == "gi_only"] <- "no_diabetes"

# 时间间隔
cohort$cancer_to_interview <- cohort$age - cohort$gi_first_dx_age
cohort$dm_to_interview <- cohort$age - cohort$dm_dx_age

# 嵌套模型附加变量清洗
cohort$hypertension_b <- as.logical(cohort$hypertension_ever)
cohort$chd_b <- as.logical(cohort$chd_ever)
cohort$stroke_b <- as.logical(cohort$stroke_ever)
cohort$phys_inactive <- ifelse(is.na(cohort$phys_active_any), NA, !as.logical(cohort$phys_active_any))
cohort$obesity <- ifelse(is.na(cohort$bmi), NA, cohort$bmi >= 30)
cohort$srh_num <- suppressWarnings(as.numeric(cohort$srh))
cohort$srh_f <- ifelse(cohort$srh_num %in% 1:5, cohort$srh_num, NA)
cohort$srh_fairpoor <- ifelse(is.na(cohort$srh_f), NA, cohort$srh_f >= 4)
cohort$smoking_cf <- ifelse(is.na(cohort$smoking_3cat), NA, cohort$smoking_3cat %in% c("current", "former"))
cohort$uninsured <- ifelse(
  is.na(cohort$insurance_type), NA,
  as.character(cohort$insurance_type) == "uninsured"
)
cohort$rx_cost_barrier <- ifelse(
  is.na(cohort$cost_barrier_rx), NA,
  as.logical(cohort$cost_barrier_rx)
)

make_factor <- function(x, ref) {
  x <- as.character(x)
  x[x == ""] <- NA
  x <- factor(x)  # 显式无序，规避 parquet 元数据可能带入的 ordered factor
  if (ref %in% levels(x)) x <- stats::relevel(x, ref = ref)
  x
}
cohort$trajectory_6cat <- make_factor(cohort$trajectory_6cat, "gi_only")
cohort$trajectory_5cat_rev <- make_factor(cohort$trajectory_5cat_rev, "no_diabetes")
cohort$sex <- make_factor(cohort$sex, "1")
cohort$race <- make_factor(cohort$race, "1")
cohort$region <- make_factor(cohort$region, "1")
cohort$smoking_3cat <- make_factor(cohort$smoking_3cat, "current")
cohort$education_4cat <- make_factor(cohort$education_4cat, "college_grad")
cohort$poverty_3cat <- make_factor(cohort$poverty_3cat, "2_0_to_3_99")          # 旧（错码）
cohort$poverty_3cat_rep <- make_factor(cohort$poverty_3cat_rep, "2_0_to_3_99") # 修复
# 修复收入缺失保留为显式类别（与原投稿口径一致，避免 17.6% 样本损失；嵌套模型另行 complete-case）
cohort$poverty_3cat_rep <- factor(cohort$poverty_3cat_rep, exclude = NULL)
levels(cohort$poverty_3cat_rep)[is.na(levels(cohort$poverty_3cat_rep))] <- "missing"
cohort$poverty_3cat_rep <- stats::relevel(cohort$poverty_3cat_rep, ref = "2_0_to_3_99")
cohort$survey_year <- make_factor(cohort$survey_year, "1997")
cohort$insurance_type <- make_factor(cohort$insurance_type, "private_only")
cohort$alcohol_status <- make_factor(cohort$alcohol_status, "never")
cohort$srh_f <- factor(cohort$srh_f)

# zero-time 敏感性时间列（G 节）
cohort$time_z05 <- ifelse(cohort$time_months == 0, 0.5, cohort$time_months)
cohort$time_z1  <- ifelse(cohort$time_months == 0, 1.0, cohort$time_months)
cohort$time_z3  <- ifelse(cohort$time_months == 0, 3.0, cohort$time_months)

SITE_FLAGS <- c("colon_flag", "esoph_flag", "gallbladder_flag", "liver_flag",
                "pancreas_flag", "rectum_flag", "stomach_flag")
for (v in SITE_FLAGS) cohort[[v]] <- as.logical(cohort[[v]])

# Domain-correct: 先在全部 ELIGSTAT=1 pooled Sample Adults 上建 design 再 subset。
# GI-only 建库时 75.9% 层为单 PSU 层, lonely.psu=remove 会丢失其方差贡献。
source(file.path("scripts", "rev1_domain_design.R"))
cohort_all <- cohort
flag_pos <- !is.na(cohort$time_months) & cohort$time_months > 0 &
  !is.na(cohort$design_psu) & !is.na(cohort$design_strata) & !is.na(cohort$sa_wgt_pool)
cohort <- cohort[flag_pos, , drop = FALSE]
combined <- rev1_attach_fullsample(cohort_all)
design_pop <- rev1_domain_population_design(combined)
cat("全样本 design: rows =", nrow(combined), " df =", survey::degf(design_pop), "\n")
flag_dom <- rev1_flag_dom(combined)
cat("domain rows (time>0):", sum(flag_dom), "\n")
design_full <- subset(design_pop, flag_dom)
design_allgi <- subset(design_pop, combined$gi_any %in% TRUE & !is.na(combined$time_months))

extract_row <- function(fit, term, label, model_name) {
  s <- summary(fit)$coefficients
  if (!term %in% rownames(s)) return(NULL)
  r <- s[term, ]
  data.frame(model = model_name, term = label, logHR = unname(r["coef"]),
             HR = exp(unname(r["coef"])), robust_se = unname(r["robust se"]),
             ci_lo = exp(unname(r["coef"] - 1.96 * r["robust se"])),
             ci_hi = exp(unname(r["coef"] + 1.96 * r["robust se"])),
             p = unname(r["Pr(>|z|)"]), n = length(fit$y[, "status"]),
             events = sum(fit$y[, "status"]), stringsAsFactors = FALSE)
}

master <- list()

## ---------- A. Principal：合并 established pre-cancer DM（修复收入） ----------
cat("\n== A. Principal 合并模型 ==\n")
fitA <- svycoxph(
  Surv(time_months, event) ~ trajectory_5cat_rev + age + sex + race + region + bmi +
    smoking_3cat + survey_year + colon_flag + esoph_flag + gallbladder_flag + liver_flag +
    pancreas_flag + rectum_flag + stomach_flag + education_4cat + poverty_3cat_rep,
  design = design_full)
print(summary(fitA)$coefficients[paste0("trajectory_5cat_rev",
  c("established_pre_cancer_dm", "peri_diagnostic", "post_cancer_dm", "dm_order_unknown")), ])
for (tr in c("established_pre_cancer_dm", "peri_diagnostic", "post_cancer_dm", "dm_order_unknown")) {
  master[[paste0("A_", tr)]] <- extract_row(fitA, paste0("trajectory_5cat_rev", tr), tr, "A_principal_5cat")
}

# Secondary year-coverage comparison uses the same period/domain design and
# covariates as the primary model, with annual weights divided by 21 years.
combined21 <- droplevels(combined[combined$year != 2007L, , drop = FALSE])
combined21$sa_wgt_pool <- combined21$sa_wgt_pool * 22 / 21
design21 <- rev1_domain_population_design(combined21)
flag21 <- rev1_flag_dom(combined21)
fit21 <- svycoxph(formula(fitA), design = subset(design21, flag21))
year_coverage <- rev1_bind_rows(list(
  extract_row(fitA, "trajectory_5cat_revestablished_pre_cancer_dm",
              "established_pre_cancer_dm", "1997-2018 (primary)"),
  extract_row(fit21, "trajectory_5cat_revestablished_pre_cancer_dm",
              "established_pre_cancer_dm", "1997-2006, 2008-2018")
))
year_coverage$population_n <- c(nrow(combined), nrow(combined21))
year_coverage$cohort_n <- c(sum(combined$gi_any), sum(combined21$gi_any))
write_rev1_csv(year_coverage, file.path(OUT, "year_coverage_sensitivity.csv"))

## ---------- A2. 固定主模型样本的 no-income 敏感性 ----------
cat("\n== A2. No-income sensitivity on principal model sample ==\n")
design_principal <- fitA$survey.design
if (!is.null(fitA$na.action)) {
  design_principal <- design_principal[-as.integer(fitA$na.action), ]
}
fit_no_income <- svycoxph(
  Surv(time_months, event) ~ trajectory_5cat_rev + age + sex + race + region + bmi +
    smoking_3cat + survey_year + colon_flag + esoph_flag + gallbladder_flag + liver_flag +
    pancreas_flag + rectum_flag + stomach_flag + education_4cat,
  design = design_principal)
if (length(fit_no_income$y[, "status"]) != length(fitA$y[, "status"]) ||
    sum(fit_no_income$y[, "status"]) != sum(fitA$y[, "status"])) {
  stop("No-income sensitivity did not retain the fixed principal sample.", call. = FALSE)
}
master$A2_no_income <- extract_row(
  fit_no_income,
  "trajectory_5cat_revestablished_pre_cancer_dm",
  "established_pre_cancer_dm",
  "A2_no_income_sensitivity"
)
print(master$A2_no_income)

## ---------- G. zero-time 约定敏感性（S11；domain 全样本） ----------
cat("\n== G. zero-time conventions ==\n")
for (zc in c("time_z05", "time_z1", "time_z3")) {
  f <- as.formula(paste(
    "Surv(", zc, ", event) ~ trajectory_5cat_rev + age + sex + race + region + bmi +",
    "smoking_3cat + survey_year + colon_flag + esoph_flag + gallbladder_flag + liver_flag +",
    "pancreas_flag + rectum_flag + stomach_flag + education_4cat + poverty_3cat_rep"))
  fit_z <- svycoxph(f, design = design_allgi)
  master[[paste0("G_", zc)]] <- extract_row(
    fit_z, "trajectory_5cat_revestablished_pre_cancer_dm",
    "established_pre_cancer_dm", paste0("G_zero_time_", zc))
  cat(zc, "HR =", round(master[[paste0("G_", zc)]]$HR, 4),
      " n =", master[[paste0("G_", zc)]]$n,
      " events =", master[[paste0("G_", zc)]]$events, "\n")
}

## ---------- B. 次要：六分类分解 + Wald β_2-10 = β_>10 ----------
cat("\n== B. 滞后异质性 ==\n")
fitB <- svycoxph(
  Surv(time_months, event) ~ trajectory_6cat + age + sex + race + region + bmi +
    smoking_3cat + survey_year + colon_flag + esoph_flag + gallbladder_flag + liver_flag +
    pancreas_flag + rectum_flag + stomach_flag + education_4cat + poverty_3cat_rep,
  design = design_full)
b1 <- coef(fitB)["trajectory_6catdm_to_gi_2_10y"]
b2 <- coef(fitB)["trajectory_6catdm_to_gi_gt10y"]
V <- vcov(fitB)
v_diff <- V["trajectory_6catdm_to_gi_2_10y", "trajectory_6catdm_to_gi_2_10y"] +
          V["trajectory_6catdm_to_gi_gt10y", "trajectory_6catdm_to_gi_gt10y"] -
          2 * V["trajectory_6catdm_to_gi_2_10y", "trajectory_6catdm_to_gi_gt10y"]
se_diff <- sqrt(v_diff)
wald_chi2 <- (b1 - b2)^2 / v_diff
p_het <- pchisq(wald_chi2, df = 1, lower.tail = FALSE)
ratio_hr <- exp(b1 - b2)
master$wald_lag <- data.frame(model = "B_lag_heterogeneity", term = "2-10y vs >10y",
  beta_2_10 = unname(b1), beta_gt10 = unname(b2), diff_loghr = unname(b1 - b2),
  se_diff = se_diff, wald_chi2 = unname(wald_chi2), p_heterogeneity = p_het,
  ratio_of_HRs = ratio_hr, ratio_ci_lo = exp(b1 - b2 - 1.96 * se_diff),
  ratio_ci_hi = exp(b1 - b2 + 1.96 * se_diff))
cat("Wald chi2 =", wald_chi2, " P_het =", p_het, " ratio of HRs =", ratio_hr, "\n")
for (tr in c("dm_to_gi_2_10y", "dm_to_gi_gt10y", "peri_diagnostic", "gi_to_dm", "dm_order_unknown")) {
  master[[paste0("B_", tr)]] <- extract_row(fitB, paste0("trajectory_6cat", tr), tr, "B_6cat_repaired_income")
}

## ---------- C. 时间间隔结构（R1-3） ----------
cat("\n== C. 时间间隔分布 ==\n")
probs <- c(0.10, 0.25, 0.5, 0.75, 0.90)
timing_rows <- list()
for (tr in levels(cohort$trajectory_6cat)) {
  sub <- subset(design_full, trajectory_6cat == tr)
  q <- as.numeric(svyquantile(~cancer_to_interview, sub, quantiles = probs, ci = FALSE, na.rm = TRUE)[[1]][1, ])
  m <- as.numeric(svymean(~cancer_to_interview, sub, na.rm = TRUE))
  timing_rows[[length(timing_rows) + 1]] <- data.frame(
    trajectory = tr, variable = "cancer_to_interview", n = nrow(sub$variables),
    mean = m[1], p10 = q[1], p25 = q[2], median = q[3], p75 = q[4], p90 = q[5])
  if (tr != "gi_only") {
    q2 <- as.numeric(svyquantile(~dm_to_interview, sub, quantiles = probs, ci = FALSE, na.rm = TRUE)[[1]][1, ])
    m2 <- as.numeric(svymean(~dm_to_interview, sub, na.rm = TRUE))
    timing_rows[[length(timing_rows) + 1]] <- data.frame(
      trajectory = tr, variable = "dm_to_interview",
      n = sum(!is.na(sub$variables$dm_to_interview)),
      mean = m2[1], p10 = q2[1], p25 = q2[2], median = q2[3], p75 = q2[4], p90 = q2[5])
  }
}
timing <- rev1_bind_rows(timing_rows)
write_rev1_csv(timing, file.path(OUT, "timing_distributions.csv"))
print(timing)

# cancer-to-interview 调整敏感性（修复收入版）
fitC <- svycoxph(
  Surv(time_months, event) ~ trajectory_5cat_rev + cancer_to_interview + age + sex + race + region + bmi +
    smoking_3cat + survey_year + colon_flag + esoph_flag + gallbladder_flag + liver_flag +
    pancreas_flag + rectum_flag + stomach_flag + education_4cat + poverty_3cat_rep,
  design = design_full)
master$C_interval_adj <- extract_row(fitC, "trajectory_5cat_revestablished_pre_cancer_dm",
                                     "established_pre_cancer_dm", "C_plus_cancer_to_interview")

## ---------- D. T1D 代理排除 ----------
cat("\n== D. T1D 代理排除 ==\n")
excl_index <- cohort$trajectory_6cat != "gi_only" &
  !is.na(cohort$dm_dx_age) & cohort$dm_dx_age < 30
excl_table <- table(as.character(cohort$trajectory_6cat[excl_index]))
excl_by_traj <- data.frame(
  trajectory_6cat = names(excl_table),
  excluded_lt30 = as.integer(excl_table),
  stringsAsFactors = FALSE
)
print(excl_by_traj)
write_rev1_csv(excl_by_traj, file.path(OUT, "t1d_proxy_exclusions.csv"))
design_t1d <- subset(design_full, !(trajectory_6cat != "gi_only" & !is.na(dm_dx_age) & dm_dx_age < 30))
fitD <- svycoxph(
  Surv(time_months, event) ~ trajectory_5cat_rev + age + sex + race + region + bmi +
    smoking_3cat + survey_year + colon_flag + esoph_flag + gallbladder_flag + liver_flag +
    pancreas_flag + rectum_flag + stomach_flag + education_4cat + poverty_3cat_rep,
  design = design_t1d)
for (tr in c("established_pre_cancer_dm", "peri_diagnostic", "post_cancer_dm", "dm_order_unknown")) {
  master[[paste0("D_", tr)]] <- extract_row(fitD, paste0("trajectory_5cat_rev", tr), tr, "D_t1d_proxy_exclusion")
}

## ---------- E. 嵌套解释模型（固定 complete-case 样本） ----------
cat("\n== E. 嵌套模型 ==\n")
nested_vars <- c("age", "sex", "race", "region", "survey_year", SITE_FLAGS,
                 "cancer_to_interview", "bmi", "smoking_3cat", "education_4cat", "poverty_3cat_rep",
                 "hypertension_b", "chd_b", "stroke_b",
                 "phys_inactive", "alcohol_status", "insurance_type", "srh_f")
cc <- complete.cases(cohort[, c("trajectory_5cat_rev", "time_months", "event", nested_vars)])
cohort_cc <- cohort[cc, ]
cat("M3B complete-case n =", nrow(cohort_cc), " events =", sum(cohort_cc$event), "\n")
cat("暴露组构成:\n"); print(table(cohort_cc$trajectory_5cat_rev))
design_cc <- subset(design_pop, flag_dom & complete.cases(
  combined[, c("trajectory_5cat_rev", "time_months", "event", nested_vars)]))

TERM <- "trajectory_5cat_revestablished_pre_cancer_dm"
rhs_M0  <- "trajectory_5cat_rev + age + sex + race + region + survey_year + colon_flag + esoph_flag + gallbladder_flag + liver_flag + pancreas_flag + rectum_flag + stomach_flag"
rhs_M0T <- paste(rhs_M0, "+ cancer_to_interview")
rhs_M1  <- paste(rhs_M0, "+ bmi + smoking_3cat + education_4cat + poverty_3cat_rep")
rhs_M2  <- paste(rhs_M1, "+ hypertension_b + chd_b + stroke_b")
rhs_M3A <- paste(rhs_M2, "+ phys_inactive + alcohol_status + insurance_type")
rhs_M3B <- paste(rhs_M3A, "+ srh_f")
fits_nested <- list()
for (nm in c("M0", "M0T", "M1", "M2", "M3A", "M3B")) {
  f <- as.formula(paste("Surv(time_months, event) ~", get(paste0("rhs_", nm))))
  fits_nested[[nm]] <- svycoxph(f, design = design_cc)
  master[[paste0("E_", nm)]] <- extract_row(fits_nested[[nm]], TERM, "established_pre_cancer_dm", paste0("E_nested_", nm))
}
# 样本限制检查：主模型在自然样本 vs M3B-complete 样本
fitM1_natural <- fitA
cat("M1 自然样本 HR =", exp(coef(fitM1_natural)[TERM]),
    " | M1 固定样本 HR =", exp(coef(fits_nested$M1)[TERM]), "\n")
master$E_sample_check <- data.frame(model = "E_sample_restriction_check", term = "established_pre_cancer_dm",
  HR_natural_sample = exp(coef(fitM1_natural)[TERM]), HR_M3B_complete_sample = exp(coef(fits_nested$M1)[TERM]),
  n_natural = length(fitM1_natural$y[, "status"]), n_complete = nrow(cohort_cc))

## ---------- F. 标准化负担比较（R2-11 第二证据链） ----------
cat("\n== F. 标准化负担比较 ==\n")
burden_vars <- c("hypertension_b", "chd_b", "stroke_b", "obesity", "smoking_cf",
                 "phys_inactive", "srh_fairpoor", "uninsured", "rx_cost_barrier")
comparison_groups <- c("no_diabetes", "established_pre_cancer_dm")

# Model-based predictive margins standardize each binary indicator over the
# outcome-specific complete-case survey population.  The adjustment model
# includes the prespecified age, sex, race, cancer-site, and survey-year set.
# Individual-level averaging avoids sparse direct-standardization cells; the
# survey::svypredmeans variance combines Taylor model and target-distribution
# components.  Prevalence CIs use a logit-delta transform; RD and SMD CIs use
# first-order delta-method normal intervals.
cohort_b <- cohort
for (v in burden_vars) {
  cohort_b[[paste0(v, "_num")]] <- ifelse(
    is.na(cohort_b[[v]]), NA_real_, as.numeric(as.logical(cohort_b[[v]]))
  )
}
# Domain-correct: burden design 同样构建于全样本再 subset。
combined_b <- combined
for (v in burden_vars) {
  combined_b[[paste0(v, "_num")]] <- ifelse(
    is.na(combined_b[[v]]), NA_real_, as.numeric(as.logical(combined_b[[v]]))
  )
}
design_b <- subset(rev1_domain_population_design(combined_b), flag_dom)
burden_adjustment_rhs <- paste(
  "splines::ns(age, df = 4) + sex + race +",
  paste(SITE_FLAGS, collapse = " + "),
  "+ survey_year"
)

logit_delta_ci <- function(estimate, se, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  if (!is.finite(estimate) || !is.finite(se) || estimate <= 0 || estimate >= 1) {
    return(pmin(1, pmax(0, estimate + c(-1, 1) * z * se)))
  }
  se_logit <- se / (estimate * (1 - estimate))
  stats::plogis(stats::qlogis(estimate) + c(-1, 1) * z * se_logit)
}

normal_delta_ci <- function(estimate, se, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  estimate + c(-1, 1) * z * se
}

standardize_binary_burden <- function(variable) {
  outcome <- paste0(variable, "_num")
  formula_adjust <- stats::as.formula(
    paste(outcome, "~", burden_adjustment_rhs)
  )
  adjust_fit <- survey::svyglm(
    formula_adjust, design = design_b, family = stats::quasibinomial()
  )
  final_fit <- update(
    adjust_fit, . ~ . + trajectory_5cat_rev,
    design = adjust_fit$survey.design
  )
  if (!isTRUE(adjust_fit$converged) || !isTRUE(final_fit$converged) ||
      anyNA(stats::coef(final_fit))) {
    stop("Burden standardization model failed for ", variable, call. = FALSE)
  }

  margins <- survey::svypredmeans(
    adjust_fit, ~trajectory_5cat_rev,
    # Keep all factor levels. Passing a character vector makes svypredmeans()
    # replace the prediction column with a one-level character and causes an
    # invalid contrasts matrix under R 4.6.1.
    predictat = factor(
      comparison_groups, levels = levels(cohort_b$trajectory_5cat_rev)
    )
  )
  estimates <- stats::coef(margins)[comparison_groups]
  ses <- survey::SE(margins)[comparison_groups]
  if (any(!is.finite(estimates)) || any(!is.finite(ses))) {
    stop("Non-finite predictive margin for ", variable, call. = FALSE)
  }

  rd_contrast <- stats::setNames(c(-1, 1), comparison_groups)
  rd <- survey::svycontrast(margins, rd_contrast)
  smd <- survey::svycontrast(
    margins,
    quote(
      (established_pre_cancer_dm - no_diabetes) /
        sqrt(
          (established_pre_cancer_dm * (1 - established_pre_cancer_dm) +
             no_diabetes * (1 - no_diabetes)) / 2
        )
    )
  )
  rd_est <- as.numeric(stats::coef(rd))
  rd_se <- as.numeric(survey::SE(rd))
  smd_est <- as.numeric(stats::coef(smd))
  smd_se <- as.numeric(survey::SE(smd))
  p0_ci <- logit_delta_ci(estimates[["no_diabetes"]], ses[["no_diabetes"]])
  p1_ci <- logit_delta_ci(
    estimates[["established_pre_cancer_dm"]],
    ses[["established_pre_cancer_dm"]]
  )
  rd_ci <- normal_delta_ci(rd_est, rd_se)
  smd_ci <- normal_delta_ci(smd_est, smd_se)

  data.frame(
    variable = variable,
    no_diabetes = estimates[["no_diabetes"]],
    no_diabetes_se = ses[["no_diabetes"]],
    no_diabetes_ci_lo = p0_ci[1], no_diabetes_ci_hi = p0_ci[2],
    established_pre_cancer_dm = estimates[["established_pre_cancer_dm"]],
    established_pre_cancer_dm_se = ses[["established_pre_cancer_dm"]],
    established_pre_cancer_dm_ci_lo = p1_ci[1],
    established_pre_cancer_dm_ci_hi = p1_ci[2],
    diff = rd_est, diff_se = rd_se, diff_ci_lo = rd_ci[1], diff_ci_hi = rd_ci[2],
    smd = smd_est, smd_se = smd_se, smd_ci_lo = smd_ci[1], smd_ci_hi = smd_ci[2],
    model_n = nrow(adjust_fit$survey.design$variables),
    model_positive = sum(stats::model.response(stats::model.frame(adjust_fit))),
    model_rank = final_fit$rank,
    model_coefficients = length(stats::coef(final_fit)),
    design_df = survey::degf(adjust_fit$survey.design),
    model_converged = isTRUE(final_fit$converged),
    standardization_target = "outcome-specific complete-case full five-state survey population",
    adjustment = burden_adjustment_rhs,
    prevalence_ci_method = "logit-delta Taylor 95% CI",
    contrast_ci_method = "first-order delta normal 95% CI",
    stringsAsFactors = FALSE
  )
}

burden_rows <- lapply(burden_vars, standardize_binary_burden)
burden_wide <- rev1_bind_rows(burden_rows)
write_rev1_csv(burden_wide, file.path(OUT, "burden_comparison.csv"))
print(burden_wide)

# Descriptive distribution of a four-item measured cardiometabolic-condition
# count (hypertension, CHD, stroke, obesity). Diabetes itself is intentionally
# not counted because it defines the comparison groups; cholesterol is excluded
# because the harmonized ever-diagnosed construct is unavailable for 58.5% of
# v4 records and is not comparable in 2005/2013.
count_components <- c("hypertension_b_num", "chd_b_num", "stroke_b_num", "obesity_num")
# Domain-correct: 计数列同时在 domain 全样本框架上派生。
count_frames <- list(cohort_b = cohort_b, combined_c = combined_b)
for (fn in names(count_frames)) {
  fr <- count_frames[[fn]]
  fr_complete <- stats::complete.cases(fr[, count_components, drop = FALSE])
  fr$cardiometabolic_condition_count <- NA_integer_
  fr$cardiometabolic_condition_count[fr_complete] <- rowSums(
    fr[fr_complete, count_components, drop = FALSE]
  )
  fr$cardiometabolic_count_cat <- ifelse(
    is.na(fr$cardiometabolic_condition_count), NA_character_,
    ifelse(
      fr$cardiometabolic_condition_count == 0, "0",
      ifelse(fr$cardiometabolic_condition_count == 1, "1", ">=2")
    )
  )
  for (level in c("0", "1", ">=2")) {
    safe_level <- if (level == ">=2") "ge2" else level
    fr[[paste0("count_cat_", safe_level)]] <- ifelse(
      is.na(fr$cardiometabolic_count_cat), NA_real_,
      as.numeric(fr$cardiometabolic_count_cat == level)
    )
  }
  count_frames[[fn]] <- fr
}
cohort_b <- count_frames$cohort_b
combined_c <- count_frames$combined_c
design_count <- subset(rev1_domain_population_design(combined_c), flag_dom)
count_rows <- list()
for (grp in comparison_groups) {
  sub <- subset(
    design_count,
    trajectory_5cat_rev == grp & !is.na(cardiometabolic_count_cat)
  )
  for (level in c("0", "1", ">=2")) {
    safe_level <- if (level == ">=2") "ge2" else level
    stat <- svymean(
      stats::as.formula(paste0("~count_cat_", safe_level)), sub,
      na.rm = TRUE
    )
    estimate <- as.numeric(stats::coef(stat))
    se <- as.numeric(survey::SE(stat))
    ci <- logit_delta_ci(estimate, se)
    count_rows[[length(count_rows) + 1L]] <- data.frame(
      group = grp, category = level,
      unweighted_n = nrow(sub$variables),
      weighted_n = sum(stats::weights(sub, type = "sampling")),
      prevalence = estimate, se = se, ci_lo = ci[1], ci_hi = ci[2],
      components = "hypertension; CHD; stroke; obesity",
      method = "survey-weighted descriptive prevalence; logit-delta Taylor 95% CI",
      stringsAsFactors = FALSE
    )
  }
}
count_distribution <- rev1_bind_rows(count_rows)
write_rev1_csv(
  count_distribution,
  file.path(OUT, "cardiometabolic_condition_count.csv")
)
print(count_distribution)

cholesterol_missing <- is.na(cohort$cholesterol_high_ever) |
  as.character(cohort$cholesterol_high_ever) == "missing"
construct_audit <- data.frame(
  variable = "cholesterol_high_ever",
  n_total = nrow(cohort),
  n_construct_available = sum(!cholesterol_missing),
  n_construct_missing = sum(cholesterol_missing),
  missing_fraction = mean(cholesterol_missing),
  available_survey_years = paste(
    sort(unique(cohort$year[!cholesterol_missing])), collapse = ";"
  ),
  decision = "excluded_from_nested_M2_and_burden_standardization",
  reason = paste(
    "Cross-year ever-diagnosed construct unavailable for most years;",
    "2005 CHOLEST is medication use and 2013 CHLYR1 is past-12-month status."
  ),
  stringsAsFactors = FALSE
)
write_rev1_csv(
  construct_audit, file.path(OUT, "burden_construct_audit.csv")
)

## ---------- 汇总输出 ----------
master_df <- rev1_bind_rows(master[sapply(master, function(x) is.data.frame(x) && "HR" %in% names(x))])
wald_df <- master$wald_lag
sample_df <- master$E_sample_check
write_rev1_csv(master_df, file.path(OUT, "master_results.csv"))
write_rev1_csv(wald_df, file.path(OUT, "wald_lag_test.csv"))
write_rev1_csv(sample_df, file.path(OUT, "sample_restriction_check.csv"))
saveRDS(list(fitA = fitA, fit_no_income = fit_no_income, fitB = fitB, fitC = fitC, fitD = fitD, nested = fits_nested),
        file.path(OUT, "rev1_fits.rds"))
write_rev1_provenance(OUT, "rev1_p1_main", c(COHORT_PATH, POVERTY_PATH, rev1_fullsample_path()))
cat("\n完成。\n")
sink(); close(log_con)
