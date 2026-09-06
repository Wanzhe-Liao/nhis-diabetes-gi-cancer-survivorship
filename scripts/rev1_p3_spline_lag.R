# rev1_p3_spline_lag.R —— 修回 P3：探索性 survey-weighted Cox spline lag 曲线
#
# 定位：探索性 Supplementary Figure（不改变主结论）。
# 规格（GPT Pro 审定）：
#   * 仅在 established pre-cancer DM 中对 diabetes-to-cancer lag 拟合 3-knot RCS；
#   * no diabetes 为参照；peri-diagnostic / post-cancer / unknown 为独立分类系数，
#     不与 spline 曲线连接；
#   * svycoxph + 当前 v4 权重 / year-prefixed strata / PSU，协变量与主模型一致；
#   * knots = established 组加权 lag 分布的 10/50/90 百分位；
#   * 95% pointwise CI：stratified PSU bootstrap，500 次，每次重拟合完整模型；
#   * 报告 P for nonlinearity（非线性项 Wald）、spline 相关项 PH 检验、
#     3-knot vs 4-knot 形状稳定性。
# R 4.6.1 / 输出：outputs/revision_round1_v4_r461_submission/spline_lag_*.csv
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
log_con <- file(file.path(OUT, "run_log_spline.txt"), open = "wt")
sink(log_con, split = TRUE)
cat("rev1_p3_spline_lag.R |", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(rev1_runtime_info(), row.names = FALSE)

## ---------- 0. 数据准备（与 rev1_p1_main.R 完全一致） ----------
cohort <- read_rev1_parquet(COHORT_PATH)
POVERTY_PATH <- file.path(SUPPORT_OUT, "income_repair", "poverty_repaired.csv")
pov <- read.csv(POVERTY_PATH, stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN"))
cohort$year <- as.integer(cohort$year)
cohort$hhx_k <- as.character(as.integer(cohort$hhx))
cohort$fmx_k <- as.character(as.integer(cohort$fmx))
pov$hhx_k <- as.character(as.integer(pov$hhx))
pov$fmx_k <- as.character(as.integer(pov$fmx))
pov$year <- as.integer(pov$year)
cohort <- rev1_left_join(
  cohort,
  pov[, c("year", "hhx_k", "fmx_k", "poverty_3cat_rep"), drop = FALSE],
  by = c("year", "hhx_k", "fmx_k")
)

cohort$time_months <- pmin(cohort$followup_years * 12, 120)
cohort$event <- as.numeric(cohort$mortstat == 1 & cohort$followup_years <= 10)
cohort <- cohort[!is.na(cohort$time_months) & cohort$time_months > 0, ]
cohort <- cohort[!is.na(cohort$design_psu) & !is.na(cohort$design_strata) & !is.na(cohort$sa_wgt_pool), ]
cohort <- rev1_period_identifiers(cohort)

make_factor <- function(x, ref) {
  x <- as.character(x); x[x == ""] <- NA
  x <- factor(x)
  if (ref %in% levels(x)) x <- stats::relevel(x, ref = ref)
  x
}
cohort$sex <- make_factor(cohort$sex, "1")
cohort$race <- make_factor(cohort$race, "1")
cohort$region <- make_factor(cohort$region, "1")
cohort$smoking_3cat <- make_factor(cohort$smoking_3cat, "current")
cohort$education_4cat <- make_factor(cohort$education_4cat, "college_grad")
cohort$poverty_3cat_rep <- make_factor(cohort$poverty_3cat_rep, "2_0_to_3_99")
cohort$poverty_3cat_rep <- factor(cohort$poverty_3cat_rep, exclude = NULL)
levels(cohort$poverty_3cat_rep)[is.na(levels(cohort$poverty_3cat_rep))] <- "missing"
cohort$poverty_3cat_rep <- stats::relevel(cohort$poverty_3cat_rep, ref = "2_0_to_3_99")
cohort$survey_year <- make_factor(cohort$survey_year, "1997")

SITE_FLAGS <- c("colon_flag", "esoph_flag", "gallbladder_flag", "liver_flag",
                "pancreas_flag", "rectum_flag", "stomach_flag")
for (v in SITE_FLAGS) cohort[[v]] <- as.logical(cohort[[v]])

## ---------- 1. lag 暴露变量 ----------
cohort$lag_years <- cohort$gi_first_dx_age - cohort$dm_dx_age
cohort$exp_est <- as.integer(cohort$trajectory_6cat %in% c("dm_to_gi_2_10y", "dm_to_gi_gt10y"))
cohort$exp_peri <- as.integer(cohort$trajectory_6cat == "peri_diagnostic")
cohort$exp_post <- as.integer(cohort$trajectory_6cat == "gi_to_dm")
cohort$exp_unk  <- as.integer(cohort$trajectory_6cat == "dm_order_unknown")
stopifnot(all(cohort$lag_years[cohort$exp_est == 1] >= 2))
stopifnot(!anyNA(cohort$lag_years[cohort$exp_est == 1]))

COV <- c("age", "sex", "race", "region", "bmi", "smoking_3cat", "survey_year",
         SITE_FLAGS, "education_4cat", "poverty_3cat_rep")
MODEL_VARS <- c("time_months", "event", "exp_est", "exp_peri", "exp_post", "exp_unk",
                "lag_years", COV, "sa_wgt_pool", "design_strata_prefixed", "design_psu_prefixed")
cc <- complete.cases(cohort[, setdiff(MODEL_VARS, "lag_years")]) &
  (cohort$exp_est == 0 | !is.na(cohort$lag_years))
adat <- cohort[cc, ]
cat("分析样本 n =", nrow(adat), " events =", sum(adat$event),
    " | established n =", sum(adat$exp_est), "\n")
principal_result <- read.csv(file.path(OUT, "master_results.csv"), stringsAsFactors = FALSE)
principal_result <- principal_result[principal_result$model == "A_principal_5cat", ][1, ]
stopifnot(nrow(adat) == principal_result$n, sum(adat$event) == principal_result$events)

# restricted cubic spline basis（Harrell 参数化，线性尾部）
rcs_nonlin <- function(x, knots) {
  k <- length(knots)
  p <- function(z) pmax(z, 0)^3
  out <- matrix(NA_real_, length(x), k - 2L)
  for (j in seq_len(k - 2L)) {
    out[, j] <- p(x - knots[j]) -
      p(x - knots[k - 1L]) * (knots[k] - knots[j]) / (knots[k] - knots[k - 1L]) +
      p(x - knots[k]) * (knots[k - 1L] - knots[j]) / (knots[k] - knots[k - 1L])
  }
  out
}
rcs_terms <- function(x, knots, l0) {
  # 返回以 l0 为中心的 (linear, nonlinear...) 基；x 可为向量
  nl <- rcs_nonlin(x, knots)
  nl0 <- rcs_nonlin(l0, knots)
  cbind(x - l0, sweep(nl, 2L, nl0, "-"))
}

## Domain-correct: 全样本 design + subset（NCHS domain 规范）
source(file.path("scripts", "rev1_domain_design.R"))
combined <- rev1_attach_fullsample(cohort)
design_pop <- rev1_domain_population_design(combined)
flag_adat <- combined$gi_any %in% TRUE &
  complete.cases(combined[, setdiff(MODEL_VARS, "lag_years")]) &
  (combined$exp_est == 0 | !is.na(combined$lag_years))
stopifnot(sum(flag_adat) == nrow(adat))
design_est_only <- subset(design_pop, flag_adat & combined$exp_est == 1)
wq <- function(probs) as.numeric(svyquantile(~lag_years, design_est_only,
                                             quantiles = probs, ci = FALSE, na.rm = TRUE)[[1]][1, ])
KNOTS3 <- wq(c(0.10, 0.50, 0.90))
KNOTS4 <- wq(c(0.05, 0.35, 0.65, 0.95))
L0 <- KNOTS3[2]                       # 参照 lag = 加权中位数
Q_RANGE <- wq(c(0.025, 0.975))        # 曲线展示范围
cat("3-knot:", KNOTS3, " 4-knot:", KNOTS4, " L0:", L0, " range:", Q_RANGE, "\n")

build_terms <- function(dat, knots, l0) {
  est <- dat$exp_est == 1
  tr <- rcs_terms(ifelse(est, dat$lag_years, l0), knots, l0)
  for (j in seq_len(ncol(tr))) {
    dat[[paste0("lag_b", j)]] <- ifelse(est, tr[, j], 0)
  }
  dat
}
adat <- build_terms(adat, KNOTS3, L0)
N_SPLINE <- ncol(rcs_terms(L0, KNOTS3, L0))   # 2 (linear + 1 nonlinear)
EXPO_TERMS <- c("exp_est", paste0("lag_b", seq_len(N_SPLINE)))
rhs <- paste(c(EXPO_TERMS, "exp_peri", "exp_post", "exp_unk", COV), collapse = " + ")
fm <- as.formula(paste("Surv(time_months, event) ~", rhs))

combined <- build_terms(combined, KNOTS3, L0)  # 非 GI 行 exp_est 为 NA, lag_b 置 NA
design_adat <- subset(rev1_domain_population_design(combined), flag_adat)
fitS <- svycoxph(fm, design = design_adat)
stopifnot(!anyNA(coef(fitS)))
cat("\nspline 模型暴露项:\n")
print(summary(fitS)$coefficients[ c(EXPO_TERMS, "exp_peri"), ])

## ---------- 2. 点估计曲线与 peri 点 ----------
GRID <- seq(max(2, floor(Q_RANGE[1] * 2) / 2), ceiling(Q_RANGE[2] * 2) / 2, by = 0.25)
curve_loghr <- function(coefs, knots, l0, grid) {
  tr <- rcs_terms(grid, knots, l0)
  drop(coefs["exp_est"] + tr %*% coefs[paste0("lag_b", seq_len(ncol(tr)))])
}
pt_loghr <- curve_loghr(coef(fitS), KNOTS3, L0, GRID)

sS <- summary(fitS)$coefficients
peri <- data.frame(
  term = "peri_diagnostic", HR = exp(sS["exp_peri", "coef"]),
  ci_lo = exp(sS["exp_peri", "coef"] - 1.96 * sS["exp_peri", "robust se"]),
  ci_hi = exp(sS["exp_peri", "coef"] + 1.96 * sS["exp_peri", "robust se"]),
  p = sS["exp_peri", "Pr(>|z|)"])

## ---------- 3. 诊断：P for nonlinearity / PH / 4-knot 稳定性 ----------
V <- vcov(fitS)
b_nl <- coef(fitS)["lag_b2"]
wald_nl <- b_nl^2 / V["lag_b2", "lag_b2"]
p_nonlin <- pchisq(wald_nl, df = 1, lower.tail = FALSE)

# PH 诊断沿用锁定管线权威方法：svycoxph 强制为 coxph 后 cox.zph
# （使用拟合对象的归一化 case weights，避免原始权重虚增卡方）
fit_zph <- fitS
class(fit_zph) <- "coxph"
ph <- cox.zph(fit_zph)
ph_tab <- ph$table[c("exp_est", paste0("lag_b", seq_len(N_SPLINE)), "GLOBAL"), , drop = FALSE]
cat("\nPH 诊断（spline 相关项 + GLOBAL）:\n"); print(ph_tab)

adat4 <- build_terms(adat, KNOTS4, L0)
# 4-knot 基有 3 列；与 3-knot 模型变量名不兼容，单独拟合
N4 <- ncol(rcs_terms(L0, KNOTS4, L0))
EXPO4 <- c("exp_est", paste0("lag_b", seq_len(N4)))
fm4 <- as.formula(paste("Surv(time_months, event) ~",
                        paste(c(EXPO4, "exp_peri", "exp_post", "exp_unk", COV), collapse = " + ")))
combined4 <- build_terms(combined, KNOTS4, L0)
design_adat4 <- subset(rev1_domain_population_design(combined4), flag_adat)
fitS4 <- svycoxph(fm4, design = design_adat4)
stopifnot(!anyNA(coef(fitS4)))
pt4 <- curve_loghr(coef(fitS4), KNOTS4, L0, GRID)
stab_max_abs <- max(abs(pt4 - pt_loghr))
stab_cor <- cor(pt4, pt_loghr)
cat(sprintf("4-knot 稳定性: max|ΔlogHR| = %.4f, cor = %.5f\n", stab_max_abs, stab_cor))

## ---------- 4. stratified PSU bootstrap（500 次，复用 rev1_p1_bootstrap 机制） ----------
bootstrap_design <- rev1_resample_domain

boot_reps <- as.integer(Sys.getenv("REV1_BOOT_REPS", "500"))
checkpoint_path <- file.path(OUT, "boot_checkpoint_spline.rds")
log_path <- file.path(OUT, "boot_progress_spline.log")
append_log <- function(text) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text),
                                 file = log_path, append = TRUE)
base_vars <- as.data.frame(design_adat$variables)
# Domain-correct: PSU 抽样框架取全样本（含无 GI 成员的 PSU）。
sampler <- rev1_prepare_domain_sampler(base_vars)
stat_names <- c(paste0("g", seq_along(GRID)), "peri")

point_stats <- c(exp(pt_loghr), peri = exp(unname(coef(fitS)["exp_peri"])))
names(point_stats) <- stat_names

checkpoint_method <- paste("rev1_spline_lag_psu_bootstrap_3knot_domain_fullsample",
                           REV1_VARIANCE_METHOD, unname(tools::md5sum(COHORT_PATH)),
                           unname(tools::md5sum(rev1_fullsample_path())), REV1_REQUIRED_R, sep = "|")
if (file.exists(checkpoint_path)) {
  checkpoint <- readRDS(checkpoint_path)
  if (is.list(checkpoint) && identical(checkpoint$method_id, checkpoint_method) &&
      !is.null(checkpoint$boot_results)) {
    boot_results <- checkpoint$boot_results
    if (!is.null(checkpoint$rng_state)) .Random.seed <<- checkpoint$rng_state
    if (length(boot_results) < boot_reps) length(boot_results) <- boot_reps
    start_b <- max(which(!vapply(boot_results, is.null, logical(1))), 0) + 1
    append_log(sprintf("Resuming bootstrap from replicate %d/%d", start_b, boot_reps))
  } else {
    boot_results <- vector("list", boot_reps); start_b <- 1
    append_log(sprintf("Ignoring incompatible checkpoint; starting B=%d", boot_reps))
  }
} else {
  boot_results <- vector("list", boot_reps); start_b <- 1
  append_log(sprintf("Starting bootstrap with B=%d", boot_reps))
}

need <- c(EXPO_TERMS, "exp_peri")
if (boot_reps > 0 && start_b <= boot_reps) {
  for (b in start_b:boot_reps) {
    b_vars <- bootstrap_design(base_vars, sampler)
    b_fit <- tryCatch(
      survival::coxph(fm, data = b_vars, weights = sa_wgt_pool),
      error = function(e) { append_log(sprintf("rep %d/%d fit skip: %s", b, boot_reps, e$message)); NULL })
    if (is.null(b_fit) || anyNA(coef(b_fit)[need])) {
      if (!is.null(b_fit)) append_log(sprintf("rep %d/%d aliased exposure coef", b, boot_reps))
      boot_results[[b]] <- NA_real_; next
    }
    bc <- coef(b_fit)
    vals <- c(exp(curve_loghr(bc, KNOTS3, L0, GRID)), peri = exp(unname(bc["exp_peri"])))
    boot_results[[b]] <- stats::setNames(vals, stat_names)
    if (b %% 25 == 0) {
      n_ok <- sum(vapply(boot_results, function(x) is.numeric(x) && length(x) == length(stat_names) && all(is.finite(x)), logical(1)))
      append_log(sprintf("rep %d/%d n_converged=%d", b, boot_reps, n_ok))
    }
    if (b %% 50 == 0) saveRDS(list(method_id = checkpoint_method, boot_results = boot_results,
                                   rng_state = .Random.seed), checkpoint_path)
  }
}
saveRDS(list(method_id = checkpoint_method, boot_results = boot_results, rng_state = .Random.seed),
        checkpoint_path)

is_ok <- function(x) is.numeric(x) && length(x) == length(stat_names) && all(is.finite(x))
boot_mat <- do.call(rbind, lapply(boot_results[vapply(boot_results, is_ok, logical(1))],
                                  function(x) as.numeric(x[stat_names])))
colnames(boot_mat) <- stat_names
n_conv <- nrow(boot_mat)
cat("bootstrap 收敛:", n_conv, "/", boot_reps, "\n")
stopifnot(n_conv == boot_reps)
ci <- apply(boot_mat, 2L, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

curve_df <- data.frame(
  lag_years = GRID, hr = exp(pt_loghr),
  ci_lo = ci[1, paste0("g", seq_along(GRID))],
  ci_hi = ci[2, paste0("g", seq_along(GRID))],
  hr_4knot = exp(pt4), n_boot_converged = n_conv, row.names = NULL)
peri$ci_lo_boot <- ci[1, "peri"]; peri$ci_hi_boot <- ci[2, "peri"]

write_rev1_csv(curve_df, file.path(OUT, "spline_lag_curve.csv"))
write_rev1_csv(peri, file.path(OUT, "spline_lag_peri.csv"))

rug <- data.frame(lag_years = adat$lag_years[adat$exp_est == 1],
                  weight = adat$sa_wgt_pool[adat$exp_est == 1])
write_rev1_csv(rug, file.path(OUT, "spline_lag_rug.csv"))

diag <- data.frame(
  item = c("n_analysis", "events", "n_established", "events_established",
           "knot3_1", "knot3_2", "knot3_3", "reference_lag",
           "range_p2.5", "range_p97.5",
           "p_nonlinearity", "ph_exp_est", "ph_lag_b1", "ph_lag_b2", "ph_global",
           "stab4_max_abs_loghr_diff", "stab4_cor", "boot_reps", "boot_converged"),
  value = c(nrow(adat), sum(adat$event), sum(adat$exp_est), sum(adat$event[adat$exp_est == 1]),
            KNOTS3, L0, Q_RANGE, p_nonlin,
            ph_tab["exp_est", "p"], ph_tab["lag_b1", "p"], ph_tab["lag_b2", "p"], ph_tab["GLOBAL", "p"],
            stab_max_abs, stab_cor, boot_reps, n_conv))
write_rev1_csv(diag, file.path(OUT, "spline_lag_diagnostics.csv"))
write_rev1_provenance(OUT, "rev1_p3_spline_lag", c(COHORT_PATH, POVERTY_PATH, rev1_fullsample_path()))
print(diag)
print(peri)
cat("\n完成。\n")
sink(); close(log_con)
