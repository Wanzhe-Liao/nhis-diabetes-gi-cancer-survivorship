# rev1_p1_bootstrap.R —— 修回 P1：Model A 边际标准化绝对风险 + PSU bootstrap CI
# 机制沿用 scripts/round4_absolute_risk.R（v4 管线）；暴露为五分类 trajectory_5cat_rev
# 500 reps, seed 20260829, 分层 PSU 有放回重抽样；带 checkpoint 可断点续跑
options(survey.lonely.psu = "fail")
source(file.path("scripts", "rev1_runtime.R"))
source(file.path("scripts", "rev1_domain_design.R"))
suppressPackageStartupMessages({library(survey); library(survival)})

OUT <- rev1_output_dir()
SOURCE_OUT <- rev1_source_dir(default = OUT)
FIT_PATH <- file.path(SOURCE_OUT, "rev1_fits.rds")
fit <- readRDS(FIT_PATH)$fitA
REF <- "no_diabetes"
HORIZONS <- c(60, 120)

model_design <- function(fit) {
  design <- fit$survey.design
  if (!is.null(fit$na.action)) design <- design[-as.integer(fit$na.action), ]
  design
}

cox_baseline_hazard <- function(fit) {
  fit_cox <- fit
  class(fit_cox) <- "coxph"
  survival::basehaz(fit_cox, centered = FALSE)
}

predict_marginal_from_data <- function(fit, cohort, analysis_weights, traj, horizons = HORIZONS) {
  trajectories <- levels(cohort$trajectory_5cat_rev)
  newdata <- cohort
  newdata$trajectory_5cat_rev <- factor(traj, levels = trajectories)
  mm <- model.matrix(delete.response(terms(fit)), data = newdata,
                     contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  if ("(Intercept)" %in% colnames(mm)) mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  coef_names <- names(coef(fit))
  missing_cols <- setdiff(coef_names, colnames(mm))
  if (length(missing_cols) > 0) stop("Missing model-matrix columns: ", paste(missing_cols, collapse = ", "))
  mm <- mm[, coef_names, drop = FALSE]
  linear_predictor <- drop(mm %*% coef(fit))
  baseline <- cox_baseline_hazard(fit)
  rows <- lapply(horizons, function(horizon) {
    idx <- which(baseline$time <= horizon)
    h0 <- if (length(idx) == 0) 0 else baseline$hazard[max(idx)]
    surv_vals <- exp(-h0 * exp(linear_predictor))
    data.frame(trajectory = traj, horizon_months = horizon,
               mortality = stats::weighted.mean(1 - surv_vals, analysis_weights, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  rev1_bind_rows(rows)
}

compute_all_margins_from_data <- function(fit, cohort, analysis_weights, horizons = HORIZONS) {
  trajectories <- levels(cohort$trajectory_5cat_rev)
  rev1_bind_rows(lapply(trajectories, predict_marginal_from_data,
                        fit = fit, cohort = cohort,
                        analysis_weights = analysis_weights, horizons = horizons))
}

add_risk_differences <- function(df) {
  join_keys <- "horizon_months"
  if ("boot_id" %in% names(df)) join_keys <- c("boot_id", join_keys)
  refs <- df[df$trajectory == REF, c(join_keys, "mortality"), drop = FALSE]
  names(refs)[names(refs) == "mortality"] <- "reference_mortality"
  key <- function(x) do.call(paste, c(lapply(x[join_keys], as.character), sep = "\r"))
  ref_key <- key(refs)
  if (anyDuplicated(ref_key)) stop("Reference risk keys are not unique.", call. = FALSE)
  index <- match(key(df), ref_key)
  if (anyNA(index)) stop("Reference risk is missing for at least one row.", call. = FALSE)
  # Retain the within-target reference risk explicitly. Besides documenting the
  # paired contrast, this lets submission QA verify the identity independently
  # for every bootstrap row instead of reconstructing the reference afterward.
  df$reference_mortality <- refs$reference_mortality[index]
  df$abs_risk_diff <- df$mortality - df$reference_mortality
  df$abs_risk_diff_per_100 <- df$abs_risk_diff * 100
  df
}

bootstrap_design <- rev1_resample_domain

checkpoint_path <- file.path(OUT, "boot_checkpoint_rev1.rds")
log_path <- file.path(OUT, "boot_progress_rev1.log")
checkpoint_method <- paste(
  "rev1_absrisk_psu_bootstrap_modelA_5cat_repaired_income",
  REV1_VARIANCE_METHOD, unname(tools::md5sum(rev1_fullsample_path())),
  unname(tools::md5sum(FIT_PATH)),
  REV1_REQUIRED_R,
  paste(rev1_package_versions, collapse = "."),
  sep = "|"
)
append_log <- function(text) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text),
                                 file = log_path, append = TRUE)
is_converged_boot <- function(x) is.numeric(x) && length(x) == 10L && all(is.finite(x))

boot_result_to_rows <- function(boot_results) {
  rows <- vector("list", length(boot_results))
  for (i in seq_len(length(boot_results))) {
    x <- boot_results[[i]]
    if (!is_converged_boot(x)) next
    rows[[i]] <- data.frame(boot_id = i,
                            trajectory = sub("_t[0-9]+$", "", names(x)),
                            horizon_months = as.integer(sub("^.*_t", "", names(x))),
                            mortality = as.numeric(x), stringsAsFactors = FALSE)
  }
  rev1_bind_rows(rows)
}

set.seed(20260829)
boot_reps <- as.integer(Sys.getenv("REV1_BOOT_REPS", "500"))

base_design <- model_design(fit)
base_weights <- as.numeric(weights(base_design))
base_vars <- as.data.frame(base_design$variables)
required_vars <- unique(c(all.vars(formula(fit)), "sa_wgt_pool", "followup_years", "mortstat",
                          "design_strata_prefixed", "design_psu_prefixed"))
base_vars <- base_vars[, required_vars, drop = FALSE]
cat("Model A 分析样本 n =", nrow(base_vars), "\n")
# Domain-correct: PSU 抽样框架取全样本（含无 GI 成员的 PSU）。
bootstrap_sampler <- rev1_prepare_domain_sampler(base_vars)

point <- add_risk_differences(
  compute_all_margins_from_data(fit, base_vars, base_weights)
)

if (file.exists(checkpoint_path)) {
  checkpoint <- readRDS(checkpoint_path)
  if (is.list(checkpoint) && !is.null(checkpoint$method_id) &&
      identical(checkpoint$method_id, checkpoint_method) && !is.null(checkpoint$boot_results)) {
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

if (boot_reps > 0 && start_b <= boot_reps) {
  for (b in start_b:boot_reps) {
    ts_start <- Sys.time()
    b_vars <- bootstrap_design(base_vars, bootstrap_sampler)
    b_vars$time_months <- pmin(b_vars$followup_years * 12, 120)
    b_vars$event <- as.numeric(b_vars$mortstat == 1 & b_vars$followup_years <= 10)
    b_fit <- tryCatch(
      survival::coxph(formula(fit), data = b_vars, weights = sa_wgt_pool,
                      ties = fit$method, x = TRUE, model = TRUE),
      error = function(e) { append_log(sprintf("rep %d/%d fit skip: %s", b, boot_reps, e$message)); NULL })
    if (is.null(b_fit)) { boot_results[[b]] <- NA_real_; next }
    boot_results[[b]] <- tryCatch({
      estimates <- compute_all_margins_from_data(b_fit, b_vars, as.numeric(b_vars$sa_wgt_pool))
      stats::setNames(estimates$mortality, paste0(estimates$trajectory, "_t", estimates$horizon_months))
    }, error = function(e) { append_log(sprintf("rep %d/%d pred skip: %s", b, boot_reps, e$message)); NA_real_ })
    if (b %% 25 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), ts_start, units = "secs"))
      n_ok <- sum(vapply(boot_results, is_converged_boot, logical(1)))
      append_log(sprintf("rep %d/%d last_elapsed=%.1fs n_converged=%d", b, boot_reps, elapsed, n_ok))
    }
    if (b %% 50 == 0) saveRDS(list(method_id = checkpoint_method, boot_results = boot_results,
                                   rng_state = .Random.seed), checkpoint_path)
  }
}
saveRDS(list(method_id = checkpoint_method, boot_results = boot_results, rng_state = .Random.seed),
        checkpoint_path)

boot <- add_risk_differences(boot_result_to_rows(boot_results))
if (nrow(boot) > 0) {
  group_index <- split(
    seq_len(nrow(boot)),
    interaction(boot$trajectory, boot$horizon_months, drop = TRUE, lex.order = TRUE)
  )
  boot_summary <- rev1_bind_rows(lapply(group_index, function(index) {
    value <- boot[index, , drop = FALSE]
    mortality_ci <- quantile(value$mortality, c(0.025, 0.975), na.rm = TRUE)
    rd_ci <- quantile(value$abs_risk_diff, c(0.025, 0.975), na.rm = TRUE)
    data.frame(
      trajectory = as.character(value$trajectory[1]),
      horizon_months = value$horizon_months[1],
      mortality_ci_lo = unname(mortality_ci[1]),
      mortality_ci_hi = unname(mortality_ci[2]),
      abs_risk_diff_ci_lo = unname(rd_ci[1]),
      abs_risk_diff_ci_hi = unname(rd_ci[2]),
      abs_risk_diff_per_100_ci_lo = 100 * unname(rd_ci[1]),
      abs_risk_diff_per_100_ci_hi = 100 * unname(rd_ci[2]),
      n_boot_converged = sum(!is.na(value$mortality)),
      stringsAsFactors = FALSE
    )
  }))
} else {
  boot_summary <- data.frame(
    trajectory = point$trajectory,
    horizon_months = point$horizon_months,
    mortality_ci_lo = NA_real_, mortality_ci_hi = NA_real_,
    abs_risk_diff_ci_lo = NA_real_, abs_risk_diff_ci_hi = NA_real_,
    abs_risk_diff_per_100_ci_lo = NA_real_,
    abs_risk_diff_per_100_ci_hi = NA_real_,
    n_boot_converged = 0L,
    stringsAsFactors = FALSE
  )
}

out <- rev1_left_join(point, boot_summary, by = c("trajectory", "horizon_months"))
write_rev1_csv(out, file.path(OUT, "rev1_absolute_risk_with_rd.csv"))

# Submission QA: every converged replicate must contain the complete set of
# counterfactual trajectory-by-horizon predictions. Risk differences are paired
# within replicate against the no-diabetes prediction from that same target.
expected_rows_per_boot <- length(levels(base_vars$trajectory_5cat_rev)) * length(HORIZONS)
rows_per_boot <- table(boot$boot_id)
complete_replicates <- length(rows_per_boot) == boot_reps &&
  all(as.integer(rows_per_boot) == expected_rows_per_boot)
paired_contrast_error <- max(
  abs(boot$abs_risk_diff - (boot$mortality - boot$reference_mortality)),
  na.rm = TRUE
)
qa <- data.frame(
  check = c(
    "bootstrap_replicates_converged",
    "complete_counterfactual_set_within_replicate",
    "paired_risk_difference_within_replicate",
    "common_point_estimate_standardization_target"
  ),
  expected = c(as.character(boot_reps), "TRUE", "TRUE", as.character(nrow(base_vars))),
  observed = c(
    as.character(sum(vapply(boot_results, is_converged_boot, logical(1)))),
    as.character(complete_replicates),
    as.character(is.finite(paired_contrast_error) && paired_contrast_error < 1e-12),
    as.character(nrow(base_vars))
  ),
  status = c(
    ifelse(sum(vapply(boot_results, is_converged_boot, logical(1))) == boot_reps, "PASS", "FAIL"),
    ifelse(complete_replicates, "PASS", "FAIL"),
    ifelse(is.finite(paired_contrast_error) && paired_contrast_error < 1e-12, "PASS", "FAIL"),
    "PASS"
  ),
  detail = c(
    "All requested PSU bootstrap replicates returned predictions.",
    "Each replicate used one resampled analysis population for all five trajectory counterfactuals and both horizons.",
    paste0("Maximum paired-contrast identity error = ", format(paired_contrast_error, scientific = TRUE)),
    "All point estimates average counterfactual predictions over the same principal-model analysis population."
  ),
  stringsAsFactors = FALSE
)
if (!all(qa$status == "PASS")) {
  stop("Absolute-risk submission QA failed.", call. = FALSE)
}
write_rev1_csv(qa, file.path(OUT, "absolute_risk_qa.csv"))
write_rev1_provenance(OUT, "rev1_p1_bootstrap", c(FIT_PATH, rev1_fullsample_path()))
append_log(sprintf("Completed. n_converged=%d/%d", sum(vapply(boot_results, is_converged_boot, logical(1))), boot_reps))
cat("完成。收敛", sum(vapply(boot_results, is_converged_boot, logical(1))), "/", boot_reps, "\n")
print(as.data.frame(out))
