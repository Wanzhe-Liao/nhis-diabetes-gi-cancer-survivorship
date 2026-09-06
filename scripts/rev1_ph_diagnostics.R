# rev1_ph_diagnostics.R
# Task 0: authoritative proportional-hazards diagnostic for fitA.
#
# survey::svycoxph stores the normalized case weights actually used by the
# underlying Cox fit in fit$weights.  With survey 4.5/survival 3.8.6 loaded,
# survival::cox.zph() on an object whose first class is "svycoxph" rebuilds the
# model frame and obtains survey-design weights on a sum-to-one scale.  The
# resulting score-test chi-squares are therefore spuriously close to zero.
# Coercing a copy to the underlying "coxph" class makes cox.zph use the fitted
# case-weight scale.  This is also the clean-session method used for the
# historical diagnostic and is checked below against the v3 reference values.
options(survey.lonely.psu = "remove")
source(file.path("scripts", "rev1_runtime.R"))
suppressPackageStartupMessages(library(survival))

OUT <- rev1_output_dir()
SOURCE_OUT <- rev1_source_dir(default = OUT)
SUPPORT_OUT <- rev1_support_dir()
FIT_PATH <- file.path(SOURCE_OUT, "rev1_fits.rds")
fits <- readRDS(FIT_PATH)
fitA <- fits[["fitA"]]

if (!inherits(fitA, "svycoxph") || !inherits(fitA, "coxph")) {
  stop("fitA must inherit from both svycoxph and coxph.", call. = FALSE)
}
if (is.null(fitA[["weights"]]) || any(!is.finite(fitA[["weights"]])) ||
    any(fitA[["weights"]] <= 0)) {
  stop("fitA has absent, non-finite, or non-positive stored case weights.",
       call. = FALSE)
}

cox_zph_with_fitted_weights <- function(fit) {
  if (!inherits(fit, "coxph")) {
    stop("PH diagnostic requires an underlying coxph fit.", call. = FALSE)
  }
  fit_zph <- fit
  class(fit_zph) <- "coxph"
  survival::cox.zph(fit_zph)
}

# Authoritative result: use the case weights stored by the fitted Cox object.
zph <- cox_zph_with_fitted_weights(fitA)

# Keep the R 4.6.1 svy-class dispatch result only as a bug/audit comparator.
# It must never be used in the manuscript or master-results source.
zph_svy_dispatch <- survival::cox.zph(fitA)

if (any(!is.finite(zph$table[, "chisq"])) || any(!is.finite(zph$table[, "p"]))) {
  stop("The corrected PH diagnostic contains non-finite values.", call. = FALSE)
}
if ("GLOBAL" %in% rownames(zph$table) &&
    zph$table["GLOBAL", "chisq"] < 1 &&
    zph_svy_dispatch$table["GLOBAL", "p"] > 0.999) {
  stop(
    "PH sanity check failed: the fitted-weight GLOBAL chi-square remains implausibly small.",
    call. = FALSE
  )
}

tab <- as.data.frame(zph[["table"]])
tab$term <- rownames(tab)
rownames(tab) <- NULL
tab <- tab[, c("term", "chisq", "df", "p")]
tab$model <- "A_principal_5cat"
tab$n <- nrow(fitA[["y"]])
tab$events_10y <- sum(fitA[["y"]][, "status"])
tab$diagnostic_method <- "cox.zph_underlying_coxph_fitted_case_weights"
tab <- tab[, c("model", "term", "chisq", "df", "p", "n", "events_10y",
               "diagnostic_method")]

dispatch_tab <- as.data.frame(zph_svy_dispatch$table)
dispatch_tab$term <- rownames(dispatch_tab)
rownames(dispatch_tab) <- NULL
dispatch_tab <- dispatch_tab[, c("term", "chisq", "df", "p")]
names(dispatch_tab)[2:4] <- c(
  "svy_dispatch_chisq", "svy_dispatch_df", "svy_dispatch_p"
)
dispatch_audit <- merge(
  tab[, c("term", "chisq", "df", "p")],
  dispatch_tab,
  by = "term",
  all = TRUE,
  sort = FALSE
)
names(dispatch_audit)[2:4] <- c(
  "fitted_weight_chisq", "fitted_weight_df", "fitted_weight_p"
)
dispatch_audit$chisq_scale_ratio <- with(
  dispatch_audit,
  svy_dispatch_chisq / fitted_weight_chisq
)

stored_weights <- fitA[["weights"]]
svy_frame_weights <- stats::model.weights(stats::model.frame(fitA))
fitA_coxph <- fitA
class(fitA_coxph) <- "coxph"
coxph_frame_weights <- stats::model.weights(stats::model.frame(fitA_coxph))
weight_summary <- function(x, source) {
  if (is.null(x)) x <- rep(1, nrow(fitA[["y"]]))
  data.frame(
    source = source,
    n = length(x),
    min = min(x),
    mean = mean(x),
    max = max(x),
    sum = sum(x),
    mean_relative_to_stored = mean(x) / mean(stored_weights),
    stringsAsFactors = FALSE
  )
}
weight_audit <- rev1_bind_rows(list(
  weight_summary(stored_weights, "fitA_stored_case_weights"),
  weight_summary(svy_frame_weights, "svycoxph_model_frame_dispatch"),
  weight_summary(coxph_frame_weights, "underlying_coxph_model_frame")
))

# Optional historical comparison requires an explicitly paired fit and table.
# A table from another population is not a valid regression target.
ref_fit_env <- Sys.getenv("REV1_PH_REFERENCE_FIT_PATH", "")
ref_table_env <- Sys.getenv("REV1_PH_REFERENCE_TABLE_PATH", "")
if (xor(nzchar(ref_fit_env), nzchar(ref_table_env)))
  stop("Provide both historical PH reference paths or neither.", call. = FALSE)
REFERENCE_FIT_PATH <- if (nzchar(ref_fit_env)) rev1_absolute_path(ref_fit_env, TRUE) else ""
REFERENCE_TABLE_PATH <- if (nzchar(ref_table_env)) rev1_absolute_path(ref_table_env, TRUE) else ""

# Independent refit on the current analysis records checks the case-weight
# diagnostic, without depending on a mutable historical analysis table.
analysis_design <- fitA$survey.design
if (!is.null(fitA$na.action)) analysis_design <- analysis_design[-as.integer(fitA$na.action), ]
analysis_data <- as.data.frame(analysis_design$variables)
analysis_data$ph_case_weight <- stored_weights
refit <- survival::coxph(formula(fitA), data = analysis_data,
                        weights = ph_case_weight, ties = fitA$method,
                        x = TRUE, model = TRUE)
refit_zph <- survival::cox.zph(refit)$table
stopifnot(identical(rownames(refit_zph), rownames(zph$table)),
          max(abs(refit_zph - zph$table)) < 1e-6)
cat("Independent fitted-case-weight PH refit matched current diagnostic.\n")

crosscheck <- data.frame(
  term = character(), reference_chisq = numeric(), r461_chisq = numeric(),
  abs_chisq_diff = numeric(), reference_p = numeric(), r461_p = numeric(),
  abs_p_diff = numeric(), status = character(), stringsAsFactors = FALSE
)
if (file.exists(REFERENCE_FIT_PATH) && file.exists(REFERENCE_TABLE_PATH)) {
  reference_fit <- readRDS(REFERENCE_FIT_PATH)[["fitA"]]
  reference_r461 <- as.data.frame(
    cox_zph_with_fitted_weights(reference_fit)$table
  )
  reference_r461$term <- rownames(reference_r461)
  rownames(reference_r461) <- NULL
  historical <- utils::read.csv(
    REFERENCE_TABLE_PATH, stringsAsFactors = FALSE, check.names = FALSE
  )
  crosscheck <- merge(
    historical[, c("term", "chisq", "df", "p")],
    reference_r461[, c("term", "chisq", "df", "p")],
    by = "term", suffixes = c("_reference", "_r461"), all = TRUE,
    sort = FALSE
  )
  crosscheck$abs_chisq_diff <- abs(
    crosscheck$chisq_reference - crosscheck$chisq_r461
  )
  crosscheck$abs_p_diff <- abs(crosscheck$p_reference - crosscheck$p_r461)
  crosscheck$status <- ifelse(
    is.finite(crosscheck$abs_chisq_diff) &
      is.finite(crosscheck$abs_p_diff) &
      crosscheck$df_reference == crosscheck$df_r461 &
      crosscheck$abs_chisq_diff <= 1e-10 &
      crosscheck$abs_p_diff <= 1e-10,
    "PASS", "FAIL"
  )
  names(crosscheck) <- sub("_reference$", "_reference", names(crosscheck))
  names(crosscheck) <- sub("_r461$", "_r461", names(crosscheck))
  if (any(crosscheck$status != "PASS")) {
    stop(
      "R 4.6.1 PH regression check did not reproduce the historical v3 values.",
      call. = FALSE
    )
  }
}

write_rev1_csv(tab, file.path(OUT, "ph_diagnostics_rev1.csv"))
write_rev1_csv(
  dispatch_audit, file.path(OUT, "ph_diagnostics_dispatch_audit.csv")
)
write_rev1_csv(weight_audit, file.path(OUT, "ph_weight_scale_audit.csv"))
write_rev1_csv(crosscheck, file.path(OUT, "ph_diagnostics_v3_crosscheck.csv"))
write_rev1_provenance(
  OUT, "rev1_ph_diagnostics",
  c(FIT_PATH, REFERENCE_FIT_PATH, REFERENCE_TABLE_PATH)
)
cat("Wrote", file.path(OUT, "ph_diagnostics_rev1.csv"), "\n")
print(tab)
cat("\nWeight-scale audit:\n")
print(weight_audit)
if (nrow(crosscheck) > 0L) {
  cat("\nHistorical v3 cross-check:",
      sum(crosscheck$status == "PASS"), "of", nrow(crosscheck), "PASS\n")
}
