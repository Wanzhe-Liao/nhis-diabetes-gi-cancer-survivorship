# rev1_p2_stomach_causes.R
# Reviewer-requested stomach-cancer competing-event support audit.
options(survey.lonely.psu = "remove")
source(file.path("scripts", "rev1_runtime.R"))

OUT <- rev1_output_dir()
SUPPORT_OUT <- rev1_support_dir()
COHORT_PATH <- rev1_cohort_path()
POVERTY_PATH <- file.path(SUPPORT_OUT, "income_repair", "poverty_repaired.csv")

cohort_raw <- read_rev1_parquet(COHORT_PATH)
pov <- read.csv(
  POVERTY_PATH, stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)
cohort_raw$year <- as.integer(cohort_raw$year)
cohort_raw$hhx_k <- as.character(as.integer(cohort_raw$hhx))
cohort_raw$fmx_k <- as.character(as.integer(cohort_raw$fmx))
pov$year <- as.integer(pov$year)
pov$hhx_k <- as.character(as.integer(pov$hhx))
pov$fmx_k <- as.character(as.integer(pov$fmx))
cohort_raw <- rev1_left_join(
  cohort_raw,
  pov[, c("year", "hhx_k", "fmx_k", "poverty_3cat_rep"), drop = FALSE],
  by = c("year", "hhx_k", "fmx_k")
)

make_factor <- function(x, ref = NULL) {
  x <- as.character(x)
  x[x == ""] <- NA
  x <- factor(x)
  if (!is.null(ref) && ref %in% levels(x)) x <- stats::relevel(x, ref = ref)
  x
}

make_trajectory <- function(x) {
  x <- as.character(x)
  x[x %in% c("dm_to_gi_2_10y", "dm_to_gi_gt10y")] <- "established_pre_cancer_dm"
  x[x == "gi_to_dm"] <- "post_cancer_dm"
  x[x == "gi_only"] <- "no_diabetes"
  factor(
    x,
    levels = c(
      "no_diabetes", "established_pre_cancer_dm", "peri_diagnostic",
      "post_cancer_dm", "dm_order_unknown"
    )
  )
}

cohort_raw$trajectory_5cat_rev <- make_trajectory(cohort_raw$trajectory_6cat)
is_stomach_raw <- !is.na(cohort_raw$stomach_flag) & as.logical(cohort_raw$stomach_flag)
stomach_cohort_n <- sum(is_stomach_raw)

cohort <- cohort_raw
cohort$time_months <- pmin(cohort$followup_years * 12, 120)
cohort$event <- as.numeric(cohort$mortstat == 1 & cohort$followup_years <= 10)
positive_followup <- !is.na(cohort$time_months) & cohort$time_months > 0
stomach_positive_n <- sum(is_stomach_raw & positive_followup)
cohort <- cohort[positive_followup, , drop = FALSE]
design_ok <- !is.na(cohort$design_psu) & !is.na(cohort$design_strata) &
  !is.na(cohort$sa_wgt_pool)
stomach_design_n <- sum(as.logical(cohort$stomach_flag) & design_ok, na.rm = TRUE)
cohort <- cohort[design_ok, , drop = FALSE]

cohort$trajectory_5cat_rev <- make_trajectory(cohort$trajectory_6cat)
cohort$sex <- make_factor(cohort$sex, "1")
cohort$race <- make_factor(cohort$race, "1")
cohort$region <- make_factor(cohort$region, "1")
cohort$smoking_3cat <- make_factor(cohort$smoking_3cat, "current")
cohort$survey_year <- make_factor(cohort$survey_year, "1997")
cohort$education_4cat <- make_factor(cohort$education_4cat, "college_grad")
cohort$poverty_3cat_rep <- make_factor(cohort$poverty_3cat_rep, "2_0_to_3_99")
cohort$poverty_3cat_rep <- factor(cohort$poverty_3cat_rep, exclude = NULL)
levels(cohort$poverty_3cat_rep)[is.na(levels(cohort$poverty_3cat_rep))] <- "missing"
cohort$poverty_3cat_rep <- stats::relevel(
  cohort$poverty_3cat_rep, ref = "2_0_to_3_99"
)

site_flags <- c(
  "colon_flag", "esoph_flag", "gallbladder_flag", "liver_flag",
  "pancreas_flag", "rectum_flag", "stomach_flag"
)
for (variable in site_flags) cohort[[variable]] <- as.logical(cohort[[variable]])

principal_vars <- c(
  "trajectory_5cat_rev", "time_months", "event", "age", "sex", "race",
  "region", "bmi", "smoking_3cat", "survey_year", site_flags,
  "education_4cat", "poverty_3cat_rep"
)
stomach <- cohort[cohort$stomach_flag, , drop = FALSE]
stomach <- stomach[stats::complete.cases(stomach[, principal_vars, drop = FALSE]), , drop = FALSE]
stomach$trajectory_5cat_rev <- factor(
  stomach$trajectory_5cat_rev,
  levels = levels(cohort$trajectory_5cat_rev)
)

flow <- data.frame(
  step = c(
    "stomach_flag cohort", "after positive follow-up", "after design variables",
    "principal covariate sample", "10y all-cause events"
  ),
  n = c(
    stomach_cohort_n, stomach_positive_n, stomach_design_n, nrow(stomach),
    sum(stomach$event)
  ),
  stringsAsFactors = FALSE
)
write_rev1_csv(flow, file.path(OUT, "stomach_subgroup_flow.csv"))

cause_code <- suppressWarnings(as.integer(as.character(stomach$ucod_leading)))
is_dead10 <- stomach$event == 1
audit_rows <- lapply(levels(stomach$trajectory_5cat_rev), function(group) {
  in_group <- stomach$trajectory_5cat_rev == group
  deaths <- in_group & is_dead10
  broad_cancer <- deaths & !is.na(cause_code) & cause_code == 2L
  non_cancer <- deaths & !is.na(cause_code) & cause_code != 2L
  missing_cause <- deaths & is.na(cause_code)
  data.frame(
    trajectory_5cat_rev = group,
    n = sum(in_group),
    broad_cancer_deaths = sum(broad_cancer),
    non_cancer_competing_deaths = sum(non_cancer),
    total_deaths_10y = sum(deaths),
    cause_missing_deaths = sum(missing_cause),
    stringsAsFactors = FALSE
  )
})
audit <- rev1_bind_rows(audit_rows)
if (sum(audit$n) != nrow(stomach) ||
    sum(audit$total_deaths_10y) != sum(stomach$event) ||
    any(audit$broad_cancer_deaths + audit$non_cancer_competing_deaths +
          audit$cause_missing_deaths != audit$total_deaths_10y) ||
    any(audit$cause_missing_deaths != 0)) {
  stop("Stomach competing-event audit failed internal reconciliation.", call. = FALSE)
}
write_rev1_csv(audit, file.path(OUT, "stomach_competing_event_audit.csv"))
write_rev1_csv(audit, file.path(OUT, "stomach_subgroup_audit.csv"))
write_rev1_provenance(
  OUT, "rev1_p2_stomach_causes", c(COHORT_PATH, POVERTY_PATH)
)
print(flow)
print(audit)
