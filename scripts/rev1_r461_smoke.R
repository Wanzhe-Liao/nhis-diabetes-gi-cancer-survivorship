source(file.path("scripts", "rev1_runtime.R"))

cohort_path <- rev1_cohort_path()
cohort <- read_rev1_parquet(cohort_path)

stopifnot(
  nrow(cohort) == 5123L,
  ncol(cohort) >= 104L,
  identical(names(cohort)[1:8],
            c("year", "publicid", "srvy_yr", "hhx", "fmx",
              "person_key", "intv_qrt", "age")),
  identical(class(cohort$year), "integer"),
  all(c("trajectory_6cat", "followup_years", "sa_wgt_pool") %in% names(cohort))
)

cat("R 4.6.1 smoke test passed: ", nrow(cohort), " rows x ",
    ncol(cohort), " columns; clean process exit expected.\n", sep = "")
print(rev1_runtime_info(), row.names = FALSE)
