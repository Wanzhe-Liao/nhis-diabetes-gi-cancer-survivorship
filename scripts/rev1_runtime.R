# Shared runtime bootstrap for revision analyses under R 4.6.1.
# Source this file from the project root before loading analysis packages.

REV1_REQUIRED_R <- "4.6.1"
REV1_EXPECTED_PACKAGES <- c(
  nanoparquet = "0.5.1",
  survey = "4.5",
  survival = "3.8.6"
)

rev1_project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(rev1_project_root, "scripts", "rev1_runtime.R"))) {
  stop("Run revision R scripts from the project root: ", rev1_project_root,
       call. = FALSE)
}

REV1_DEFAULT_COHORT_PATH <- file.path(
  "outputs", "cohort", "analytic_cohort_v4inc2007.parquet"
)
REV1_DEFAULT_OUTPUT_DIR <- file.path(
  "outputs", "revision_round1_v4_r461_sens2007"
)
REV1_LOCKED_OUTPUT_DIR <- file.path("outputs", "revision_round1")

rev1_absolute_path <- function(path, mustWork = FALSE) {
  if (!grepl("^([A-Za-z]:[/\\\\]|/)", path)) {
    path <- file.path(rev1_project_root, path)
  }
  normalizePath(path, winslash = "/", mustWork = mustWork)
}

rev1_cohort_path <- function() {
  path <- Sys.getenv("REV1_COHORT_PATH", REV1_DEFAULT_COHORT_PATH)
  rev1_absolute_path(path, mustWork = TRUE)
}

rev1_output_dir <- function(create = TRUE) {
  path <- Sys.getenv("REV1_OUTPUT_DIR", REV1_DEFAULT_OUTPUT_DIR)
  absolute <- rev1_absolute_path(path, mustWork = FALSE)
  locked <- rev1_absolute_path(REV1_LOCKED_OUTPUT_DIR, mustWork = TRUE)
  if (identical(tolower(absolute), tolower(locked))) {
    stop(
      "The historical output directory is locked and cannot be used by the R 4.6.1 pipeline: ",
      locked,
      call. = FALSE
    )
  }
  if (create) dir.create(absolute, showWarnings = FALSE, recursive = TRUE)
  absolute
}

rev1_source_dir <- function(default = rev1_output_dir(create = FALSE)) {
  rev1_absolute_path(Sys.getenv("REV1_SOURCE_DIR", default), mustWork = TRUE)
}

rev1_support_dir <- function() {
  legacy_source <- Sys.getenv("REV1_SOURCE_DIR", REV1_DEFAULT_OUTPUT_DIR)
  rev1_absolute_path(Sys.getenv("REV1_SUPPORT_DIR", legacy_source), mustWork = TRUE)
}

if (!identical(as.character(getRversion()), REV1_REQUIRED_R)) {
  stop("Revision analyses require R ", REV1_REQUIRED_R,
       "; running ", as.character(getRversion()), ".", call. = FALSE)
}

rev1_local_library <- file.path(rev1_project_root, ".r-lib", REV1_REQUIRED_R)
if (!dir.exists(rev1_local_library)) {
  stop("Project R library is missing. Run scripts/setup_r_4_6_1.R first.",
       call. = FALSE)
}
.libPaths(unique(c(rev1_local_library, .libPaths())))

rev1_package_versions <- vapply(names(REV1_EXPECTED_PACKAGES), function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}, character(1))

bad_packages <- names(REV1_EXPECTED_PACKAGES)[
  is.na(rev1_package_versions) |
    rev1_package_versions != unname(REV1_EXPECTED_PACKAGES)
]
if (length(bad_packages) > 0L) {
  details <- paste0(
    bad_packages, " expected ", REV1_EXPECTED_PACKAGES[bad_packages],
    ", found ", rev1_package_versions[bad_packages]
  )
  stop("R package version check failed: ", paste(details, collapse = "; "),
       ". Run scripts/setup_r_4_6_1.R.", call. = FALSE)
}

read_rev1_parquet <- function(path) {
  if (!file.exists(path)) stop("Parquet input does not exist: ", path, call. = FALSE)

  out <- nanoparquet::read_parquet(path)
  schema <- nanoparquet::read_parquet_schema(path)
  integer_names <- schema$name[
    !is.na(schema$r_col) & schema$type %in% c("INT32", "INT64")
  ]

  # nanoparquet conservatively maps Parquet integers to R doubles. Restore
  # ordinary 32-bit integer columns when every observed value is representable;
  # this matches arrow's classes for this frozen cohort without risking overflow.
  for (name in intersect(integer_names, names(out))) {
    value <- out[[name]]
    finite <- value[!is.na(value)]
    safe_integer <- is.numeric(value) &&
      all(is.finite(finite)) &&
      all(finite == trunc(finite)) &&
      all(finite >= -2147483648 & finite <= .Machine$integer.max)
    if (safe_integer) storage.mode(out[[name]]) <- "integer"
  }

  as.data.frame(out, stringsAsFactors = FALSE)
}

rev1_bind_rows <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  if (length(rows) == 0L) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(columns, names(x))
    for (name in missing) x[[name]] <- NA
    x[, columns, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

rev1_left_join <- function(x, y, by) {
  if (!all(by %in% names(x)) || !all(by %in% names(y))) {
    stop("Join keys are absent from one of the inputs: ", paste(by, collapse = ", "),
         call. = FALSE)
  }
  x_key <- do.call(paste, c(lapply(x[by], as.character), sep = "\r"))
  y_key <- do.call(paste, c(lapply(y[by], as.character), sep = "\r"))
  if (anyDuplicated(y_key)) {
    stop("Right-hand join keys are not unique.", call. = FALSE)
  }
  index <- match(x_key, y_key)
  for (name in setdiff(names(y), by)) x[[name]] <- y[[name]][index]
  x
}

write_rev1_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "", quote = TRUE)
}

write_rev1_provenance <- function(output_dir, script_name, inputs = character()) {
  inputs <- unique(inputs[file.exists(inputs)])
  input_rows <- if (length(inputs) > 0L) {
    data.frame(
      input = vapply(inputs, rev1_absolute_path, character(1), mustWork = TRUE),
      bytes = unname(file.info(inputs)$size),
      md5 = unname(tools::md5sum(inputs)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(input = character(), bytes = numeric(), md5 = character())
  }
  runtime <- rev1_runtime_info()
  runtime$script <- script_name
  runtime$run_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  write_rev1_csv(runtime, file.path(output_dir, paste0(script_name, "_runtime.csv")))
  write_rev1_csv(input_rows, file.path(output_dir, paste0(script_name, "_inputs.csv")))
  invisible(list(runtime = runtime, inputs = input_rows))
}

rev1_runtime_info <- function() {
  data.frame(
    component = c("R", names(REV1_EXPECTED_PACKAGES)),
    version = c(REV1_REQUIRED_R, rev1_package_versions),
    stringsAsFactors = FALSE
  )
}
