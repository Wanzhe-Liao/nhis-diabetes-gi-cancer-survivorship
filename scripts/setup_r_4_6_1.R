# One-time setup for the project-local R 4.6.1 package library.

required_r <- "4.6.1"
if (!identical(as.character(getRversion()), required_r)) {
  stop("Use R ", required_r, " to run this setup; running ",
       as.character(getRversion()), ".", call. = FALSE)
}
if (!file.exists(file.path("scripts", "setup_r_4_6_1.R"))) {
  stop("Run this script from the project root.", call. = FALSE)
}

expected <- c(
  nanoparquet = "0.5.1",
  survey = "4.5",
  survival = "3.8.6"
)
local_library <- file.path(getwd(), ".r-lib", required_r)
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(local_library, .libPaths())))

installed_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}
current <- vapply(names(expected), installed_version, character(1))
to_install <- names(expected)[is.na(current) | current != unname(expected)]
if (length(to_install) > 0L) {
  utils::install.packages(
    to_install,
    lib = local_library,
    repos = "https://cloud.r-project.org",
    type = "binary"
  )
}

current <- vapply(names(expected), installed_version, character(1))
bad <- names(expected)[is.na(current) | current != unname(expected)]
if (length(bad) > 0L) {
  stop(
    "Unable to establish the pinned package set: ",
    paste0(bad, " expected ", expected[bad], ", found ", current[bad],
           collapse = "; "),
    call. = FALSE
  )
}

cat("R 4.6.1 environment is ready.\n")
print(data.frame(package = names(expected), version = current, row.names = NULL))
