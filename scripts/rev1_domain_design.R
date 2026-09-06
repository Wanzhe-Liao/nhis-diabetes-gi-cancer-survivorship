# rev1_domain_design.R —— NCHS domain/subpopulation 规范的共享 design 构建器。
# 在全部 ELIGSTAT=1 pooled Sample Adults 上建立 survey design, 再 subset 到
# GI survivor domain; 避免先裁剪数据导致 PSU/strata 设计信息丢失
# (GI-only 建库时 75.9% 层仅剩单个 PSU, lonely.psu=remove 使旧 SE 系统性低估)。
# 由 rev1_runtime.R 之后 source。

rev1_fullsample_path <- function() {
  rev1_absolute_path(
    Sys.getenv("REV1_FULLSAMPLE_PATH",
               file.path("outputs", "cohort", "analytic_cohort_v4inc2007_domain_fullsample.parquet")),
    mustWork = TRUE)
}

# NCHS pooling guidance: retain common variance units within each design
# period; make units distinct only across redesigns (2018 Survey Description,
# Appendix IV). Recompute from raw design fields without rewriting cohorts.
REV1_VARIANCE_METHOD <- "nhis_design_periods_1997_2006_2016_raowu_nminus1"
rev1_period_identifiers <- function(data) {
  year <- as.integer(data$year)
  if (anyNA(year) || any(year < 1997L | year > 2018L))
    stop("Unsupported NHIS survey year.", call. = FALSE)
  period <- ifelse(year <= 2005L, "1997-2005",
                   ifelse(year <= 2015L, "2006-2015", "2016-2018"))
  fmt <- function(x) {
    v <- suppressWarnings(as.numeric(as.character(x)))
    if (anyNA(v) || any(!is.finite(v)) || any(v != trunc(v)))
      stop("Missing or noninteger public-use design identifier.", call. = FALSE)
    as.character(as.integer(v))
  }
  data$design_period <- period
  data$design_strata_prefixed <- paste(period, fmt(data$design_strata), sep = ".")
  data$design_psu_prefixed <- paste(data$design_strata_prefixed, fmt(data$design_psu), sep = ".")
  data
}

# cohort_all: 已完成全部派生(含 zero-time 行)的 GI cohort 数据框。
# 返回 combined: 646,201 行, 设计列来自 skeleton, 分析列按 publicid 匹配;
# 非 GI 行分析列为 NA。
rev1_attach_fullsample <- function(cohort_all, fullsample_path = rev1_fullsample_path()) {
  skeleton <- rev1_period_identifiers(read_rev1_parquet(fullsample_path))
  skeleton$year <- as.integer(skeleton$year)
  skeleton$gi_any <- as.logical(skeleton$gi_any)
  overlap <- intersect(names(cohort_all), names(skeleton))
  part <- cohort_all[, setdiff(names(cohort_all), overlap), drop = FALSE]
  idx <- match(skeleton$publicid, cohort_all$publicid)
  if (sum(!is.na(idx)) != nrow(cohort_all)) {
    stop("Full-sample skeleton does not cover every cohort row.", call. = FALSE)
  }
  combined <- cbind(skeleton, part[idx, , drop = FALSE])
  rownames(combined) <- NULL
  combined
}

rev1_domain_population_design <- function(combined) {
  combined <- rev1_period_identifiers(combined)
  counts <- table(unique(combined[c("design_strata_prefixed", "design_psu_prefixed")])$design_strata_prefixed)
  if (any(counts < 2L)) stop("Full-population stratum has fewer than two PSUs.", call. = FALSE)
  survey::svydesign(ids = ~design_psu_prefixed, strata = ~design_strata_prefixed,
                    weights = ~sa_wgt_pool, nest = TRUE, data = combined)
}

# 常用 domain 标记（在 combined 上求值）
rev1_flag_dom <- function(combined) {
  combined$gi_any %in% TRUE & !is.na(combined$time_months) &
    combined$time_months > 0 & !is.na(combined$design_psu) &
    !is.na(combined$design_strata) & !is.na(combined$sa_wgt_pool)
}

# Domain-correct PSU bootstrap 采样器：抽样框架为全样本全部 strata/PSU
# （含无 GI 成员的 PSU），PSU 内只实体化 GI 分析行——估计量分布与
# “全样本重抽样后 subset” 完全一致，但避免每次实体化 64 万行。
rev1_prepare_domain_sampler <- function(base_vars, full_frame = NULL) {
  if (is.null(full_frame)) full_frame <- read_rev1_parquet(rev1_fullsample_path())
  full_frame <- rev1_period_identifiers(full_frame)
  if (!all(base_vars$design_psu_prefixed %in% full_frame$design_psu_prefixed))
    stop("Analysis PSU identifiers do not match the full period design.", call. = FALSE)
  gi_rows_by_psu <- split(seq_len(nrow(base_vars)),
                          as.character(base_vars$design_psu_prefixed), drop = TRUE)
  stratum_psu <- split(as.character(full_frame$design_psu_prefixed),
                       as.character(full_frame$design_strata_prefixed), drop = TRUE)
  lapply(names(stratum_psu), function(st) {
    psu_names <- unique(stratum_psu[[st]])
    if (length(psu_names) < 2L) stop("Cannot resample a singleton full-population stratum.", call. = FALSE)
    list(stratum = st,
         psu_indices = stats::setNames(
           lapply(psu_names, function(p) {
             r <- gi_rows_by_psu[[p]]
             if (is.null(r)) integer(0) else r
           }), psu_names))
  })
}

# Rao-Wu rescaled bootstrap: draw n_h-1 PSUs and multiply each occurrence's
# original weight by n_h/(n_h-1). Apply multiplicity as a replicate weight,
# keeping each original record once (as in survey::subbootweights). Duplicating
# records changes Efron tie handling. Empty GI PSUs remain in the frame.
rev1_resample_domain <- function(base_vars, sampler) {
  indices <- scales <- list()
  k <- 0L
  for (info in sampler) {
    n <- length(info$psu_indices)
    if (n < 2L) stop("Rao-Wu resampling requires at least two PSUs per stratum.", call. = FALSE)
    multiplicity <- tabulate(sample.int(n, n - 1L, replace = TRUE), nbins = n)
    for (j in which(multiplicity > 0L)) {
      k <- k + 1L
      idx <- info$psu_indices[[j]]
      indices[[k]] <- idx
      scales[[k]] <- rep(multiplicity[j] * n / (n - 1), length(idx))
    }
  }
  out <- base_vars[unlist(indices, use.names = FALSE), , drop = FALSE]
  out$sa_wgt_pool <- out$sa_wgt_pool * unlist(scales, use.names = FALSE)
  rownames(out) <- NULL
  out
}
