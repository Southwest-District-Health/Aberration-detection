# surveillance_functions.R --------------------------------------------------
# Shared functions for the Health District 3 surveillance pipeline:
# name canonicalization, weekly series building, Farrington Flexible and
# GLR-NB detection, automatic method selection, and result extraction.
#
# generate_cache.R sources this file. Nothing here needs Shiny, and nothing
# here is deployed anywhere: it all runs locally.
# ---------------------------------------------------------------------------

data.table::setDTthreads(1)   # keeps memory use predictable on large exports

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)     # use as lubridate::year/isoyear/isoweek/floor_date
  library(surveillance)
  library(MASS)           # glm.nb, used to auto-detect seasonality per condition
  library(ggplot2)
  library(dplyr)
})

# ---- Default data files (RDS) ----
DEFAULT_CASES_PATH <- "data.rds"
DEFAULT_POP_PATH   <- "pop.rds"

# ---- Helper: safe date ----
to_Date <- function(x) {
  if (inherits(x, "Date")) return(x)
  as.Date(x)
}

# ==========================================
# Canonicalization (normalize Condition names)
# Merges known duplicate/legacy labels (e.g. pre/post-2018
# Salmonellosis labels, renamed TB reporting category) so they
# are treated as one condition. Anything not listed here passes
# through unchanged (trimmed/whitespace-collapsed), so EVERY
# condition present in the data is usable -- this is not a
# curated allow-list anymore.
# ==========================================
CANONICAL_MAP <- c(
  "campylobacteriosis" = "Campylobacteriosis",
  "campylobacter" = "Campylobacteriosis",
  "campy" = "Campylobacteriosis",
  
  "cryptosporidiosis" = "Cryptosporidiosis",
  "crypto" = "Cryptosporidiosis",
  "cryptosporidium" = "Cryptosporidiosis",
  
  "giardiasis" = "Giardiasis",
  "giardia" = "Giardiasis",
  
  "norovirus" = "Noroviruses",
  "noroviruses" = "Noroviruses",
  "norovirus infection" = "Noroviruses",
  "noro" = "Noroviruses",
  
  "pertussis" = "Pertussis",
  "whooping cough" = "Pertussis",
  
  "rabies pep" = "Rabies PEP",
  "rabies post-exposure prophylaxis" = "Rabies PEP",
  "rabies, pep" = "Rabies PEP",
  "pep rabies" = "Rabies PEP",
  
  "rabies, animal" = "Rabies, animal",
  "rabies (animal)" = "Rabies, animal",
  "animal rabies" = "Rabies, animal",
  
  "respiratory syncytial virus (rsv)" = "Respiratory syncytial virus (RSV)",
  "respiratory syncytial virus" = "Respiratory syncytial virus (RSV)",
  "rsv" = "Respiratory syncytial virus (RSV)",
  
  "shiga toxin-producing escherichia coli (stec)" = "Shiga toxin-producing Escherichia coli (STEC)",
  "stec" = "Shiga toxin-producing Escherichia coli (STEC)",
  "e. coli stec" = "Shiga toxin-producing Escherichia coli (STEC)",
  "shiga toxin producing e. coli" = "Shiga toxin-producing Escherichia coli (STEC)",
  "shiga toxin-producing e. coli" = "Shiga toxin-producing Escherichia coli (STEC)",
  
  # ---- legacy / renamed reporting categories seen in the export ----
  "salmonellosis (excl s. typhi and s. paratyphi)" = "Salmonellosis (excl S. Typhi and S. Paratyphi)",
  "salmonellosis - prior to 2018" = "Salmonellosis (excl S. Typhi and S. Paratyphi)",
  "salmonellosis 2018 (excl paratyphoid and typhoid)" = "Salmonellosis (excl S. Typhi and S. Paratyphi)",
  
  "tuberculosis" = "Tuberculosis",
  "tuberculosis (2020 rvct)" = "Tuberculosis",
  
  "streptococcus pneumoniae, invasive disease (ipd)" = "Streptococcus pneumoniae, invasive disease (IPD)",
  "strep pneumoniae, invasive" = "Streptococcus pneumoniae, invasive disease (IPD)"
)

canonicalize_organism <- function(x_raw) {
  x <- tolower(trimws(x_raw))
  x <- gsub("\\s+", " ", x)
  x <- gsub("[\\.;:]+$", "", x)
  
  mapped <- CANONICAL_MAP[x]
  out <- ifelse(is.na(mapped), trimws(x_raw), mapped)
  out
}

# ==========================================
# Conditions excluded from the app entirely (not surveillance
# targets for this dashboard, not comparable weekly counts, or
# not relevant for aberration detection).
# ==========================================
EXCLUDED_CONDITIONS <- c(
  "2019-nCoV",
  "Severe reactions to any immunization",
  "Phenylketonuria (PKU) - positive",
  "Congenital hypothyroidism",
  "Influenza",
  "Influenza Outbreak",
  "Extraordinary occurrence of illness"
)


# ==========================================
# Core pipeline functions
# ==========================================
make_weekly_series_with_pop <- function(dt, org, pop_dt,
                                        week_start_day = 1,
                                        ir_mult = 1e6) {
  x <- dt[organism == org]
  if (nrow(x) == 0) stop(sprintf("No rows found for organism = '%s'", org))
  
  x[, week_start := lubridate::floor_date(as.Date(report_date), unit = "week", week_start = week_start_day)]
  wk <- x[, .(cases = .N), by = .(week_start)][order(week_start)]
  
  all_weeks <- data.table(week_start = seq(min(wk$week_start), max(wk$week_start), by = "week"))
  out <- merge(all_weeks, wk, by = "week_start", all.x = TRUE)
  out[is.na(cases), cases := 0L]
  
  out[, year := lubridate::isoyear(week_start)]
  pop <- copy(pop_dt)[, .(year = as.integer(year), population = as.numeric(population))]
  out <- merge(out, pop, by = "year", all.x = TRUE)
  
  if (anyNA(out$population)) {
    miss <- unique(out[is.na(population), year])
    stop("Missing population for years: ", paste(miss, collapse = ", "))
  }
  
  out[, IR := (cases / population) * ir_mult]
  out[, organism := org]
  setcolorder(out, c("organism", "week_start", "year", "cases", "population", "IR"))
  out[]
}

# ==========================================
# Pad a condition's weekly series with zero-case weeks out to a
# shared, dataset-wide ending week. Without this, "the most recent
# 2 weeks" for a condition whose last case was reported months ago
# would be those old weeks, not the dataset's actual most recent
# weeks, which is misleading for the flagged-conditions check.
# ==========================================
extend_weekly_to_global_end <- function(weekly_dt, global_end_week, pop_dt, ir_mult = 1e6) {
  weekly_dt <- copy(weekly_dt)[order(week_start)]
  cur_max <- max(weekly_dt$week_start)
  if (cur_max >= global_end_week) return(weekly_dt)
  
  new_weeks <- seq(cur_max + 7, global_end_week, by = "week")
  if (length(new_weeks) == 0) return(weekly_dt)
  
  org <- weekly_dt$organism[1]
  pad <- data.table(week_start = new_weeks, cases = 0L)
  pad[, year := lubridate::isoyear(week_start)]
  pop <- copy(pop_dt)[, .(year = as.integer(year), population = as.numeric(population))]
  pad <- merge(pad, pop, by = "year", all.x = TRUE)
  pad[, IR := (cases / population) * ir_mult]
  pad[, organism := org]
  setcolorder(pad, names(weekly_dt))
  
  rbind(weekly_dt, pad)[order(week_start)]
}

make_sts_weekly_with_pop <- function(weekly_dt) {
  weekly_dt <- copy(weekly_dt)[order(week_start)]
  if (nrow(weekly_dt) == 0) stop("weekly_dt has 0 rows.")
  
  start_year <- lubridate::isoyear(weekly_dt$week_start[1])
  start_week <- lubridate::isoweek(weekly_dt$week_start[1])
  
  y_obs <- matrix(weekly_dt$cases,      ncol = 1)
  y_pop <- matrix(weekly_dt$population, ncol = 1)
  
  colnames(y_obs) <- "series1"
  colnames(y_pop) <- "series1"
  
  surveillance::sts(
    observed   = y_obs,
    start      = c(start_year, start_week),
    frequency  = 52,
    population = y_pop
  )
}

# ==========================================
# Shared evaluation window
#
# Both methods need b years of history before they can judge a week, and
# Farrington also needs a margin either side of the matching calendar week,
# so the first testable week is (b + 1) years in. Using ONE definition for
# both methods keeps the baseline period, the alarm counts and the charts
# comparable across conditions no matter which method was chosen.
# ==========================================
test_range <- function(sts_obj, b = 5) {
  n    <- length(surveillance::observed(sts_obj))
  freq <- frequency(sts_obj)
  start_test <- (b + 1) * freq + 1
  if (start_test > n) {
    stop("Not enough history: need more than ", b + 1, " years of weeks, have ",
         round(n / freq, 1), " years.")
  }
  start_test:n
}

run_farrington <- function(sts_obj, alpha = 0.05, b = 5, w = 3) {
  rng <- test_range(sts_obj, b = b)
  
  res <- surveillance::farringtonFlexible(
    sts_obj,
    control = list(
      range = rng,
      b = b,
      w = w,
      reweight = TRUE,
      weightsThreshold = 3,
      alpha = alpha,
      trend = TRUE,
      pThresholdTrend = 1,
      limit54 = c(5, 4),
      powertrans = "2/3",
      thresholdMethod = "delta",
      populationOffset = TRUE
    )
  )
  
  list(res = res, range = rng)
}

# ==========================================
# GLR-NB (negative binomial generalized likelihood ratio)
# Better fit than Farrington Flexible for diseases with NO
# seasonal pattern and/or low, sparse counts (e.g. Rabies PEP),
# because it does not try to fit a seasonal baseline and does
# not require a multi-year reference window per calendar week.
# mu0 = list(S = 0, trend = FALSE) tells it to fit a FLAT
# (non-seasonal) negative-binomial baseline using the b years
# of history before the test range, instead of week-matched
# historical windows like Farrington uses.
#
# Note on what this returns: GLR-NB has no single-week threshold
# the way Farrington does. It accumulates evidence week to week and
# alarms when the accumulated statistic crosses c.ARL, so a week can
# be flagged on the strength of the weeks before it.
#
# ret = "cases" would return the count that would have triggered an
# alarm that week, but surveillance computes it brute force and it is
# far too slow across every condition. So we use ret = "value" (the
# GLR statistic, which is fast) and separately return the fitted
# baseline rate, mu0. The statistic is NOT a case threshold and must
# never be charted as one; the baseline is the honest comparison line.
# ==========================================
run_glrnb <- function(sts_obj, b = 5, c_ARL = 5) {
  rng <- test_range(sts_obj, b = b)
  
  res <- surveillance::glrnb(
    sts_obj,
    control = list(
      range     = rng,
      c.ARL     = c_ARL,        # detection threshold (lower = more sensitive/more alarms)
      mu0       = list(S = 0, trend = FALSE),  # flat, non-seasonal NB baseline
      change    = "intercept",  # detect a step increase in the rate
      dir       = "inc",        # only flag increases (not decreases)
      ret       = "value"       # fast; see note above
    )
  )
  
  # The fitted in-control mean for each tested week. surveillance stores it
  # back into the control list once fitted; guard in case that changes.
  mu0 <- tryCatch({
    m <- res@control$mu0
    if (is.list(m)) m <- m$mu0
    m <- as.numeric(m)
    if (length(m) == 1L) m <- rep(m, length(rng))
    if (length(m) == length(rng)) m else NULL
  }, error = function(e) NULL)
  
  if (is.null(mu0)) {
    warning("GLR-NB baseline (mu0) not available; the chart will show no ",
            "baseline line for this condition.", call. = FALSE)
  }
  
  list(res = res, range = rng, baseline = mu0)
}

# ==========================================
# Automatic method selection ("Auto" mode)
#
# For each condition, fit two negative-binomial models on its
# weekly series:
#   flat model:     cases ~ offset(log(population))
#   seasonal model: cases ~ sin(2*pi*woy) + cos(2*pi*woy) + offset(log(population))
# and compare them with a likelihood-ratio test. If adding the
# seasonal (harmonic) terms significantly improves fit, the
# condition has a detectable seasonal pattern -> Farrington
# Flexible. Otherwise -> GLR-NB (flat baseline).
#
# This replaces a hand-curated list of "seasonal" diseases, which
# doesn't scale once every condition in the data is selectable.
#
# If there isn't enough history to fit either model reliably
# (too few years or too few total cases), returns "insufficient"
# so the app can fall back to showing counts/IR without alarms
# rather than producing an unreliable or crashing model fit.
# ==========================================
detect_method <- function(weekly_dt, b = 5, min_years = 2, min_cases = 24) {
  n_years     <- length(unique(weekly_dt$year))
  total_cases <- sum(weekly_dt$cases, na.rm = TRUE)
  
  if (n_years < min_years || total_cases < min_cases) {
    return("insufficient")
  }
  
  # Neither method can judge a week until (b + 1) years of history exist
  # (see test_range()). Without this check a short series produced a
  # descending, meaningless evaluation range instead of a clean skip.
  if (nrow(weekly_dt) <= (b + 1) * 52) {
    return("insufficient")
  }
  
  wk <- copy(weekly_dt)
  woy <- (lubridate::isoweek(wk$week_start) - 1) / 52
  wk[, sin1 := sin(2 * pi * woy)]
  wk[, cos1 := cos(2 * pi * woy)]
  
  fit_flat <- tryCatch(
    MASS::glm.nb(cases ~ offset(log(population)), data = wk),
    error = function(e) NULL, warning = function(w) NULL
  )
  fit_season <- tryCatch(
    MASS::glm.nb(cases ~ sin1 + cos1 + offset(log(population)), data = wk),
    error = function(e) NULL, warning = function(w) NULL
  )
  
  # If either model fails to fit (e.g. too sparse/degenerate), default
  # to the more conservative, less data-hungry choice: GLR-NB.
  if (is.null(fit_flat) || is.null(fit_season)) return("glrnb")
  
  lr <- tryCatch(anova(fit_flat, fit_season), error = function(e) NULL)
  pval <- tryCatch(lr[["Pr(Chi)"]][2], error = function(e) NA_real_)
  
  if (!is.na(pval) && pval < 0.05) "farrington" else "glrnb"
}

# baseline: optional vector of fitted in-control means, one per week in rng
#           (GLR-NB). is_threshold = FALSE means the algorithm's upperbound is
#           NOT a case threshold, so it is dropped rather than charted.
extract_results_full_ir <- function(sts_obj, res_obj, weekly_dt, rng, ir_mult = 1e6,
                                    baseline = NULL, is_threshold = TRUE) {
  weekly_dt <- copy(weekly_dt)[order(week_start)]
  n <- length(surveillance::observed(sts_obj))
  stopifnot(nrow(weekly_dt) == n)
  
  obs_full <- as.integer(surveillance::observed(sts_obj))
  ub_raw   <- as.numeric(surveillance::upperbound(res_obj))
  al_raw   <- as.integer(surveillance::alarms(res_obj))
  
  ub_full <- rep(NA_real_, n)
  al_full <- rep(NA_integer_, n)
  
  if (length(ub_raw) == n) {
    ub_full <- ub_raw
  } else if (length(ub_raw) == length(rng)) {
    ub_full[rng] <- ub_raw
  } else {
    stop(sprintf("upperbound length unexpected: %d (n=%d, range=%d)",
                 length(ub_raw), n, length(rng)))
  }
  
  if (length(al_raw) == n) {
    al_full <- al_raw
  } else if (length(al_raw) == length(rng)) {
    al_full[rng] <- al_raw
  } else {
    stop(sprintf("alarms length unexpected: %d (n=%d, range=%d)",
                 length(al_raw), n, length(rng)))
  }
  
  pop <- weekly_dt$population
  
  if (!is_threshold) ub_full <- rep(NA_real_, n)
  
  base_full <- rep(NA_real_, n)
  if (!is.null(baseline) && length(baseline) == length(rng)) {
    base_full[rng] <- as.numeric(baseline)
  }
  
  data.table(
    organism      = weekly_dt$organism,
    week_start    = weekly_dt$week_start,
    cases         = obs_full,
    population    = pop,
    IR            = (obs_full / pop) * ir_mult,
    upper_cases   = ub_full,
    upper_IR      = (ub_full   / pop) * ir_mult,
    baseline_cases = base_full,
    baseline_IR   = (base_full / pop) * ir_mult,
    alarm         = al_full
  )
}

plot_farr_ir <- function(dt, title) {
  dt <- data.table::as.data.table(dt)
  if (!"alarm" %in% names(dt)) stop("Column 'alarm' not found in dt.")
  dt[, alarm := as.integer(alarm)]
  
  ggplot(dt, aes(x = week_start)) +
    geom_line(aes(y = IR), linewidth = 0.6) +
    geom_line(aes(y = upper_IR), linewidth = 0.6, linetype = "dashed", color = "pink3") +
    geom_point(data = dt[alarm == 1L], aes(y = IR), color = "red", size = 2) +
    labs(
      title = title,
      x = "Week start",
      y = "Incidence rate (per 1,000,000)",
      subtitle = "IR (solid), Farrington IR threshold (dashed), alarms (red)"
    ) +
    theme_minimal()
}

plot_farr_ir_mark_all <- function(dt, title) {
  dt <- data.table::as.data.table(dt)
  if (!"alarm" %in% names(dt)) stop("Column 'alarm' not found in dt.")
  dt[, alarm := as.integer(alarm)]
  
  xmin <- min(dt$week_start, na.rm = TRUE)
  xmax <- max(dt$week_start, na.rm = TRUE)
  month_breaks <- seq(lubridate::floor_date(xmin, unit = "month"),
                      lubridate::floor_date(xmax, unit = "month"),
                      by = "1 month")
  
  ggplot(dt, aes(x = week_start)) +
    geom_line(aes(y = IR), linewidth = 0.6, color = "steelblue3") +
    geom_point(aes(y = IR), color = "steelblue4", size = 1.8, alpha = 0.9) +
    geom_point(data = dt[alarm == 1L], aes(y = IR), color = "red3", size = 3.2) +
    scale_x_date(breaks = month_breaks, date_labels = "%b\n%Y",
                 expand = expansion(mult = c(0.01, 0.02))) +
    labs(
      title = title,
      x = "",
      y = "Weekly Incidence rate",
      subtitle = "Observed IR (blue); alarms (red)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 9))
}

# ==========================================
# Data loaders (RDS-based; robust to column names)
# ==========================================
load_population <- function(path = DEFAULT_POP_PATH) {
  if (!file.exists(path)) stop("Population file not found: ", path, call. = FALSE)
  p <- readRDS(path)
  if (!is.data.table(p)) setDT(p)
  setnames(p, old = names(p), new = tolower(names(p)))
  if (!"population" %in% names(p) && "pop" %in% names(p)) {
    setnames(p, "pop", "population")
  }
  req_cols <- c("year", "population")
  miss <- setdiff(req_cols, names(p))
  if (length(miss)) stop("Population data missing columns: ", paste(miss, collapse = ", "))
  p[, year := as.integer(year)]
  p[, population := as.numeric(population)]
  p[]
}

load_cases <- function(path = DEFAULT_CASES_PATH, filter_juris = TRUE) {
  if (!file.exists(path)) stop("Cases file not found: ", path, call. = FALSE)
  df <- readRDS(path)
  if (!is.data.frame(df)) stop("Cases RDS must contain a data.frame or data.table")
  setDT(df)
  setnames(df, old = names(df), new = tolower(names(df)))
  exp_cols <- c("investigation id","condition","report date","jurisdiction")
  missing <- setdiff(exp_cols, names(df))
  if (length(missing)) stop("Cases data missing columns: ", paste(missing, collapse = ", "))
  
  if (isTRUE(filter_juris)) {
    df <- df[jurisdiction == "Health District 3"]
  }
  df <- df %>% dplyr::rename(
    patient_id  = "investigation id",
    organism    = "condition",
    report_date = "report date"
  )
  df[, report_date := to_Date(report_date)]
  df[, organism_std := canonicalize_organism(organism)]
  df <- df[!(organism_std %in% EXCLUDED_CONDITIONS)]
  df[]
}
