# build_surveillance.R ------------------------------------------------------
# Turns the cached detection results into index.html: a standalone
# surveillance dashboard that needs no R server.
#
# Author: Lekshmi Rita-Venugopal, Southwest District Health
#
# Order each time the data updates:
#   1. Run_this_first.R        (data.xlsx -> data.rds)
#   2. generate_cache.R        (runs detection, writes cache/)
#   3. source("build_surveillance.R", echo = TRUE)
#
# Needs in the project folder:
#   cache/all_results.rds, cache/method_overview.rds,
#   cache/flagged.rds, cache/meta.rds
#   surveillance_template.html
# ---------------------------------------------------------------------------

PROJECT_DIR <- "~/Outbreak"
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

TEMPLATE  <- "surveillance_template.html"
OUT_HTML  <- "index.html"
CACHE_DIR <- "cache"

# Weeks before this are dropped from the page. The models still use the full
# history; this only limits what the chart can show.
FIRST_YEAR <- 2018

# --- Load the cache --------------------------------------------------------
need <- file.path(CACHE_DIR, c("all_results.rds", "method_overview.rds",
                               "flagged.rds", "meta.rds"))
missing <- need[!file.exists(need)]
if (length(missing)) {
  stop("Missing cache file(s): ", paste(basename(missing), collapse = ", "),
       ". Run generate_cache.R first.", call. = FALSE)
}
if (!file.exists(TEMPLATE)) stop(TEMPLATE, " not found in ", getwd(), call. = FALSE)

all_results <- readRDS(file.path(CACHE_DIR, "all_results.rds"))
overview    <- as.data.table(readRDS(file.path(CACHE_DIR, "method_overview.rds")))
flagged     <- as.data.table(readRDS(file.path(CACHE_DIR, "flagged.rds")))
meta        <- readRDS(file.path(CACHE_DIR, "meta.rds"))

cache_age <- as.numeric(difftime(Sys.time(), meta$generated_at, units = "days"))
if (cache_age > 8) {
  warning("The cache is ", round(cache_age), " days old. Rerun generate_cache.R ",
          "before sharing this page.", call. = FALSE)
}

# --- One compact series per condition --------------------------------------
# Weeks are evenly spaced, so the page gets a start date and the values only.
series_for <- function(dt) {
  dt <- as.data.table(dt)[order(week_start)]
  dt <- dt[data.table::year(week_start) >= FIRST_YEAR]
  if (!nrow(dt)) return(NULL)

  gaps <- unique(as.integer(diff(dt$week_start)))
  even <- length(gaps) == 0 || all(gaps == 7)

  out <- list(
    start = format(min(dt$week_start)),
    cases = I(as.integer(dt$cases)),
    ir    = I(round(as.numeric(dt$IR), 2)),
    upper = I(round(as.numeric(dt$upper_IR), 2)),
    base  = I(round(as.numeric(dt$baseline_IR), 2)),
    alarm = I(as.integer(dt$alarm)),
    train = I(as.integer(!dt$in_test))          # 1 = baseline period
  )
  if (!even) out$weeks <- I(format(dt$week_start))   # fallback if a week is missing
  out
}

series <- lapply(all_results, series_for)
series <- series[!vapply(series, is.null, logical(1))]

# --- Consistency check ------------------------------------------------------
# A flagged week should sit at or above its own threshold. If it doesn't, the
# threshold and the alarms are not measuring the same thing. That happens, for
# example, when glrnb() is called with ret = "value" (which returns the GLR
# test statistic) instead of ret = "cases".
check_consistency <- function() {
  bad <- rbindlist(lapply(names(all_results), function(org) {
    dt <- as.data.table(all_results[[org]])
    # Only meaningful where the upperbound really is a per-week threshold,
    # which is Farrington. GLR-NB alarms build up over several weeks, so a
    # flagged week is not required to clear any single-week line.
    x  <- dt[!is.na(alarm) & alarm == 1L & !is.na(upper_IR) & IR < upper_IR]
    if (!nrow(x)) return(NULL)
    data.table(org = org, method = dt$method[1], n = nrow(x))
  }))

  if (!nrow(bad)) {
    message("Check: every flagged week is at or above its threshold.")
    return(invisible(NULL))
  }

  warning("\n", sum(bad$n), " flagged week(s) across ", nrow(bad),
          " condition(s) fall BELOW their own threshold. The threshold and the ",
          "alarms disagree, so the chart will look wrong.\n",
          "These should all be Farrington Flexible conditions; if a GLR-NB ",
          "condition appears here, its GLR statistic is being stored as a ",
          "threshold. Check extract_results_full_ir() is called with ",
          "is_threshold = FALSE for GLR-NB, then rerun generate_cache.R.\n",
          paste0("  ", bad$org, " (", bad$method, "): ", bad$n, collapse = "\n"),
          call. = FALSE)
  invisible(bad)
}
check_consistency()

# --- Per-condition summary the page lists and searches ---------------------
summary_for <- function(org) {
  dt <- as.data.table(all_results[[org]])[order(week_start)]
  ev <- dt[in_test == TRUE]
  last <- dt[nrow(dt)]
  alarms <- dt[!is.na(alarm) & alarm == 1L]

  list(
    id        = org,
    method    = overview[Organism == org, Recommended_Method][1],
    cases     = sum(dt$cases, na.rm = TRUE),
    alarms    = nrow(alarms),
    lastAlarm = if (nrow(alarms)) format(max(alarms$week_start)) else NA_character_,
    lastCases = as.integer(last$cases),
    lastIR    = round(as.numeric(last$IR), 2),
    lastUpper = if (is.na(last$upper_IR)) NA_real_ else round(as.numeric(last$upper_IR), 2),
    lastBase  = if (is.na(last$baseline_IR)) NA_real_ else round(as.numeric(last$baseline_IR), 2),
    watched   = nrow(ev) > 0
  )
}

conditions <- lapply(names(series), summary_for)

# --- Recent signals, straight from the cache -------------------------------
signals <- if (nrow(flagged)) {
  lapply(seq_len(nrow(flagged)), function(i) {
    r <- flagged[i]
    list(id = r$Organism, week = r$Week, cases = as.integer(r$Cases),
         ir = as.numeric(r$IR), upper = as.numeric(r$Threshold_IR),
         method = r$Method)
  })
} else list()

# --- Build the page --------------------------------------------------------
p <- meta$params
payload <- list(
  dataThrough  = as.character(meta$global_end_week),
  generatedAt  = format(meta$generated_at, "%b %d, %Y"),
  irMult       = p$ir_mult,
  params       = list(alpha = p$alpha, b = p$b, w = p$w, c_ARL = p$c_ARL,
                      dedup = p$dedup_days, jurisdiction = p$jurisdiction),
  conditions   = conditions,
  signals      = signals,
  series       = series
)

json <- toJSON(payload, auto_unbox = TRUE, na = "null", null = "null", digits = NA)

tpl   <- paste(readLines(TEMPLATE, warn = FALSE), collapse = "\n")
parts <- strsplit(tpl, "/*__DATA__*/", fixed = TRUE)[[1]]
if (length(parts) != 2) {
  stop("Could not find the /*__DATA__*/ placeholder in ", TEMPLATE, call. = FALSE)
}
writeLines(paste0(parts[1], json, parts[2]), OUT_HTML, useBytes = TRUE)

# --- Report ----------------------------------------------------------------
message("\nWrote ", normalizePath(OUT_HTML), " (",
        round(file.size(OUT_HTML) / 1024), " KB)")
message("Data through week starting ", meta$global_end_week,
        " | cache generated ", format(meta$generated_at, "%Y-%m-%d"),
        " (", round(cache_age, 1), " days ago)")
message("Conditions on the page: ", length(series),
        " | signals in the last 2 weeks: ", length(signals))

if (length(signals)) {
  message("\nFlagged:")
  for (s in signals) message("  ", s$id, "  week of ", s$week,
                             "  ", s$cases, " cases")
}
