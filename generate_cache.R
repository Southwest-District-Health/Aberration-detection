

# ==========================================
# generate_cache.R
# ==========================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(surveillance)
  library(MASS)
  library(dplyr)
  library(readxl)
})

# ==========================================
# Make sure we're running from the project folder
# (the folder containing app.R, generate_cache.R, pop.rds)
#
# This works whether you ran source("generate_cache.R") in RStudio,
# clicked "Run" with the file open in the editor, or ran
# Rscript generate_cache.R from a terminal in that folder.
# ==========================================
PROJECT_DIR <- tryCatch({
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    rstudioapi::getActiveDocumentContext()$path
  } else {
    NA_character_
  }
}, error = function(e) NA_character_)

if (is.na(PROJECT_DIR) || !nzchar(PROJECT_DIR)) {
  # Case 2: Rscript from terminal
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file_arg) > 0) PROJECT_DIR <- file_arg
}

if (!is.na(PROJECT_DIR) && nzchar(PROJECT_DIR)) {
  setwd(dirname(normalizePath(PROJECT_DIR)))
}

# Final check: app.R must be in the current working directory.
# If none of the above worked (e.g. you ran source() without the file
# open as an active tab), set PROJECT_DIR manually below and uncomment:
# setwd("C:/Users/lekshmi.venugopal/OneDrive - Southwest District Health/Documents/Outbreak")

if (!file.exists("surveillance_functions.R") && !file.exists("app.R")) {
  stop(
    "Could not find surveillance_functions.R or app.R in the working directory (",
    getwd(), ").\n",
    "Fix: either (a) open generate_cache.R as a file tab in RStudio and use the ",
    "'Source' button (not the console), or (b) uncomment the setwd() line near ",
    "the top of this script and set it to your project folder, e.g.\n",
    "  setwd(\"C:/Users/lekshmi.venugopal/OneDrive - Southwest District Health/Documents/Outbreak\")"
  )
}

cat("Working directory set to:", getwd(), "\n\n")

# ==========================================
# DATA FILE PATHS  <-- update these if files move
# ==========================================
DATA_FOLDER  <- getwd()  # data files live alongside this script by default

# Cases file: can be .xlsx (raw export) OR .rds (already converted)
CASES_FILE   <- file.path(DATA_FOLDER, "data.xlsx")

# Population file: must be .rds with columns 'year' and 'population'
# (or 'pop'). This is the pop.rds already in your project folder.
POP_FILE     <- "pop.rds"

# Where to write the app's cached results (relative to this script)
CACHE_DIR    <- "cache"

# Sheet name in the Excel file (change if different)
CASES_SHEET  <- 1   # use 1 for the first sheet, or a quoted sheet name

# ==========================================
# Source shared functions from app.R
# (CANONICAL_MAP, EXCLUDED_CONDITIONS, pipeline functions, etc.)
# Shiny is NOT loaded -- the ui/server blocks are parsed but not run.
# ==========================================
FUNCTIONS_FILE <- if (file.exists("surveillance_functions.R")) {
  "surveillance_functions.R"
} else {
  "app.R"   # older layout: functions still live inside the Shiny app file
}
suppressMessages(
  sys.source(FUNCTIONS_FILE, envir = environment(), chdir = TRUE)
)

# ---- Sanity check: make sure app.R actually defines everything we need ----
# If this fails, the app.R in this folder is an older/different version than
# expected. Re-save the latest app.R into this project folder and try again.
required_objects <- c(
  "EXCLUDED_CONDITIONS", "CANONICAL_MAP", "canonicalize_organism",
  "to_Date", "make_weekly_series_with_pop", "extend_weekly_to_global_end",
  "make_sts_weekly_with_pop", "run_farrington", "run_glrnb",
  "detect_method", "extract_results_full_ir"
)
missing_objects <- required_objects[!sapply(required_objects, exists, envir = environment())]
if (length(missing_objects) > 0) {
  stop(
    "app.R in this folder is missing required object(s): ",
    paste(missing_objects, collapse = ", "), ".\n",
    "This usually means the app.R file here is an older version. ",
    "Re-save the latest app.R into this project folder (", getwd(), ") and try again."
  )
}
cat("Shared functions loaded from", FUNCTIONS_FILE, "\n")

# ==========================================
# Detection parameters -- must match app.R defaults
# so cached results align with what the app shows
# ==========================================
CACHE_PARAMS <- list(
  b            = 5,
  w            = 3,
  alpha        = 0.05,
  c_ARL        = 5,
  ir_mult      = 1e6,   # incidence rate per 1,000,000 population
  dedup_days   = 30,
  jurisdiction = "Health District 3"
)

dir.create(CACHE_DIR, showWarnings = FALSE)
cat(sprintf("Parameters: b=%d  w=%d  alpha=%.3f  c_ARL=%.1f  ir_mult=%.0f\n",
            CACHE_PARAMS$b, CACHE_PARAMS$w, CACHE_PARAMS$alpha,
            CACHE_PARAMS$c_ARL, CACHE_PARAMS$ir_mult))

# ==========================================
# Load data
# ==========================================
cat("\nLoading population data from:", POP_FILE, "\n")
if (!file.exists(POP_FILE)) stop("Population file not found: ", POP_FILE)
pop <- readRDS(POP_FILE)
if (!is.data.table(pop)) setDT(pop)
setnames(pop, names(pop), tolower(names(pop)))
if (!"population" %in% names(pop) && "pop" %in% names(pop)) setnames(pop, "pop", "population")
pop[, year := as.integer(year)]
pop[, population := as.numeric(population)]
cat(sprintf("  Population loaded: %d years (%d to %d)\n",
            nrow(pop), min(pop$year), max(pop$year)))

cat("Loading cases from:", CASES_FILE, "\n")
if (!file.exists(CASES_FILE)) stop("Cases file not found: ", CASES_FILE)

ext <- tolower(tools::file_ext(CASES_FILE))
if (ext == "rds") {
  df <- readRDS(CASES_FILE)
} else if (ext %in% c("xlsx", "xls")) {
  cat("  Reading Excel file (this may take a moment for large files)...\n")
  df <- as.data.table(readxl::read_excel(CASES_FILE, sheet = CASES_SHEET))
} else {
  stop("Unrecognised file type: ", ext, ". Use .xlsx or .rds.")
}

setDT(df)
setnames(df, names(df), tolower(names(df)))

# Check required columns
required <- c("investigation id", "condition", "report date", "jurisdiction")
missing  <- setdiff(required, names(df))
if (length(missing) > 0) stop("Missing columns in cases file: ", paste(missing, collapse = ", "))

df <- df %>% dplyr::rename(
  patient_id  = "investigation id",
  organism    = "condition",
  report_date = "report date"
)
df[, report_date  := to_Date(report_date)]
df[, organism_std := canonicalize_organism(organism)]
df <- df[jurisdiction == CACHE_PARAMS$jurisdiction]
df <- df[!(organism_std %in% EXCLUDED_CONDITIONS)]

cat(sprintf("  Cases loaded: %d rows for %s\n", nrow(df), CACHE_PARAMS$jurisdiction))

# Episode dedup
setorder(df, patient_id, organism_std, report_date)
df[, lag_date := shift(report_date), by = .(patient_id, organism_std)]
df[, day_diff := as.integer(report_date - lag_date)]
df[, new_ep   := is.na(lag_date) | day_diff > CACHE_PARAMS$dedup_days]
df[, ep_id    := cumsum(new_ep), by = .(patient_id, organism_std)]
df <- df[, .SD[1], by = .(patient_id, organism_std, ep_id)]

orgs            <- sort(unique(df$organism_std))
global_end_week <- lubridate::floor_date(max(df$report_date, na.rm = TRUE),
                                         unit = "week", week_start = 1)

cat(sprintf("  After dedup: %d records, %d conditions\n", nrow(df), length(orgs)))
cat(sprintf("  Dataset end week: %s\n\n", as.character(global_end_week)))

# Also save as .rds alongside the project so the app can load it quickly
# (the app's DEFAULT_CASES_PATH / DEFAULT_POP_PATH point to data.rds / pop.rds
#  in the project root, which shinyapps.io will pick up from the bundle)
cat("Saving data.rds for app to use at runtime...\n")
saveRDS(df, "data.rds", compress = "xz")
cat("  Saved data.rds\n\n")

# ==========================================
# Run detection for every condition
# ==========================================
cat(sprintf("Running detection models for %d conditions...\n", length(orgs)))

all_results  <- list()
method_log   <- list()
flagged_rows <- list()

process_one_condition <- function(org, df, pop, global_end_week, params) {
  x <- df[organism_std == org]
  x[, organism := organism_std]
  
  weekly_dt     <- make_weekly_series_with_pop(x, org, pop, ir_mult = params$ir_mult)
  weekly_dt_ext <- extend_weekly_to_global_end(weekly_dt, global_end_week,
                                               pop, ir_mult = params$ir_mult)
  
  method       <- detect_method(weekly_dt, b = params$b)
  n_cases      <- nrow(x)
  method_label <- switch(method,
                         farrington   = "Farrington Flexible",
                         glrnb        = "GLR-NB",
                         insufficient = "Insufficient data",
                         method
  )
  
  if (method == "insufficient") {
    result_dt <- copy(weekly_dt_ext)[order(week_start)]
    result_dt[, `:=`(upper_cases = NA_real_, upper_IR = NA_real_,
                     baseline_cases = NA_real_, baseline_IR = NA_real_,
                     alarm = NA_integer_, method = "insufficient",
                     in_test = FALSE, row_status = "insufficient_data")]
    cat(sprintf("insufficient (%d cases)\n", n_cases))
    return(list(org = org, n_cases = n_cases, method = method_label,
                flagged = NULL, result_dt = result_dt))
  }
  
  sts_obj <- make_sts_weekly_with_pop(weekly_dt_ext)
  ff <- if (method == "glrnb") {
    run_glrnb(sts_obj, b = params$b, c_ARL = params$c_ARL)
  } else {
    run_farrington(sts_obj, alpha = params$alpha, b = params$b, w = params$w)
  }
  
  # GLR-NB's upperbound is the GLR statistic, not a case threshold, so it is
  # dropped; its fitted baseline is carried through instead.
  result_dt <- extract_results_full_ir(sts_obj, ff$res, weekly_dt_ext,
                                       ff$range, ir_mult = params$ir_mult,
                                       baseline = ff$baseline,
                                       is_threshold = (method != "glrnb"))
  result_dt[, method := method]
  in_test <- seq_len(nrow(result_dt)) %in% ff$range
  result_dt[, in_test    := in_test]
  result_dt[, row_status := fifelse(!in_test, "training",
                                    fifelse(is.na(alarm), "eval_na", "eval_ok"))]
  result_dt <- result_dt[order(week_start)]
  
  n_alarms <- sum(result_dt$alarm == 1L, na.rm = TRUE)
  
  # Flagged in last 2 weeks?
  recent  <- result_dt[week_start %in% c(global_end_week, global_end_week - 7)]
  flagged <- recent[!is.na(alarm) & alarm == 1L]
  flagged_entry <- if (nrow(flagged) > 0) {
    data.frame(
      Organism     = org,
      Week         = format(flagged$week_start, "%Y-%m-%d"),
      Cases        = flagged$cases,
      IR           = round(flagged$IR, 2),
      Threshold_IR = round(flagged$upper_IR, 2),
      Baseline_IR  = round(flagged$baseline_IR, 2),
      Method       = method_label,
      stringsAsFactors = FALSE
    )
  } else NULL
  
  flag_str <- if (!is.null(flagged_entry)) sprintf(" *** FLAGGED (%d week(s))", nrow(flagged_entry)) else ""
  cat(sprintf("%s  [%d alarms]%s\n", method_label, n_alarms, flag_str))
  
  list(org = org, n_cases = n_cases, method = method_label,
       flagged = flagged_entry, result_dt = result_dt)
}

for (i in seq_along(orgs)) {
  org <- orgs[i]
  cat(sprintf("  [%d/%d] %-55s", i, length(orgs), substr(org, 1, 55)))
  
  result <- tryCatch(
    process_one_condition(org, df, pop, global_end_week, CACHE_PARAMS),
    error = function(e) {
      cat(sprintf("ERROR: %s\n", conditionMessage(e)))
      list(org = org, n_cases = NA_integer_, method = "error",
           flagged = NULL, result_dt = NULL)
    }
  )
  
  if (!is.null(result$result_dt)) all_results[[org]] <- result$result_dt
  method_log[[org]] <- result
  if (!is.null(result$flagged)) flagged_rows[[org]] <- result$flagged
}

# ==========================================
# Save cache files
# ==========================================
cat("\nWriting cache files...\n")

saveRDS(all_results, file.path(CACHE_DIR, "all_results.rds"), compress = "xz")
cat(sprintf("  cache/all_results.rds   -- %d conditions\n", length(all_results)))

method_overview_df <- do.call(rbind, lapply(method_log, function(m) {
  data.frame(Organism = m$org, Recommended_Method = m$method,
             Cases_all_years = coalesce(m$n_cases, NA_integer_),
             stringsAsFactors = FALSE)
}))
saveRDS(method_overview_df, file.path(CACHE_DIR, "method_overview.rds"), compress = "xz")
cat(sprintf("  cache/method_overview.rds -- %d rows\n", nrow(method_overview_df)))

flagged_df <- if (length(flagged_rows) > 0) do.call(rbind, flagged_rows) else
  data.frame(Organism = character(0), Week = character(0), Cases = integer(0),
             IR = numeric(0), Threshold_IR = numeric(0), Method = character(0),
             stringsAsFactors = FALSE)
saveRDS(flagged_df, file.path(CACHE_DIR, "flagged.rds"), compress = "xz")
cat(sprintf("  cache/flagged.rds         -- %d flagged condition-weeks\n", nrow(flagged_df)))

cache_meta <- list(
  generated_at    = Sys.time(),
  global_end_week = as.character(global_end_week),
  params          = CACHE_PARAMS,
  n_conditions    = length(orgs),
  orgs            = orgs
)
saveRDS(cache_meta, file.path(CACHE_DIR, "meta.rds"), compress = "xz")
cat("  cache/meta.rds\n")

# ==========================================
# Summary
# ==========================================
n_flagged <- nrow(flagged_df)
n_ff      <- sum(sapply(method_log, function(m) m$method == "Farrington Flexible"), na.rm = TRUE)
n_glrnb   <- sum(sapply(method_log, function(m) m$method == "GLR-NB"), na.rm = TRUE)
n_insuf   <- sum(sapply(method_log, function(m) m$method == "Insufficient data"), na.rm = TRUE)

cat(sprintf("
Done.
  %d conditions cached
  %d Farrington Flexible  |  %d GLR-NB  |  %d insufficient data
  %d condition-week(s) flagged in last 2 weeks
  Generated: %s

Next step:
  Run build_surveillance.R to rebuild index.html from this cache.
",
            length(all_results), n_ff, n_glrnb, n_insuf,
            n_flagged,
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
))


rsconnect::deployApp()