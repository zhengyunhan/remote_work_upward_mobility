#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(WeightIt)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(arrow)
})

set.seed(42)

# =========================================================
# Paths
# =========================================================
d <- Sys.getenv("d")
stopifnot(nchar(d) > 0)

data_parquet <- file.path(d, "data", "job_entry_data.parquet")

dep_var <- "weighting_result"
out_dir <- file.path(d, "data", "modeling_outputs", dep_var)
prep_dir <- file.path(out_dir, "analysis_data")
dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Data
# =========================================================
df <- arrow::read_parquet(data_parquet)

covars <- c(
  "f_prob", "ethnicity_predicted", "prev_role_k50", "prev_seniority",
  "prev_wage", "prev_pagerank_z", "prev_parent_headcount", "prev_naics2",
  "prev_metro_area", "entry_year", "country"
)

df_use <- df %>%
  select(
    row_id,
    eligible_sharp,
    all_of(covars)
  ) %>%
  tidyr::drop_na(eligible_sharp, all_of(covars)) %>%
  mutate(
    ethnicity_predicted = factor(ethnicity_predicted),
    prev_role_k50       = factor(prev_role_k50),
    prev_seniority      = factor(prev_seniority),
    prev_naics2         = factor(prev_naics2),
    prev_metro_area     = factor(prev_metro_area),
    entry_year          = factor(entry_year)
  )

req_vars <- c(
  "row_id",
  "country",
  "eligible_sharp",
  "f_prob", "ethnicity_predicted", "prev_role_k50",
  "prev_seniority", "prev_wage", "prev_pagerank_z", "prev_parent_headcount",
  "prev_naics2", "prev_metro_area", "entry_year"
)
stopifnot(all(req_vars %in% names(df_use)))

d0 <- df_use %>%
  select(all_of(req_vars)) %>%
  tidyr::drop_na(
    eligible_sharp, f_prob, ethnicity_predicted, prev_role_k50,
    prev_seniority, prev_wage, prev_pagerank_z, prev_parent_headcount,
    prev_naics2, prev_metro_area, entry_year
  ) %>%
  mutate(
    log_prev_wage             = log(prev_wage + 1),
    log_prev_parent_headcount = log(prev_parent_headcount + 1),
    prev_seniority = factor(
      prev_seniority,
      levels = 1:7,
      labels = c("1 Entry","2 Junior","3 Associate","4 Manager","5 Director","6 Executive","7 SrExecutive"),
      ordered = TRUE
    ),
    prev_naics2         = factor(prev_naics2),
    ethnicity_predicted = factor(ethnicity_predicted),
    prev_role_k50       = factor(prev_role_k50),
    prev_metro_area     = factor(prev_metro_area),
    entry_year          = factor(entry_year)
  ) %>%
  droplevels()

cat("Sample size used for weighting:", nrow(d0), "\n")


ps_formula <- eligible_sharp ~
  f_prob + ethnicity_predicted + prev_role_k50 + prev_seniority +
  log_prev_wage + prev_pagerank_z + log_prev_parent_headcount +
  prev_naics2 + entry_year

w_out <- weightit(ps_formula, data = d0, method = "ps", estimand = "ATT")

weights_out <- d0 %>%
  transmute(
    row_id = row_id,
    w = w_out$weights
  )


readr::write_csv(
  weights_out,
  file.path(prep_dir, "d0_US_clean.csv.gz")
)

cat("Saved weighting file:\n")
cat(file.path(prep_dir, "d0_US_clean.csv.gz"), "\n")