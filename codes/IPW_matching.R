#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(WeightIt)
  library(fixest)
  library(cobalt)
  library(pROC)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(arrow)
  library(tibble)
  library(stringr)
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
dir.create(file.path(out_dir, "first_stage"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "love_plots"), showWarnings = FALSE, recursive = TRUE)

prep_dir <- file.path(out_dir, "analysis_data")
dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Data
# =========================================================
df <- arrow::read_parquet(data_parquet)

covars <- c(
  "f_prob","ethnicity_predicted","prev_role_k50","prev_seniority",
  "prev_wage","prev_pagerank_z","prev_parent_headcount","prev_naics2",
  "prev_metro_area","entry_year","country"
)

df_use <- df %>%
  select(
    row_id,
    log_wage_growth,
    promotion_at_entry,
    cross_metro_move,
    curr_parent_headcount,
    curr_pagerank_z,
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
  "f_prob","ethnicity_predicted","prev_role_k50",
  "prev_seniority","prev_wage","prev_pagerank_z","prev_parent_headcount",
  "prev_naics2","prev_metro_area","entry_year",
  "log_wage_growth","promotion_at_entry","cross_metro_move",
  "curr_parent_headcount","curr_pagerank_z"
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
    prev_naics2 = factor(prev_naics2),
    ethnicity_predicted = factor(ethnicity_predicted),
    prev_role_k50 = factor(prev_role_k50),
    prev_metro_area = factor(prev_metro_area),
    entry_year = factor(entry_year)
  ) %>%
  droplevels()

cat("US sample size:", nrow(d0), "\n")

ps_formula <- eligible_sharp ~
  f_prob + ethnicity_predicted + prev_role_k50 + prev_seniority +
  log_prev_wage + prev_pagerank_z + log_prev_parent_headcount +
  prev_naics2 + entry_year

pretty_var <- function(x) {
  x <- gsub("^prev_role_k50_",               "Previous role: ",                    x)
  x <- gsub("^ethnicity_predicted_",         "Ethnicity: ",                        x)
  x <- gsub("^prev_naics2_",                 "NAICS-2: ",                          x)
  x <- gsub("^prev_metro_area_",             "Previous metro area: ",              x)
  x <- gsub("^prev_seniority_",              "Previous seniority: ",               x)
  x <- gsub("^entry_year_",                  "Entry year: ",                       x)
  x <- gsub("^prev_pagerank_z$",             "Previous PageRank (z)",              x)
  x <- gsub("^log_prev_parent_headcount$",   "Log(1 + prev parent headcount)",     x)
  x <- gsub("^log_prev_wage$",               "Log(1 + previous wage)",             x)
  x <- gsub("^prev_parent_headcount$",       "Prev parent headcount (raw)",        x)
  x <- gsub("^prev_wage$",                   "Previous wage (raw)",                x)
  x <- gsub("^f_prob$",                      "Probability of being female",        x)
  x <- gsub("^prop.score$",                  "Propensity score",                   x)
  x <- gsub("_", " ", x); x
}

save_balance_and_love <- function(wobj, dd, scope_label, out_dir, pretty_var_fun) {
  bt <- bal.tab(
    wobj,
    un = TRUE,
    stats = "mean.diffs",
    abs = TRUE,
    s.d.denom = "treated",
    binary = "std",
    m.threshold = 0.1
  )

  bal_df <- tibble::as_tibble(bt$Balance, rownames = "covariate")

  w <- wobj$weights
  g1 <- dd$eligible_sharp == 1
  ess_overall <- (sum(w)^2) / sum(w^2)
  ess_treated <- (sum(w[g1])^2) / sum(w[g1]^2)
  ess_control <- (sum(w[!g1])^2) / sum(w[!g1]^2)

  max_abs_smd_un  <- max(bal_df$Diff.Un,  na.rm = TRUE)
  max_abs_smd_adj <- max(bal_df$Diff.Adj, na.rm = TRUE)
  all_balanced    <- max_abs_smd_adj <= 0.1

  summary_row <- tibble(
    covariate       = ".summary",
    Diff.Un         = NA_real_, Diff.Adj = NA_real_,
    V.Ratio.Un      = NA_real_, V.Ratio.Adj = NA_real_,
    max_abs_smd_un  = max_abs_smd_un,
    max_abs_smd_adj = max_abs_smd_adj,
    all_balanced    = all_balanced,
    ess_overall     = ess_overall,
    ess_treated     = ess_treated,
    ess_control     = ess_control,
    n_obs           = nrow(dd)
  )

  readr::write_csv(
    bind_rows(summary_row, bal_df),
    file.path(out_dir, "first_stage", sprintf("balance_%s.csv", scope_label))
  )

  p <- love.plot(
    wobj,
    stats = "mean.diffs",
    abs = TRUE,
    thresholds = c(m = .1),
    s.d.denom = "treated",
    binary = "std",
    var.order = "unadjusted",
    line = TRUE,
    sample.names = c("Unweighted","Weighted"),
    title = sprintf("Covariate Balance — %s", scope_label)
  ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position="none",
      axis.text.y = element_text(size=11, hjust=1),
      axis.text.x = element_text(size=12),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face="bold", size=18, hjust=0.5)
    ) +
    labs(x = "Absolute Standardized Mean Difference", y = NULL) +
    scale_y_discrete(labels = pretty_var_fun)

  n_vars <- nrow(bal_df); h_in <- max(8, 0.22 * n_vars)
  ggsave(
    file.path(out_dir, "love_plots", sprintf("love_%s.png", scope_label)),
    p, width = 8.5, height = h_in, dpi = 600
  )
  ggsave(
    file.path(out_dir, "love_plots", sprintf("love_%s.tiff", scope_label)),
    p, width = 8.5, height = h_in, dpi = 600, compression = "lzw"
  )

  list(
    ess_overall = ess_overall, ess_treated = ess_treated, ess_control = ess_control,
    max_abs_smd_un = max_abs_smd_un, max_abs_smd_adj = max_abs_smd_adj,
    all_balanced = all_balanced
  )
}

w_out_all <- weightit(ps_formula, data = d0, method = "ps", estimand = "ATT")
d0$w  <- w_out_all$weights
d0$ps <- w_out_all$ps

all_stats <- save_balance_and_love(w_out_all, d0, "All", out_dir, pretty_var)

ps_glm <- glm(ps_formula, data = d0, family = binomial())
p_hat  <- as.numeric(predict(ps_glm, type = "response"))
y      <- as.numeric(d0$eligible_sharp)
w      <- d0$w
auc_val <- tryCatch(as.numeric(pROC::auc(y, p_hat)), error = function(e) NA_real_)
brier   <- mean((y - p_hat)^2, na.rm = TRUE)
brier_w <- sum(w * (y - p_hat)^2, na.rm = TRUE) / sum(w, na.rm = TRUE)

write_csv(
  tibble(
    scope = "All",
    auc = auc_val, brier = brier, brier_weighted = brier_w,
    ess_overall = all_stats$ess_overall, ess_treated = all_stats$ess_treated, ess_control = all_stats$ess_control,
    max_abs_smd_un = all_stats$max_abs_smd_un, max_abs_smd_adj = all_stats$max_abs_smd_adj,
    all_balanced = all_stats$all_balanced,
    n_obs = nrow(d0)
  ),
  file.path(out_dir, "first_stage", "first_stage_fit_summary_All.csv")
)

saveRDS(d0, file.path(prep_dir, "d0_All_clean.rds"))
readr::write_csv(d0, file.path(prep_dir, "d0_All_clean.csv.gz"))

anchor <- "eligible_sharp"

ind_remote_tbl <- d0 %>%
  dplyr::group_by(prev_naics2) %>%
  dplyr::summarise(
    remote_intensity = mean(.data[[anchor]], na.rm = TRUE),
    n_in_industry    = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    pr = dplyr::percent_rank(remote_intensity),
    ind_remote_tercile = dplyr::case_when(
      pr <= 1/3 ~ "L",
      pr <= 2/3 ~ "M",
      TRUE      ~ "H"
    ),
    ind_remote_tercile = factor(ind_remote_tercile, levels = c("L","M","H"))
  ) %>%
  dplyr::select(-pr)

readr::write_csv(
  ind_remote_tbl,
  file.path(prep_dir, "industry_remote_intensity_by_prev_naics2.csv")
)

d0_ind <- d0 %>%
  dplyr::left_join(ind_remote_tbl, by = "prev_naics2") %>%
  droplevels()

saveRDS(d0_ind, file.path(prep_dir, "d0_All_with_industry_terciles.rds"))
readr::write_csv(d0_ind, file.path(prep_dir, "d0_All_with_industry_terciles.csv.gz"))

terciles <- c("L","M","H")

for (tlev in terciles) {
  dd <- d0_ind %>% filter(ind_remote_tercile == tlev) %>% droplevels()
  scope <- sprintf("All_indTercile_%s", tlev)

  cat(sprintf("\n[Weighting] Tercile = %s | n = %s\n", tlev, nrow(dd)))

  if (nrow(dd) < 2 || length(unique(dd$eligible_sharp)) < 2) {
    warning(sprintf("Skip %s: not enough observations or only one treatment level.", scope))
    next
  }

  w_out_t <- weightit(ps_formula, data = dd, method = "ps", estimand = "ATT")
  dd$w <- w_out_t$weights
  dd$ps <- w_out_t$ps

  stats_t <- save_balance_and_love(w_out_t, dd, scope, out_dir, pretty_var)

  ps_glm_t <- glm(ps_formula, data = dd, family = binomial())
  p_hat_t  <- as.numeric(predict(ps_glm_t, type = "response"))
  y_t      <- as.numeric(dd$eligible_sharp)
  w_t      <- dd$w
  auc_t <- tryCatch(as.numeric(pROC::auc(y_t, p_hat_t)), error = function(e) NA_real_)
  brier_t  <- mean((y_t - p_hat_t)^2, na.rm = TRUE)
  brierw_t <- sum(w_t * (y_t - p_hat_t)^2, na.rm = TRUE) / sum(w_t, na.rm = TRUE)

  readr::write_csv(
    tibble(
      scope = scope,
      auc = auc_t, brier = brier_t, brier_weighted = brierw_t,
      ess_overall = stats_t$ess_overall, ess_treated = stats_t$ess_treated, ess_control = stats_t$ess_control,
      max_abs_smd_un = stats_t$max_abs_smd_un, max_abs_smd_adj = stats_t$max_abs_smd_adj,
      all_balanced = stats_t$all_balanced,
      n_obs = nrow(dd)
    ),
    file.path(out_dir, "first_stage", sprintf("first_stage_fit_summary_%s.csv", scope))
  )

  saveRDS(dd, file.path(prep_dir, sprintf("d0_%s.rds", scope)))
  readr::write_csv(dd, file.path(prep_dir, sprintf("d0_%s.csv.gz", scope)))

  readr::write_csv(
    dd %>% select(row_id, ind_remote_tercile, w),
    file.path(prep_dir, sprintf("weights_%s.csv", scope))
  )

  cat(sprintf("Saved outputs for %s (weights, love plots, balance, fit summary)\n", scope))
}

cat(
  "\nAll tercile-specific runs completed. Outputs under:\n",
  out_dir, " and ", prep_dir, "\n", sep = ""
)
