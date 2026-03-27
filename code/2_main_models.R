#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(fixest)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(arrow)
  library(stringr)
  library(tibble)
  library(purrr)
  library(forcats)
})

set.seed(42)

# =========================================================
# Paths
# =========================================================
d <- Sys.getenv("d")
stopifnot(nchar(d) > 0)

data_parquet <- file.path(d, "data", "job_entry_data.parquet")
stopifnot(file.exists(data_parquet))

w_root   <- file.path(d, "data", "modeling_outputs", "weighting_result")
prep_dir <- file.path(w_root, "analysis_data")
stopifnot(dir.exists(prep_dir))

out_root <- file.path(d, "data", "modeling_outputs", "analysis_effects")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# Progress helpers
# =========================================================
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logi <- function(...) { cat(sprintf("[%s] ", ts()), sprintf(...), "\n", sep=""); flush.console() }
tic <- function() proc.time()[["elapsed"]]
toc <- function(t0) proc.time()[["elapsed"]] - t0

# =========================================================
# Load + rebuild the analysis base
# =========================================================
logi("[Step 1] Reading base parquet: %s", data_parquet)
stopifnot(file.exists(data_parquet))
t0 <- tic()
df0 <- arrow::read_parquet(data_parquet)
logi("[Step 1] Done reading parquet. n=%s cols=%s (%.1fs)", nrow(df0), ncol(df0), toc(t0))

need_vars <- c(
  "row_id",
  "eligible_sharp","f_prob","ethnicity_predicted","prev_role_k50","prev_seniority",
  "prev_wage","prev_pagerank_z","prev_parent_headcount",
  "curr_parent_headcount","curr_pagerank_z",
  "prev_naics2","prev_metro_area","entry_year","country","job_category",
  "log_wage_growth","promotion_at_entry","cross_metro_move"
)
stopifnot(all(need_vars %in% names(df0)))

t0 <- tic()
df_use <- df0 %>%
  select(all_of(need_vars)) %>%
  mutate(
    prev_seniority_num = suppressWarnings(as.integer(prev_seniority)),

    managerial_role = ifelse(!is.na(prev_seniority_num) & prev_seniority_num >= 4 & prev_seniority_num <= 7,
                             "Managerial","NonManagerial"),
    managerial_role = factor(managerial_role, levels = c("NonManagerial","Managerial")),

    engineer_role   = ifelse(job_category == "Engineer","Engineer","NonEngineer"),
    engineer_role   = factor(engineer_role, levels = c("NonEngineer","Engineer")),

    log_prev_wag = log1p(prev_wage),
    log_prev_parent_headcount = log1p(prev_parent_headcount),

    prev_seniority    = factor(prev_seniority,
                               levels = 1:7,
                               labels = c("1 Entry","2 Junior","3 Associate",
                                          "4 Manager","5 Director","6 Executive","7 SrExecutive"),
                               ordered = TRUE),
    prev_naics2       = factor(prev_naics2),
    ethnicity_predicted = factor(ethnicity_predicted),
    prev_role_k50     = factor(prev_role_k50),
    prev_metro_area   = factor(prev_metro_area),
    entry_year        = factor(entry_year)
  ) %>%
  droplevels()
logi("[Step 1b] Done df_use (%.1fs)", toc(t0))

# =========================================================
# Read saved weights
# =========================================================
logi("[Step 2] Reading saved weight file ...")

f_us <- file.path(prep_dir, "d0_US_clean.csv.gz")
stopifnot(file.exists(f_us))

t0 <- tic()
d_us <- readr::read_csv(f_us, show_col_types = FALSE)
logi("[Step 2] Loaded US weights file (%.1fs)", toc(t0))

keep_cols <- function(dd) {
  ks <- intersect(c("row_id","w"), names(dd))
  dd %>% select(all_of(ks))
}

w_us <- keep_cols(d_us) %>% mutate(scope = "US")

weights_all <- w_us
stopifnot(all(c("row_id","w","scope") %in% names(weights_all)))

make_scope_data <- function(scope_key) {
  wa <- weights_all %>% filter(scope == scope_key) %>% select(row_id, w)
  df_use %>% inner_join(wa, by = "row_id") %>% droplevels()
}

scopes <- c("US")

logi("[Step 2b] Joining weights into df_use  ...")
t0 <- tic()
data_scoped <- setNames(lapply(scopes, make_scope_data), scopes)
logi("[Step 2b] Done joining scopes (%.1fs).", toc(t0))
for (sc in scopes) logi("  - n(%s) = %s", sc, nrow(data_scoped[[sc]]))

# =========================================================
# Helpers / constants
# =========================================================
CTRL_full <- c("f_prob","log_prev_wag","prev_pagerank_z","log_prev_parent_headcount")
CLUSTER_full <- "~ prev_metro_area + entry_year"

p_stars <- function(p) {
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.1]   <- "."
  out[!is.na(p) & p < 0.05]  <- "*"
  out[!is.na(p) & p < 0.01]  <- "**"
  out[!is.na(p) & p < 0.001] <- "***"
  out
}

lhs_vector <- function(dd, outcome, log_headcount = FALSE) {
  if (log_headcount) return(log(dd$curr_parent_headcount + 1))
  dd[[outcome]]
}

has_variation <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) < 2) return(FALSE)
  stats::var(x) > 0
}

build_formula <- function(outcome, rhs_ctrl, drop_FE = NULL,
                          log_headcount = FALSE) {
  fe_terms <- c("prev_seniority","ethnicity_predicted","prev_role_k50",
                "prev_naics2","prev_metro_area","entry_year")
  if (!is.null(drop_FE)) fe_terms <- setdiff(fe_terms, drop_FE)
  fe_part <- paste(fe_terms, collapse = " + ")

  lhs <- outcome
  if (log_headcount) lhs <- "log(curr_parent_headcount + 1)"

  as.formula(paste0(lhs, " ~ eligible_sharp + ",
                    paste(rhs_ctrl, collapse = " + "),
                    " | ", fe_part))
}

run_model_once <- function(dd, outcome, rhs_ctrl, drop_FE = NULL,
                           is_binary = FALSE, log_headcount = FALSE,
                           return_fit = FALSE) {

  out_var <- if (log_headcount) "curr_parent_headcount" else outcome

  vars_needed <- unique(c(
    "eligible_sharp", "w",
    rhs_ctrl,
    "prev_seniority","ethnicity_predicted","prev_role_k50",
    "prev_naics2","prev_metro_area","entry_year",
    out_var
  ))

  dd_clean <- dd %>%
    select(all_of(vars_needed)) %>%
    tidyr::drop_na()

  if (nrow(dd_clean) < 2 || length(unique(dd_clean$eligible_sharp)) < 2) {
    out_tbl <- tibble(
      coef_eligible = NA_real_,
      se_eligible   = NA_real_,
      ci_low_95     = NA_real_,
      ci_high_95    = NA_real_,
      p_value       = NA_real_,
      p_stars       = "",
      adj_r2        = NA_real_,
      n             = nrow(dd_clean)
    )
    return(if (return_fit) list(summary = out_tbl, fit = NULL) else out_tbl)
  }

  y <- lhs_vector(dd_clean, outcome, log_headcount = log_headcount)
  if (!has_variation(y)) {
    out_tbl <- tibble(
      coef_eligible = NA_real_,
      se_eligible   = NA_real_,
      ci_low_95     = NA_real_,
      ci_high_95    = NA_real_,
      p_value       = NA_real_,
      p_stars       = "",
      adj_r2        = NA_real_,
      n             = length(y)
    )
    return(if (return_fit) list(summary = out_tbl, fit = NULL) else out_tbl)
  }

  fml <- build_formula(outcome, rhs_ctrl, drop_FE = drop_FE,
                     log_headcount = log_headcount)

  fit <- feols(fml, data = dd_clean, weights = ~ w, vcov = as.formula(CLUSTER_full))

  ct   <- coeftable(fit)
  est  <- as.numeric(unname(coef(fit)["eligible_sharp"]))
  se   <- as.numeric(unname(ct["eligible_sharp","Std. Error"]))
  pval <- as.numeric(unname(ct["eligible_sharp","Pr(>|t|)"]))

  ci95   <- confint(fit, "eligible_sharp", level = 0.95)
  ci_low <- as.numeric(ci95[1, 1, drop = TRUE])
  ci_hi  <- as.numeric(ci95[1, 2, drop = TRUE])

  ar2 <- suppressWarnings(as.numeric(fitstat(fit, "ar2")))
  n   <- as.numeric(nobs(fit))

  out_tbl <- tibble(
    coef_eligible = est,
    se_eligible   = se,
    ci_low_95     = ci_low,
    ci_high_95    = ci_hi,
    p_value       = pval,
    p_stars       = p_stars(pval),
    adj_r2        = ar2,
    n             = n
  )

  if (return_fit) list(summary = out_tbl, fit = fit) else out_tbl
}

# =========================================================
# Export model results and specification metadata as CSV files
# =========================================================
dict_overall <- tibble::tribble(
  ~term, ~term_label,
  "eligible_sharp", "Remote-eligible",
  "f_prob", "Prob. female",
  "log_prev_wag", "log(Previous wage + 1)",
  "prev_pagerank_z", "Previous employer prestige (z)",
  "log_prev_parent_headcount", "log(Previous employer headcount + 1)"
)

write_multiscope_support_csvs <- function(outcome, fits_named, out_dir_outcome) {
  fits_named <- fits_named[vapply(fits_named, function(x) inherits(x, "fixest"), logical(1))]
  if (length(fits_named) == 0) {
    logi("[Step 4] No valid fixest fits for outcome=%s; skip multiscope-support CSVs.", outcome)
    return(invisible(NULL))
  }

  keep_terms <- dict_overall$term

  coef_long <- bind_rows(lapply(names(fits_named), function(sc) {
    fit <- fits_named[[sc]]
    ct  <- as.data.frame(coeftable(fit))
    ct$term <- rownames(ct)
    rownames(ct) <- NULL

    ct %>%
      transmute(
        outcome = outcome,
        scope   = sc,
        term,
        estimate  = .data[["Estimate"]],
        std_error = .data[["Std. Error"]],
        p_value   = .data[["Pr(>|t|)"]]
      ) %>%
      filter(term %in% keep_terms) %>%
      left_join(dict_overall, by = "term")
  }))

  out_coef_csv <- file.path(out_dir_outcome, sprintf("SI_main_multiscope_%s_coef_long.csv", outcome))
  write_csv(coef_long, out_coef_csv)

  stats_df <- bind_rows(lapply(names(fits_named), function(sc) {
    fit <- fits_named[[sc]]
    tibble(
      outcome = outcome,
      scope   = sc,
      N       = as.numeric(nobs(fit)),
      adj_r2  = suppressWarnings(as.numeric(fitstat(fit, "ar2")))
    )
  }))
  out_stats_csv <- file.path(out_dir_outcome, sprintf("SI_main_multiscope_%s_model_stats.csv", outcome))
  write_csv(stats_df, out_stats_csv)

  spec_info <- tibble(
    outcome = outcome,
    controls = paste(CTRL_full, collapse = " + "),
    fixed_effects = "prev_role_k50 + prev_seniority + ethnicity_predicted + prev_naics2 + prev_metro_area + entry_year",
    weights = "IPW (ATT), pre-computed",
    cluster = CLUSTER_full
  )
  out_spec_csv <- file.path(out_dir_outcome, sprintf("SI_main_multiscope_%s_spec_info.csv", outcome))
  write_csv(spec_info, out_spec_csv)

  logi("[Step 4] Saved multiscope-support CSVs (NO tex):")
  logi("  - %s", out_coef_csv)
  logi("  - %s", out_stats_csv)
  logi("  - %s", out_spec_csv)

  invisible(list(coef=out_coef_csv, stats=out_stats_csv, spec=out_spec_csv))
}

# =========================================================
# Main estimation routine:
# (1) estimate the overall model
# (2) estimate heterogeneity analyses
# =========================================================
run_effects_for_outcome <- function(scope_key, dd_scope, outcome,
                                    run_hetero = TRUE,
                                    is_binary = FALSE,
                                    log_headcount = FALSE) {

  out_dir_outcome <- file.path(out_root, outcome)
  dir.create(out_dir_outcome, recursive = TRUE, showWarnings = FALSE)
  logi("Estimating models for scope=%s, outcome=%s, heterogeneity=%s", scope_key, outcome, run_hetero)

  # ---------- Overall ----------
  pack_overall <- run_model_once(dd_scope, outcome, CTRL_full, NULL,
                                 is_binary = is_binary,
                                 log_headcount = log_headcount,
                                 return_fit = TRUE)

  res_overall <- pack_overall$summary %>%
    mutate(scope = scope_key, spec = "overall")

  write_csv(res_overall, file.path(out_dir_outcome, sprintf("overall_%s.csv", scope_key)))
  fit_overall <- pack_overall$fit

  # ---------- Heterogeneity analyses ----------
  if (!isTRUE(run_hetero)) return(fit_overall)

  # 2) Gender
  dd_g <- dd_scope %>%
    mutate(female_ind = ifelse(is.na(f_prob), NA, as.integer(f_prob >= 0.5))) %>%
    filter(!is.na(female_ind))
  rhs_gender <- setdiff(CTRL_full, "f_prob")

  res_gender <- list()
  for (g in c(0,1)) {
    dsub <- dd_g %>% filter(female_ind == g)
    if (nrow(dsub) >= 50 && length(unique(dsub$eligible_sharp)) >= 2) {
      rg <- run_model_once(dsub, outcome, rhs_gender, NULL,
                           is_binary = is_binary,
                           log_headcount = log_headcount) %>%
        mutate(scope=scope_key, spec=paste0("gender_", ifelse(g==1,"Female","Male")))
      res_gender[[length(res_gender)+1]] <- rg
    }
  }
  if (length(res_gender)>0) {
    write_csv(bind_rows(res_gender),
              file.path(out_dir_outcome, sprintf("hetero_gender_%s.csv", scope_key)))
  }

  # 3) Seniority
  res_sen <- list()
  for (lev in levels(dd_scope$prev_seniority)) {
    dsub <- dd_scope %>% filter(prev_seniority == lev)
    if (nrow(dsub)>=50 && length(unique(dsub$eligible_sharp))>=2) {
      rs <- run_model_once(dsub, outcome, CTRL_full, drop_FE="prev_seniority",
                           is_binary = is_binary,
                           log_headcount = log_headcount) %>%
        mutate(scope=scope_key, spec=paste0("seniority_", lev))
      res_sen[[length(res_sen)+1]] <- rs
    }
  }
  if (length(res_sen)>0) {
    write_csv(bind_rows(res_sen),
              file.path(out_dir_outcome, sprintf("hetero_seniority_%s.csv", scope_key)))
  }

  # 4) Wage terciles
  dd_w <- dd_scope %>%
    filter(!is.na(prev_wage)) %>%
    mutate(wage_tercile = ntile(prev_wage, 3),
           wage_tercile = factor(wage_tercile, levels=1:3, labels=c("L","M","H")))
  rhs_wage <- setdiff(CTRL_full, "log_prev_wag")

  res_wage <- list()
  for (wt in c("L","M","H")) {
    dsub <- dd_w %>% filter(wage_tercile==wt)
    if (nrow(dsub)>=50 && length(unique(dsub$eligible_sharp))>=2) {
      rw <- run_model_once(dsub, outcome, rhs_wage, NULL,
                           is_binary = is_binary,
                           log_headcount = log_headcount) %>%
        mutate(scope=scope_key, spec=paste0("wageTercile_",wt))
      res_wage[[length(res_wage)+1]] <- rw
    }
  }
  if (length(res_wage)>0) {
    write_csv(bind_rows(res_wage),
              file.path(out_dir_outcome, sprintf("hetero_wageTerciles_%s.csv", scope_key)))
  }

  # 5) Managerial
  dd_mgr <- dd_scope %>%
    filter(!is.na(prev_seniority_num)) %>%
    mutate(managerial_role = ifelse(prev_seniority_num >= 4 & prev_seniority_num <= 7,
                                    "Managerial","NonManagerial"),
           managerial_role = factor(managerial_role,
                                    levels = c("NonManagerial","Managerial")))
  res_mgr <- list()
  for (grp in levels(dd_mgr$managerial_role)) {
    dsub <- dd_mgr %>% filter(managerial_role == grp)
    if (nrow(dsub) >= 50 && length(unique(dsub$eligible_sharp)) >= 2) {
      rmgr <- run_model_once(dsub, outcome, CTRL_full, NULL,
                             is_binary = is_binary,
                             log_headcount = log_headcount) %>%
        mutate(scope=scope_key, spec=paste0("managerialRole_", grp))
      res_mgr[[length(res_mgr)+1]] <- rmgr
    }
  }
  if (length(res_mgr)>0) {
    write_csv(bind_rows(res_mgr),
              file.path(out_dir_outcome, sprintf("hetero_managerialRole_%s.csv", scope_key)))
  }

  # 6) Engineer
  dd_eng <- dd_scope %>% filter(!is.na(engineer_role))
  res_eng <- list()
  for (grp in levels(dd_eng$engineer_role)) {
    dsub <- dd_eng %>% filter(engineer_role == grp)
    if (nrow(dsub) >= 50 && length(unique(dsub$eligible_sharp)) >= 2) {
      reng <- run_model_once(dsub, outcome, CTRL_full, NULL,
                             is_binary = is_binary,
                             log_headcount = log_headcount) %>%
        mutate(scope=scope_key, spec=paste0("engineerRole_", grp))
      res_eng[[length(res_eng)+1]] <- reng
    }
  }
  if (length(res_eng)>0) {
    write_csv(bind_rows(res_eng),
              file.path(out_dir_outcome, sprintf("hetero_engineerRole_%s.csv", scope_key)))
  }

  fit_overall
}

# =========================================================
# Outcomes
# =========================================================
outcomes <- tribble(
  ~outcome,               ~is_binary, ~log_headcount,
  "log_wage_growth",      FALSE,      FALSE,
  "promotion_at_entry",   TRUE,       FALSE,
  "cross_metro_move",     TRUE,       FALSE,
  "curr_parent_headcount",FALSE,      TRUE,
  "curr_pagerank_z",      FALSE,      FALSE
)

# =========================================================
# Run (US only)
# =========================================================
logi("[Step 3] Running models (US only) ...")

for (i in 1:nrow(outcomes)) {
  oc <- outcomes$outcome[i]
  logi("===============================================")
  logi("[Outcome] %s", oc)
  logi("===============================================")

  fits_named <- list()

  sc <- "US"
  dd_sc <- data_scoped[[sc]]
  if (nrow(dd_sc) < 50 || length(unique(dd_sc$eligible_sharp)) < 2) {
    logi("  [Skip] %s: insufficient data.", sc)
    next
  }

  fit_overall <- run_effects_for_outcome(
    scope_key = sc,
    dd_scope  = dd_sc,
    outcome   = oc,
    run_hetero = TRUE,
    is_binary = outcomes$is_binary[i],
    log_headcount = outcomes$log_headcount[i]
  )

  fits_named[[sc]] <- fit_overall

  logi("[Step 4] Writing multiscope-support CSVs ...")
  out_dir_outcome <- file.path(out_root, oc)
  dir.create(out_dir_outcome, recursive = TRUE, showWarnings = FALSE)
  write_multiscope_support_csvs(oc, fits_named, out_dir_outcome)
}

logi("[Done] Results saved under: %s", out_root)
