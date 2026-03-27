#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(knitr)
})

# =========================================================
# Paths
#   - Set d = project root (e.g., "$d/modeling_outputs/...")
# =========================================================
d <- Sys.getenv("d")
stopifnot(nchar(d) > 0)

out_root <- file.path(d, "data", "modeling_outputs", "analysis_effects")
stopifnot(dir.exists(out_root))

latex_dir <- file.path(out_root, "_SI_tables")
dir.create(latex_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# Labels
# =========================================================
scope_order <- c("US", "US_indTercile_L", "US_indTercile_M", "US_indTercile_H")
scope_label <- c(
  US = "US",
  US_indTercile_L = "Low-remote",
  US_indTercile_M = "Mid-remote",
  US_indTercile_H = "High-remote"
)

outcome_label <- function(x) {
  dplyr::recode(
    x,
    "log_wage_growth"       = "Log Wage Growth",
    "promotion_at_entry"    = "Promotion at Entry",
    "curr_pagerank_z"       = "Employer Prestige Score",
    "curr_parent_headcount" = "Log Employer Size",
    "cross_metro_move"      = "Cross-Metro Job Change",
    .default = x
  )
}

outcome_slug <- function(x) {
  x2 <- outcome_label(x)
  x2 <- str_to_lower(x2)
  x2 <- str_replace_all(x2, "[^a-z0-9]+", "-")
  x2 <- str_replace_all(x2, "^-|-$", "")
  x2
}

fmt_fixed <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

fmt_sci <- function(x, digits = 1) {
  ifelse(is.na(x), "", formatC(x, format = "e", digits = digits))
}

p_stars <- function(p) {
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.1]   <- "."
  out[!is.na(p) & p < 0.05]  <- "*"
  out[!is.na(p) & p < 0.01]  <- "**"
  out[!is.na(p) & p < 0.001] <- "***"
  out
}

tiny_cutoff <- 0.0005

fmt_value_global <- function(x, digits_main = 3, digits_sci = 1) {
  ifelse(
    is.na(x), "",
    ifelse(abs(x) < tiny_cutoff,
           fmt_sci(x, digits_sci),
           fmt_fixed(x, digits_main))
  )
}

fmt_cell_global <- function(est, se, p, digits_main = 3, digits_sci = 1) {
  ok <- !(is.na(est) | is.na(se))
  out <- rep("", length(est))
  out[ok] <- paste0(
    fmt_value_global(est[ok], digits_main, digits_sci),
    p_stars(p[ok]),
    " (",
    fmt_value_global(se[ok], digits_main, digits_sci),
    ")"
  )
  out
}

wrap_tabular_resizebox <- function(tex) {
  tex <- stringr::str_replace(
    tex,
    stringr::fixed("\\begin{tabular}"),
    "\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}"
  )
  tex <- stringr::str_replace(
    tex,
    stringr::fixed("\\end{tabular}"),
    "\\end{tabular}\n}"
  )
  tex
}

# =========================================================
# Build LaTeX table for one outcome
# =========================================================
build_one_table <- function(outcome) {
  od <- file.path(out_root, outcome)
  if (!dir.exists(od)) return(NULL)

  f_coef <- file.path(od, sprintf("SI_main_multiscope_%s_coef_long.csv", outcome))
  f_stat <- file.path(od, sprintf("SI_main_multiscope_%s_model_stats.csv", outcome))
  f_spec <- file.path(od, sprintf("SI_main_multiscope_%s_spec_info.csv", outcome))

  if (!file.exists(f_coef) || !file.exists(f_stat) || !file.exists(f_spec)) {
    message("Skip (missing multiscope-support CSVs): ", outcome)
    return(NULL)
  }

  coef_long <- readr::read_csv(f_coef, show_col_types = FALSE)
  stats_df  <- readr::read_csv(f_stat, show_col_types = FALSE)

  if ("term_label" %in% names(coef_long)) {
    coef_long <- coef_long %>%
      mutate(
        term_label = ifelse(is.na(term_label), term, term_label),
        term_label = str_replace_all(term_label, fixed("Previous employer prestige (z)"), "Previous employer prestige"),
        term_label = str_replace_all(term_label, fixed("Employer prestige (z)"), "Previous employer prestige"),
        term_label = str_replace_all(term_label, fixed("Employer prestige"), "Previous employer prestige"),
        term_label = str_replace_all(term_label, fixed("Previous employer prestige"), "Previous employer prestige"),
        term_label = str_replace_all(term_label, fixed("log(Previous wage + 1)"), "Log previous wage"),
        term_label = str_replace_all(term_label, fixed("log(Previous employer headcount + 1)"), "Log previous employer size"),
        term_label = str_replace_all(term_label, fixed("Prob. female"), "Prob. female")
      )
  } else {
    coef_long <- coef_long %>% mutate(term_label = term)
  }

  coef_long <- coef_long %>%
    mutate(
      scope = factor(scope, levels = scope_order),
      cell  = fmt_cell_global(estimate, std_error, p_value, digits_main = 3, digits_sci = 0)
    ) %>%
    filter(!is.na(scope)) %>%
    arrange(scope)

  tab_main <- coef_long %>%
    select(term_label, scope, cell) %>%
    distinct() %>%
    mutate(scope_lbl = scope_label[as.character(scope)]) %>%
    select(term_label, scope_lbl, cell) %>%
    tidyr::pivot_wider(names_from = scope_lbl, values_from = cell)

  stats_w <- stats_df %>%
    mutate(scope = factor(scope, levels = scope_order)) %>%
    filter(!is.na(scope)) %>%
    mutate(scope_lbl = scope_label[as.character(scope)]) %>%
    select(scope_lbl, N, adj_r2) %>%
    tidyr::pivot_longer(cols = c(N, adj_r2), names_to = "stat", values_to = "value") %>%
    mutate(
      stat_label = dplyr::recode(stat, N = "N", adj_r2 = "Adj. $R^2$"),
      value_fmt  = ifelse(stat == "N", as.character(as.integer(value)), fmt_fixed(value, 3))
    ) %>%
    select(stat_label, scope_lbl, value_fmt) %>%
    tidyr::pivot_wider(names_from = scope_lbl, values_from = value_fmt)

  tab_all <- bind_rows(
    tab_main,
    tibble::tibble(term_label = "\\midrule", US = "", `Low-remote` = "", `Mid-remote` = "", `High-remote` = ""),
    stats_w %>% rename(term_label = stat_label)
  )

  want_cols <- c("term_label", "US", "Low-remote", "Mid-remote", "High-remote")
  for (cc in setdiff(want_cols, names(tab_all))) tab_all[[cc]] <- ""
  tab_all <- tab_all %>% select(all_of(want_cols))

  cap <- sprintf(
    "Main regression results by industry remote-eligibility terciles (outcome: %s).",
    outcome_label(outcome)
  )
  lab <- sprintf("tab:si-main-%s", outcome_slug(outcome))

  notes <- c(
    "The unit of observation is a job-transition record.",
    "All specifications are weighted using inverse probability weights.",
    "All regressions include fixed effects for previous role, seniority, ethnicity, industry, previous metropolitan area, and entry year.",
    "Standard errors are clustered by previous metropolitan area and entry year.",
    "To accommodate zero values, wage and employer size are log-transformed using log(x + 1).",
    "Entries report coefficient estimates with standard errors in parentheses.",
    "Significance levels are $^{***}p<0.001$, $^{**}p<0.01$, $^{*}p<0.05$, and $^{.}p<0.1$."
  )

  tex <- knitr::kable(
    tab_all %>% mutate(term_label = ifelse(term_label == "\\midrule", "", term_label)),
    format   = "latex",
    booktabs = TRUE,
    escape   = FALSE,
    col.names = c("", "All", "Low-remote", "Mid-remote", "High-remote"),
    caption  = cap,
    align    = c("l", "c", "c", "c", "c")
  ) %>% as.character()

  tex <- stringr::str_replace(tex, "\nN &", "\n\\\\midrule\nN &")
  tex <- wrap_tabular_resizebox(tex)
  tex <- stringr::str_replace(tex, stringr::fixed("\\begin{table}"), "\\begin{table}[!htbp]")
  tex <- stringr::str_replace(
    tex,
    stringr::fixed(sprintf("\\caption{%s}", cap)),
    sprintf("\\caption{%s}\\label{%s}", cap, lab)
  )

  note_tex <- paste0(
    "\n\\begin{flushleft}\n",
    "\\footnotesize\n",
    "\\textit{Notes:} ", paste(notes, collapse = " "),
    "\n\\end{flushleft}\n"
  )

  tex <- stringr::str_replace(
    tex,
    stringr::fixed("\\end{table}"),
    paste0(note_tex, "\\end{table}")
  )

  stopifnot(stringr::str_detect(tex, stringr::fixed("\\resizebox{\\textwidth}{!}{%")))
  stopifnot(stringr::str_detect(tex, stringr::fixed(sprintf("\\label{%s}", lab))))
  stopifnot(stringr::str_detect(tex, stringr::fixed("\\begin{flushleft}")))
  stopifnot(stringr::str_detect(tex, stringr::fixed("\\end{table}")))

  out_tex <- file.path(latex_dir, sprintf("SI_multiscope_%s.tex", outcome))
  writeLines(tex, out_tex)
  message("Wrote: ", out_tex)
  message("  Label: ", lab)

  out_tex
}

# =========================================================
# Run
# =========================================================
keep_outcomes <- c("log_wage_growth","promotion_at_entry","curr_pagerank_z","curr_parent_headcount","cross_metro_move")
outcome_dirs <- list.dirs(out_root, recursive = FALSE, full.names = FALSE)
outcomes <- intersect(outcome_dirs[outcome_dirs != "_SI_tables"], keep_outcomes)

written <- purrr::map(outcomes, build_one_table)
written <- written[!vapply(written, is.null, logical(1))]

index_file <- file.path(latex_dir, "SI_tables_index.txt")
writeLines(c("Generated LaTeX tables:", unlist(written)), index_file)

message("Done. Index: ", index_file)
message("LaTeX output folder: ", latex_dir)
