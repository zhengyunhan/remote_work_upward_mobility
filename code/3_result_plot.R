#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(scales)
})

# =========================================================
# Paths
# =========================================================
d <- Sys.getenv("d")
stopifnot(nchar(d) > 0)

effect_root <- file.path(d, "data", "modeling_outputs", "analysis_effects")
stopifnot(dir.exists(effect_root))

figure_dir <- file.path(effect_root, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

scope <- "US"

# =========================================================
# Outcomes and labels
# =========================================================
outcomes <- c(
  "log_wage_growth",
  "promotion_at_entry",
  "cross_metro_move",
  "curr_parent_headcount",
  "curr_pagerank_z"
)

nice_outcome <- c(
  "log_wage_growth"       = "Log wage growth",
  "promotion_at_entry"    = "Upward seniority move",
  "cross_metro_move"      = "Cross-metro job change",
  "curr_parent_headcount" = "Log employer size",
  "curr_pagerank_z"       = "Employer prestige score"
)

outcome_order <- c(
  "Log wage growth",
  "Upward seniority move",
  "Employer prestige score",
  "Log employer size",
  "Cross-metro job change"
)

group_order <- c(
  "Overall",
  "Sex",
  "Managerial role",
  "Engineer role",
  "Previous wage tercile"
)

level_order <- c(
  "Overall",
  "Female", "Male",
  "Non-managerial", "Managerial",
  "Non-engineer", "Engineer",
  "Wage: Low", "Wage: Middle", "Wage: High"
)

# =========================================================
# Readers
# =========================================================
read_overall_outcome <- function(outcome_name, scope_name) {
  f <- file.path(effect_root, outcome_name, sprintf("overall_%s.csv", scope_name))
  if (!file.exists(f)) return(NULL)

  df0 <- read_csv(f, show_col_types = FALSE)
  required_cols <- c("coef_eligible", "ci_low_95", "ci_high_95", "p_value")
  stopifnot(all(required_cols %in% names(df0)))

  df0 %>%
    filter(!is.na(coef_eligible), !is.na(ci_low_95), !is.na(ci_high_95)) %>%
    slice(1) %>%
    mutate(
      outcome = outcome_name,
      group = "Overall",
      level = "Overall",
      spec = "overall"
    )
}

read_heterogeneity_outcome <- function(outcome_name, scope_name) {
  outcome_dir <- file.path(effect_root, outcome_name)

  files <- c(
    gender     = file.path(outcome_dir, sprintf("hetero_gender_%s.csv", scope_name)),
    wage       = file.path(outcome_dir, sprintf("hetero_wageTerciles_%s.csv", scope_name)),
    managerial = file.path(outcome_dir, sprintf("hetero_managerialRole_%s.csv", scope_name)),
    engineer   = file.path(outcome_dir, sprintf("hetero_engineerRole_%s.csv", scope_name))
  )

  existing_files <- files[file.exists(files)]
  if (length(existing_files) == 0) return(NULL)

  bind_rows(lapply(names(existing_files), function(tag) {
    read_csv(existing_files[[tag]], show_col_types = FALSE) %>%
      mutate(.hetero = tag, outcome = outcome_name)
  }))
}

# =========================================================
# Load data
# =========================================================
overall_df <- bind_rows(lapply(outcomes, read_overall_outcome, scope_name = scope))
heterogeneity_df <- bind_rows(lapply(outcomes, read_heterogeneity_outcome, scope_name = scope))
df <- bind_rows(overall_df, heterogeneity_df)

if (is.null(df) || nrow(df) == 0) {
  stop("No valid model result files were found under: ", effect_root)
}

required_cols <- c("outcome", "coef_eligible", "ci_low_95", "ci_high_95", "p_value")
stopifnot(all(required_cols %in% names(df)))

if (!"group" %in% names(df)) df$group <- NA_character_
if (!"level" %in% names(df)) df$level <- NA_character_
if (!"spec" %in% names(df)) df$spec <- NA_character_

# =========================================================
# Standardize labels
# =========================================================
df <- df %>%
  filter(!is.na(coef_eligible), !is.na(ci_low_95), !is.na(ci_high_95)) %>%
  mutate(
    outcome_lab = recode(outcome, !!!nice_outcome),
    outcome_lab = factor(outcome_lab, levels = outcome_order),

    group = case_when(
      !is.na(group) & nzchar(group)        ~ group,
      str_starts(spec, "gender_")          ~ "Sex",
      str_starts(spec, "wageTercile_")     ~ "Previous wage tercile",
      str_starts(spec, "managerialRole_")  ~ "Managerial role",
      str_starts(spec, "engineerRole_")    ~ "Engineer role",
      TRUE                                 ~ "Other"
    ),

    level = case_when(
      !is.na(level) & nzchar(level)        ~ level,
      str_starts(spec, "gender_")          ~ str_remove(spec, "^gender_"),
      str_starts(spec, "wageTercile_")     ~ str_remove(spec, "^wageTercile_"),
      str_starts(spec, "managerialRole_")  ~ str_remove(spec, "^managerialRole_"),
      str_starts(spec, "engineerRole_")    ~ str_remove(spec, "^engineerRole_"),
      TRUE                                 ~ spec
    ),

    level = case_when(
      group == "Sex" & level %in% c("Female", "female", "F") ~ "Female",
      group == "Sex" & level %in% c("Male", "male", "M")     ~ "Male",

      group == "Previous wage tercile" & level == "L"        ~ "Low",
      group == "Previous wage tercile" & level == "M"        ~ "Middle",
      group == "Previous wage tercile" & level == "H"        ~ "High",

      group == "Managerial role" &
        level %in% c("Managerial", "1", "Y", "Yes", "yes")   ~ "Managerial",
      group == "Managerial role" &
        level %in% c("NonManagerial", "0", "N", "No", "no")  ~ "Non-managerial",

      group == "Engineer role" &
        level %in% c("Engineer", "1", "Y", "Yes", "yes")     ~ "Engineer",
      group == "Engineer role" &
        level %in% c("NonEngineer", "0", "N", "No", "no")    ~ "Non-engineer",

      TRUE ~ level
    ),

    label = case_when(
      group == "Overall"               ~ "Overall",
      group == "Sex"                   ~ level,
      group == "Managerial role"       ~ level,
      group == "Engineer role"         ~ level,
      group == "Previous wage tercile" ~ paste0("Wage: ", level),
      TRUE                             ~ level
    )
  ) %>%
  filter(group %in% group_order) %>%
  mutate(
    group = factor(group, levels = group_order),
    group_facet = fct_recode(
      group,
      "Previous wage\ntercile" = "Previous wage tercile"
    ),
    level_fac = factor(label, levels = level_order)
  )

# =========================================================
# Plot variables
# =========================================================
df <- df %>%
  mutate(
    coef_plot = coef_eligible,
    lo_plot   = ci_low_95,
    hi_plot   = ci_high_95
  )

vline_at <- 0
x_label <- "Coefficient estimate (95% confidence interval)"

# =========================================================
# Significance labels
# =========================================================
df <- df %>%
  mutate(
    sig_star = case_when(
      !is.na(p_value) & p_value < 0.001 ~ "***",
      !is.na(p_value) & p_value < 0.01  ~ "**",
      !is.na(p_value) & p_value < 0.05  ~ "*",
      !is.na(p_value) & p_value < 0.1   ~ ".",
      TRUE                              ~ ""
    ),
    sign = case_when(
      coef_plot > 0 ~ "Positive",
      coef_plot < 0 ~ "Negative",
      TRUE          ~ "Zero"
    ),
    value_label = paste0(
      case_when(
        outcome %in% c("log_wage_growth", "promotion_at_entry", "cross_metro_move") ~
          formatC(coef_plot, format = "f", digits = 3),
        TRUE ~
          formatC(coef_plot, format = "f", digits = 2)
      ),
      sig_star
    )
  )

# =========================================================
# Axis limits
# =========================================================
df_limits_base <- df %>%
  group_by(outcome, outcome_lab) %>%
  summarise(
    min_val = min(pmin(lo_plot, hi_plot, vline_at), na.rm = TRUE),
    max_val = max(pmax(lo_plot, hi_plot, vline_at), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    range_val = max_val - min_val,
    range_val = if_else(range_val == 0, 1, range_val),
    x_min     = min_val - 0.35 * range_val,
    x_max     = max_val + 0.35 * range_val
  ) %>%
  select(outcome_lab, x_min, x_max)

df_limits <- df %>%
  distinct(group_facet, outcome_lab, level_fac) %>%
  left_join(df_limits_base, by = "outcome_lab")

# =========================================================
# Plot
# =========================================================
p <- ggplot(df, aes(x = coef_plot, y = level_fac)) +
  geom_blank(data = df_limits, aes(x = x_min, y = level_fac)) +
  geom_blank(data = df_limits, aes(x = x_max, y = level_fac)) +
  geom_vline(xintercept = vline_at, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(
    aes(xmin = lo_plot, xmax = hi_plot, color = sign),
    height = 0.15,
    linewidth = 0.6
  ) +
  geom_point(
    aes(color = sign),
    size = 2.4
  ) +
  geom_text(
    aes(label = value_label),
    nudge_y = 0.28,
    size = 4.2,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "Positive" = "#4C72B0",
      "Negative" = "#C44E52",
      "Zero"     = "grey50"
    ),
    guide = "none"
  ) +
  labs(x = x_label, y = NULL) +
  scale_x_continuous(
    expand = c(0, 0),
    labels = label_number(accuracy = 0.01, trim = TRUE),
    guide = guide_axis(check.overlap = TRUE)
  ) +
  facet_grid(
    rows = vars(group_facet),
    cols = vars(outcome_lab),
    scales = "free",
    space = "free_y",
    switch = "y",
    labeller = labeller(outcome_lab = label_wrap_gen(width = 16))
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(colour = "black"),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 15),
    strip.text.x = element_text(face = "bold", size = 14),
    strip.text.y.left = element_text(
      face = "bold",
      size = 14,
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      margin = margin(r = 6)
    ),
    strip.background = element_rect(fill = "grey90", colour = NA),
    axis.ticks.y = element_blank(),
    strip.placement = "outside",
    panel.spacing.x = grid::unit(2.0, "lines"),
    plot.margin = margin(5.5, 25, 15, 5.5)
  )

print(p)

# =========================================================
# Save figure
# =========================================================
out_file_png <- file.path(figure_dir, sprintf("heterogeneity_grid_%s.png", scope))
out_file_pdf <- file.path(figure_dir, sprintf("heterogeneity_grid_%s.pdf", scope))

ggsave(
  filename = out_file_png,
  plot = p,
  width = 15,
  height = 9.5,
  units = "in",
  dpi = 300
)

ggsave(
  filename = out_file_pdf,
  plot = p,
  width = 15,
  height = 9.5,
  units = "in"
)

message("Saved figure to: ", out_file_png)
message("Saved figure to: ", out_file_pdf)