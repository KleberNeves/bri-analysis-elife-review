### ### ### ### ### ### ### ### ### ### ### ###

# tables.r
# This script generates formatted tables for reporting.

### ### ### ### ### ### ### ### ### ### ### ###

options(warn = -1)

p_load(rcompanion)

# Tables ------------------------------------------------------------------

## Replication Rates -------------------------------------------------------

df_replication_rate <- read_tsv(paste0("output/", results_path, "/Replication Rate Aggregate.tsv"), show_col_types = FALSE) |> 
  mutate(MetricShortName = case_when(
    MetricShortName == "Same-sign significance (p<0.05)" ~ "FEMA is significant and has same signal as original",
    MetricShortName == "Original in replication 95% PI" ~ "Original estimate within PI of REMA",
    MetricShortName == "Replication in original 95% CI" & Level == "Set of replications" ~ "REMA estimate within CI of Original",
    MetricShortName == "≥50% subjectively replicated" ~ "Subjective assessment majority vote (with ties as success)",
    MetricShortName == "≥50% significant" ~ "t-test majority vote (with ties as success)",
    MetricShortName == "Subjectively replicated (Indiv)" ~ "Subjective assessment, individual",
    MetricShortName == "Replication in original CI (Indiv)" & Level == "Individual replication" ~ "Replication estimate within CI of Original",
    MetricShortName == "Same-sign significance (Indiv)" ~ "Replication is significant and has same sign as original",
    .default = MetricShortName
  )) |>
  filter(
    ## By experiment
    MetricShortName == "FEMA is significant and has same signal as original" |
      MetricShortName == "Original estimate within PI of REMA" |
      MetricShortName == "REMA estimate within CI of Original" |
      MetricShortName == "Subjective assessment majority vote (with ties as success)" |
      MetricShortName == "t-test majority vote (with ties as success)" |
      ## By replication
      MetricShortName == "Subjective assessment, individual" |
      MetricShortName == "Replication estimate within CI of Original" |
      MetricShortName == "Replication is significant and has same sign as original"
  )


### By experiment -----------------------------------------------------------
desired_order_exp <- c(
  "Original estimate within PI of REMA",
  "REMA estimate within CI of Original",
  "FEMA is significant and has same signal as original",
  "t-test majority vote (with ties as success)",
  "Subjective assessment majority vote (with ties as success)"
)

df_by_experiment <- df_replication_rate |>
  filter(
    ## By experiment
    MetricShortName == "FEMA is significant and has same signal as original" |
      MetricShortName == "Original estimate within PI of REMA" |
      MetricShortName == "REMA estimate within CI of Original" |
      MetricShortName == "Subjective assessment majority vote (with ties as success)" |
      MetricShortName == "t-test majority vote (with ties as success)"
  ) |>
  mutate(Value = paste0(round(Value * 100, 2), "%")) |>
  mutate(MetricShortName = factor(MetricShortName, levels = desired_order_exp)) |>
  arrange(MetricShortName)



### By replication ----------------------------------------------------------
desired_order_rep <- c(
  "Replication estimate within CI of Original",
  "Replication is significant and has same sign as original",
  "Subjective assessment, individual"
)

df_by_replication <- df_replication_rate |>
  filter(
    ## By replication
    MetricShortName == "Subjective assessment, individual" |
      MetricShortName == "Replication estimate within CI of Original" |
      MetricShortName == "Replication is significant and has same sign as original"
  ) |>
  mutate(Value = paste0(round(Value * 100, 2), "%")) |>
  mutate(MetricShortName = factor(MetricShortName, levels = desired_order_rep)) |>
  arrange(MetricShortName)

### Functions ----

create_tbl_by_method_pcr <- function(analysis_type, distribution) {
  # Ensure that distribution is in lowercase
  distribution <- tolower(distribution)

  # Build wide tables and avoid relying on column order
  exp_wide <- df_by_experiment |>
    filter(Inclusion_Set == analysis_type) |>
    filter(Method != "ALL_ALTPCR" & Method != "ALL_ALTMTT" & Method != "ALL_ALTPCR_ALTMTT" & Method != "ALTMTT") |>
    filter(MA_Dist == distribution) |>
    mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
    select(MetricShortName, Method, Value_Denominator) |>
    pivot_wider(names_from = Method, values_from = Value_Denominator)

  rep_wide <- df_by_replication |>
    filter(Method != "ALL_ALTPCR" & Method != "ALL_ALTMTT" & Method != "ALL_ALTPCR_ALTMTT" & Method != "ALTMTT") |>
    filter(Inclusion_Set == analysis_type) |>
    filter(MA_Dist == distribution) |>
    mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
    select(MetricShortName, Method, Value_Denominator) |>
    pivot_wider(names_from = Method, values_from = Value_Denominator)

  # For bigexp distribution, add missing rows for commensurability with other tables
  if (distribution == "bigexp") {
    # Get data from t distribution for missing metrics
    exp_wide_t <- df_by_experiment |>
      filter(Inclusion_Set == analysis_type) |>
      filter(Method != "ALL_ALTPCR" & Method != "ALL_ALTMTT" & Method != "ALL_ALTPCR_ALTMTT" & Method != "ALTMTT") |>
      filter(MA_Dist == "t") |>
      mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
      select(MetricShortName, Method, Value_Denominator) |>
      pivot_wider(names_from = Method, values_from = Value_Denominator)
    
    rep_wide_t <- df_by_replication |>
      filter(Method != "ALL_ALTPCR" & Method != "ALL_ALTMTT" & Method != "ALL_ALTPCR_ALTMTT" & Method != "ALTMTT") |>
      filter(Inclusion_Set == analysis_type) |>
      filter(MA_Dist == "t") |>
      mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
      select(MetricShortName, Method, Value_Denominator) |>
      pivot_wider(names_from = Method, values_from = Value_Denominator)
    
    # Define all expected metrics for experiment level (in desired order)
    expected_exp_metrics <- c(
      "Original estimate within PI of REMA",
      "REMA estimate within CI of Original",
      "FEMA is significant and has same signal as original",
      "t-test majority vote (with ties as success)",
      "Subjective assessment majority vote (with ties as success)"
    )
    
    # Define all expected metrics for replication level (in desired order)
    expected_rep_metrics <- c(
      "Replication estimate within CI of Original",
      "Replication is significant and has same sign as original",
      "Subjective assessment, individual"
    )
    
    # Get available columns (methods)
    available_methods <- setdiff(colnames(exp_wide), "MetricShortName")
    
    # Add missing experiment-level metrics as NA (don't use t distribution for experiment level)
    for (metric in expected_exp_metrics) {
      if (!(metric %in% exp_wide$MetricShortName)) {
        new_row <- tibble(MetricShortName = metric)
        for (m in available_methods) new_row[[m]] <- "NA"
        exp_wide <- bind_rows(exp_wide, new_row)
      }
    }
    
    # Add missing replication-level metrics (use t distribution values)
    available_rep_methods <- setdiff(colnames(rep_wide), "MetricShortName")
    for (metric in expected_rep_metrics) {
      if (!(metric %in% rep_wide$MetricShortName)) {
        # Check if metric is available in t distribution
        if (metric %in% rep_wide_t$MetricShortName) {
          new_row <- rep_wide_t |> filter(MetricShortName == metric) |> select(MetricShortName, any_of(available_rep_methods))
          rep_wide <- bind_rows(rep_wide, new_row)
        } else {
          new_row <- tibble(MetricShortName = metric)
          for (m in available_rep_methods) new_row[[m]] <- "NA"
          rep_wide <- bind_rows(rep_wide, new_row)
        }
      }
    }
    
    # Reorder exp_wide according to expected_exp_metrics order
    exp_wide <- exp_wide |>
      mutate(order = match(MetricShortName, expected_exp_metrics)) |>
      arrange(order) |>
      select(-order)
    
    # Reorder rep_wide according to expected_rep_metrics order
    rep_wide <- rep_wide |>
      mutate(order = match(MetricShortName, expected_rep_metrics)) |>
      arrange(order) |>
      select(-order)
  }

  # Desired method order and present columns
  desired_methods <- c("ALL_PCR_MTT", "MTT", "PCR", "ALTPCR", "EPM")
  present_methods_exp <- desired_methods[desired_methods %in% colnames(exp_wide)]
  present_methods_rep <- desired_methods[desired_methods %in% colnames(rep_wide)]
  present_methods <- desired_methods[desired_methods %in% unique(c(present_methods_exp, present_methods_rep))]

  # Reorder columns
  if (length(present_methods_exp) > 0) {
    exp_wide <- exp_wide |>
      select(MetricShortName, all_of(present_methods_exp))
  }
  if (length(present_methods_rep) > 0) {
    rep_wide <- rep_wide |>
      select(MetricShortName, all_of(present_methods_rep))
  }

  # Header rows carrying column names (ensures no NA/"" in names)
  header_exp <- tibble(MetricShortName = "By experiment")
  for (m in present_methods) header_exp[[m]] <- ifelse(m == "ALL_PCR_MTT", "ALL", m)
  header_rep <- tibble(MetricShortName = "By replication")
  for (m in present_methods) header_rep[[m]] <- ifelse(m == "ALL_PCR_MTT", "ALL", m)

  # Combine in desired order
  tbl_by_method_all_exps <- bind_rows(
    header_exp,
    exp_wide,
    header_rep,
    rep_wide
  ) |>
    mutate(
      MetricShortName = case_when(
        MetricShortName == "FEMA is significant and has same signal as original" ~ "Same-sign significance (p<0.05)",
        MetricShortName == "Original estimate within PI of REMA" ~ "Original in replication’s 95% PI",
        MetricShortName == "REMA estimate within CI of Original" ~ "Replication in original 95% CI",
        MetricShortName == "Subjective assessment majority vote (with ties as success)" ~ "≥50% subjectively replicated",
        MetricShortName == "t-test majority vote (with ties as success)" ~ "≥50% replications significant",
        MetricShortName == "Subjective assessment, individual" ~ "Subjectively replicated",
        MetricShortName == "Replication estimate within CI of Original" ~ "Replication in original 95% CI",
        MetricShortName == "Replication is significant and has same sign as original" ~ "Same-sign significance (p<0.05)",
        TRUE ~ MetricShortName
      )
    ) |>
    select(-any_of(c("ALTPCR")))

  # Adjust column names and build the flextable
  colnames(tbl_by_method_all_exps) <- as.character(tbl_by_method_all_exps[1, ])

  footer_text <- paste0("Same-sign significance is based on a fixed meta-analysis estimate, while effect size comparisons are based on random-effects meta-analysis. All statistical tests use the ", toupper(distribution), " distribution. PI, prediction interval; CI, confidence interval. For more information on criteria, see https://osf.io/9rnuj.")

  tbl_by_method_all_exps <- tbl_by_method_all_exps |>
    slice(-1) |>
    flextable() |>
    bold(j = 2) |>
    bold(i = 6) |>
    hline(i = 5:6) |>
    add_footer_lines(value = as_paragraph(footer_text)) |>
    add_header_lines(values = paste0("Table S - Replication Rates - ", analysis_type, " - dist ", toupper(distribution))) |>
    bold(i = 2, part = "header") |>
    set_table_properties(layout = "autofit")

  return(tbl_by_method_all_exps)
}

create_tbl_by_method_altpcr <- function(analysis_type, distribution) {
  # Ensure that distribution is in lowercase
  distribution <- tolower(distribution)

  # Build wide tables and avoid relying on column order
  exp_wide <- df_by_experiment |>
    filter(Inclusion_Set == analysis_type) |>
    filter(Method != "ALL" & Method != "ALL_ALTMTT" & Method != "ALL_ALTPCR_ALTMTT" & Method != "ALTMTT") |>
    filter(MA_Dist == distribution) |>
    mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
    select(MetricShortName, Method, Value_Denominator) |>
    pivot_wider(names_from = Method, values_from = Value_Denominator)

  rep_wide <- df_by_replication |>
    filter(Method != "ALL" & Method != "ALL_ALTMTT" & Method != "ALL_ALTPCR_ALTMTT" & Method != "ALTMTT") |>
    filter(Inclusion_Set == analysis_type) |>
    filter(MA_Dist == distribution) |>
    mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
    select(MetricShortName, Method, Value_Denominator) |>
    pivot_wider(names_from = Method, values_from = Value_Denominator)

  # Desired method order and present columns
  desired_methods <- c("ALL_ALTPCR_MTT", "MTT", "PCR", "ALTPCR", "EPM")
  present_methods_exp <- desired_methods[desired_methods %in% colnames(exp_wide)]
  present_methods_rep <- desired_methods[desired_methods %in% colnames(rep_wide)]
  present_methods <- desired_methods[desired_methods %in% unique(c(present_methods_exp, present_methods_rep))]

  # Reorder columns
  if (length(present_methods_exp) > 0) {
    exp_wide <- exp_wide |>
      select(MetricShortName, all_of(present_methods_exp))
  }
  if (length(present_methods_rep) > 0) {
    rep_wide <- rep_wide |>
      select(MetricShortName, all_of(present_methods_rep))
  }

  # Header rows carrying column names (ensures no NA/"" in names)
  header_exp <- tibble(MetricShortName = "By experiment")
  for (m in present_methods) header_exp[[m]] <- ifelse(m == "ALL_ALTPCR_MTT", "ALL", m)
  header_rep <- tibble(MetricShortName = "By replication")
  for (m in present_methods) header_rep[[m]] <- ifelse(m == "ALL_ALTPCR_MTT", "ALL", m)

  # Combine in desired order
  tbl_by_method_all_exps <- bind_rows(
    header_exp,
    exp_wide,
    header_rep,
    rep_wide
  ) |>
    mutate(
      MetricShortName = case_when(
        MetricShortName == "FEMA is significant and has same signal as original" ~ "Same-sign significance (p<0.05)",
        MetricShortName == "Original estimate within PI of REMA" ~ "Original in replication’s 95% PI",
        MetricShortName == "REMA estimate within CI of Original" ~ "Replication in original 95% CI",
        MetricShortName == "Subjective assessment majority vote (with ties as success)" ~ "≥50% subjectively replicated",
        MetricShortName == "t-test majority vote (with ties as success)" ~ "≥50% replications significant",
        MetricShortName == "Subjective assessment, individual" ~ "Subjectively replicated",
        MetricShortName == "Replication estimate within CI of Original" ~ "Replication in original 95% CI",
        MetricShortName == "Replication is significant and has same sign as original" ~ "Same-sign significance (p<0.05)",
        TRUE ~ MetricShortName
      )
    ) |>
    select(-any_of(c("PCR")))

  # Adjust column names
  colnames(tbl_by_method_all_exps) <- as.character(tbl_by_method_all_exps[1, ])

  # Build the flextable
  footer_text <- paste0("Same-sign significance is based on a fixed meta-analysis estimate, while effect size comparisons are based on random-effects meta-analysis. All statistical tests use the ", toupper(distribution), " distribution. PI, prediction interval; CI, confidence interval. For more information on criteria, see https://osf.io/9rnuj.")

  tbl_df <- tbl_by_method_all_exps |> slice(-1)
  
  ft <- flextable(tbl_df)
  
  # column 2 may not exist if only 1 method column survived filtering
  if (ncol(tbl_df) >= 2) {
    ft <- bold(ft, j = 2)
  }
  
  # only apply row-based styling if those rows exist
  if (nrow(tbl_df) >= 6) {
    ft <- bold(ft, i = 6)
  }
  if (nrow(tbl_df) >= 6) {
    ft <- hline(ft, i = 5:6)
  } else if (nrow(tbl_df) >= 5) {
    ft <- hline(ft, i = 5)
  }
  
  ft <- ft |>
    add_footer_lines(value = as_paragraph(footer_text)) |>
    add_header_lines(values = paste0("Table S - Replication Rates ...")) |>
    bold(i = 2, part = "header") |>
    set_table_properties(layout = "autofit")
  
  return(ft)
  
}

### PCR ---------------------------------------------------------------------

# Setting up
# List of analysis methods
methods <- c(
  "primary",
  "all_exps_lab_units",
  "at_least_2_reps",
  "included_by_lab",
  "only_3_reps",
  "only_80_power_a_posteriori_T",
  "only_80_power_a_posteriori_Z",
  "only_80_power_a_posteriori_KNHA"
)

# List of distributions to check
distributions <- c("t", "z", "knha", "bigexp")

# Loop through each method and distribution
for (method in methods) {
  for (distribution in distributions) {
    # Create a regex pattern to check if the analysis exists in list_of_analyses
    analysis_pattern <- paste0(".*", method, " ", distribution, "$")

    # If the analysis is present, execute the function and store the result
    if (any(grepl(analysis_pattern, list_of_analyses))) {
      assign(
        paste0("tbl_by_method_", method, "_", distribution),
        create_tbl_by_method_pcr(method, distribution)
      )
    }
  }
}

#### t -----------------------------------------------------------------------

# Creating a list of existing tables for saving
existing_tables_t <- list()

# Define expected table names
table_names_t <- c(
  "tbl_by_method_primary_t",
  "tbl_by_method_all_exps_lab_units_t",
  "tbl_by_method_at_least_2_reps_t",
  "tbl_by_method_included_by_lab_t",
  "tbl_by_method_only_3_reps_t",
  "tbl_by_method_only_80_power_a_posteriori_T_t"
)

# Check if each table exists before adding it to the list
for (table in table_names_t) {
  if (exists(table)) {
    existing_tables_t <- append(existing_tables_t, list(get(table)))

    # Create separated tables

    if (table == "tbl_by_method_primary_t") {
      save_tbl(
        list(get(table)),
        paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table 1.docx")
      )
      cat("\n### Table 1 generated! ###\n")
    }

    if (table == "tbl_by_method_only_80_power_a_posteriori_T_t") {
      save_tbl(
        list(get(table)),
        paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S8.docx")
      )
      cat("\n### Table S8 generated! ###\n")
    }

    if (table == "tbl_by_method_all_exps_lab_units_t") {
      save_tbl(
        list(get(table)),
        paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S9.docx")
      )
      cat("\n### Table S9 generated! ###\n")
    }

    if (table == "tbl_by_method_included_by_lab_t") {
      save_tbl(
        list(get(table)),
        paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S10.docx")
      )
      cat("\n### Table S10 generated! ###\n")
    }

    if (table == "tbl_by_method_at_least_2_reps_t") {
      save_tbl(
        list(get(table)),
        paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S11.docx")
      )
      cat("\n### Table S11 generated! ###\n")
    }

    if (table == "tbl_by_method_only_3_reps_t") {
      save_tbl(
        list(get(table)),
        paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S12.docx")
      )
      cat("\n### Table S12 generated! ###\n")
    }
  }
}

# Save tables only if at least one was created
if (length(existing_tables_t) > 0) {
  save_tbl(
    existing_tables_t,
    paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table - Replication Rates - PCR - t.docx")
  )
  cat("\n### Table - Replication Rates - PCR - t generated! ###\n")
} else {
  message("No replication rates table were generated and saved.")
}

#### z -----------------------------------------------------------------------

# Creating a list of existing tables for saving
existing_tables_z <- list()

# Define expected table names
table_names_z <- c(
  "tbl_by_method_primary_z",
  "tbl_by_method_all_exps_lab_units_z",
  "tbl_by_method_at_least_2_reps_z",
  "tbl_by_method_included_by_lab_z",
  "tbl_by_method_only_3_reps_z",
  "tbl_by_method_only_80_power_a_posteriori_Z_z"
)

# Check if each table exists before adding it to the list
for (table in table_names_z) {
  if (exists(table)) {
    existing_tables_z <- append(existing_tables_z, list(get(table)))

    if (table == "tbl_by_method_primary_z") {
      save_tbl(
        list(get(table)),
        paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S13.docx")
      )
      cat("\n### Table S13 generated! ###\n")
    }
  }
}

# Save tables only if at least one was created
if (length(existing_tables_z) > 0) {
  save_tbl(
    existing_tables_z,
    paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table - Replication Rates - PCR - z.docx")
  )
  cat("\n### Table - Replication Rates - PCR - z generated! ###\n")
} else {
  message("No tables were generated and saved.")
}

####  khna -----------------------------------------------------------------------

# Creating a list of existing tables for saving
existing_tables_knha <- list()

# Define expected table names
table_names_knha <- c(
  "tbl_by_method_primary_knha",
  "tbl_by_method_all_exps_lab_units_knha",
  "tbl_by_method_at_least_2_reps_knha",
  "tbl_by_method_included_by_lab_knha",
  "tbl_by_method_only_3_reps_knha",
  "tbl_by_method_only_80_power_a_posteriori_KNHA_knha"
)

# Check if each table exists before adding it to the list
for (table in table_names_knha) {
  if (exists(table)) {
    existing_tables_knha <- append(existing_tables_knha, list(get(table)))

    if (table == "tbl_by_method_primary_knha") {
      save_tbl(
        list(get(table)),
        paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S15.docx")
      )
      cat("\n### Table S15 generated! ###\n")
    }
  }
}

# Save tables only if at least one was created
if (length(existing_tables_knha) > 0) {
  save_tbl(
    existing_tables_knha,
    paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table - Replication Rates - PCR - knha.docx")
  )
  cat("\n### Table - Replication Rates - PCR - knha generated! ###\n")
} else {
  message("No tables were generated and saved.")
}

#### bigexp -----------------------------------------------------------------------

# Creating a list of existing tables for saving
existing_tables_bigexp <- list()

# Define expected table names
table_names_bigexp <- c(
  "tbl_by_method_primary_bigexp",
  "tbl_by_method_all_exps_lab_units_bigexp",
  "tbl_by_method_at_least_2_reps_bigexp",
  "tbl_by_method_included_by_lab_bigexp",
  "tbl_by_method_only_3_reps_bigexp"
)

# Check if each table exists before adding it to the list
for (table in table_names_bigexp) {
  if (exists(table)) {
    existing_tables_bigexp <- append(existing_tables_bigexp, list(get(table)))

    if (table == "tbl_by_method_primary_bigexp") {
      save_tbl(
        list(get(table)),
        paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S17.docx")
      )
      cat("\n### Table S17 generated! ###\n")
    }


  }
}

# Save tables only if at least one was created
if (length(existing_tables_bigexp) > 0) {
  save_tbl(
    existing_tables_bigexp,
    paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table - Replication Rates - PCR - bigexp.docx")
  )
  cat("\n### Table - Replication Rates - PCR - bigexp generated! ###\n")
} else {
  message("No bigexp tables were generated and saved.")
}

### ALTPCR ---------------------------------------------------------------------

# Loop through each method and distribution
for (method in methods) {
  for (distribution in distributions) {
    # Create a regex pattern to check if the analysis exists in list_of_analyses
    analysis_pattern <- paste0(".*", method, " ", distribution, "$")

    # If the analysis is present, execute the function and store the result
    if (any(grepl(analysis_pattern, list_of_analyses))) {
      assign(
        paste0("tbl_by_method_", method, "_", distribution, "_ALTPCR"),
        create_tbl_by_method_altpcr(method, distribution)
      )
    }
  }
}

#### t -----------------------------------------------------------------------

# Creating a list of existing ALTPCR tables for t-distribution
existing_tables_t_altpcr <- list()

# Define expected table names
table_names_t_altpcr <- c(
  "tbl_by_method_primary_t_ALTPCR",
  "tbl_by_method_all_exps_lab_units_t_ALTPCR",
  "tbl_by_method_at_least_2_reps_t_ALTPCR",
  "tbl_by_method_included_by_lab_t_ALTPCR",
  "tbl_by_method_only_3_reps_t_ALTPCR",
  "tbl_by_method_only_80_power_a_posteriori_T_t_ALTPCR"
)

# Check if each table exists before adding it to the list
for (table in table_names_t_altpcr) {
  if (exists(table)) {
    existing_tables_t_altpcr <- append(existing_tables_t_altpcr, list(get(table)))
  }
}

# Save tables only if at least one was created
if (length(existing_tables_t_altpcr) > 0) {
  save_tbl(
    existing_tables_t_altpcr,
    paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table - Replication Rates - ALTPCR - t.docx")
  )
  cat("\n### Table - Replication Rates - ALTPCR - t generated! ###\n")
} else {
  message("No ALTPCR tables for t-distribution were generated and saved.")
}

#### Z -----------------------------------------------------------------------

# Creating a list of existing ALTPCR tables for z-distribution
existing_tables_z_altpcr <- list()

# Define expected table names
table_names_z_altpcr <- c(
  "tbl_by_method_primary_z_ALTPCR",
  "tbl_by_method_all_exps_lab_units_z_ALTPCR",
  "tbl_by_method_at_least_2_reps_z_ALTPCR",
  "tbl_by_method_included_by_lab_z_ALTPCR",
  "tbl_by_method_only_3_reps_z_ALTPCR",
  "tbl_by_method_only_80_power_a_posteriori_Z_z_ALTPCR"
)

# Check if each table exists before adding it to the list
for (table in table_names_z_altpcr) {
  if (exists(table)) {
    existing_tables_z_altpcr <- append(existing_tables_z_altpcr, list(get(table)))
  }
}

# Save tables only if at least one was created
if (length(existing_tables_z_altpcr) > 0) {
  save_tbl(
    existing_tables_z_altpcr,
    paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table - Replication Rates - ALTPCR - z.docx")
  )
  cat("\n### Table - Replication Rates - ALTPCR - z generated! ###\n")
} else {
  message("No ALTPCR tables for z-distribution were generated and saved.")
}

#### knha -----------------------------------------------------------------------

# Creating a list of existing ALTPCR tables for z-distribution
existing_tables_knha_altpcr <- list()

# Define expected table names
table_names_knha_altpcr <- c(
  "tbl_by_method_primary_knha_ALTPCR",
  "tbl_by_method_all_exps_lab_units_knha_ALTPCR",
  "tbl_by_method_at_least_2_reps_knha_ALTPCR",
  "tbl_by_method_included_by_lab_knha_ALTPCR",
  "tbl_by_method_only_3_reps_knha_ALTPCR",
  "tbl_by_method_only_80_power_a_posteriori_KNHA_knha_ALTPCR"
)

# Check if each table exists before adding it to the list
for (table in table_names_knha_altpcr) {
  if (exists(table)) {
    existing_tables_knha_altpcr <- append(existing_tables_knha_altpcr, list(get(table)))
  }
}

# Save tables only if at least one was created
if (length(existing_tables_knha_altpcr) > 0) {
  save_tbl(
    existing_tables_knha_altpcr,
    paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table - Replication Rates - ALTPCR - knha.docx")
  )
  cat("\n### Table - Replication Rates - ALTPCR - knha generated! ###\n")
} else {
  message("No ALTPCR tables for knha-distribution were generated and saved.")
}

#### bigexp -----------------------------------------------------------------------

# Creating a list of existing ALTPCR tables for bigexp
existing_tables_bigexp_altpcr <- list()

# Define expected table names
table_names_bigexp_altpcr <- c(
  "tbl_by_method_primary_bigexp_ALTPCR",
  "tbl_by_method_all_exps_lab_units_bigexp_ALTPCR",
  "tbl_by_method_at_least_2_reps_bigexp_ALTPCR",
  "tbl_by_method_included_by_lab_bigexp_ALTPCR",
  "tbl_by_method_only_3_reps_bigexp_ALTPCR"
)

# Check if each table exists before adding it to the list
for (table in table_names_bigexp_altpcr) {
  if (exists(table)) {
    existing_tables_bigexp_altpcr <- append(existing_tables_bigexp_altpcr, list(get(table)))
  }
}

# Save tables only if at least one was created
if (length(existing_tables_bigexp_altpcr) > 0) {
  save_tbl(
    existing_tables_bigexp_altpcr,
    paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table - Replication Rates - ALTPCR - bigexp.docx")
  )
  cat("\n### Table - Replication Rates - ALTPCR - bigexp generated! ###\n")
} else {
  message("No ALTPCR tables for bigexp were generated and saved.")
}

## Table 2 -----------------------------------------------------
tbl_2 <- df_by_experiment |>
  filter(Method == "ALL_PCR_MTT") |>
  filter(MA_Dist == "t") |>
  mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
  select(MetricShortName, Method, Value_Denominator, Inclusion_Set) |>
  pivot_wider(names_from = c(Inclusion_Set, Method), values_from = Value_Denominator) |>
  add_row(
    data.frame(
      MetricShortName = "By experiment"
    ),
    .before = 0
  ) |>
  add_row(
    df_by_replication |>
      filter(Method == "ALL_PCR_MTT") |>
      filter(MA_Dist == "t") |>
      mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
      select(MetricShortName, Method, Value_Denominator, Inclusion_Set) |>
      pivot_wider(names_from = c(Inclusion_Set, Method), values_from = Value_Denominator) |>
      add_row(
        data.frame(
          MetricShortName = "By replication"
        ),
        .before = 0
      )
  ) |>
  mutate(
    MetricShortName = case_when(
      MetricShortName == "FEMA is significant and has same signal as original" ~ "Same-sign significance (p<0.05)",
      MetricShortName == "Original estimate within PI of REMA" ~ "Original in replication’s 95% PI",
      MetricShortName == "REMA estimate within CI of Original" ~ "Replication in original 95% CI",
      MetricShortName == "Subjective assessment majority vote (with ties as success)" ~ "≥50% subjectively replicated",
      MetricShortName == "t-test majority vote (with ties as success)" ~ "≥50% replications significant",
      MetricShortName == "Subjective assessment, individual" ~ "Subjectively replicated",
      MetricShortName == "Replication estimate within CI of Original" ~ "Replication in original 95% CI",
      MetricShortName == "Replication is significant and has same sign as original" ~ "Same-sign significance (p<0.05)",
      TRUE ~ MetricShortName
    )
  ) |>
  slice(1:7, 8, 9, 10)

ordered_cols_tbl_2 <- tibble(
  colname = colnames(tbl_2),
  priority = case_when(
    str_starts(colname, "MetricShortName") ~ 1,
    str_starts(colname, "primary") ~ 2,
    str_starts(colname, "included_by_lab") ~ 3,
    str_starts(colname, "all_exps_lab_units") ~ 4,
    str_starts(colname, "at_least_2_reps") ~ 5,
    str_starts(colname, "only_3_reps") ~ 6,
    str_starts(colname, "only_80_power_a_posteriori") ~ 7,
    TRUE ~ 8
  )
) |>
  arrange(priority) |>
  pull(colname)

tbl_2 <- tbl_2 |>
  select(any_of(ordered_cols_tbl_2)) |>
  dplyr::rename(any_of(c(
    `Primary` = "primary_ALL_PCR_MTT",
    `Lab's choice` = "included_by_lab_ALL_PCR_MTT",
    `All Exps` = "all_exps_lab_units_ALL_PCR_MTT",
    `≥ 2 copies` = "at_least_2_reps_ALL_PCR_MTT",
    `3 copies` = "only_3_reps_ALL_PCR_MTT",
    `≥80% power` = "only_80_power_a_posteriori_T_ALL_PCR_MTT"
  )))

footer_text_tbl_2 <- "Replication rates for the primary and secondary analyses. Effect size comparisons are based on random-effects meta-analysis, while same-sign significance is based on a fixed meta-analysis estimate. The prediction interval criterion only uses experiments with more than one replication. Subsets for secondary analyses include all experiments judged valid by the replicating lab (Lab’s Choice), all concluded experiments (All Exps), both of which use the experimental unit as defined by the lab, only experiments with at least 2 and 3 replications, and only experiments with ≥ 80% a posteriori power. All statistical tests use the t distribution. PI, prediction interval; CI, confidence interval. For more information on replication criteria, see https://osf.io/9rnuj."

tbl_2 <- tbl_2 |>
  slice(-1) |>
  flextable() |>
  bold(j = 2) |>
  bold(i = 6) |>
  hline(i = 5:6) |>
  set_header_labels(MetricShortName = "By experiment") |>
  add_footer_lines(value = as_paragraph(footer_text_tbl_2)) |>
  add_header_lines(values = paste0("Table 2 - Replication rates for different analysis sets.")) |>
  bold(i = 2, part = "header") |>
  set_table_properties(layout = "autofit")

save_tbl(
  tbl_2,
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table 2.docx")
)
cat("\n### Table 2 generated! ###\n")

## Table 3 -----------------------------------------------------------------

df_assessment_by_experiment <- read_tsv(paste0("output/", results_path, "/primary t/Replication Assessment by Experiment.tsv"), show_col_types = FALSE) |> # You can get it from output folder: primary t result
  mutate(category = case_when(
    str_detect(EXP, "ALTPCR") == TRUE ~ "ALTPCR",
    str_detect(EXP, "ALTMTT") == TRUE ~ "ALTMTT",
    str_detect(EXP, "EPM") == TRUE ~ "EPM",
    str_detect(EXP, "PCR") == TRUE ~ "PCR",
    str_detect(EXP, "^MTT") == TRUE ~ "MTT"
  )) |>
  relocate(
    category,
    .after = EXP
  ) |>
  # Linear scale!
  mutate(
    corrected_sign_original_es = exp(corrected_sign_original_es),
    corrected_sign_replication_es = exp(corrected_sign_replication_es),
    mean_abs_diff_reps_linear = exp(mean_abs_diff_reps),
    mean_abs_diff_reps_orig_linear = exp(mean_abs_diff_reps_orig)
  )


df_success_by_experiment <- read_tsv(paste0("output/", results_path, "/primary t/Replication Success by Experiment.tsv"), show_col_types = FALSE) |>
  mutate(category = case_when(
    str_detect(EXP, "ALTPCR") == TRUE ~ "ALTPCR",
    str_detect(EXP, "ALTMTT") == TRUE ~ "ALTMTT",
    str_detect(EXP, "EPM") == TRUE ~ "EPM",
    str_detect(EXP, "PCR") == TRUE ~ "PCR",
    str_detect(EXP, "^MTT") == TRUE ~ "MTT"
  )) |>
  relocate(
    category,
    .after = EXP
  )

get_round_digits <- function(x) {
  ifelse(abs(x) < 10, 2, 1)
}

round_digits <- 2

### Effect size (original)  -------------------------------------------------

effect_size_original_all <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  summarise(
    median_range_median_corrected_sign_original_es = paste0(
      sprintf("%.2f", median(corrected_sign_original_es)),
      "\n(",
      sprintf("%.2f", min(corrected_sign_original_es)),
      ", ",
      round(max(corrected_sign_original_es), get_round_digits(max(corrected_sign_original_es))),
      ")"
    )
  )

effect_size_original_exp <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  group_by(category) |>
  summarise(
    median_range_median_corrected_sign_original_es = paste0(
      round(median(corrected_sign_original_es), get_round_digits(median(corrected_sign_original_es))),
      "\n(",
      sprintf("%.2f", min(corrected_sign_original_es)),
      ", ",
      round(max(corrected_sign_original_es), get_round_digits(max(corrected_sign_original_es))),
      ")"
    )
  )

effect_size_original_all[[1]]
effect_size_original_epm <- effect_size_original_exp$median_range_median_corrected_sign_original_es[effect_size_original_exp$category == "EPM"]
effect_size_original_pcr <- effect_size_original_exp$median_range_median_corrected_sign_original_es[effect_size_original_exp$category == "PCR"]
effect_size_original_mtt <- effect_size_original_exp$median_range_median_corrected_sign_original_es[effect_size_original_exp$category == "MTT"]


### Effect size (replication)  ----------------------------------------------

effect_size_replication_all <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  summarise(
    median_range_median_corrected_sign_replication_es = paste0(
      sprintf("%.2f", median(corrected_sign_replication_es)),
      "\n(",
      sprintf("%.2f", min(corrected_sign_replication_es)),
      ", ",
      round(max(corrected_sign_replication_es), get_round_digits(max(corrected_sign_replication_es))),
      ")"
    )
  )

effect_size_replication_exp <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  group_by(category) |>
  summarise(
    median_range_median_corrected_sign_replication_es = paste0(
      format(round(median(corrected_sign_replication_es), round_digits), nsmall = 2),
      "\n(",
      format(round(min(corrected_sign_replication_es), round_digits), nsmall = 2),
      ", ",
      format(round(max(corrected_sign_replication_es), round_digits), nsmall = 2),
      ")"
    )
  )

effect_size_replication_all[[1]]
effect_size_replication_epm <- effect_size_replication_exp$median_range_median_corrected_sign_replication_es[effect_size_replication_exp$category == "EPM"]
effect_size_replication_pcr <- effect_size_replication_exp$median_range_median_corrected_sign_replication_es[effect_size_replication_exp$category == "PCR"]
effect_size_replication_mtt <- effect_size_replication_exp$median_range_median_corrected_sign_replication_es[effect_size_replication_exp$category == "MTT"]


### Effect size ratio -------------------------------------------------------

effect_size_ratio_log_diff_all <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  summarise(
    median_range_median_es_ratio = paste0(
      sprintf("%.2f", median(es_ratio)),
      "\n(",
      round(min(es_ratio), 2),
      ", ",
      round(max(es_ratio), get_round_digits(max(es_ratio))),
      ")"
    )
  )

effect_size_ratio_log_diff_exp <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  group_by(category) |>
  summarise(
    median_range_median_es_ratio = paste0(
      round(median(es_ratio), 2),
      "\n(",
      sprintf("%.2f", min(es_ratio)),
      ", ",
      round(max(es_ratio), get_round_digits(max(es_ratio))),
      ")"
    )
  )

effect_size_ratio_log_diff_all[[1]]
effect_size_ratio_log_diff_epm <- effect_size_ratio_log_diff_exp$median_range_median_es_ratio[effect_size_ratio_log_diff_exp$category == "EPM"]
effect_size_ratio_log_diff_pcr <- effect_size_ratio_log_diff_exp$median_range_median_es_ratio[effect_size_ratio_log_diff_exp$category == "PCR"]
effect_size_ratio_log_diff_mtt <- effect_size_ratio_log_diff_exp$median_range_median_es_ratio[effect_size_ratio_log_diff_exp$category == "MTT"]


### Coefficient of variation (original) -------------------------------------

coefficient_variation_original_all <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  summarise(
    median_range_median_original_cv = paste0(
      sprintf("%.2f", median(original_cv, na.rm = TRUE)),
      "\n(",
      round(min(original_cv, na.rm = TRUE), 2),
      ", ",
      sprintf("%.2f", max(original_cv, na.rm = TRUE)),
      ")"
    )
  )

coefficient_variation_original_exp <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  group_by(category) |>
  summarise(
    median_range_median_original_cv = paste0(
      round(median(original_cv, na.rm = TRUE), 2),
      "\n(",
      sprintf("%.2f", min(original_cv, na.rm = TRUE)),
      ", ",
      sprintf("%.2f", max(original_cv, na.rm = TRUE)),
      ")"
    )
  )

coefficient_variation_original_all[[1]]
coefficient_variation_original_epm <- coefficient_variation_original_exp$median_range_median_original_cv[coefficient_variation_original_exp$category == "EPM"]
coefficient_variation_original_pcr <- coefficient_variation_original_exp$median_range_median_original_cv[coefficient_variation_original_exp$category == "PCR"]
coefficient_variation_original_mtt <- coefficient_variation_original_exp$median_range_median_original_cv[coefficient_variation_original_exp$category == "MTT"]


### Coefficient of variation (replication) ---------------------------------

coefficient_variation_replication_all <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","ALTPCR","MTT")) |>
  summarise(
    median_range_median_replication_cv = paste0(
      sprintf("%.2f", median(replication_cv, na.rm = TRUE)),
      "\n(",
      sprintf("%.2f", min(replication_cv, na.rm = TRUE)),
      ", ",
      round(max(replication_cv, na.rm = TRUE), 2),
      ")"
    )
  )

coefficient_variation_replication_exp <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","ALTPCR","MTT")) |>
  group_by(category) |>
  summarise(
    median_range_median_replication_cv = paste0(
      round(median(replication_cv, na.rm = TRUE), 2),
      "\n(",
      format(round(min(replication_cv, na.rm = TRUE), 2), nsmall = 2),
      ", ",
      round(max(replication_cv, na.rm = TRUE), 2),
      ")"
    )
  )

coefficient_variation_replication_all[[1]]
coefficient_variation_replication_epm <- coefficient_variation_replication_exp$median_range_median_replication_cv[coefficient_variation_replication_exp$category == "EPM"]
coefficient_variation_replication_altpcr <- coefficient_variation_replication_exp$median_range_median_replication_cv[coefficient_variation_replication_exp$category == "ALTPCR"]
coefficient_variation_replication_mtt <- coefficient_variation_replication_exp$median_range_median_replication_cv[coefficient_variation_replication_exp$category == "MTT"]



### Coefficient of Variation (ratio) ----------------------------------------

df_assessment_by_experiment <- df_assessment_by_experiment |>
  mutate(cv_ratio = original_cv / replication_cv)

coefficient_variation_ratio_all <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","ALTPCR","MTT")) |>
  summarise(
    median_range_median_cv_ratio = paste0(
      sprintf("%.2f", median(cv_ratio, na.rm = TRUE)),
      "\n(",
      round(min(cv_ratio, na.rm = TRUE), round_digits),
      ", ",
      round(max(cv_ratio, na.rm = TRUE), round_digits),
      ")"
    )
  )

coefficient_variation_ratio_exp <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","ALTPCR","MTT")) |>
  group_by(category) |>
  summarise(
    median_range_median_cv_ratio = paste0(
      sprintf("%.2f", median(cv_ratio, na.rm = TRUE)),
      "\n(",
      round(min(cv_ratio, na.rm = TRUE), round_digits),
      ", ",
      round(max(cv_ratio, na.rm = TRUE), round_digits),
      ")"
    )
  )

coefficient_variation_ratio_all[[1]]
coefficient_variation_ratio_epm <- coefficient_variation_ratio_exp$median_range_median_cv_ratio[coefficient_variation_ratio_exp$category == "EPM"]
coefficient_variation_ratio_pcr <- coefficient_variation_ratio_exp$median_range_median_cv_ratio[coefficient_variation_ratio_exp$category == "ALTPCR"]
coefficient_variation_ratio_mtt <- coefficient_variation_ratio_exp$median_range_median_cv_ratio[coefficient_variation_ratio_exp$category == "MTT"]


### Mean effect size difference (original vs. individual replications) ---------

mean_effect_size_diff_original_vs_individual_all <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  summarise(
    median_range_median_mean_abs_diff_reps_orig = paste0(
      sprintf("%.2f", median(mean_abs_diff_reps_orig_linear, na.rm = TRUE)),
      "\n(",
      round(min(mean_abs_diff_reps_orig_linear, na.rm = TRUE), round_digits),
      ", ",
      round(max(mean_abs_diff_reps_orig_linear, na.rm = TRUE), round_digits),
      ")"
    )
  )

mean_effect_size_diff_original_vs_individual_exp <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  group_by(category) |>
  summarise(
    median_range_median_mean_abs_diff_reps_orig = paste0(
      sprintf("%.2f", median(mean_abs_diff_reps_orig_linear, na.rm = TRUE)),
      "\n(",
      round(min(mean_abs_diff_reps_orig_linear, na.rm = TRUE), round_digits),
      ", ",
      round(max(mean_abs_diff_reps_orig_linear, na.rm = TRUE), round_digits),
      ")"
    )
  )

mean_effect_size_diff_original_vs_individual_all <- mean_effect_size_diff_original_vs_individual_all[[1]]
mean_effect_size_diff_original_vs_individual_epm <- mean_effect_size_diff_original_vs_individual_exp$median_range_median_mean_abs_diff_reps_orig[mean_effect_size_diff_original_vs_individual_exp$category == "EPM"]
mean_effect_size_diff_original_vs_individual_pcr <- mean_effect_size_diff_original_vs_individual_exp$median_range_median_mean_abs_diff_reps_orig[mean_effect_size_diff_original_vs_individual_exp$category == "PCR"]
mean_effect_size_diff_original_vs_individual_mtt <- mean_effect_size_diff_original_vs_individual_exp$median_range_median_mean_abs_diff_reps_orig[mean_effect_size_diff_original_vs_individual_exp$category == "MTT"]


### Mean effect size difference (between individual replications) ---------

mean_effect_size_diff_between_individual_all <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  summarise(
    median_range_median_mean_abs_diff_reps = paste0(
      sprintf("%.2f", median(mean_abs_diff_reps_linear, na.rm = TRUE)),
      "\n(",
      round(min(mean_abs_diff_reps_linear, na.rm = TRUE), round_digits),
      ", ",
      round(max(mean_abs_diff_reps_linear, na.rm = TRUE), 1),
      ")"
    )
  )

mean_effect_size_diff_between_individual_exp <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  group_by(category) |>
  summarise(
    median_range_median_mean_abs_diff_reps = paste0(
      sprintf("%.2f", median(mean_abs_diff_reps_linear, na.rm = TRUE)),
      "\n(",
      sprintf("%.2f", min(mean_abs_diff_reps_linear, na.rm = TRUE)),
      ", ",
      round(max(mean_abs_diff_reps_linear, na.rm = T), get_round_digits(max(mean_abs_diff_reps_linear, na.rm = T))),
      ")"
    )
  )

mean_effect_size_diff_between_individual_all <- mean_effect_size_diff_between_individual_all[[1]]
mean_effect_size_diff_between_individual_epm <- mean_effect_size_diff_between_individual_exp$median_range_median_mean_abs_diff_reps[mean_effect_size_diff_between_individual_exp$category == "EPM"]
mean_effect_size_diff_between_individual_pcr <- mean_effect_size_diff_between_individual_exp$median_range_median_mean_abs_diff_reps[mean_effect_size_diff_between_individual_exp$category == "PCR"]
mean_effect_size_diff_between_individual_mtt <- mean_effect_size_diff_between_individual_exp$median_range_median_mean_abs_diff_reps[mean_effect_size_diff_between_individual_exp$category == "MTT"]


### Ratio of mean absolute differences ----------------------------------------


df_assessment_by_experiment <- df_assessment_by_experiment |>
  mutate(
    mad_ratio = ifelse(mean_abs_diff_reps == 0, NA, exp(mean_abs_diff_reps_orig - mean_abs_diff_reps))
  )

ratio_mean_absolute_diff_all <- df_assessment_by_experiment |>
  filter(category != "ALTPCR") |>
  summarise(
    median_range_median_mad_ratio = paste0(
      sprintf("%.2f", median(mad_ratio, na.rm = TRUE)),
      "\n(",
      format(round(min(mad_ratio, na.rm = TRUE), 2), nsmall = 2),
      ", ",
      round(max(mad_ratio, na.rm = TRUE), 1),
      ")"
    )
  )

ratio_mean_absolute_diff_exp <- df_assessment_by_experiment |>
  filter(category != "ALTPCR") |>
  group_by(category) |>
  summarise(
    median_range_median_mad_ratio = paste0(
      sprintf("%.2f", median(mad_ratio, na.rm = TRUE)),
      "\n(",
      format(round(min(mad_ratio, na.rm = TRUE), 2), nsmall = 2),
      ", ",
      round(max(mad_ratio, na.rm = T), get_round_digits(max(mad_ratio, na.rm = T))),
      ")"
    )
  )

ratio_mean_absolute_diff_all <- ratio_mean_absolute_diff_all[[1]]
ratio_mean_absolute_diff_epm <- ratio_mean_absolute_diff_exp$median_range_median_mad_ratio[ratio_mean_absolute_diff_exp$category == "EPM"]
ratio_mean_absolute_diff_pcr <- ratio_mean_absolute_diff_exp$median_range_median_mad_ratio[ratio_mean_absolute_diff_exp$category == "PCR"]
ratio_mean_absolute_diff_mtt <- ratio_mean_absolute_diff_exp$median_range_median_mad_ratio[ratio_mean_absolute_diff_exp$category == "MTT"]

### Signal error (% of significant) ----------------------------------------

signal_error_significant_all <- df_success_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  summarise(
    signal_error_percentage = paste0(
      sum(SignalErrorAll, na.rm = TRUE),
      "/",
      sum(REP_FEMA_SSS_Orig, na.rm = TRUE),
      " (",
      sprintf("%.0f%%", (sum(SignalErrorAll, na.rm = TRUE) / sum(REP_FEMA_SSS_Orig, na.rm = TRUE)) * 100),
      ")"
    )
  )

signal_error_significant_exp <- df_success_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  group_by(category) |>
  summarise(
    signal_error_percentage = paste0(
      sum(SignalErrorAll, na.rm = TRUE),
      "/",
      sum(REP_FEMA_SSS_Orig, na.rm = TRUE),
      " (",
      sprintf("%.0f%%", (sum(SignalErrorAll, na.rm = TRUE) / sum(REP_FEMA_SSS_Orig, na.rm = TRUE)) * 100),
      ")"
    )
  )

signal_error_significant_all <- signal_error_significant_all$signal_error_percentage

signal_error_significant_epm <- signal_error_significant_exp$signal_error_percentage[signal_error_significant_exp$category == "EPM"]
signal_error_significant_pcr <- signal_error_significant_exp$signal_error_percentage[signal_error_significant_exp$category == "PCR"]
signal_error_significant_mtt <- signal_error_significant_exp$signal_error_percentage[signal_error_significant_exp$category == "MTT"]

### Signal error (% of total) ----------------------------------------

signal_error_total_all <- df_success_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  summarise(
    signal_error_percentage = paste0(
      sum(SignalErrorAll, na.rm = TRUE),
      "/",
      n(),
      " (",
      sprintf("%.0f%%", (sum(SignalErrorAll, na.rm = TRUE) / n()) * 100),
      ")"
    )
  )

signal_error_total_exp <- df_success_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  group_by(category) |>
  summarise(
    signal_error_percentage = paste0(
      sum(SignalErrorAll, na.rm = TRUE),
      "/",
      n(),
      " (",
      sprintf("%.0f%%", (sum(SignalErrorAll, na.rm = TRUE) / n()) * 100),
      ")"
    )
  )

signal_error_total_all <- signal_error_total_all$signal_error_percentage

signal_error_total_epm <- signal_error_total_exp$signal_error_percentage[signal_error_total_exp$category == "EPM"]
signal_error_total_pcr <- signal_error_total_exp$signal_error_percentage[signal_error_total_exp$category == "PCR"]
signal_error_total_mtt <- signal_error_total_exp$signal_error_percentage[signal_error_total_exp$category == "MTT"]




### Opposite signal ----------------------------------------------------------

opposite_signal_total_all <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  summarise(
    opposite_signal_percentage = paste0(
      sum((original_es > 0 & replication_es < 0) | (original_es < 0 & replication_es > 0), na.rm = TRUE),
      "/",
      n(),
      " (",
      sprintf("%.0f%%", (sum((original_es > 0 & replication_es < 0) | (original_es < 0 & replication_es > 0), na.rm = TRUE) / n()) * 100),
      ")"
    )
  )

opposite_signal_exp <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  group_by(category) |>
  summarise(
    opposite_signal_percentage = paste0(
      sum((original_es > 0 & replication_es < 0) | (original_es < 0 & replication_es > 0), na.rm = TRUE),
      "/",
      n(),
      " (",
      sprintf("%.0f%%", (sum((original_es > 0 & replication_es < 0) | (original_es < 0 & replication_es > 0), na.rm = TRUE) / n()) * 100),
      ")"
    )
  )

opposite_signal_total_all <- opposite_signal_total_all$opposite_signal_percentage

opposite_signal_epm <- opposite_signal_exp$opposite_signal_percentage[opposite_signal_exp$category == "EPM"]
opposite_signal_pcr <- opposite_signal_exp$opposite_signal_percentage[opposite_signal_exp$category == "PCR"]
opposite_signal_mtt <- opposite_signal_exp$opposite_signal_percentage[opposite_signal_exp$category == "MTT"]

### Merging ----------------------------------------------------------------
# Build dynamic column header names based on current number of experiments
n_all <- df_assessment_by_experiment |>
  filter(category %in% c("EPM","PCR","MTT")) |>
  summarise(n = n_distinct(EXP)) |>
  pull(n)

n_mtt <- df_assessment_by_experiment |>
  filter(category == "MTT") |>
  summarise(n = n_distinct(EXP)) |>
  pull(n)

n_pcr <- df_assessment_by_experiment |>
  filter(category == "PCR") |>
  summarise(n = n_distinct(EXP)) |>
  pull(n)

n_epm <- df_assessment_by_experiment |>
  filter(category == "EPM") |>
  summarise(n = n_distinct(EXP)) |>
  pull(n)

col_headers_tbl3 <- c(
  "By experiment",
  paste0("All (n=", n_all, ")"),
  paste0("MTT (n=", n_mtt, ")"),
  paste0("PCR (n=", n_pcr, ")"),
  paste0("EPM (n=", n_epm, ")")
)

tbl_3 <- tribble(
  ~`By experiment`, ~`All`, ~`MTT`, ~`PCR`, ~`EPM`, 

  # Dados da primeira parte (By experiment)
  "Effect size (original)",
  effect_size_original_all[[1]],
  effect_size_original_mtt,
  effect_size_original_pcr,
  effect_size_original_epm,
  "Effect size (replication)",
  effect_size_replication_all[[1]],
  effect_size_replication_mtt,
  effect_size_replication_pcr,
  effect_size_replication_epm,
  "Effect size ratio (from log difference)",
  effect_size_ratio_log_diff_all[[1]],
  effect_size_ratio_log_diff_mtt,
  effect_size_ratio_log_diff_pcr,
  effect_size_ratio_log_diff_epm,
  "Coefficient of variation (original)",
  coefficient_variation_original_all[[1]],
  coefficient_variation_original_mtt,
  coefficient_variation_original_pcr,
  coefficient_variation_original_epm,
  "Coefficient of variation (replication)",
  coefficient_variation_replication_all[[1]],
  coefficient_variation_replication_mtt,
  coefficient_variation_replication_altpcr,
  coefficient_variation_replication_epm,
  "Coefficient of variation ratio",
  coefficient_variation_ratio_all[[1]],
  coefficient_variation_ratio_mtt,
  coefficient_variation_ratio_pcr,
  coefficient_variation_ratio_epm,
  "Mean effect size difference (original vs. individual replications)",
  mean_effect_size_diff_original_vs_individual_all[[1]],
  mean_effect_size_diff_original_vs_individual_mtt,
  mean_effect_size_diff_original_vs_individual_pcr,
  mean_effect_size_diff_original_vs_individual_epm,
  "Mean effect size difference (between individual replications)",
  mean_effect_size_diff_between_individual_all[[1]],
  mean_effect_size_diff_between_individual_mtt,
  mean_effect_size_diff_between_individual_pcr,
  mean_effect_size_diff_between_individual_epm,
  "Ratio of mean absolute differences (paired)",
  ratio_mean_absolute_diff_all[[1]],
  ratio_mean_absolute_diff_mtt,
  ratio_mean_absolute_diff_pcr,
  ratio_mean_absolute_diff_epm,

  # "Sign error (% of significant)",
  # signal_error_significant_all[[1]],
  # signal_error_significant_mtt,
  # signal_error_significant_pcr,
  # signal_error_significant_epm,

  "Sign error (% of total)",
  signal_error_total_all[[1]],
  signal_error_total_mtt,
  signal_error_total_pcr,
  signal_error_total_epm,
  "Opposite sign (total)",
  opposite_signal_total_all[[1]],
  opposite_signal_mtt,
  opposite_signal_pcr,
  opposite_signal_epm
)

colnames(tbl_3) <- col_headers_tbl3

tbl_3 <- tbl_3 |>
  flextable() |>
  add_footer_lines("Median (Range)") |>
  flextable::footnote(
    i = c(5, 6), j = 4,
    value = as_paragraph("ALTPCR")
  ) |>
  bold(j = 2) |>
  align(
    align = "center", part = "all",
    j = 2:5
  ) |>
  set_table_properties(layout = "autofit") |>
  add_name("Table - Additional features of replications (Linear Scale).") |>
  bold(i = 2, part = "header")

### Saving ----
save_tbl(
  list(tbl_3),
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table 3.docx")
)

cat("\n### Table 3 generated! ###\n")

## Table S1 ----------------------------------------------------------------

df_registered <- read_excel("other-data/intermediate_steps.xlsx", 1) |>
  clean_names() |>
  var_labels(
    geographic_region = "Geographic Region",
    state = "State"
  )

df_included <- read_excel("other-data/intermediate_steps.xlsx", 2) |>
  clean_names() |>
  var_labels(
    geographic_region = "Geographic Region",
    state = "State"
  )

df_concluded <- read_excel("other-data/intermediate_steps.xlsx", 3) |>
  clean_names() |>
  var_labels(
    geographic_region = "Geographic Region",
    state = "State"
  )

tbl_registered <- df_registered |>
  select(-author) |>
  tbl_summary(
    sort = everything() ~ "frequency"
  )

tbl_included <- df_included |>
  select(-lab, -author) |>
  tbl_summary(
    sort = everything() ~ "frequency"
  )

tbl_concluded <- df_concluded |>
  select(geographic_region, state) |>
  tbl_summary(
    sort = everything() ~ "frequency"
  )


# By
df_all <- bind_rows(
  df_registered |> mutate(source = "Registered"),
  df_included |> mutate(source = "Included"),
  df_concluded |> mutate(source = "Concluded")
) |>
  sjlabelled::var_labels(
    geographic_region = "Geographic Region",
    state = "Region"
  )

tbl_s1 <- df_all |>
  select(geographic_region, state, source) |>
  mutate(source = factor(source,
    levels = c("Registered", "Included", "Concluded"),
    labels = c("Registered", "Included", "Contributed Data")
  )) |>
  tbl_summary(
    by = source,
    label = c(
      geographic_region ~ "Geographic Region",
      state ~ "Region"
    ),
    sort = everything() ~ "frequency"
  ) |>
  modify_header(label = "") |>
  add_name("Table - Geographic Distribution of Laboratories.") |>
  add_footer_lines("Numbers represent number (% of total) of labs in each region/state of Brazil for registered laboratories (those registering to participate in the project), included ones (those that received experiments to replicate) and those that concluded at least one experiment and contributed data to the project.")

### Saving ----
save_tbl(
  list(tbl_s1),
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S1.docx")
)

cat("\n### Table S1 generated! ###\n")


## Table S2 ----------------------------------------------------------------
theme_gtsummary_language(language = "en", big.mark = "")

df_authorship_list <- read_excel("other-data/authorship_list.xlsx") |> clean_names()

df_authorship_list <- df_authorship_list |>
  dplyr::distinct(author, .keep_all = T) |> 
  mutate(year_first_paper = as.numeric(year_first_paper)) |>
  mutate(highest_degree = factor(highest_degree, labels = c("Doctorate", "High School", "Undergraduate Degree", "Master's Degree"))) |>
  mutate(ruf2023 = as.numeric(ruf2023)) |>
  mutate(geographic_region = if_else(geographic_region == "?", NA, geographic_region)) |>
  mutate(geographic_region = factor(geographic_region, labels = c("Center-West", "Northeast", "Southeast", "South"))) |>
  mutate(today_date = ymd(today())) |>
  mutate(years_since_first_paper = year(today_date) - year_first_paper) |>
  mutate(years_since_first_paper = if_else(is.na(years_since_first_paper), 0, years_since_first_paper)) |>
  mutate(researchers_no_papers = if_else(n_papers == 0, "Yes", "No")) |>
  mutate(
    authorship_data_collection = if_else(authorship_data_collection == "NA", "No", authorship_data_collection)
  ) |>
  var_labels(
    highest_degree = "Highest Degree",
    highest_degree_year_conclusion = "Year of Conclusion of Highest Degree",
    year_first_paper = "Year of First Paper",
    years_since_first_paper = "Years Since First Paper",
    h_index = "H-Index",
    n_papers = "Number of Papers",
    researchers_no_papers = "Researchers Without Publications",
    n_citations = "Number of Citations",
    geographic_region = "Geographic Location",
    ruf2023 = "Folha Ranking - 2023"
  )

tbl_s2_only_data_collection <- df_authorship_list |>
  dplyr::filter(authorship_data_collection != "No" & authorship_protocol_design == "No") |>
  select(
    highest_degree,
    highest_degree_year_conclusion,
    years_since_first_paper,
    h_index,
    n_papers,
    researchers_no_papers,
    n_citations,
    ruf2023
  ) |>
  tbl_summary(
    digits = c(
      highest_degree_year_conclusion ~ 0,
      h_index ~ 0,
      years_since_first_paper ~ 0,
      all_categorical() ~ c(0, 1),
      n_papers ~ 0
    ),
    label = c(highest_degree ~ "Highest Degree"),
    missing_text = "Missing",
    missing = "no",
    statistic = c(all_continuous() ~ "{median} ({min}, {max})"),
    type = c(h_index ~ "continuous"),
    sort = list(
      highest_degree ~ "frequency"
    )
  )

tbl_s2_only_protocol <- df_authorship_list |>
  filter(authorship_data_collection == "No" & authorship_protocol_design != "No") |>
  select(
    highest_degree,
    highest_degree_year_conclusion,
    years_since_first_paper,
    h_index,
    n_papers,
    researchers_no_papers,
    n_citations,
    ruf2023
  ) |>
  tbl_summary(
    digits = c(
      highest_degree_year_conclusion ~ 0,
      h_index ~ 0,
      all_categorical() ~ c(0, 1),
      n_papers ~ 0
    ),
    label = c(highest_degree ~ "Highest Degree"),
    missing_text = "Missing",
    missing = "no",
    statistic = c(all_continuous() ~ "{median} ({min}, {max})"),
    type = c(h_index ~ "continuous"),
    sort = list( # geographic_region ~ "frequency",
      highest_degree ~ "frequency"
    )
  )

tbl_s2_both <- df_authorship_list |>
  filter(authorship_data_collection != "No" & authorship_protocol_design != "No") |>
  select(
    highest_degree,
    highest_degree_year_conclusion,
    years_since_first_paper,
    h_index,
    n_papers,
    researchers_no_papers,
    n_citations,
    ruf2023
  ) |>
  tbl_summary(
    digits = c(
      highest_degree_year_conclusion ~ 0,
      h_index ~ 0,
      all_categorical() ~ c(0, 1),
      n_papers ~ 0
    ),
    label = c(highest_degree ~ "Highest Degree"),
    missing_text = "Missing",
    missing = "no",
    statistic = c(all_continuous() ~ "{median} ({min}, {max})"),
    type = c(h_index ~ "continuous"),
    sort = list( # geographic_region ~ "frequency",
      highest_degree ~ "frequency"
    )
  )

tbl_s2_all <- df_authorship_list |>
  filter(!(authorship_data_collection == "No" & authorship_protocol_design == "No")) |>
  select(
    highest_degree,
    highest_degree_year_conclusion,
    years_since_first_paper,
    h_index,
    n_papers,
    researchers_no_papers,
    n_citations,
    ruf2023
  ) |>
  tbl_summary(
    digits = c(
      highest_degree_year_conclusion ~ 0,
      h_index ~ 0,
      all_categorical() ~ c(0, 1),
      n_papers ~ 0
    ),
    label = c(highest_degree ~ "Highest Degree"),
    missing_text = "Missing",
    missing = "no",
    statistic = c(all_continuous() ~ "{median} ({min}, {max})"),
    type = c(h_index ~ "continuous"),
    sort = list( # geographic_region ~ "frequency",
      highest_degree ~ "frequency"
    )
  )

tbl_s2 <- tbl_merge(list(tbl_s2_all, tbl_s2_only_protocol, tbl_s2_only_data_collection, tbl_s2_both),
  tab_spanner = c("**All**", "**Only Protocol**", "**Only Data Collection**", "**Both**")
) |>
  add_name("Table – Researchers Involved in Protocol Design and Data Collection.") |>
  set_table_properties(layout = "autofit")

### Saving ----
save_tbl(
  list(tbl_s2),
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S2.docx")
)

cat("\n### Table S2 generated! ###\n")

## Table S3 ----------------------------------------------------------------

df_predictor_data <- read_excel("other-data/predictor-data.xlsx") |>
  clean_names() |>
  select(!c(doi, figure_table))

df_extracted_data <- read_tsv("other-data/Original Experiments Statistical Summaries.tsv", show_col_types = FALSE) |>
  clean_names() |>
  select(!c(original_article_doi, selected_figure_table, method))

df_tbl_s3 <- df_predictor_data |> left_join(df_extracted_data, by = "exp", suffix = c("_predictor_data", "_extracted_data"))

df_predictor_data <- df_predictor_data |>
  mutate(
    n_control = str_replace_all(n_control, "NA", NA_character_),
    n_control = as.numeric(n_control),
    n_treated = str_replace_all(n_treated, "NA", NA_character_),
    n_treated = as.numeric(n_treated),
    normalized_impact_factor = str_replace_all(normalized_impact_factor, "-", NA_character_),
    normalized_impact_factor = as.numeric(normalized_impact_factor),
    position_in_the_folha_ranking_current = str_replace_all(position_in_the_folha_ranking_current, "not applicable", NA_character_),
    position_in_the_folha_ranking_current = as.numeric(position_in_the_folha_ranking_current),
    institution_type = stringr::str_to_sentence(institution_type),
    animal_model_species = fct_relabel(animal_model_species, ~ str_to_sentence(.)),
    biological_model = fct_relabel(biological_model, ~ str_to_sentence(.)),
    randomization = str_squish(randomization),
    randomization = if_else(randomization == "yes(any)", "yes (any)", randomization),
    randomization = if_else(randomization == "yes(experiment)", "yes (experiment)", randomization),
    randomization = fct_relabel(randomization, ~ str_to_sentence(.)),
    blinded_conduct_of_experiment = fct_relabel(blinded_conduct_of_experiment, ~ str_to_sentence(.)),
    blinding = fct_relabel(blinding, ~ str_to_sentence(.)),
    sample_size_calculation = fct_relabel(sample_size_calculation, ~ str_to_sentence(.)),
    inclusion_exclusion_criteria = fct_relabel(inclusion_exclusion_criteria, ~ str_to_sentence(.)),
    geographical_region = fct_recode(geographical_region,
      "Central-West" = "CO",
      "Central-West" = "CW",
      "Northeast" = "NE",
      "South" = "S",
      "Southeast" = "SE"
    ),
    `Hedges's g` = hedges_g(cohen_s_d, assumed_control_sample_size + assumed_treated_sample_size),
    cohen_s_d = abs(cohen_s_d),
    `Hedges's g` = abs(`Hedges's g`)
  ) |>
  var_labels(
    biological_model = "Biological model",
    animal_model_species = "Species of animal model",
    n_control = "Number of control group subjects",
    n_treated = "Number of treated group subjects",
    cohen_s_d = "Cohen's d effect size",
    p_value_calculated = "P-value of statistical test",
    original_p_value = "Original P-value",
    quality_of_reporting_index_from_gaps_in_original_protocol = "Quality of reporting index (gaps in original protocol)",
    randomization = "Randomization method",
    blinded_conduct_of_experiment = "Blinded conduct of the experiment",
    blinding = "Blinding method used",
    sample_size_calculation = "Sample size calculation",
    inclusion_exclusion_criteria = "Inclusion and exclusion criteria",
    year_of_online_publication = "Year publication",
    number_of_citations_first_2_years = "Number of citations (2 years)",
    normalized_impact_factor = "Area-normalized journal impact factor",
    last_author_s_academic_age = "Last author's academic age",
    number_of_papers_5_years_before_publishing_the_article = "Number of papers published in the 5 years prior to this article",
    position_in_the_folha_ranking_current = "Current position in Folha university ranking",
    geographical_region = "Geographical region",
    institution_type = "Type of institution"
  )

theme_gtsummary_language(language = "en", big.mark = "")

tbl_s3_1 <- df_tbl_s3 |>
  mutate(sample_size = (control_sample_size_for_d_fs + treated_sample_size_for_d_fs) / 2) |>
  mutate(hedgesg = hedges_g(cohen_s_d, control_sample_size_for_d_fs + treated_sample_size_for_d_fs)) |>
  mutate(hedgesg = abs(hedgesg)) |> 
  mutate(coef_var = pooled_sd / ((control_mean + treated_mean) / 2)) |>
  select(
    biological_model,
    sample_size,
    hedgesg,
    coef_var,
    animal_model_species,
    original_p_value_extracted_data,
    p_value_calculated
  ) |>
  mutate(
    original_p_value_extracted_data = str_squish(original_p_value_extracted_data),
    original_p_value_extracted_data = str_to_lower(original_p_value_extracted_data),
    original_p_value_extracted_data = case_match(
      original_p_value_extracted_data,
      "p<0.01" ~ "p < 0.01",
      "p < 0,01" ~ "p < 0.01",
      "<0.01" ~ "p < 0.01",
      "p=0.04" ~ "p < 0.05",
      "p = 0.04" ~ "p < 0.05",
      "p<0.001" ~ "p < 0.001",
      "<0,001" ~ "p < 0.001",
      "p<0.0001" ~ "p < 0.0001",
      "p<0,0001" ~ "p < 0.0001",
      "0,001" ~ "p < 0.001",
      "p<0.05" ~ "p < 0.05",
      "p<0,05" ~ "p < 0.05",
      "p < 0,05" ~ "p < 0.05",
      "0,05" ~ "p < 0.05",
      "<0.05" ~ "p < 0.05",
      "<0,05" ~ "p < 0.05",
      "p<0,01" ~ "p < 0.01",
      "p<o.001" ~ "p < 0.001",
      .default = original_p_value_extracted_data
    )
  ) |>
  mutate(
    biological_model = if_else(biological_model == "animal", animal_model_species, biological_model),
    biological_model = base::droplevels(as.factor(biological_model)),
    biological_model = factor(biological_model, labels = c("Cell line", "Mouse", "Primary Culture", "Rat"))
  ) |>
  tbl_summary(
    include = -animal_model_species,
    missing = "no",
    digits = c(p_value_calculated ~ function(x) format(x, digits = 6, scientific = TRUE)),
    type = all_continuous() ~ "continuous",
    sort = all_categorical() ~ "frequency",
    label = c(
      biological_model ~ "Biological Model",
      original_p_value_extracted_data ~ "Reported p value",
      hedgesg ~ "Effect size (Hedges's g)",
      sample_size ~ "Sample size (per group)",
      coef_var ~ "Coefficient of variation (original)",
      p_value_calculated ~ "Calculated p value"
    ),
    statistic = all_continuous() ~ c("{median} ({min}, {max})")
  ) |>
  bold_labels()

tbl_s3_2 <- df_predictor_data |>
  select(
    year_of_online_publication,
    number_of_citations_first_2_years,
    normalized_impact_factor,
    last_author_s_academic_age,
    number_of_papers_5_years_before_publishing_the_article,
    position_in_the_folha_ranking_current,
    geographical_region,
    institution_type
  ) |>
  tbl_summary(
    digits = c(year_of_online_publication ~ 0),
    sort = c(
      geographical_region ~ "frequency",
      institution_type ~ "frequency"
    ),
    label = c(
      geographical_region ~ "Geographical region",
      institution_type ~ "Type of institution"
    ),
    type = all_continuous() ~ "continuous",
    statistic = all_continuous() ~ c("{median} ({min}, {max})"),
    missing = "no"
  ) |>
  bold_labels()


tbl_s3_3 <- df_predictor_data |>
  filter(str_detect(exp, "EPM")) |>
  mutate(inclusion_exclusion_criteria = if_else(inclusion_exclusion_criteria == "Yes (any)" | inclusion_exclusion_criteria == "Yes (experiment)", "Yes", inclusion_exclusion_criteria)) |>
  mutate(sample_size_calculation = if_else(sample_size_calculation == "Yes (any)" | sample_size_calculation == "Yes (experiment)", "Yes", sample_size_calculation)) |>
  mutate(randomization = if_else(randomization == "Yes (any)" | randomization == "Yes (experiment)", "Yes", randomization)) |>
  mutate(blinding = if_else(blinding != "No/not reported", "Yes", blinding)) |>
  select(
    randomization,
    blinding,
    sample_size_calculation,
    inclusion_exclusion_criteria
  ) |>
  tbl_summary(
    missing = "no",
    digits = all_categorical() ~ c(0, 0),
    value = c(
      blinding ~ "Yes",
      randomization ~ "Yes",
      sample_size_calculation ~ "Yes",
      inclusion_exclusion_criteria ~ "Yes"
    ),
    label = c(
      blinding ~ "Blinded/automated outcome assessment",
      randomization ~ "Randomization",
      inclusion_exclusion_criteria ~ "Inclusion/exclusion criteria",
      sample_size_calculation ~ "Sample size calculation"
    ),
    statistic = all_continuous() ~ c("{median} ({min}, {max})")
  ) |>
  bold_labels()

tbl_s3 <- tbl_stack(
  list(
    tbl_s3_1,
    tbl_s3_3,
    tbl_s3_2
  ),
  group_header = c("Experiment (n = 60)", "Bias control (EPM experiments, n=20)", "Article, Corresponding author, Institution (n=60)")
) |>
  add_name("Table – Features of selected experiments, bias control, articles, authors and institutions") |>
  set_table_properties(layout = "autofit")

### Saving ----
save_tbl(
  list(tbl_s3),
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S3.docx")
)

cat("\n### Table S3 generated! ###\n")


## Table S4 ----------------------------------------------------------------

# From inclusions sets

df_inclusion_sets <- read_excel("other-data/inclusion_sets.xlsx")

all_rep_total <- df_inclusion_sets |>
  group_by(
    EXP
  ) |>
  summarise(
    number_replications = sum(analysis_all_exps_lab_units == "INCLUDE", na.rm = T)
  ) |>
  mutate(
    category = case_when(
      str_detect(EXP, "EPM") ~ "EPM",
      str_detect(EXP, "PCR") ~ "PCR",
      str_detect(EXP, "MTT") ~ "MTT"
    )
  ) |>
  group_by(number_replications) |>
  summarise(
    n()
  ) |>
  mutate(
    category = "Total"
  ) |>
  pivot_wider(
    names_from = "number_replications",
    values_from = "n()"
  )


validated_rep_total <- df_inclusion_sets |>
  group_by(
    EXP
  ) |>
  summarise(
    number_replications = sum(analysis_primary == "INCLUDE", na.rm = T)
  ) |>
  mutate(
    category = case_when(
      str_detect(EXP, "EPM") ~ "EPM",
      str_detect(EXP, "PCR") ~ "PCR",
      str_detect(EXP, "MTT") ~ "MTT"
    )
  ) |>
  group_by(number_replications) |>
  summarise(
    n()
  ) |>
  mutate(
    category = "Total"
  ) |>
  pivot_wider(
    names_from = "number_replications",
    values_from = "n()"
  )

tbl_s4 <- df_inclusion_sets |>
  group_by(
    EXP
  ) |>
  summarise(
    number_replications = sum(analysis_all_exps_lab_units == "INCLUDE", na.rm = T)
  ) |>
  mutate(
    category = case_when(
      str_detect(EXP, "EPM") ~ "EPM",
      str_detect(EXP, "PCR") ~ "PCR",
      str_detect(EXP, "MTT") ~ "MTT"
    )
  ) |>
  group_by(number_replications, category) |>
  summarise(
    n()
  ) |>
  arrange(category) |>
  pivot_wider(
    names_from = "number_replications",
    values_from = "n()"
  ) |>
  mutate(across(everything(), ~ replace_na(., 0))) |>
  add_row(all_rep_total) |>
  select(category, `3`, `2`, `1`, `0`) |>
  mutate(
    Total = c(20, 20, 20, 60),
  ) |>
  mutate(across(everything(), as.character)) |>
  add_row(
    tibble(
      category = "Validated replications",
      `3` = "",
      `2` = "",
      `1` = "",
      `0` = "",
      Total = ""
    )
  ) |>
  add_row(
    tibble(
      category = "# of replications",
      `3` = "3",
      `2` = "2",
      `1` = "1",
      `0` = "0",
      Total = "Total"
    )
  ) |>
  add_row(
    df_inclusion_sets |>
      group_by(
        EXP
      ) |>
      summarise(
        number_replications = sum(analysis_primary == "INCLUDE", na.rm = T)
      ) |>
      mutate(
        category = case_when(
          str_detect(EXP, "EPM") ~ "EPM",
          str_detect(EXP, "PCR") ~ "PCR",
          str_detect(EXP, "MTT") ~ "MTT"
        )
      ) |>
      group_by(number_replications, category) |>
      summarise(
        n()
      ) |>
      arrange(category) |>
      pivot_wider(
        names_from = "number_replications",
        values_from = "n()"
      ) |>
      mutate(across(everything(), ~ replace_na(., 0))) |>
      add_row(validated_rep_total) |>
      select(category, `3`, `2`, `1`, `0`) |>
      mutate(
        Total = c(20, 20, 20, 60),
      ) |>
      mutate(across(everything(), as.character))
  )

tbl_s4 <- rbind(
  tibble(
    category = "All replications",
    `3` = "",
    `2` = "",
    `1` = "",
    `0` = "",
    Total = ""
  ),
  tibble(
    category = "# of replications",
    `3` = "3",
    `2` = "2",
    `1` = "1",
    `0` = "0",
    Total = "Total"
  ),
  tbl_s4
)

tbl_s4 <- tbl_s4 |>
  flextable() |>
  hline(i = 1:2) |>
  hline(i = 6:7) |>
  bold(i = 1:2) |>
  bold(j = 1) |>
  bold(i = 7) |>
  delete_rows(i = 1, part = "header") |>
  add_footer_lines("Number of replications per experiment before (top) and after (bottom) the validation process. Cells display the number of experiments for each method with three, two, one or zero independent replications in these two samples.") |>
  add_name("Table - Number of Replications per Experiment") |>
  set_table_properties(layout = "autofit")

### Saving ----
save_tbl(
  list(tbl_s4),
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S4.docx")
)

cat("\n### Table S4 generated! ###\n")

## Table S5 ----------------------------------------------------------------

df_experiments_numbers_overview_wide <- read_tsv(paste0("output/", results_path, "/Experiment Numbers Overview - wide.tsv"), show_col_types = FALSE)

form1 <- read_excel("./other-data/self-assessment/Dificuldades Experimentais (Responses) - manual edit.xlsx")

df_inclusion_sets$`Validation - Decision`[inclusion_sets$`Validation - Decision` == "Low sample size"] <- "Excluded"

inclusion_sets_short <- df_inclusion_sets |>
  filter(UNIT == "BRI") |>
  select(c("EXP", "LAB", "Validation - Mean score", "Validation - Decision", "Validation - Excluded by"))

form1_val <- full_join(form1, inclusion_sets_short, by = c("LAB", "EXP")) |>
  mutate(method = str_extract(EXP, "[A-Za-z]+")) |>
  mutate(Validation_new = if_else(condition = `Validation - Decision` == "Included",
    true = `Validation - Decision`,
    false = `Validation - Excluded by`
  )) |>
  filter(!is.na(decision_agree)) |>
  filter(!is.na(Validation_new)) |>
  clean_names() |>
  select(exp, decision_agree, validation_new)

form1_val

# Build Table S5 data directly from inclusion_sets to avoid category mismatches
# Use inclusion_sets as the authoritative source for validation decisions

# Get experiment counts per validation category (only Done/BRI experiments)
experiment_counts <- df_inclusion_sets |>
  filter(UNIT == "BRI") |>
  filter(done == "Yes") |>
  mutate(
    # Use is.na(Validation - Excluded by) to determine Included, matching plots.R logic
    validation_decision = case_when(
      is.na(`Validation - Excluded by`) ~ "Included",
      TRUE ~ `Validation - Excluded by`
    )
  ) |>
  count(validation_decision, name = "n_experiments")

# Get lab agreement data: join form1 survey with inclusion_sets
# Note: form1 has decision_agree column from survey responses
form1_with_validation <- full_join(
  form1, 
  df_inclusion_sets |> 
    filter(UNIT == "BRI") |>
    filter(done == "Yes") |>
    select(EXP, LAB, `Validation - Decision`, `Validation - Excluded by`),
  by = c("LAB", "EXP")
) |>
  filter(!is.na(decision_agree)) |>
  mutate(
    # Use is.na(Validation - Excluded by) to determine Included, matching plots.R logic
    validation_decision = case_when(
      is.na(`Validation - Excluded by`) ~ "Included",
      TRUE ~ `Validation - Excluded by`
    )
  ) |>
  filter(!is.na(validation_decision))

# Count agreements and total responses per category
lab_agreement_summary <- form1_with_validation |>
  group_by(validation_decision) |>
  summarise(
    lab_agree = sum(decision_agree == "Sim", na.rm = TRUE),
    total_responses = n()
  )

# Join experiment counts with lab agreement
df_s5_data <- experiment_counts |>
  left_join(lab_agreement_summary, by = "validation_decision") |>
  rename(n = n_experiments) |>
  arrange(desc(n))


# Add total row (sum only Done experiments)
df_s5_data <- df_s5_data |>
  add_row(
    validation_decision = "Total",
    n = sum(df_s5_data$n, na.rm = TRUE),
    lab_agree = sum(df_s5_data$lab_agree, na.rm = TRUE),
    total_responses = sum(df_s5_data$total_responses, na.rm = TRUE)
  )

# Format output - percentage is calculated as agreement among respondents
df_experiments_numbers_overview_wide_join <- df_s5_data |>
  mutate(
    percent = lab_agree / total_responses * 100,
    `Lab Agreement` = paste0(lab_agree, " (", round(percent, 0), "%)")
  ) |>
  mutate(`Lab Agreement` = if_else(is.na(lab_agree) | is.na(percent), "", `Lab Agreement`)) |>
  select(validation_decision, n, `Lab Agreement`) |>
  rename(`Validation decision` = validation_decision, `# of cases` = n)


tbl_s5 <- df_experiments_numbers_overview_wide_join |>
  flextable() |>
  add_footer_lines("Table shows numbers for the different decisions reached at the validation process, as well as the rate of agreement of labs with each of them and a summary description of each category. For more information on the validation process, see https://osf.io/e3fjg.") |>
  add_name("Table - Agreement with validation decisions.") |>
  set_table_properties(layout = "autofit")

### Saving ----
save_tbl(
  list(tbl_s5),
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S5.docx")
)

cat("\n### Table S5 generated! ###\n")

## Table S7 -----------------------------------------------------
df_final_sub_avaliation <- read_excel("other-data/self-assessment/Avaliacao_subjetiva_final.xlsx", sheet = 4) |>
  clean_names() |>
  select(replicou, categoria1, categoria2, categoria3) |>
  mutate(categoria1 = str_to_sentence(categoria1)) |>
  mutate(categoria2 = str_to_sentence(categoria2)) |>
  mutate(categoria3 = str_to_sentence(categoria3)) |>
  mutate(categoria1 = case_match(
    categoria1,
    "Efeito oposto" ~ "Effect in the opposite direction",
    "Tamanho de efeito diferente" ~ "Different effect size",
    "Baixo tamanho amostral" ~ "Low sample size",
    "Alta variância" ~ "High variance",
    "Ausência de efeito" ~ "No effect",
    "Ausência de significância estatística" ~ "Not statistically significant",
    "Método não replicado" ~ "Failure to replicate method",
    "Resultado diferente" ~ "Different result (not specified)",
    "Tamanho de efeito comparável" ~ "Similar effect size",
    "Resultado significativo" ~ "Statistically significant",
    "Efeito na mesma direção" ~ "Effect in the same direction",
    "Resultado semelhante" ~ "Similar results (not specified)",
    "Sem resposta" ~ "No rationale"
  )) |>
  mutate(categoria2 = case_match(
    categoria2,
    "Efeito oposto" ~ "Effect in the opposite direction",
    "Tamanho de efeito diferente" ~ "Different effect size",
    "Baixo tamanho amostral" ~ "Low sample size",
    "Alta variância" ~ "High variance",
    "Ausência de efeito" ~ "No effect",
    "Ausência de significância estatística" ~ "Not statistically significant",
    "Método não replicado" ~ "Failure to replicate method",
    "Resultado diferente" ~ "Different result (not specified)",
    "Tamanho de efeito comparável" ~ "Similar effect size",
    "Resultado significativo" ~ "Statistically significant",
    "Efeito na mesma direção" ~ "Effect in the same direction",
    "Resultado semelhante" ~ "Similar results (not specified)",
    "Sem resposta" ~ "No rationale"
  )) |>
  mutate(categoria3 = case_match(
    categoria3,
    "Efeito oposto" ~ "Effect in the opposite direction",
    "Tamanho de efeito diferente" ~ "Different effect size",
    "Baixo tamanho amostral" ~ "Low sample size",
    "Alta variância" ~ "High variance",
    "Ausência de efeito" ~ "No effect",
    "Ausência de significância estatística" ~ "Not statistically significant",
    "Método não replicado" ~ "Failure to replicate method",
    "Resultado diferente" ~ "Different result (not specified)",
    "Tamanho de efeito comparável" ~ "Similar effect size",
    "Resultado significativo" ~ "Statistically significant",
    "Efeito na mesma direção" ~ "Effect in the same direction",
    "Resultado semelhante" ~ "Similar results (not specified)",
    "Sem resposta" ~ "No rationale"
  )) |>
  mutate(replicou = if_else(replicou == "Nao", "No", "Yes"))

n_successful <- df_final_sub_avaliation |>
  filter(replicou == "Yes") |>
  nrow()

n_unsuccessful <- df_final_sub_avaliation |>
  filter(replicou == "No") |>
  nrow()

tbl_s7 <- tibble(
  Category = paste0("Successful Replications (n = ", n_successful, ")"),
  `N (%)` = NA
) |>
  add_row(
    df_final_sub_avaliation |>
      filter(replicou == "Yes") |>
      select(-replicou) |>
      pivot_longer(cols = everything(), values_to = "category") |>
      count(category) |>
      pivot_wider(names_from = category, values_from = n, values_fill = 0) |>
      pivot_longer(cols = everything(), values_to = "category") |>
      filter(name != "NA") |>
      arrange(desc(category)) |>
      mutate(percent = if (n_successful > 0) category / n_successful * 100 else NA_real_) |>
      mutate(category = paste0(category, " (", round(percent, 1), "%)")) |>
      rename(
        Category = name,
        `N (%)` = category
      ) |>
      select(-percent)
  ) |>
  add_row(
    Category = paste0("Unsuccessful Replications (n = ", n_unsuccessful, ")")
  ) |>
  add_row(
    df_final_sub_avaliation |>
      filter(replicou == "No") |>
      select(-replicou) |>
      pivot_longer(cols = everything(), values_to = "category") |>
      count(category) |>
      pivot_wider(names_from = category, values_from = n, values_fill = 0) |>
      pivot_longer(cols = everything(), values_to = "category") |>
      filter(name != "NA") |>
      arrange(desc(category)) |>
      mutate(percent = if (n_unsuccessful > 0) category / n_unsuccessful * 100 else NA_real_) |>
      mutate(category = paste0(category, " (", round(percent, 1), "%)")) |>
      rename(
        Category = name,
        `N (%)` = category
      ) |>
      select(-percent) |> 
      slice(1:6, 9, 8, 7)
  ) |>
  flextable() |>
  add_name("Table - Rationale for subjective assessment of replication.") |>
  bold(i = 2, part = "header") |>
  bold(i = c(1, 7), part = "body") |>
  set_table_properties(layout = "autofit")

### Saving ----
save_tbl(
  list(tbl_s7),
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S7.docx")
)

cat("\n### Table S7 generated! ###\n")


## Table S14 -----------------------------------------------------
tbl_s14 <- df_by_experiment |>
  filter(Method == "ALL_PCR_MTT") |>
  filter(MA_Dist == "z") |>
  mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
  select(MetricShortName, Method, Value_Denominator, Inclusion_Set) |>
  pivot_wider(names_from = c(Inclusion_Set, Method), values_from = Value_Denominator) |>
  add_row(
    data.frame(
      MetricShortName = "By experiment"
    ),
    .before = 0
  ) |>
  add_row(
    df_by_replication |>
      filter(Method == "ALL_PCR_MTT") |>
      filter(MA_Dist == "z") |>
      mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
      select(MetricShortName, Method, Value_Denominator, Inclusion_Set) |>
      pivot_wider(names_from = c(Inclusion_Set, Method), values_from = Value_Denominator) |>
      add_row(
        data.frame(
          MetricShortName = "By replication"
        ),
        .before = 0
      )
  ) |>
  mutate(
    MetricShortName = case_when(
      MetricShortName == "FEMA is significant and has same signal as original" ~ "Same-sign significance (p<0.05)",
      MetricShortName == "Original estimate within PI of REMA" ~ "Original in replication’s 95% PI",
      MetricShortName == "REMA estimate within CI of Original" ~ "Replication in original 95% CI",
      MetricShortName == "Subjective assessment majority vote (with ties as success)" ~ "≥50% subjectively replicated",
      MetricShortName == "t-test majority vote (with ties as success)" ~ "≥50% replications significant",
      MetricShortName == "Subjective assessment, individual" ~ "Subjectively replicated",
      MetricShortName == "Replication estimate within CI of Original" ~ "Replication in original 95% CI",
      MetricShortName == "Replication is significant and has same sign as original" ~ "Same-sign significance (p<0.05)",
      TRUE ~ MetricShortName
    )
  ) |>
  slice(1:7, 8, 9, 10)

ordered_cols_tbl_s14 <- tibble(
  colname = colnames(tbl_s14),
  priority = case_when(
    str_starts(colname, "MetricShortName") ~ 1,
    str_starts(colname, "primary") ~ 2,
    str_starts(colname, "included_by_lab") ~ 3,
    str_starts(colname, "all_exps_lab_units") ~ 4,
    str_starts(colname, "at_least_2_reps") ~ 5,
    str_starts(colname, "only_3_reps") ~ 6,
    str_starts(colname, "only_80_power_a_posteriori") ~ 7,
    TRUE ~ 8
  )
) |>
  arrange(priority) |>
  pull(colname)

tbl_s14 <- tbl_s14 |>
  select(any_of(ordered_cols_tbl_s14)) |>
  dplyr::rename(any_of(c(
    `Primary` = "primary_ALL_PCR_MTT",
    `Lab's choice` = "included_by_lab_ALL_PCR_MTT",
    `All Exps` = "all_exps_lab_units_ALL_PCR_MTT",
    `≥ 2 copies` = "at_least_2_reps_ALL_PCR_MTT",
    `3 copies` = "only_3_reps_ALL_PCR_MTT",
    `≥80% power` = "only_80_power_a_posteriori_Z_ALL_PCR_MTT"
  )))

footer_text_tbl_s14 <- "Replication rates for the primary and secondary analyses. Results are the same as in Table 2, but using a z distribution for statistical tests and prediction intervals derived from meta-analyses (as the distribution was not specified in the primary analysis). Same-sign significance is based on a fixed meta-analysis estimate, while effect size comparisons are based on random-effects meta-analysis. Subsets for secondary analyses include all experiments judged valid by the replicating lab (Lab’s Choice), all concluded experiments (using the lab unit as defined by the replicating lab), only experiments with at least 2 and 3 copies, and only experiments with ≥ 80% a posteriori power (using a z distribution, which leads to a different subset than that included in Table 2 and Table S8. PI, prediction interval, CI, confidence interval. For more information on replication criteria, see https://osf.io/9rnuj."

tbl_s14 <- tbl_s14 |>
  slice(-1) |>
  flextable() |>
  bold(j = 2) |>
  bold(i = 6) |>
  hline(i = 5:6) |>
  set_header_labels(MetricShortName = "By experiment") |>
  add_footer_lines(value = as_paragraph(footer_text_tbl_s14)) |>
  add_header_lines(values = paste0("Table S14 - Replication rates for different analysis sets using the z distribution")) |>
  bold(i = 2, part = "header") |>
  set_table_properties(layout = "autofit")

### Saving -----

save_tbl(
  tbl_s14,
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S14.docx")
)
cat("\n### Table S14 generated! ###\n")

## Table S16 -----------------------------------------------------
tbl_s16 <- df_by_experiment |>
  filter(Method == "ALL_PCR_MTT") |>
  filter(MA_Dist == "knha") |>
  mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
  select(MetricShortName, Method, Value_Denominator, Inclusion_Set) |>
  pivot_wider(names_from = c(Inclusion_Set, Method), values_from = Value_Denominator) |>
  add_row(
    data.frame(
      MetricShortName = "By experiment"
    ),
    .before = 0
  ) |>
  add_row(
    df_by_replication |>
      filter(Method == "ALL_PCR_MTT") |>
      filter(MA_Dist == "knha") |>
      mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
      select(MetricShortName, Method, Value_Denominator, Inclusion_Set) |>
      pivot_wider(names_from = c(Inclusion_Set, Method), values_from = Value_Denominator) |>
      add_row(
        data.frame(
          MetricShortName = "By replication"
        ),
        .before = 0
      )
  ) |>
  mutate(
    MetricShortName = case_when(
      MetricShortName == "FEMA is significant and has same signal as original" ~ "Same-sign significance (p<0.05)",
      MetricShortName == "Original estimate within PI of REMA" ~ "Original in replication’s 95% PI",
      MetricShortName == "REMA estimate within CI of Original" ~ "Replication in original 95% CI",
      MetricShortName == "Subjective assessment majority vote (with ties as success)" ~ "≥50% subjectively replicated",
      MetricShortName == "t-test majority vote (with ties as success)" ~ "≥50% replications significant",
      MetricShortName == "Subjective assessment, individual" ~ "Subjectively replicated",
      MetricShortName == "Replication estimate within CI of Original" ~ "Replication in original 95% CI",
      MetricShortName == "Replication is significant and has same sign as original" ~ "Same-sign significance (p<0.05)",
      TRUE ~ MetricShortName
    )
  ) |>
  slice(1:7, 8, 9, 10)

ordered_cols_tbl_s16 <- tibble(
  colname = colnames(tbl_s16),
  priority = case_when(
    str_starts(colname, "MetricShortName") ~ 1,
    str_starts(colname, "primary") ~ 2,
    str_starts(colname, "included_by_lab") ~ 3,
    str_starts(colname, "all_exps_lab_units") ~ 4,
    str_starts(colname, "at_least_2_reps") ~ 5,
    str_starts(colname, "only_3_reps") ~ 6,
    str_starts(colname, "only_80_power_a_posteriori") ~ 7,
    TRUE ~ 8
  )
) |>
  arrange(priority) |>
  pull(colname)

tbl_s16 <- tbl_s16 |>
  select(any_of(ordered_cols_tbl_s16)) |>
  dplyr::rename(any_of(c(
    `Primary` = "primary_ALL_PCR_MTT",
    `Lab's choice` = "included_by_lab_ALL_PCR_MTT",
    `All Exps` = "all_exps_lab_units_ALL_PCR_MTT",
    `≥ 2 copies` = "at_least_2_reps_ALL_PCR_MTT",
    `3 copies` = "only_3_reps_ALL_PCR_MTT",
    `≥80% power` = "only_80_power_a_posteriori_KNHA_ALL_PCR_MTT"
  )))

footer_text_tbl_s16 <- "Replication rates for the primary and secondary analyses. Results are the same as in Table 2, but using a knha distribution for statistical tests and prediction intervals derived from meta-analyses (as the distribution was not specified in the primary analysis). Same-sign significance is based on a fixed meta-analysis estimate, while effect size comparisons are based on random-effects meta-analysis. Subsets for secondary analyses include all experiments judged valid by the replicating lab (Lab’s Choice), all concluded experiments (using the lab unit as defined by the replicating lab), only experiments with at least 2 and 3 copies, and only experiments with ≥ 80% a posteriori power (using a z distribution, which leads to a different subset than that included in Table 2 and Table S8. PI, prediction interval, CI, confidence interval. For more information on replication criteria, see https://osf.io/9rnuj."

tbl_s16 <- tbl_s16 |>
  slice(-1) |>
  flextable() |>
  bold(j = 2) |>
  bold(i = 6) |>
  hline(i = 5:6) |>
  set_header_labels(MetricShortName = "By experiment") |>
  add_footer_lines(value = as_paragraph(footer_text_tbl_s16)) |>
  add_header_lines(values = paste0("Table S16 - Replication rates for different analysis sets using the knha distribution")) |>
  bold(i = 2, part = "header") |>
  set_table_properties(layout = "autofit")

### Saving -----

save_tbl(
  tbl_s16,
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S16.docx")
)
cat("\n### Table S16 generated! ###\n")

## Table S18 (bigexp equivalent of Table 2) -----------------------------------------------------
# Build experiment-level data from bigexp
tbl_s18_exp <- df_by_experiment |>
  filter(Method == "ALL_PCR_MTT") |>
  filter(MA_Dist == "bigexp") |>
  mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
  select(MetricShortName, Method, Value_Denominator, Inclusion_Set) |>
  pivot_wider(names_from = c(Inclusion_Set, Method), values_from = Value_Denominator)

# Build replication-level data from bigexp
tbl_s18_rep <- df_by_replication |>
  filter(Method == "ALL_PCR_MTT") |>
  filter(MA_Dist == "bigexp") |>
  mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
  select(MetricShortName, Method, Value_Denominator, Inclusion_Set) |>
  pivot_wider(names_from = c(Inclusion_Set, Method), values_from = Value_Denominator)

# Get replication-level data from t distribution for missing metrics
tbl_s18_rep_t <- df_by_replication |>
  filter(Method == "ALL_PCR_MTT") |>
  filter(MA_Dist == "t") |>
  mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
  select(MetricShortName, Method, Value_Denominator, Inclusion_Set) |>
  pivot_wider(names_from = c(Inclusion_Set, Method), values_from = Value_Denominator)

# Define expected metrics and their sources
expected_exp_metrics <- c(
  "Original estimate within PI of REMA",
  "REMA estimate within CI of Original",
  "FEMA is significant and has same signal as original",
  "t-test majority vote (with ties as success)",
  "Subjective assessment majority vote (with ties as success)"
)

expected_rep_metrics <- c(
  "Replication estimate within CI of Original",
  "Replication is significant and has same sign as original",
  "Subjective assessment, individual"
)

# Get available columns (analysis sets)
available_cols_exp <- setdiff(colnames(tbl_s18_exp), "MetricShortName")
available_cols_rep <- setdiff(colnames(tbl_s18_rep), "MetricShortName")

# Add missing experiment-level metrics as NA
for (metric in expected_exp_metrics) {
  if (!(metric %in% tbl_s18_exp$MetricShortName)) {
    new_row <- tibble(MetricShortName = metric)
    for (col in available_cols_exp) new_row[[col]] <- "NA"
    tbl_s18_exp <- bind_rows(tbl_s18_exp, new_row)
  }
}

# Add missing replication-level metrics from t distribution
for (metric in expected_rep_metrics) {
  if (!(metric %in% tbl_s18_rep$MetricShortName)) {
    if (metric %in% tbl_s18_rep_t$MetricShortName) {
      new_row <- tbl_s18_rep_t |> 
        filter(MetricShortName == metric) |> 
        select(MetricShortName, any_of(available_cols_rep))
      tbl_s18_rep <- bind_rows(tbl_s18_rep, new_row)
    } else {
      new_row <- tibble(MetricShortName = metric)
      for (col in available_cols_rep) new_row[[col]] <- "NA"
      tbl_s18_rep <- bind_rows(tbl_s18_rep, new_row)
    }
  }
}

# Reorder according to expected order
tbl_s18_exp <- tbl_s18_exp |>
  mutate(order = match(MetricShortName, expected_exp_metrics)) |>
  arrange(order) |>
  select(-order)

tbl_s18_rep <- tbl_s18_rep |>
  mutate(order = match(MetricShortName, expected_rep_metrics)) |>
  arrange(order) |>
  select(-order)

# Combine experiment and replication level data
tbl_s18 <- tbl_s18_exp |>
  add_row(
    data.frame(MetricShortName = "By experiment"),
    .before = 0
  ) |>
  add_row(
    tbl_s18_rep |>
      add_row(
        data.frame(MetricShortName = "By replication"),
        .before = 0
      )
  ) |>
  mutate(
    MetricShortName = case_when(
      MetricShortName == "FEMA is significant and has same signal as original" ~ "Same-sign significance (p<0.05)",
      MetricShortName == "Original estimate within PI of REMA" ~ "Original in replication's 95% PI",
      MetricShortName == "REMA estimate within CI of Original" ~ "Replication in original 95% CI",
      MetricShortName == "Subjective assessment majority vote (with ties as success)" ~ "≥50% subjectively replicated",
      MetricShortName == "t-test majority vote (with ties as success)" ~ "≥50% replications significant",
      MetricShortName == "Subjective assessment, individual" ~ "Subjectively replicated",
      MetricShortName == "Replication estimate within CI of Original" ~ "Replication in original 95% CI",
      MetricShortName == "Replication is significant and has same sign as original" ~ "Same-sign significance (p<0.05)",
      TRUE ~ MetricShortName
    )
  )

# Order columns - keep only bigexp power column
ordered_cols_tbl_s18 <- tibble(
  colname = colnames(tbl_s18),
  priority = case_when(
    str_starts(colname, "MetricShortName") ~ 1,
    str_starts(colname, "primary") ~ 2,
    str_starts(colname, "included_by_lab") ~ 3,
    str_starts(colname, "all_exps_lab_units") ~ 4,
    str_starts(colname, "at_least_2_reps") ~ 5,
    str_starts(colname, "only_3_reps") ~ 6,
    colname == "only_80_power_a_posteriori_bigexp_ALL_PCR_MTT" ~ 7,  # Keep only bigexp power (lowercase)
    str_starts(colname, "only_80_power") ~ 99,  # Exclude all other power columns
    TRUE ~ 8
  )
) |>
  filter(priority < 99) |>
  arrange(priority) |>
  pull(colname)

tbl_s18 <- tbl_s18 |>
  select(all_of(ordered_cols_tbl_s18))

# Rename columns to human-readable names (only if they exist)
if ("primary_ALL_PCR_MTT" %in% colnames(tbl_s18)) {
  tbl_s18 <- tbl_s18 |> dplyr::rename(`Primary` = primary_ALL_PCR_MTT)
}
if ("included_by_lab_ALL_PCR_MTT" %in% colnames(tbl_s18)) {
  tbl_s18 <- tbl_s18 |> dplyr::rename(`Lab's choice` = included_by_lab_ALL_PCR_MTT)
}
if ("all_exps_lab_units_ALL_PCR_MTT" %in% colnames(tbl_s18)) {
  tbl_s18 <- tbl_s18 |> dplyr::rename(`All Exps` = all_exps_lab_units_ALL_PCR_MTT)
}
if ("at_least_2_reps_ALL_PCR_MTT" %in% colnames(tbl_s18)) {
  tbl_s18 <- tbl_s18 |> dplyr::rename(`≥ 2 copies` = at_least_2_reps_ALL_PCR_MTT)
}
if ("only_3_reps_ALL_PCR_MTT" %in% colnames(tbl_s18)) {
  tbl_s18 <- tbl_s18 |> dplyr::rename(`3 copies` = only_3_reps_ALL_PCR_MTT)
}
# Rename bigexp power column with shorter title (lowercase bigexp)
if ("only_80_power_a_posteriori_bigexp_ALL_PCR_MTT" %in% colnames(tbl_s18)) {
  tbl_s18 <- tbl_s18 |> dplyr::rename(`>80% power` = only_80_power_a_posteriori_bigexp_ALL_PCR_MTT)
}

footer_text_tbl_s18 <- "Replication rates for the primary and secondary analyses using the bigexp distribution (experiment-level analysis). Effect size comparisons are based on random-effects meta-analysis, while same-sign significance is based on a fixed meta-analysis estimate. Metrics marked 'NA' are not available in the bigexp analysis. Replication-level metrics use values from the t distribution analysis (Table 2) for commensurability. Subsets for secondary analyses include all experiments judged valid by the replicating lab (Lab's Choice), all concluded experiments (All Exps), both using the experimental unit as defined by the lab, and only experiments with at least 2 and 3 replications. PI, prediction interval; CI, confidence interval. For more information on replication criteria, see https://osf.io/9rnuj."

tbl_s18 <- tbl_s18 |>
  slice(-1) |>
  flextable() |>
  bold(j = 2) |>
  bold(i = 6) |>
  hline(i = 5:6) |>
  set_header_labels(MetricShortName = "By experiment") |>
  add_footer_lines(value = as_paragraph(footer_text_tbl_s18)) |>
  add_header_lines(values = paste0("Table S18 - Replication rates for different analysis sets using the bigexp distribution")) |>
  bold(i = 2, part = "header") |>
  set_table_properties(layout = "autofit")

### Saving -----

save_tbl(
  tbl_s18,
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S18.docx")
)
cat("\n### Table S18 generated! ###\n")


## Table S20 -----------------------------------------------------------------
tbl_s20 <- df_by_experiment |>
  filter(Method == "PCR" | Method == "ALTPCR") |>
  filter(MA_Dist == "t") |>
  filter(Inclusion_Set == "all_exps_lab_units" | Inclusion_Set == "primary") |>
  mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
  select(MetricShortName, Method, Value_Denominator, Inclusion_Set, MA_Dist) |>
  pivot_wider(names_from = c(Inclusion_Set, Method, MA_Dist), values_from = Value_Denominator) |>
  add_row(
    data.frame(
      MetricShortName = "By experiment"
    ),
    .before = 0
  ) |>
  add_row(
    df_by_replication |>
      filter(Method == "PCR" | Method == "ALTPCR") |>
      filter(MA_Dist == "t") |>
      filter(Inclusion_Set == "all_exps_lab_units" | Inclusion_Set == "primary") |>
      mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
      select(MetricShortName, Method, Value_Denominator, Inclusion_Set, MA_Dist) |>
      pivot_wider(names_from = c(Inclusion_Set, Method, MA_Dist), values_from = Value_Denominator) |>
      add_row(
        data.frame(
          MetricShortName = "By replication"
        ),
        .before = 0
      )
  ) |>
  mutate(
    MetricShortName = case_when(
      MetricShortName == "FEMA is significant and has same signal as original" ~ "Same-sign significance (p<0.05)",
      MetricShortName == "Original estimate within PI of REMA" ~ "Original in replication’s 95% PI",
      MetricShortName == "REMA estimate within CI of Original" ~ "Replication in original 95% CI",
      MetricShortName == "Subjective assessment majority vote (with ties as success)" ~ "≥50% subjectively replicated",
      MetricShortName == "t-test majority vote (with ties as success)" ~ "≥50% replications significant",
      MetricShortName == "Subjective assessment, individual" ~ "Subjectively replicated",
      MetricShortName == "Replication estimate within CI of Original" ~ "Replication in original 95% CI",
      MetricShortName == "Replication is significant and has same sign as original" ~ "Same-sign significance (p<0.05)",
      TRUE ~ MetricShortName
    )
  ) |>
  slice(1:7, 8, 9, 10)

ordered_cols_tbl_s20 <- tibble(
  colname = colnames(tbl_s20),
  priority = case_when(
    str_starts(colname, "MetricShortName") ~ 1,
    str_starts(colname, "primary") ~ 2,
    str_starts(colname, "all_exps_lab_units") ~ 3,
    TRUE ~ 4
  )
) |>
  arrange(priority) |>
  pull(colname)

tbl_s20 <- tbl_s20 |>
  select(any_of(ordered_cols_tbl_s20)) |>
  dplyr::rename(any_of(c(
    `Primary (Log values)` = "primary_PCR_t",
    `Primary (Linear values)` = "primary_ALTPCR_t",
    `All experiments (Log values)` = "all_exps_lab_units_PCR_t",
    `All experiments (Linear values)` = "all_exps_lab_units_ALTPCR_t"
  )))

footer_text_tbl_s20 <- "Log values columns show results of analysis aggregating ΔCt values for real-time PCR or relative expression in log2 scale for conventional PCR, as in the primary analysis. Linear values columns show results aggregating linearized values (usually 2-ΔΔCt) for real-time PCR and relative expression in linear scale for conventional PCR, as originally planned in the protocols. Same-sign significance is based on a fixed meta-analysis estimate, while effect size comparisons are based on random-effects meta-analysis. All statistical tests use the t distribution. PI, prediction interval, CI, confidence interval. For more information on replication criteria, see the registered replication protocol at https://osf.io/9rnuj."

tbl_s20 <- tbl_s20 |>
  slice(-1) |>
  flextable() |>
  bold(j = 2) |>
  bold(i = 6) |>
  hline(i = 5:6) |>
  set_header_labels(MetricShortName = "By experiment") |>
  add_footer_lines(value = as_paragraph(footer_text_tbl_s20)) |>
  add_header_lines(values = paste0("Table S20 - Replication rates for PCR experiments using different aggregation methods. ")) |>
  bold(i = 2, part = "header") |>
  set_table_properties(layout = "autofit")

### Saving ----
save_tbl(
  tbl_s20,
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S20.docx")
)
cat("\n### Table S20 generated! ###\n")

## Table S19 -----------------------------------------------------------------
tbl_s19 <- df_by_experiment |>
  filter(Method == "MTT" | Method == "ALTMTT") |>
  filter(MA_Dist == "t") |>
  filter(Inclusion_Set == "all_exps_lab_units" | Inclusion_Set == "primary") |>
  mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
  select(MetricShortName, Method, Value_Denominator, Inclusion_Set, MA_Dist) |>
  pivot_wider(names_from = c(Inclusion_Set, Method, MA_Dist), values_from = Value_Denominator) |>
  add_row(
    data.frame(
      MetricShortName = "By experiment"
    ),
    .before = 0
  ) |>
  add_row(
    df_by_replication |>
      filter(Method == "MTT" | Method == "ALTMTT") |>
      filter(MA_Dist == "t") |>
      filter(Inclusion_Set == "all_exps_lab_units" | Inclusion_Set == "primary") |>
      mutate(Value_Denominator = paste0(round(successful), "/", N, " (", Value, ")")) |>
      select(MetricShortName, Method, Value_Denominator, Inclusion_Set, MA_Dist) |>
      pivot_wider(names_from = c(Inclusion_Set, Method, MA_Dist), values_from = Value_Denominator) |>
      add_row(
        data.frame(
          MetricShortName = "By replication"
        ),
        .before = 0
      )
  ) |>
  mutate(
    MetricShortName = case_when(
      MetricShortName == "FEMA is significant and has same signal as original" ~ "Same-sign significance (p<0.05)",
      MetricShortName == "Original estimate within PI of REMA" ~ "Original in replication’s 95% PI",
      MetricShortName == "REMA estimate within CI of Original" ~ "Replication in original 95% CI",
      MetricShortName == "Subjective assessment majority vote (with ties as success)" ~ "≥50% subjectively replicated",
      MetricShortName == "t-test majority vote (with ties as success)" ~ "≥50% replications significant",
      MetricShortName == "Subjective assessment, individual" ~ "Subjectively replicated",
      MetricShortName == "Replication estimate within CI of Original" ~ "Replication in original 95% CI",
      MetricShortName == "Replication is significant and has same sign as original" ~ "Same-sign significance (p<0.05)",
      TRUE ~ MetricShortName
    )
  ) |>
  slice(1:7, 8, 9, 10)

ordered_cols_tbl_s19 <- tibble(
  colname = colnames(tbl_s19),
  priority = case_when(
    str_starts(colname, "MetricShortName") ~ 1,
    str_starts(colname, "primary") ~ 2,
    str_starts(colname, "all_exps_lab_units") ~ 3,
    TRUE ~ 4
  )
) |>
  arrange(priority) |>
  pull(colname)

tbl_s19 <- tbl_s19 |>
  select(any_of(ordered_cols_tbl_s19)) |>
  dplyr::rename(any_of(c(
    `Primary (Paired)` = "primary_MTT_t",
    `Primary (Original)` = "primary_ALTMTT_t",
    `All experiments (Paired)` = "all_exps_lab_units_MTT_t",
    `All experiments (Original)` = "all_exps_lab_units_ALTMTT_t"
  )))

### Coefficient of variation row for MTT replications ----
# Load replication-level assessment data for primary and all experiments
df_rep_primary_mtt <- read_tsv(paste0("output/", results_path, "/primary t/Replication Assessment by Replication.tsv"), show_col_types = FALSE) |>
  filter(str_detect(EXP, "^MTT")) |>
  mutate(Method = "MTT")

df_rep_primary_altmtt <- read_tsv(paste0("output/", results_path, "/primary t/Replication Assessment by Replication.tsv"), show_col_types = FALSE) |>
  filter(str_detect(EXP, "^ALTMTT")) |>
  mutate(Method = "ALTMTT")

allexp_t_path <- paste0("output/", results_path, "/all_exps_lab_units t/Replication Assessment by Replication.tsv")
if (file.exists(allexp_t_path)) {
  df_rep_allexp_mtt <- read_tsv(allexp_t_path, show_col_types = FALSE) |>
    filter(str_detect(EXP, "^MTT")) |>
    mutate(Method = "MTT")

  df_rep_allexp_altmtt <- read_tsv(allexp_t_path, show_col_types = FALSE) |>
    filter(str_detect(EXP, "^ALTMTT")) |>
    mutate(Method = "ALTMTT")
} else {
  df_rep_allexp_mtt <- df_rep_primary_mtt[0, ]
  df_rep_allexp_altmtt <- df_rep_primary_altmtt[0, ]
}

# Calculate median (range) for each combination
cv_primary_paired <- df_rep_primary_mtt |>
  summarise(
    cv_summary = paste0(
      sprintf("%.3f", median(replication_cv, na.rm = TRUE)),
      " (",
      sprintf("%.3f", min(replication_cv, na.rm = TRUE)),
      ", ",
      sprintf("%.3f", max(replication_cv, na.rm = TRUE)),
      ")"
    )
  ) |>
  pull(cv_summary)

cv_primary_original <- df_rep_primary_altmtt |>
  summarise(
    cv_summary = paste0(
      sprintf("%.3f", median(replication_cv, na.rm = TRUE)),
      " (",
      sprintf("%.3f", min(replication_cv, na.rm = TRUE)),
      ", ",
      sprintf("%.3f", max(replication_cv, na.rm = TRUE)),
      ")"
    )
  ) |>
  pull(cv_summary)

cv_allexp_paired <- df_rep_allexp_mtt |>
  summarise(
    cv_summary = paste0(
      sprintf("%.3f", median(replication_cv, na.rm = TRUE)),
      " (",
      sprintf("%.3f", min(replication_cv, na.rm = TRUE)),
      ", ",
      sprintf("%.3f", max(replication_cv, na.rm = TRUE)),
      ")"
    )
  ) |>
  pull(cv_summary)

cv_allexp_original <- df_rep_allexp_altmtt |>
  summarise(
    cv_summary = paste0(
      sprintf("%.3f", median(replication_cv, na.rm = TRUE)),
      " (",
      sprintf("%.3f", min(replication_cv, na.rm = TRUE)),
      ", ",
      sprintf("%.3f", max(replication_cv, na.rm = TRUE)),
      ")"
    )
  ) |>
  pull(cv_summary)

# Add CV row to the table data
cv_row <- tibble(
  MetricShortName = "Coefficient of variation (replications)",
  `Primary (Paired)` = cv_primary_paired,
  `Primary (Original)` = cv_primary_original,
  `All experiments (Paired)` = cv_allexp_paired,
  `All experiments (Original)` = cv_allexp_original
)

tbl_s19 <- tbl_s19 |>
  add_row(cv_row)

footer_text_tbl_s19 <- "Paired columns show results of MTT analyzed with paired comparisons (MTT), while Original columns show the alternative non-paired analysis (ALTMTT). Same-sign significance is based on a fixed meta-analysis estimate, while effect size comparisons are based on random-effects meta-analysis. All statistical tests use the t distribution. PI, prediction interval; CI, confidence interval. Coefficient of variation is reported as median (range). For more information on replication criteria, see https://osf.io/9rnuj."

tbl_s19 <- tbl_s19 |>
  slice(-1) |>
  flextable() |>
  bold(j = 2) |>
  bold(i = 6) |>
  hline(i = 5:6) |>
  set_header_labels(MetricShortName = "By experiment") |>
  add_footer_lines(value = as_paragraph(footer_text_tbl_s19)) |>
  add_header_lines(values = paste0("Table S19 - Replication rates for MTT experiments using paired vs original (non-paired) analyses. ")) |>
  bold(i = 2, part = "header") |>
  set_table_properties(layout = "autofit")

### Saving ----
save_tbl(
  tbl_s19,
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S19.docx")
)
cat("\n### Table S19 generated! ###\n")


## Table S24 ---------------------------------------------------------------
df_coord_difficulties <- read_excel("other-data/Coordinating team assessment of difficulties.xlsx") |>
  clean_names() |>
  select(categoria, dificuldade, soma_um)

tbl_s24 <- df_coord_difficulties |>
  group_by(categoria) |>
  reframe(
    dificuldade = dificuldade,
    soma = round(soma_um / 6, 1),
    total_soma = sum(soma_um)
  ) |>
  arrange(desc(total_soma)) |>
  select(-total_soma) |>
  slice_max(soma, n = 3, by = categoria) |>
  mutate(
    categoria = case_when(
      categoria == "Dificuldades relacionadas à organização dos laboratórios" ~ "Large-scale project management",
      categoria == "Dificuldades relacionadas a práticas e formação de pesquisadores" ~ "Researcher practices and expertise",
      categoria == "Dificuldades relacionadas à literatura" ~ "Scientific literature issues",
      categoria == "Dificuldades intrínsecas" ~ "Intrinsic experimental difficulties",
      categoria == "Dificuldades estruturais" ~ "Logistics and infrastructure"
    )
  ) |>
  mutate(dificuldade = case_when(
    dificuldade == "Comunicação deficiente com a equipe coordenadora" ~ "Communication difficulties between the coordinating team and participating labs.",
    dificuldade == "Não antecipação de perdas experimentais no planejamento da análise" ~ "Failure to anticipate loss of labs and experiments when planning analyses.",
    dificuldade == "Pouca atenção a orientações/manuais escritos" ~ "Lack of attention to written guidelines/manuals by participant labs.",
    dificuldade == "Ausência de terminologia comum para definir unidade experimental" ~ "Absence of standardized terminologies to define experimental units and other aspects.",
    dificuldade == "Pouca atenção a protocolos na elaboração" ~ "Lack of attention to detail when preparing protocols.",
    dificuldade == "Pouca atenção no gerenciamento de dados" ~ "Poor data management.",
    dificuldade == "Metodologia mal descrita" ~ "Poorly described methods in original articles.",
    dificuldade == "Unidade experimental mal descrita" ~ "Poorly described experimental unit in original articles.",
    dificuldade == "Variabilidade possivelmente subestimada" ~ "Likely underestimation of biological variability in original articles.",
    dificuldade == "Modelos com funcionamento diferente do previsto" ~ "Models behaving differently than expected.",
    dificuldade == "Dilema entre \"replicar metodologia\" vs. \"replicar modelo\"" ~ "Tension between replicating method as described and replicating model features.",
    dificuldade == "Dificuldades de pré-registro sem testes prévios com o modelo" ~ "Difficulties in pre-registration without prior testing of the model.",
    dificuldade == "Logística: fornecedores" ~ "Issues with materials suppliers, including delays and communication problems.",
    dificuldade == "Dificuldade em obter licenças" ~ "Difficulties in obtaining licenses with Brazilian regulatory agencies.",
    dificuldade == "Dificuldade com fornecimento de biotérios" ~ "Limitations in animal suppliers and facilities.",
    TRUE ~ dificuldade
  )) |>
  rename(Category = categoria) |>
  rename(Difficulty = dificuldade) |>
  rename(Rating = soma) |>
  flextable() |>
  add_name("Table S24 – Project difficulties - Coordinating team assessment") |>
  set_table_properties(layout = "autofit")

### Saving ----
save_tbl(
  list(tbl_s24),
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S24.docx")
)

cat("\n### Table S24 generated! ###\n")

## Table S21 ---------------------------------------------------------------
# Survey participant demographics

### Data loading and cleaning ----
df_participants <- read.csv2("other-data/survey-data/joined_participant_information_anonymous.csv") |>
  mutate(
    # Fix gender
    Gender = if_else(str_detect(Gender_other, "non-binary"), "Non-binary", Gender, Gender),
    # Convert numeric fields
    Institutional_age = as.numeric(Institutional_age),
    General_rep_rate = as.numeric(General_rep_rate),
    # Consolidate positions
    Institutional_position = case_when(
      Institutional_position_other == "Research Analyst" ~ "Post-doc/Research associate",
      Institutional_position %in% c("Post-doctoral researcher", "Research associate") ~ "Post-doc/Research associate",
      TRUE ~ Institutional_position
    ),
    # Fix academic title
    Highest_academic_title = if_else(Highest_academic_title_other == "Pos Doc", "PhD", Highest_academic_title, Highest_academic_title),
    # Shorten country name
    Country_of_birth = str_replace(Country_of_birth, "United Kingdom of Great Britain and Northern Ireland", "United Kingdom")
  ) |>
  rename(Practical_knowledge = Practical_knowledge_technique, Theoretical_knowledge = Theoretical_knowledge_technique) |>
  # Add geographical region
  mutate(Geographical_location = case_when(
    Country_of_living == "Brazil" ~ "Brazil",
    Country_of_living %in% c("Germany", "United Kingdom of Great Britain and Northern Ireland", "Sweden", 
                              "Netherlands", "Austria", "Switzerland", "France", "Belgium", "Finland") ~ "Europe",
    Country_of_living == "Australia" ~ "Oceania",
    Country_of_living == "Canada" ~ "North America",
    TRUE ~ "Not specified"
  ))

### Summary functions ----
summarize_vars <- function(df, group_col = NULL) {
  num_vars <- c("Age", "Years_in_research", "General_rep_rate")
  
  # Numerical: median (range)
  num_summary <- df |>
    summarise(across(all_of(num_vars), list(
      stat = ~sprintf("%.0f (%.0f - %.0f)", median(.x, na.rm = TRUE), min(.x, na.rm = TRUE), max(.x, na.rm = TRUE))
    ))) |>
    pivot_longer(everything(), names_to = "Variable", values_to = "Stat") |>
    mutate(Variable = str_remove(Variable, "_stat"), Label = Variable)
  
  # Categorical: n (%)
  cat_vars <- c("Gender", "BRI_member", "Geographical_location", "Highest_academic_title", 
                "Institutional_position", "Theoretical_knowledge", "Practical_knowledge", 
                "Stats_knowledge", "Research_area1")
  
  cat_summary <- map_df(cat_vars, function(v) {
    df |>
      count(Label = as.character(.data[[v]])) |>
      filter(!is.na(Label)) |>
      mutate(Variable = v, Stat = sprintf("%d (%d%%)", n, round(100 * n / sum(n)))) |>
      select(Variable, Label, Stat)
  })
  
  bind_rows(num_summary, cat_summary)
}

### Generate summaries by technique ----
techniques <- c("Combined", "MTT", "PCR", "EPM")

all_summaries <- map_df(techniques, function(tech) {
  df_tech <- if (tech == "Combined") df_participants else filter(df_participants, Technique == tech)
  summarize_vars(df_tech) |> mutate(Technique = tech)
})

### Format table ----
# Define display order for labels
label_order <- tribble(
  ~Variable, ~Label, ~order,
  "BRI_member", "Yes", 1, "BRI_member", "No", 2,
  "Geographical_location", "Brazil", 1, "Geographical_location", "Europe", 2,
  "Geographical_location", "North America", 3, "Geographical_location", "Oceania", 4,
  "Highest_academic_title", "PhD", 1, "Highest_academic_title", "Master's", 2,
  "Institutional_position", "Professor", 1, "Institutional_position", "Post-doc/Research associate", 2,
  "Institutional_position", "PhD student", 3, "Institutional_position", "Master's student", 4,
  "Institutional_position", "Technician", 5, "Institutional_position", "Other (please specify)", 6
)

# Variable display order
var_order <- c("Age", "Gender", "BRI_member", "Geographical_location", "Years_in_research",
               "Highest_academic_title", "Institutional_position", "Research_area1",
               "Theoretical_knowledge", "Practical_knowledge", "Stats_knowledge", "General_rep_rate")

tbl_s21_data <- all_summaries |>
  left_join(label_order, by = c("Variable", "Label")) |>
  arrange(factor(Variable, levels = var_order), order) |>
  select(-order) |>
  pivot_wider(names_from = Technique, values_from = Stat, values_fill = "0 (0%)") |>
  # Rename Research_area1 for display
  mutate(Variable = if_else(Variable == "Research_area1", "Research_area", Variable)) |>
  # Collapse rows by variable
  group_by(Variable) |>
  summarise(
    Labels = str_c(Label, collapse = "\n"),
    across(all_of(techniques), ~ str_c(.x, collapse = "\n")),
    .groups = "drop"
  ) |>
  # Format labels
  mutate(
    is_numeric = Variable %in% c("Age", "Years_in_research", "General_rep_rate"),
    Labels = if_else(is_numeric, Variable, paste0(Variable, "\n", Labels)),
    across(all_of(techniques), ~ if_else(is_numeric, .x, paste0("\n", .x)))
  ) |>
  arrange(factor(Variable, levels = var_order)) |>
  select(-Variable, -is_numeric) |>
  mutate(Labels = str_replace_all(Labels, "_", " ")) |>
  rename("Combined\nn=70" = Combined, "MTT\nn=22" = MTT, "PCR\nn=18" = PCR, "EPM\nn=30" = EPM)

### Create flextable ----
tbl_s21 <- tbl_s21_data |>
  flextable() |>
  add_header_lines(values = "Table S21 – Survey participant demographics") |>
  bold(i = 1, part = "header") |>
  set_table_properties(layout = "autofit")

### Saving ----
save_tbl(
  list(tbl_s21),
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S21.docx")
)

cat("\n### Table S21 generated! ###\n")

## Text-Cited Numbers ------------------------------------------------------

### General ------------------------------------------------------------------

df_inclusion_sets <- df_inclusion_sets |> clean_names()

# New word document
doc <- read_docx()

# a) Count and % of experiments with log_es_ratio > 0 (PCR analysis)
result_a <- df_assessment_by_experiment |>
  filter(!str_detect(EXP, "ALTPCR") & !str_detect(EXP, "ALTMTT")) |>
  summarise(
    count = sum(log_es_ratio > 0, na.rm = TRUE),
    total = n(),
    percent = mean(log_es_ratio > 0, na.rm = TRUE) * 100
  )

doc <- doc |>
  body_add(fpar(ftext("a) Count and % of experiments with log_es_ratio > 0 (PCR Analysis)", prop = fp_text(bold = TRUE)))) |>
  body_add_par(paste("Count:", result_a$count, "Total:", result_a$total, "Percent:", sprintf("%.2f%%", result_a$percent)))

# b) Count and % of experiments where original_cv > replication_cv (ALTPCR analysis)
result_b <- df_assessment_by_experiment |>
  filter(str_detect(EXP, "^(EPM|MTT|ALTPCR)")) |>
  summarise(
    count = sum(original_cv < replication_cv, na.rm = TRUE),
    percent = mean(original_cv < replication_cv, na.rm = TRUE) * 100
  )

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("b) Count and % of experiments where original_cv > replication_cv (ALTPCR Analysis)", prop = fp_text(bold = TRUE)))) |>
  body_add_par(paste("Count:", result_b$count, "| Percent:", sprintf("%.2f%%", result_b$percent)))

# c) Wilcoxon signed-rank test: original_cv vs replication_cv (using ALTPCR)
result_c <- df_assessment_by_experiment |>
  filter(str_detect(EXP, "^(EPM|MTT|ALTPCR)")) |>
  summarise(
    p_value = wilcox.test(original_cv, replication_cv, paired = TRUE)$p.value,
    statistic = wilcoxonZ(original_cv, replication_cv, paired = TRUE)
  )

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("c) Wilcoxon Signed-Rank Test: original_cv vs replication_cv (ALTPCR)", prop = fp_text(bold = TRUE)))) |>
  body_add_par(paste("P-value:", result_c$p_value)) |> 
  body_add_par(paste("Statistic (z):", result_c$statistic))

# d) Wilcoxon test: original_cv vs replication_cv by technique
result_d <- df_assessment_by_experiment |>
  mutate(
    technique = case_when(
      str_detect(EXP, "EPM") ~ "EPM",
      str_detect(EXP, "^(MTT)") ~ "MTT",
      str_detect(EXP, "^(PCR)(?!ALT)") ~ "PCR",
      str_detect(EXP, "^(ALTPCR)") ~ "ALTPCR"
    )
  ) |>
  group_by(technique) |>
  summarise(
    p_value = wilcox.test(original_cv, replication_cv, paired = TRUE, na.action = na.omit)$p.value,
    statistic = wilcoxonZ(original_cv, replication_cv, paired = TRUE)
    
  )

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("d) Wilcoxon Signed-Rank Test: original_cv vs replication_cv by Technique", prop = fp_text(bold = TRUE))))
  
for (i in 1:nrow(result_d)) {
  doc <- doc |>
    body_add_par(paste("Technique:", result_d$technique[i], "P-value:", result_d$p_value[i],  "Statistic (z):", result_d$statistic[i])) 
}

# e) Wilcoxon signed-rank test: mean_abs_diff_reps vs mean_abs_diff_reps_orig (excluding ALTPCR)
result_e <- df_assessment_by_experiment |>
  filter(!str_detect(EXP, "ALTPCR")) |>
  summarise(
    p_value = wilcox.test(mean_abs_diff_reps, mean_abs_diff_reps_orig, paired = TRUE, na.action = na.omit)$p.value,
    statistic = wilcoxonZ(mean_abs_diff_reps, mean_abs_diff_reps_orig, paired = TRUE)
  )

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("e) Wilcoxon Signed-Rank Test: mean_abs_diff_reps vs mean_abs_diff_reps_orig", prop = fp_text(bold = TRUE)))) |>
  body_add_par(paste("P-value:", result_e$p_value)) |> 
  body_add_par(paste("Statistic (z):", result_e$statistic))

# f) Wilcoxon test: mean_abs_diff_reps vs mean_abs_diff_reps_orig by technique
result_f <- df_assessment_by_experiment |>
  mutate(
    technique = case_when(
      str_detect(EXP, "EPM") ~ "EPM",
      str_detect(EXP, "^(MTT)") ~ "MTT",
      str_detect(EXP, "^(PCR)(?!ALT)") ~ "PCR",
      str_detect(EXP, "^(ALTPCR)") ~ "ALTPCR"
    )
  ) |>
  group_by(technique) |>
  summarise(
    p_value = wilcox.test(mean_abs_diff_reps, mean_abs_diff_reps_orig, paired = TRUE, na.action = na.omit)$p.value,
    statistic = wilcoxonZ(mean_abs_diff_reps, mean_abs_diff_reps_orig, paired = TRUE)
  )

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("f) Wilcoxon Signed-Rank Test: mean_abs_diff_reps vs mean_abs_diff_reps_orig by Technique", prop = fp_text(bold = TRUE)))) 
  
for (i in 1:nrow(result_f)) {
  doc <- doc |>
    body_add_par(paste("Technique:", result_f$technique[i], "P-value:", result_f$p_value[i],  "Statistic (z):", result_f$statistic[i]))
}

# g) Count and % of experiments with Validation-mean score == 1
result_g <- df_inclusion_sets |>
  filter(!is.na(validation_mean_score)) |>
  filter(unit != "LABUNIT") |>
  mutate(
    technique = case_when(
      str_detect(exp, "EPM") ~ "EPM",
      str_detect(exp, "^(MTT)") ~ "MTT",
      str_detect(exp, "^(PCR)(?!ALT)") ~ "PCR",
      str_detect(exp, "^(ALTPCR)") ~ "ALTPCR"
    )
  ) |>
  group_by(technique) |>
  summarise(
    total_count = sum(validation_mean_score == 1),
    total_percent = mean(validation_mean_score == 1) * 100
  )

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("g) Count and % of experiments with Validation-mean score == 1", prop = fp_text(bold = TRUE))))
  
for (i in 1:nrow(result_g)) {
  doc <- doc |>
    body_add_par(paste(
      "Technique:", result_g$technique[i],
      "Total Count:", result_g$total_count[i],
      "Total Percent:", sprintf("%.2f%%", result_g$total_percent[i])
    ))
}

# h) Mean ± SD of Validation-mean score (all experiments)
result_h <- df_inclusion_sets |>
  filter(unit != "LABUNIT") |>
  summarise(
    mean = mean(validation_mean_score, na.rm = TRUE),
    sd = sd(validation_mean_score, na.rm = TRUE)
  )

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("h) Mean ± SD of Validation-mean score (all experiments)", prop = fp_text(bold = TRUE)))) |>
  body_add_par(paste("Mean:", sprintf("%.2f", result_h$mean), "SD:", sprintf("%.2f", result_h$sd)))

# h-2) Mean ± SD of Validation-mean score (excluding rediscussed experiments)
result_h2 <- df_inclusion_sets |>
  filter(unit != "LABUNIT") |>
  filter(validation_discussed != "E-mail" | is.na(validation_discussed)) |>
  summarise(
    mean = mean(validation_mean_score, na.rm = TRUE),
    sd = sd(validation_mean_score, na.rm = TRUE)
  )

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("h-2) Mean ± SD of Validation-mean score (excluding rediscussed experiments)", prop = fp_text(bold = TRUE)))) |>
  body_add_par(paste("Mean:", sprintf("%.2f", result_h2$mean), "SD:", sprintf("%.2f", result_h2$sd)))

# i) Mean ± SD of Validation-mean score (Primary analysis only)
result_i <- df_inclusion_sets |>
  filter(analysis_primary == "INCLUDE") |>
  filter(unit != "LABUNIT") |>
  summarise(
    mean = mean(validation_mean_score, na.rm = TRUE),
    sd = sd(validation_mean_score, na.rm = TRUE)
  )

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("i) Mean ± SD of Validation-mean score (Primary Analysis Only)", prop = fp_text(bold = TRUE)))) |>
  body_add_par(paste("Mean:", sprintf("%.2f", result_i$mean), "SD:", sprintf("%.2f", result_i$sd)))

# j) Mean ± SD of changes_from_preregistration (all experiments)
result_j <- df_inclusion_sets |>
  filter(unit != "LABUNIT") |>
  summarise(
    mean = mean(q1_changes_from_preregistration, na.rm = TRUE),
    sd = sd(q1_changes_from_preregistration, na.rm = TRUE)
  )

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("j) Mean ± SD of changes_from_preregistration (all experiments)", prop = fp_text(bold = TRUE)))) |>
  body_add_par(paste("Mean:", sprintf("%.2f", result_j$mean), "SD:", sprintf("%.2f", result_j$sd)))

# j-2) Mean ± SD of changes_from_preregistration (excluding rediscussed experiments)
result_j2 <- df_inclusion_sets |>
  filter(unit != "LABUNIT") |>
  filter(validation_discussed != "E-mail" | is.na(validation_discussed)) |>
  summarise(
    mean = mean(q1_changes_from_preregistration, na.rm = TRUE),
    sd = sd(q1_changes_from_preregistration, na.rm = TRUE)
  )

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("j-2) Mean ± SD of changes_from_preregistration (excluding rediscussed experiments)", prop = fp_text(bold = TRUE)))) |>
  body_add_par(paste("Mean:", sprintf("%.2f", result_j2$mean), "SD:", sprintf("%.2f", result_j2$sd)))

# k) Mean ± SD of changes_from_preregistration (Primary analysis only)
result_k <- df_inclusion_sets |>
  filter(unit != "LABUNIT") |>
  filter(analysis_primary == "INCLUDE") |>
  summarise(
    mean = mean(q1_changes_from_preregistration, na.rm = TRUE),
    sd = sd(q1_changes_from_preregistration, na.rm = TRUE)
  )

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("k) Mean ± SD of changes_from_preregistration (Primary Analysis Only)", prop = fp_text(bold = TRUE)))) |>
  body_add_par(paste("Mean:", sprintf("%.2f", result_k$mean), "SD:", sprintf("%.2f", result_k$sd)))

### Replicated by most criteria -----------------------------------------------------

tbl_most_criteria <- read_tsv(paste0("output/", results_path, "/primary t/Replication Rate Summary.tsv"), show_col_types = FALSE) |> 
  filter(Metric == "REP_Most_Criteria_WithTies") |> 
  pivot_longer(cols = -c(Metric, MetricLongName, MetricShortName)) |> 
  mutate(
    Method = name |> str_remove("^(Value|N|successful)_"),
    What = name |> str_extract("^(Value|N|successful)_") |> str_remove("_")
  ) |> 
  select(-name) |>  
  pivot_wider(id_cols = c(Metric, MetricShortName, Method), names_from = What, values_from = value)

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("l) Most criteria replicated", prop = fp_text(bold = TRUE)))) |>
  body_add_par(paste0(tbl_most_criteria$successful[tbl_most_criteria$Method == "ALL_PCR_MTT"], 
                     " (", round(tbl_most_criteria$Value[tbl_most_criteria$Method == "ALL_PCR_MTT"],2)*100, "%) ",
                     "experiments were replicated by at least half of the applicable criteria (primary, t, ALL_PCR_MTT)."
                                          )
               )

### Specification curve range -----------------------------------------------------
replication_rate_list <- read_tsv(paste0("output/", results_path, "/Replication Rate Aggregate - without Method.tsv"), show_col_types = FALSE) |> 
  rename(repro_rate = Value) |> 
  pull(repro_rate)

doc <- doc |>
  body_add_par(run_linebreak()) |>
  body_add(fpar(ftext("m) Specification curve range", prop = fp_text(bold = TRUE)))) |> 
  body_add_par(paste0(
    "Minimum: ", min(replication_rate_list), " | ",
    "Maximum: ", max(replication_rate_list), "."
  )
  )

### Fleiss Kappa ------------------------------------------------------------

predictor_correlations_exp <- read_tsv(paste0("output/", results_path, "/primary t/predictors/by experiment/Table for Predictor Correlations.tsv")) |>
  select(
    `Original in replication PI`,
    `Replication in original CI`,
    `Same sense significance`,
    `Voting (with ties)`,
    `Subjective (with ties)`
  )

data_exp <- predictor_correlations_exp |>
  mutate(across(everything(), ~ as.integer(. == TRUE)))


predictor_correlations_rep <- read_tsv(paste0("output/", results_path, "/primary t/predictors/by replication/Table for Predictor Correlations.tsv")) |>
  filter(EXP != "(Protocol Team)") |>
  filter(EXP != "(Data Collection Team)") |>
  select(
    `Replication in original CI`,
    `Same sense significance`,
    `Subjective`
  )

data_rep <- predictor_correlations_rep |>
  mutate(across(everything(), ~ as.integer(. == TRUE)))

fleiss_kappa_exp <- kappam.fleiss(as.matrix(data_exp |>
  drop_na()))

fleiss_kappa_rep <- kappam.fleiss(as.matrix(data_rep |>
  drop_na()))

# Recalculate p-values from z-statistic to avoid numeric underflow
# The kappam.fleiss$p.value can underflow to 0 when the z-statistic is very large,
# because R's double precision cannot represent values smaller than ~1e-308.
# Using pnorm() directly with the z-statistic avoids this issue.
p_value_exp <- 2 * pnorm(abs(fleiss_kappa_exp$statistic), lower.tail = FALSE)
p_value_rep <- 2 * pnorm(abs(fleiss_kappa_rep$statistic), lower.tail = FALSE)

# Create a data frame with the information
data_summary <- data.frame(
  `Criteria` = c("Experiment Level", "Replication Level"),
  `Subjects` = c(fleiss_kappa_exp$subjects, fleiss_kappa_rep$subjects),
  `Raters` = c(fleiss_kappa_exp$raters, fleiss_kappa_rep$raters),
  `Kappa` = c(round(fleiss_kappa_exp$value, 2), round(fleiss_kappa_rep$value, 2)),
  `p-value` = c(sprintf("%.2e", p_value_exp), sprintf("%.2e", p_value_rep))
)

# Create a flextable table
fleiss_kappa_table <- flextable(data_summary) |>
  set_header_labels(
    Criteria = "Replication Outcome Criteria",
    Subjects = "#",
    Raters = "Number of Criterias",
    Kappa = "Fleiss' Kappa",
    `p-value` = "p-value"
  ) |>
  autofit() |>
  bold(part = "header") |>
  align(align = "center", part = "all") |>
  add_name("Table - Fleiss' kappa for agreement between replication measures in primary analysis (t distribution).") |>
  flextable::footnote(
    i = 2, j = 3,
    value = as_paragraph(
      c(
        "Experiment-level criteria: Original in replication PI, Replication in original CI, Same sense significance, Voting (with ties), Subjective (with ties).
        Replication-level criteria: Replication in original CI, Same sense significance, Subjective."
      )
    ),
    ref_symbols = c("a"),
    part = "header"
  )


### Experiments per lab -----------------------------------------------------

tbl_n_exp_lab_primary <- df_inclusion_sets |>
  mutate(
    analysis_primary = if_else(analysis_primary == "INCLUDE", TRUE, FALSE)
  ) |>
  group_by(lab) |>
  summarise(
    n = sum(analysis_primary, na.rm = T)
  ) |>
  group_by(n) |>
  summarise(
    number_labs = n()
  ) |>
  pivot_wider(
    names_from = "n",
    values_from = "number_labs"
  ) |>
  mutate(
    Total = 62
  ) |>
  mutate(across(everything(), as.character)) |>
  mutate(
    `# of experiments` = "# of laboratories"
  ) |>
  relocate(`# of experiments`, .before = `0`) |>
  flextable() |>
  add_name("Table - Number of experiments (included in the primary analysis) per laboratory.") |>
  bold(i = 2, part = "header") |>
  bold(part = "header") |>
  set_table_properties(layout = "autofit")


tbl_n_exp_lab_all <- df_inclusion_sets |>
  mutate(
    analysis_all_exps = if_else(analysis_all_exps_lab_units == "INCLUDE", TRUE, FALSE)
  ) |>
  group_by(lab) |>
  summarise(
    n = sum(analysis_all_exps, na.rm = T)
  ) |>
  group_by(n) |>
  summarise(
    number_labs = n()
  ) |>
  pivot_wider(
    names_from = "n",
    values_from = "number_labs"
  ) |>
  mutate(
    Total = 62
  ) |>
  mutate(across(everything(), as.character)) |>
  mutate(
    `# of experiments` = "# of laboratories"
  ) |>
  relocate(`# of experiments`, .before = `0`) |>
  flextable() |>
  add_name("Table - Number of experiments (all experiments) per laboratory.") |>
  bold(i = 2, part = "header") |>
  bold(part = "header") |>
  set_table_properties(layout = "autofit")

### Saving ------------------------------------------------------------------

# Task 4: Supplementary table with sample size comparison (normal t-test vs log ratio of means)
df_ss_comp <- read_tsv("./other-data/Original Experiments Statistical Summaries.tsv", show_col_types = FALSE) |>
  select(
    EXP,
    `Estimated Sample Size for 95% power`,
    `Estimated Sample Size for 95% power (Log Ratio)`
  ) |>
  rename(
    `Experiment Code` = EXP,
    `Sample Size (Normal t-test, 95% Power)` = `Estimated Sample Size for 95% power`,
    `Sample Size (Log Ratio of Means, 95% Power)` = `Estimated Sample Size for 95% power (Log Ratio)`
  ) |>
  mutate(
    `Sample Size (Normal t-test, 95% Power)` = round(`Sample Size (Normal t-test, 95% Power)`, 2),
    `Sample Size (Log Ratio of Means, 95% Power)` = round(`Sample Size (Log Ratio of Means, 95% Power)`, 2)
  )

tbl_sample_size_comparison <- df_ss_comp |>
  flextable() |>
  set_caption("Supplementary Table: Sample size comparison for 95% power using normal t-tests vs. log ratio of means") |>
  autofit()

dir.create(paste0("output/", results_path, "/_manuscript figures and tables/tables"), showWarnings = FALSE, recursive = TRUE)
save_as_docx(tbl_sample_size_comparison, path = paste0("output/", results_path, "/_manuscript figures and tables/tables/Supplementary Table - Sample Size Comparison.docx"))

doc <- doc |>
  body_add_break(pos = "after") |>
  body_add_flextable(fleiss_kappa_table) |>
  body_add_par(run_linebreak()) |>
  body_add_par(run_linebreak()) |>
  body_add_par(run_linebreak()) |>
  body_add_flextable(tbl_n_exp_lab_primary) |>
  body_add_par(run_linebreak()) |>
  body_add_par(run_linebreak()) |>
  body_add_par(run_linebreak()) |>
  body_add_flextable(tbl_n_exp_lab_all) |>
  body_add_par(run_linebreak()) |>
  body_add_par(run_linebreak()) |>
  body_add_par(run_linebreak()) |>
  body_add_flextable(tbl_sample_size_comparison)

print(doc, target = paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Document - Text-Cited Numbers.docx"))

cat("\n### Text-Cited Numbers Document generated! ###\n")

options(warn = 0)
