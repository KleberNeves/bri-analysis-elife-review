### ### ### ### ### ### ### ### ### ### ### ###

# plots.r
# This script contains all the plotting functions and aesthetics definitions.

### ### ### ### ### ### ### ### ### ### ### ###

#### Definitions for aesthetics, plot and color defaults ####

bri_color_scale <- brewer.pal(9, "Blues")

bri_color <- c(
  main = "#4292C6",
  low = "#CC3D3D",
  mid = "#ACACAC",
  high = "#4292C6",
  none = "#333333",
  dark = "#111111",
  light = "#C6DBEF",
  emphasis = "#082857",
  qa1 = "#4292C7",
  qb1 = "#001847",
  qa2 = "#CC3D3D",
  qb2 = "#982020",
  q1 = "#4292C7",
  q2 = "#57C746",
  q3 = "#C79645",
  q4 = "#C74180",
  q5 = "#5C4474",
  q6 = "#333333",
  epm = "#2f52a3",
  mtt = "#1cb075",
  pcr = "#d95f02",
  epm_dark = "#003267",
  mtt_dark = "#005700",
  pcr_dark = "#573600"
)

bri_3color <- unname(c(bri_color["q1"], bri_color["q2"], bri_color["q3"]))
bri_6color <- unname(c(bri_color["q1"], bri_color["q2"], bri_color["q3"], bri_color["q4"], bri_color["q5"], bri_color["q6"]))

update_geom_defaults("point", list(
  color = bri_color[["main"]],
  size = 2
))

update_geom_defaults("bar", list(
  fill = bri_color[["main"]]
))

update_geom_defaults("col", list(
  fill = bri_color[["main"]]
))

bri_theme <- list(
  theme_bw(),
  theme(
    panel.border = element_rect(fill = NA, linewidth = 0.7, color = "black"),
    legend.position = "bottom",
    axis.line = element_blank(),
    axis.text = element_text(size = 11),
    strip.background = element_rect(fill = bri_color[["main"]]),
    strip.text = element_text(colour = "white")
  )
)

scale_y_continuous_col <- function(...) ggplot2:::scale_y_continuous(..., expand = c(0, 0))

#### Functions ####

bri_ggsave <- function(
    filename,
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = NA,
    height = NA,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = "white",
    create.dir = FALSE,
    ...) {
  ggsave(
    filename,
    plot,
    device,
    path,
    scale,
    width,
    height,
    units,
    dpi,
    limitsize,
    bg,
    create.dir,
    ...
  )
}

# Runs the correlation calculation, individual plots and then plots a correlation table of predictors/outcomes
plot_cortable <- function(FDATA, vars1, vars2, title, fn, show_label, make_individual_plots) {
  combs <- expand_grid(c1 = vars1, c2 = vars2)

  test_results <- pmap_dfr(combs, function(c1, c2) {
    run_correlation(FDATA, c1, c2, dirname(fn), make_plot = make_individual_plots)
  })

  vars1 <- ifelse(vars1 == "Log ES Ratio", "-log ES Ratio", vars1)
  vars2 <- ifelse(vars2 == "Log ES Ratio", "-log ES Ratio", vars2)

  test_results <- test_results |>
    filter(col1 != col2) |>
    mutate(
      # Change sign for logES ratio so that it intuitively matches the others
      rho = ifelse(
        col1 == "Log ES Ratio" | col2 == "Log ES Ratio",
        -rho, rho
      ),
      pearson_r = ifelse(
        col1 == "Log ES Ratio" | col2 == "Log ES Ratio",
        -pearson_r, pearson_r
      ),
      col1 = ifelse(col1 == "Log ES Ratio", "-log ES Ratio", col1),
      col2 = ifelse(col2 == "Log ES Ratio", "-log ES Ratio", col2),

      # Create formatted labels
      label = paste0(round(rho, 2), "\n(n = ", n, ")"),
      label1 = paste0(round(rho, 2), "\n"),
      label2 = paste0("\n(n = ", n, ")"),
      capped.p_rho = ifelse(p.value < 0.001, 0.001, p.value),
      label_p = paste0(round(pearson_r, 2), "\n(n = ", n, ")"),
      label1_p = paste0(round(pearson_r, 2), "\n"),
      label2_p = paste0("\n(n = ", n, ")"),
      capped.p_pearson = ifelse(pearson_p < 0.001, 0.001, pearson_p),

      # Change to factor to force the ordering
      col1 = factor(col1, levels = vars1),
      col2 = factor(col2, levels = vars2),

      # Change axis labels
      col1 = col1 |>
        fct_recode(
          # Predictors
          "Animal" = "Animal Model",
          "Coeff. Variation" = "Original CV",
          "Survey\n(Replication)" = "Survey - Replication",
          "Survey\n(Effect Size)" = "Survey - Effect Size",
          "Survey\n(Challenge)" = "Survey - Challenge",
          "Protocol Deviation\n(Lab)" = "Protocol deviation (Lab)",
          "Protocol Deviation\n(Committee)" = "Protocol deviation (Committee)",
          "Academic Age\n(Protocol Team)" = "Academic age (Protocol Team)",
          "Academic Age\n(Data Collection Team)" = "Academic age (Data Collection Team)",
          "Publications\n(Protocol Team)" = "Publications (Protocol Team)",
          "Publications\n(Data Collection Team)" = "Publications (Data Collection Team)",
          "University Ranking" = "Institution Ranking",

          # Outcomes
          "Original in replication's 95% PI" = "Original in replication PI",
          "Replication in original 95% CI" = "Replication in original CI",
          "Same-sign significance (p<0.05)" = "Same sense significance",
          "≥50% replications significant" = "Voting (with ties)",
          "≥50% subjectively replicated" = "Subjective (with ties)",
          "Subjectively replicated" = "Subjective"
        ),
      col2 = col2 |>
        fct_recode(
          # Predictors
          "Animal" = "Animal Model",
          "Coeff. Variation" = "Original CV",
          "Survey\n(Replication)" = "Survey - Replication",
          "Survey\n(Effect Size)" = "Survey - Effect Size",
          "Survey\n(Challenge)" = "Survey - Challenge",
          "Protocol Deviation\n(Lab)" = "Protocol deviation (Lab)",
          "Protocol Deviation\n(Committee)" = "Protocol deviation (Committee)",
          "Academic Age\n(Protocol Team)" = "Academic age (Protocol Team)",
          "Academic Age\n(Data Collection Team)" = "Academic age (Data Collection Team)",
          "Publications\n(Protocol Team)" = "Publications (Protocol Team)",
          "Publications\n(Data Collection Team)" = "Publications (Data Collection Team)",
          "University Ranking" = "Institution Ranking",

          # Outcomes
          "Original in replication's 95% PI" = "Original in replication PI",
          "Replication in original 95% CI" = "Replication in original CI",
          "Same-sign significance (p<0.05)" = "Same sense significance",
          "≥50% replications significant" = "Voting (with ties)",
          "≥50% subjectively replicated" = "Subjective (with ties)",
          "Subjectively replicated" = "Subjective"
        )
    )

  p <- ggplot(test_results) +
    aes(
      x = fct_relevel(col1, "Coeff. Variation", after = 6),
      y = fct_rev(fct_relevel(
        col2,
        "t Value",
        "-log ES Ratio",
        "Replication Effect Size",
        "SMD Effect Size",
        "Original in replication's 95% PI"
      )),
      fill = sign(rho) * -log10(p.value),
      label = label
    ) +
    geom_tile() +
    labs(x = "", y = "", title = title, fill = "-log10(p-value)") +
    scale_fill_gradient2(
      low = bri_color[["low"]],
      mid = bri_color[["mid"]],
      high = bri_color[["high"]],
      na.value = scales::alpha("white", 0)
    ) +
    bri_theme +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      title = element_text(size = 9)
    )

  p_p <- ggplot(test_results) +
    aes(
      x = col1, y = fct_rev(fct_relevel(
        col2, "t Value",
        "-log ES Ratio",
        "Replication Effect Size",
        "SMD Effect Size",
        "Original in replication's 95% PI"
      )),
      fill = sign(pearson_r) * -log10(pearson_p),
      label = label_p
    ) +
    geom_tile() +
    labs(x = "", y = "", title = title, fill = "-log10(p-value)") +
    scale_fill_gradient2(
      low = bri_color[["low"]],
      mid = bri_color[["mid"]],
      high = bri_color[["high"]],
      na.value = scales::alpha("white", 0)
    ) +
    bri_theme +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      title = element_text(size = 9)
    )

  if (show_label) {
    p <- p +
      geom_text(
        aes(label = label1),
        color = bri_color[["dark"]], size = 3.0
      ) +
      geom_text(
        aes(label = label2),
        color = bri_color[["dark"]], size = 2.5
      )

    p_p <- p_p +
      geom_text(
        aes(label = label1_p),
        color = bri_color[["dark"]], size = 3.0
      ) +
      geom_text(
        aes(label = label2_p),
        color = bri_color[["dark"]], size = 2.5
      )
  }

  p <- p +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(panel.border = element_blank())

  write_tsv(test_results |> select(-label1, -label2, -label1_p, -label2_p), paste0("./", fn, " cordata.tsv"))

  bri_ggsave(paste0("./", fn, " - spearman.png"),
    plot = p,
    width = max(1 + length(vars1) * 0.5, 8), height = max(1 + length(vars2) * 0.5, 5)
  )

  bri_ggsave(paste0("./", fn, " - pearson.png"),
    plot = p_p,
    width = max(1 + length(vars1) * 0.5, 8), height = max(1 + length(vars2) * 0.5, 5)
  )

  p
}

# An alternative plot_cortable function to combine predictors and outcomes into a single figure.
plot_cortable_alternative <- function(FDATA, FDATA_EXP, FDATA_REP, vars1, vars2, title, fn, show_label, make_individual_plots) {
  global_vars <- setdiff(colnames(FDATA_EXP), "EXP")
  individual_vars <- setdiff(colnames(FDATA_REP), c("EXP", "LAB"))

  # Compute correlations between global (EXP) variables
  global_grouped <- expand_grid(
    c1 = global_vars,
    c2 = global_vars
  ) |>
    filter(c1 != c2)

  # Compute correlations between individual (REP) variables
  within_individual <- expand_grid(
    c1 = individual_vars,
    c2 = individual_vars
  ) |>
    filter(c1 != c2)

  # Compute correlations between global (EXP) and individual (REP) variables
  between <- expand_grid(
    c1 = global_vars,
    c2 = individual_vars
  )

  # Combine results
  results <- bind_rows(global_grouped, within_individual, between)

  test_results <- pmap_dfr(results, function(c1, c2) {
    run_correlation_alternative(FDATA, c1, c2, dirname(fn), make_plot = F, global_vars = setdiff(colnames(FDATA_EXP), "EXP"))
  })

  test_results <- test_results |>
    rowwise() |>
    mutate(
      idx1 = match(col1, vars1),
      idx2 = match(col2, vars1),
      
      # Swap if in upper triangle (idx1 > idx2) to enforce lower triangle (idx1 <= idx2)
      col1_new = ifelse(idx1 > idx2, col2, col1),
      col2_new = ifelse(idx1 > idx2, col1, col2),
      
      col1 = col1_new,
      col2 = col2_new
    ) |>
    ungroup() |>
    select(-idx1, -idx2, -col1_new, -col2_new) |>
    distinct(col1, col2, .keep_all = TRUE) |>
    filter(col1 != col2) |>
    mutate(
      # Create formatted labels
      label = paste0(round(rho, 2), "\n(n = ", n, ")"),
      # label1 = paste0(round(rho, 2), "\n"),
      # label2 = paste0("\n(n = ", n, ")"),
      label1 = if_else(is.na(rho), "", paste0(round(rho, 2), "\n")),
      label2 = if_else(is.na(rho), "", paste0("\n(n = ", n, ")")),
      capped.p_rho = ifelse(p.value < 0.001, 0.001, p.value),
      label_p = paste0(round(pearson_r, 2), "\n(n = ", n, ")"),
      label1_p = paste0(round(pearson_r, 2), "\n"),
      label2_p = paste0("\n(n = ", n, ")"),
      capped.p_pearson = ifelse(pearson_p < 0.001, 0.001, pearson_p),

      # Change to factor to force the ordering
      col1 = factor(col1, levels = vars1),
      col2 = factor(col2, levels = vars2),

      # Change axis labels
      col1 = col1 |>
        fct_recode(
          # Predictors
          "Animal" = "Animal Model",
          "Coeff. Variation" = "Original CV",
          "Survey\n(Replication)" = "Survey - Replication",
          "Survey\n(Effect Size)" = "Survey - Effect Size",
          "Survey\n(Challenge)" = "Survey - Challenge",
          "Protocol Deviation\n(Lab)" = "Protocol deviation (Lab)",
          "Protocol Deviation\n(Committee)" = "Protocol deviation (Committee)",
          "Academic Age\n(Protocol Team)" = "Academic age (Protocol Team)",
          "Academic Age\n(Data Collection Team)" = "Academic age (Data Collection Team)",
          "Publications\n(Protocol Team)" = "Publications (Protocol Team)",
          "Publications\n(Data Collection Team)" = "Publications (Data Collection Team)"
        ),
      col2 = col2 |>
        fct_recode(
          # Predictors
          "Animal" = "Animal Model",
          "Coeff. Variation" = "Original CV",
          "Survey\n(Replication)" = "Survey - Replication",
          "Survey\n(Effect Size)" = "Survey - Effect Size",
          "Survey\n(Challenge)" = "Survey - Challenge",
          "Protocol Deviation\n(Lab)" = "Protocol deviation (Lab)",
          "Protocol Deviation\n(Committee)" = "Protocol deviation (Committee)",
          "Academic Age\n(Protocol Team)" = "Academic age (Protocol Team)",
          "Academic Age\n(Data Collection Team)" = "Academic age (Data Collection Team)",
          "Publications\n(Protocol Team)" = "Publications (Protocol Team)",
          "Publications\n(Data Collection Team)" = "Publications (Data Collection Team)"
        )
    )

  # Prepare data with separate fill variables for positive and negative correlations
  plot_data <- test_results |>
    filter(col2 != "Animal") |>
    distinct(pearson_r, .keep_all = T) |>
    mutate(
      fill_negative = ifelse(rho < 0, -log10(capped.p_rho), NA),
      fill_positive = ifelse(rho > 0, -log10(capped.p_rho), NA)
    )

  p <- ggplot(plot_data) +
    aes(x = col1, y = fct_rev(col2), label = label) +
    
    # First scale: Negative correlations (grey to red)
    geom_tile(aes(fill = fill_negative), data = subset(plot_data, rho < 0)) +
    scale_fill_gradient(
      low = bri_color[["mid"]],
      high = bri_color[["low"]],
      na.value = "transparent",
      name = "Negative correlation\n(p-value)",
      breaks = c(1, 2, 3),
      labels = c(expression(10^-1), expression(10^-2), expression(10^-3)),
      guide = guide_colorbar(order = 1, barwidth = unit(0.6, "cm"), barheight = unit(3, "cm"))
    ) +
    
    # Reset scale for second layer
    ggnewscale::new_scale_fill() +
    
    # Second scale: Positive correlations (grey to blue)
    geom_tile(aes(fill = fill_positive), data = subset(plot_data, rho > 0)) +
    scale_fill_gradient(
      low = bri_color[["mid"]],
      high = bri_color[["high"]],
      na.value = "transparent",
      name = "Positive correlation\n(p-value)",
      breaks = c(1, 2, 3),
      labels = c(expression(10^-1), expression(10^-2), expression(10^-3)),
      guide = guide_colorbar(order = 2, barwidth = unit(0.6, "cm"), barheight = unit(3, "cm"))
    ) +
    
    labs(x = "", y = "", title = title) +
    bri_theme +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      title = element_text(size = 9)
    )

  # Prepare data for Pearson correlation with separate fill variables
  plot_data_p <- test_results |>
    mutate(
      fill_negative = ifelse(pearson_r < 0, -log10(capped.p_pearson), NA),
      fill_positive = ifelse(pearson_r > 0, -log10(capped.p_pearson), NA)
    )

  p_p <- ggplot(plot_data_p) +
    aes(
      x = col1, y = fct_rev(fct_relevel(
        col2, "t Value",
        "-log ES Ratio",
        "Replication Effect Size",
        "SMD Effect Size",
        "Original in replication's 95% PI"
      )),
      label = label_p
    ) +
    
    # First scale: Negative correlations (grey to red)
    geom_tile(aes(fill = fill_negative), data = subset(plot_data_p, pearson_r < 0)) +
    scale_fill_gradient(
      low = bri_color[["mid"]],
      high = bri_color[["low"]],
      na.value = "transparent",
      name = "Negative correlation\n(p-value)",
      breaks = c(1, 2, 3),
      labels = c(expression(10^-1), expression(10^-2), expression(10^-3)),
      guide = guide_colorbar(order = 1, barwidth = unit(0.6, "cm"), barheight = unit(3, "cm"))
    ) +
    
    # Reset scale for second layer
    ggnewscale::new_scale_fill() +
    
    # Second scale: Positive correlations (grey to blue)
    geom_tile(aes(fill = fill_positive), data = subset(plot_data_p, pearson_r > 0)) +
    scale_fill_gradient(
      low = bri_color[["mid"]],
      high = bri_color[["high"]],
      na.value = "transparent",
      name = "Positive correlation\n(p-value)",
      breaks = c(1, 2, 3),
      labels = c(expression(10^-1), expression(10^-2), expression(10^-3)),
      guide = guide_colorbar(order = 2, barwidth = unit(0.6, "cm"), barheight = unit(3, "cm"))
    ) +
    
    labs(x = "", y = "", title = title) +
    bri_theme +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      title = element_text(size = 9)
    )

  if (show_label) {
    p <- p +
      geom_text(
        aes(label = label1),
        color = bri_color[["dark"]], size = 3.0
      ) +
      geom_text(
        aes(label = label2),
        color = bri_color[["dark"]], size = 2.5
      )

    p_p <- p_p +
      geom_text(
        aes(label = label1_p),
        color = bri_color[["dark"]], size = 3.0
      ) +
      geom_text(
        aes(label = label2_p),
        color = bri_color[["dark"]], size = 2.5
      )
  }

  # Final aesthetics changes
  p <- p +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(panel.border = element_blank()) +
    theme(legend.position = "right", legend.key.height = unit(0.8, "cm"))

  # Saving
  write_tsv(test_results |> select(-label1, -label2, -label1_p, -label2_p), paste0("./", fn, " cordata.tsv"))

  bri_ggsave(paste0("./", fn, " - spearman.png"),
    plot = p,
    width = max(1 + length(vars1) * 0.5, 8), height = max(1 + length(vars2) * 0.5, 5)
  )

  bri_ggsave(paste0("./", fn, " - pearson.png"),
    plot = p_p,
    width = max(1 + length(vars1) * 0.5, 8), height = max(1 + length(vars2) * 0.5, 5)
  )

  p
}

# Plots clusters within predictor x rep outcome combined figures
plot_cortable_cluster <- function(FDATA, vars1, vars2, title, fn, show_label, make_individual_plots) {
  combs <- expand_grid(c1 = vars1, c2 = vars2)

  test_results <- pmap_dfr(combs, function(c1, c2) {
    run_correlation(FDATA, c1, c2, dirname(fn), make_plot = make_individual_plots)
  })

  vars1 <- ifelse(vars1 == "Log ES Ratio", "-log ES Ratio", vars1)
  vars2 <- ifelse(vars2 == "Log ES Ratio", "-log ES Ratio", vars2)

  test_results <- test_results |>
    filter(col1 != col2) |>
    mutate(
      # Change sign for logES ratio so that it intuitively matches the others
      rho = ifelse(
        col1 == "Log ES Ratio" | col2 == "Log ES Ratio",
        -rho, rho
      ),
      pearson_r = ifelse(
        col1 == "Log ES Ratio" | col2 == "Log ES Ratio",
        -pearson_r, pearson_r
      ),
      col1 = ifelse(col1 == "Log ES Ratio", "-log ES Ratio", col1),
      col2 = ifelse(col2 == "Log ES Ratio", "-log ES Ratio", col2),

      # Create formatted labels
      label = paste0(round(rho, 2), "\n(n = ", n, ")"),
      label1 = paste0(round(rho, 2), "\n"),
      label2 = paste0("\n(n = ", n, ")"),
      capped.p_rho = ifelse(p.value < 0.001, 0.001, p.value),
      minus_p_rho = 1 - p.value,
      label_p = paste0(round(pearson_r, 2), "\n(n = ", n, ")"),
      label1_p = paste0(round(pearson_r, 2), "\n"),
      label2_p = paste0("\n(n = ", n, ")"),
      capped.p_pearson = ifelse(pearson_p < 0.001, 0.001, pearson_p),
      minus_p_pearson = 1 - pearson_p,

      # Change to factor to force the ordering
      col1 = factor(col1, levels = vars1),
      col2 = factor(col2, levels = vars2),

      # Change axis labels
      col1 = col1 |>
        fct_recode(
          # Predictors
          "Animal" = "Animal Model",
          "Coeff. Variation" = "Original CV",
          "Survey\n(Replication)" = "Survey - Replication",
          "Survey\n(Effect Size)" = "Survey - Effect Size",
          "Survey\n(Challenge)" = "Survey - Challenge",
          "Protocol Deviation\n(Lab)" = "Protocol deviation (Lab)",
          "Protocol Deviation\n(Committee)" = "Protocol deviation (Committee)",
          "Academic Age\n(Protocol Team)" = "Academic age (Protocol Team)",
          "Academic Age\n(Data Collection Team)" = "Academic age (Data Collection Team)",
          "Publications\n(Protocol Team)" = "Publications (Protocol Team)",
          "Publications\n(Data Collection Team)" = "Publications (Data Collection Team)",
          "University Ranking" = "Institution Ranking",
          "Protocol" = "Reporting",

          # Outcomes
          "Original in replication's 95% PI" = "Original in replication PI",
          "Replication in original 95% CI" = "Replication in original CI",
          "Same-sign significance (p<0.05)" = "Same sense significance",
          "≥50% replications significant" = "Voting (with ties)",
          "≥50% subjectively replicated" = "Subjective (with ties)",
          "Subjectively replicated" = "Subjective"
        ),
      col2 = col2 |>
        fct_recode(
          # Predictors
          "Animal" = "Animal Model",
          "Coeff. Variation" = "Original CV",
          "Survey\n(Replication)" = "Survey - Replication",
          "Survey\n(Effect Size)" = "Survey - Effect Size",
          "Survey\n(Challenge)" = "Survey - Challenge",
          "Protocol Deviation\n(Lab)" = "Protocol deviation (Lab)",
          "Protocol Deviation\n(Committee)" = "Protocol deviation (Committee)",
          "Academic Age\n(Protocol Team)" = "Academic age (Protocol Team)",
          "Academic Age\n(Data Collection Team)" = "Academic age (Data Collection Team)",
          "Publications\n(Protocol Team)" = "Publications (Protocol Team)",
          "Publications\n(Data Collection Team)" = "Publications (Data Collection Team)",
          "University Ranking" = "Institution Ranking",
          "Protocol" = "Reporting",

          # Outcomes
          "Original in replication's 95% PI" = "Original in replication PI",
          "Replication in original 95% CI" = "Replication in original CI",
          "Same-sign significance (p<0.05)" = "Same sense significance",
          "≥50% replications significant" = "Voting (with ties)",
          "≥50% subjectively replicated" = "Subjective (with ties)",
          "Subjectively replicated" = "Subjective"
        ),
    )

  # Plots spearman correlation
  plot_data <- test_results |>
    mutate(group_x = case_when(
      col1 == "Animal" ~ "Model/Method",
      col1 == "EPM" ~ "Model/Method",
      col1 == "MTT" ~ "Model/Method",
      col1 == "PCR" ~ "Model/Method",
      col1 == "Hedges's g" ~ "Statistics",
      col1 == "p-value" ~ "Statistics",
      col1 == "Coeff. Variation" ~ "Statistics",
      col1 == "Protocol" ~ "Reporting",
      col1 == "Bias Control" ~ "Reporting",
      col1 == "Impact Factor" ~ "Article",
      col1 == "Year" ~ "Article",
      col1 == "Citations" ~ "Article",
      col1 == "Academic Age" ~ "Author",
      col1 == "Publications" ~ "Author",
      col1 == "University Ranking" ~ "Lab",
      col1 == "Survey\n(Replication)" ~ "Prediction",
      col1 == "Survey\n(Effect Size)" ~ "Prediction",
      col1 == "Survey\n(Challenge)" ~ "Prediction",
      col1 == "Protocol Deviation\n(Lab)" ~ "Protocol Deviation",
      col1 == "Protocol Deviation\n(Committee)" ~ "Protocol Deviation",
      col1 == "Academic Age\n(Protocol Team)" ~ "Replication Team",
      col1 == "Publications\n(Protocol Team)" ~ "Replication Team",
      col1 == "Academic Age\n(Data Collection Team)" ~ "Replication Team",
      col1 == "Publications\n(Data Collection Team)" ~ "Replication Team",
      col1 == "University Ranking" ~ "Lab",
      .default = "Statistics"
    )) |>
    drop_na()
  
  # Guard against empty data after dropping NAs
  if (nrow(plot_data) == 0) {
    message("Warning: No valid data for correlation plot after dropping NAs. Skipping plot generation.")
    return(NULL)
  }
  
  plot_data <- plot_data |>
    mutate(group_x = fct_relevel(
      group_x,
      "Model/Method",
      "Statistics",
      "Reporting",
      "Article",
      "Author",
      "Protocol Deviation",
      "Replication Team",
      "Lab",
      "Prediction"
    )) |>
    # Brutal force ordering: Prepend index to string to force alphabetical order in ggplot
    mutate(
      desired_order = factor(col1, levels = c(
        "Animal", "EPM", "MTT", "PCR",
        "Hedges's g", "p-value", "Coeff. Variation",
        "Protocol", "Bias Control",
        "Impact Factor", "Year", "Citations",
        "Academic Age", "Publications",
        "Protocol Deviation\n(Lab)", "Protocol Deviation\n(Committee)",
        "Academic Age\n(Protocol Team)", "Publications\n(Protocol Team)",
        "Academic Age\n(Data Collection Team)", "Publications\n(Data Collection Team)",
        "University Ranking",
        "Survey\n(Replication)",
        "Survey\n(Effect Size)",
        "Survey\n(Challenge)"
      )),
      col1_idx = as.integer(desired_order),
      # Fallback for unexpected items: place at end
      col1_idx = ifelse(is.na(col1_idx), 99, col1_idx),
      col1 = paste0(sprintf("%02d_", col1_idx), col1)
    ) |>
    mutate(
      # Create separate fill variables for positive and negative correlations using -log10(p)
      fill_negative = ifelse(rho < 0, -log10(capped.p_rho), NA_real_),
      fill_positive = ifelse(rho > 0, -log10(capped.p_rho), NA_real_)
    )

  compute_scale_params <- function(values) {
    max_val <- suppressWarnings(max(values, na.rm = TRUE))
    if (!is.finite(max_val) || max_val <= 0) {
      list(limits = c(0, 0.1), breaks = c(0, 0.1))
    } else {
      breaks <- scales::pretty_breaks(n = 3)(c(0, max_val))
      breaks <- unique(breaks[breaks >= 0 & breaks <= max_val])
      if (length(breaks) == 0) {
        breaks <- c(0, max_val)
      }
      list(limits = c(0, max_val), breaks = breaks)
    }
  }

  neg_scale <- compute_scale_params(plot_data$fill_negative)
  pos_scale <- compute_scale_params(plot_data$fill_positive)

  p <- ggplot(plot_data) +
    aes(
      x = col1,
      y = fct_rev(fct_relevel(
        col2,
        "t Value",
        "-log ES Ratio",
        "Replication Effect Size",
        "SMD Effect Size",
        "Original in replication's 95% PI"
      )),
      label = label
    ) +
    
    # First scale: Negative correlations (grey to red)
    geom_tile(aes(fill = fill_negative), data = subset(plot_data, rho < 0)) +
    scale_fill_gradient(
      low = bri_color[["mid"]],
      high = bri_color[["low"]],
      na.value = bri_color[["mid"]],
      name = "Negative correlation\n(p-value)",
      breaks = c(0, 1, 2, 3),
      labels = c(expression(1), expression(10^-1), expression(10^-2), expression(10^-3)),
      limits = c(0, 3),
      guide = guide_colorbar(order = 1, barwidth = unit(0.6, "cm"), barheight = unit(2.5, "cm"), direction = "vertical")
    ) +
    
    # Reset scale for second layer
    ggnewscale::new_scale_fill() +
    
    # Second scale: Positive correlations (grey to blue)
    geom_tile(aes(fill = fill_positive), data = subset(plot_data, rho > 0)) +
    scale_fill_gradient(
      low = bri_color[["mid"]],
      high = bri_color[["high"]],
      na.value = bri_color[["mid"]],
      name = "Positive correlation\n(p-value)",
      breaks = c(0, 1, 2, 3),
      labels = c(expression(1), expression(10^-1), expression(10^-2), expression(10^-3)),
      limits = c(0, 3),
      guide = guide_colorbar(order = 2, barwidth = unit(0.6, "cm"), barheight = unit(2.5, "cm"), direction = "vertical")
    ) +
    
    # Third scale: Zero correlations (grey background)
    ggnewscale::new_scale_fill() +
    geom_tile(fill = bri_color[["mid"]], data = subset(plot_data, rho == 0)) +
    
    facet_grid(cols = ggplot2::vars(group_x), scales = "free_x", space = "free_x") +
    labs(x = "", y = "") +
    bri_theme +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      title = element_text(size = 9),
      strip.background = element_rect(fill = "transparent"),
      strip.text = element_text(color = "black"), # center: 8
      panel.spacing.x = unit(0, "line"),
      panel.border = element_rect(color = "#333", fill = NA, size = 0.5),
      legend.position = c(1.25, 0.5),
      legend.box = "vertical"
    )

  # Plots Pearson correlation
  p_p <- test_results |>
    mutate(group_x = case_when(
      col1 == "Animal" ~ "Model/Method",
      col1 == "EPM" ~ "Model/Method",
      col1 == "MTT" ~ "Model/Method",
      col1 == "PCR" ~ "Model/Method",
      col1 == "Hedges's g" ~ "Statistics",
      col1 == "p-value" ~ "Statistics",
      col1 == "Coeff. Variation" ~ "Statistics",
      col1 == "Protocol" ~ "Reporting",
      col1 == "Bias Control" ~ "Reporting",
      col1 == "Impact Factor" ~ "Article",
      col1 == "Year" ~ "Article",
      col1 == "Citations" ~ "Article",
      col1 == "Academic Age" ~ "Author",
      col1 == "Publications" ~ "Author",
      col1 == "University Ranking" ~ "Lab",
      col1 == "Survey - Replication" ~ "Prediction",
      col1 == "Survey - Effect Size" ~ "Prediction",
      col1 == "Survey - Challenge" ~ "Prediction",
      col1 == "Protocol Deviation\n(Lab)" ~ "Protocol Deviation",
      col1 == "Protocol Deviation\n(Committee)" ~ "Protocol Deviation",
      col1 == "Academic Age\n(Protocol Team)" ~ "Replication Team",
      col1 == "Publications\n(Protocol Team)" ~ "Replication Team",
      col1 == "Academic Age\n(Data Collection Team)" ~ "Replication Team",
      col1 == "Publications\n(Data Collection Team)" ~ "Replication Team",
      col1 == "University Ranking" ~ "Lab",
      .default = "Statistics"
    )) |>
    mutate(col1 = col1 |>
      fct_recode(
        "Institution\nRanking" = "Institution Ranking",
        "Survey\n(Replication)" = "Survey - Replication",
        "Survey\n(Effect Size)" = "Survey - Effect Size",
        "Survey\n(Challenge)" = "Survey - Challenge",
      )) |>
    drop_na() |>
    mutate(group_x = fct_relevel(
      group_x,
      "Model/Method",
      "Statistics",
      "Reporting",
      "Article",
      "Author",
      "Protocol Deviation",
      "Replication Team",
      "Center",
      # "Inst.",
      # "Lab",
      "Prediction"
    )) |>
    # Force correct order within each facet
    mutate(col1 = fct_relevel(
      col1,
      "Animal", "EPM", "MTT", "PCR",
      "Hedges's g", "p-value", "Coeff. Variation",
      "Protocol", "Bias Control",
      "Impact Factor", "Year", "Citations",
      "Academic Age", "Publications",
      "Protocol Deviation\n(Lab)", "Protocol Deviation\n(Committee)",
      "Academic Age\n(Protocol Team)", "Publications\n(Protocol Team)",
      "Academic Age\n(Data Collection Team)", "Publications\n(Data Collection Team)",
      "Institution\nRanking",
      "Survey\n(Replication)", 
      "Survey\n(Effect Size)", 
      "Survey\n(Challenge)"
    )) |>
    ggplot() +
    aes(
      x = col1,
      y = fct_rev(fct_relevel(
        col2,
        "t Value",
        "-log ES Ratio",
        "Replication Effect Size",
        "SMD Effect Size",
        "Original in replication's 95% PI"
      )),
      fill = sign(pearson_r) * -log10(pearson_p),
      label = label_p
    ) +
    geom_tile() +
    facet_grid(cols = ggplot2::vars(group_x), scales = "free_x", space = "free_x") +
    labs(x = "", y = "", fill = "-log10(p-value)") +
    scale_fill_gradient2(
      low = bri_color[["low"]],
      mid = bri_color[["mid"]],
      high = bri_color[["high"]],
      na.value = scales::alpha("white", 0)
    ) +
    bri_theme +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      title = element_text(size = 9),
      strip.background = element_rect(fill = "transparent"),
      strip.text = element_text(color = "black"), # center: 8
      panel.spacing.x = unit(0, "line"),
      panel.border = element_rect(color = "#333", fill = NA, size = 0.5),
      plot.margin = margin(t = 5, r = 150, b = 5, l = 5, unit = "pt")
    )

  if (show_label) {
    p <- p +
      geom_text(
        aes(label = label1),
        color = bri_color[["dark"]], size = 3.0
      ) +
      geom_text(
        aes(label = label2),
        color = bri_color[["dark"]], size = 2.5
      )

    p_p <- p_p +
      geom_text(
        aes(label = label1_p),
        color = bri_color[["dark"]], size = 3.0
      ) +
      geom_text(
        aes(label = label2_p),
        color = bri_color[["dark"]], size = 2.5
      )
  }

  p <- p +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0), labels = function(x) sub("^[0-9]+_", "", x))

  p_p <- p_p +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0))

  write_tsv(test_results |> select(-label1, -label2, -label1_p, -label2_p), paste0("./", fn, " cordata.tsv"))

  bri_ggsave(paste0("./", fn, " - spearman.png"),
    plot = p,
    width = max(1 + length(vars1) * 0.5, 8), height = max(1 + length(vars2) * 0.5, 5)
  )

  bri_ggsave(paste0("./", fn, " - pearson.png"),
    plot = p_p,
    width = max(1 + length(vars1) * 0.5, 8), height = max(1 + length(vars2) * 0.5, 5)
  )

  p
}

# Plots combined figure for presentation
plot_combined_cortable <- function(ct_exp, ct_rep, replication_datapath) {
  ct_exp <- ct_exp + theme(legend.position = "none") +
    labs(title = "Correlation between experiment features and replication success")

  ct_rep <- ct_rep + theme(
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    plot.margin = unit(c(0, 7.5, 3, 0), "cm")
  ) +
    labs(title = "Correlation between replication features and replication success")

  p <- plot_grid(
    plotlist = list(ct_exp, ct_rep),
    ncol = 1,
    labels = "AUTO",
    rel_heights = c(0.95, 1)
  ) +
    theme(plot.margin = unit(c(0, 0, -3, 0), "cm"))

  bri_ggsave(paste0(replication_datapath, "/predictors/cortable - predictors x replication - combined spearman.png"), plot = p, width = 10, height = 9)
}

# Plots combined figure for presentation
plot_combined_cortable_rep_outcome <- function(ct_exp, ct_rep, replication_datapath, half = TRUE) {
  ct_exp <- ct_exp +
    labs(title = "Correlation between replication criteria (experiment)") +
    scale_y_discrete(limits = c(
      "t Value",
      "-log ES Ratio",
      "≥50% subjectively replicated",
      "≥50% replications significant",
      "Same-sign significance (p<0.05)",
      "Replication in original 95% CI",
      "Original in replication's 95% PI"
    )) +
    scale_x_discrete(limits = c(
      "Original in replication's 95% PI",
      "Replication in original 95% CI",
      "Same-sign significance (p<0.05)",
      "≥50% replications significant",
      "≥50% subjectively replicated",
      "-log ES Ratio",
      "t Value"
    )) +
    theme(legend.position = "none")

  ct_rep <- ct_rep +
    labs(title = "Correlation between replication criteria (individual replication)") +
    scale_y_discrete(limits = c(
      "t Value",
      "-log ES Ratio",
      "Subjectively replicated",
      "Same-sign significance (p<0.05)",
      "Replication in original 95% CI"
    )) +
    scale_x_discrete(limits = c(
      "Replication in original 95% CI",
      "Same-sign significance (p<0.05)",
      "Subjectively replicated",
      "-log ES Ratio",
      "t Value"
    )) +
    theme(legend.position = "right")

  if (half == TRUE) {
    ct_exp <- ct_exp +
      scale_y_discrete(
        limits = c(
          "t Value",
          "-log ES Ratio",
          "≥50% subjectively replicated",
          "≥50% replications significant",
          "Same-sign significance (p<0.05)",
          "Replication in original 95% CI"
        ),
        expand = c(0, 0)
      ) +
      scale_x_discrete(
        limits = c(
          "Original in replication's 95% PI",
          "Replication in original 95% CI",
          "Same-sign significance (p<0.05)",
          "≥50% replications significant",
          "≥50% subjectively replicated",
          "-log ES Ratio"
        ),
        expand = c(0, 0)
      )

    ct_exp$data <- ct_exp$data |>
      mutate(col1 = fct_relevel(
        col1,
        "Original in replication's 95% PI",
        "Replication in original 95% CI",
        "Same-sign significance (p<0.05)",
        "≥50% replications significant",
        "≥50% subjectively replicated",
        "-log ES Ratio"
      )) |>
      arrange(col1) |>
      distinct(rho, .keep_all = TRUE)

    ct_rep <- ct_rep +
      scale_y_discrete(
        limits = c(
          "t Value",
          "-log ES Ratio",
          "Subjectively replicated",
          "Same-sign significance (p<0.05)"
        ),
        expand = c(0, 0)
      ) +
      scale_x_discrete(
        limits = c(
          "Replication in original 95% CI",
          "Same-sign significance (p<0.05)",
          "Subjectively replicated",
          "-log ES Ratio"
        ),
        expand = c(0, 0)
      )

    ct_rep$data <- ct_rep$data |>
      mutate(col1 = fct_relevel(
        col1,
        "Replication in original 95% CI",
        "Same-sign significance (p<0.05)",
        "Subjectively replicated",
        "-log ES Ratio"
      )) |>
      arrange(col1) |>
      distinct(rho, .keep_all = TRUE)


    ct_rep <- ct_rep +
      scale_fill_gradient2(
        low = bri_color[["low"]],
        mid = bri_color[["mid"]],
        high = bri_color[["high"]],
        na.value = bri_color[["none"]],
        breaks = c(0, 3, 5, 7),
        labels = c(expression(1), expression(10^-3), expression(10^-5), expression(10^-7))
      ) +
      labs(fill = "p-value")



    p <- plot_grid(
      plotlist = list(ct_exp + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")), ct_rep + theme(legend.position = "right", plot.margin = unit(c(0, 0, 2, 0), "cm"))),
      ncol = 1,
      labels = "AUTO",
      rel_heights = c(1, 1)
    ) +
      theme(
        plot.margin = unit(c(0.1, 0.2, -2, 0), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  } else {
    p <- plot_grid(
      plotlist = list(ct_exp + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")), ct_rep + theme(legend.position = "right", plot.margin = unit(c(0, 0, 2, 0), "cm"))),
      ncol = 1,
      labels = "AUTO",
      rel_heights = c(1, 1)
    ) +
      theme(plot.margin = unit(c(0.1, 0.2, -2, 0), "cm"))
  }

  bri_ggsave(paste0(replication_datapath, "/predictors/cortable - replication - combined spearman.png"), plot = p, width = 6, height = 8)
}

# Plots combined figure for presentation
plot_combined_cortable_predictor <- function(ct_exp, ct_rep, replication_datapath, half = FALSE) {
  ct_exp <- ct_exp +
    labs(title = "Correlation between experiment-level predictors") +
    theme(legend.position = "none")

  ct_rep <- ct_rep +
    labs(title = "Correlation between replication-level predictors") +
    theme(legend.position = "right", legend.key.height = unit(0.8, "cm")) +
    scale_fill_gradient2(
      low = bri_color[["low"]],
      mid = bri_color[["mid"]],
      high = bri_color[["high"]],
      na.value = scales::alpha("white", 0),
      breaks = c(-3, -2, -1, 0, 1, 2, 3),
      limits = c(-3, 3),
      labels = c(expression(10^-3), expression(10^-2), expression(10^-1), expression(1), expression(10^-1), expression(10^-2), expression(10^-3))
    )
  ct_exp$data <- ct_exp$data |>
    mutate(
      label2 = if_else(is.na(rho), "", label2),
      label1 = if_else(is.na(rho), "", label1), # NA
      capped.p.value = ifelse(p.value < 0.001, 0.001, p.value),
      fill_value = sign(rho) * -log10(capped.p.value) # Update fill value
    )

  ct_rep$data <- ct_rep$data |>
    mutate(
      capped.p.value = ifelse(p.value < 0.001, 0.001, p.value),
      fill_value = sign(rho) * -log10(capped.p.value) # Update fill value
    )

  ct_exp <- ct_exp +
    aes(fill = fill_value)

  ct_rep <- ct_rep +
    aes(fill = fill_value) +
    labs(fill = "p-value")

  if (half == TRUE) {
    ct_exp$data <- ct_exp$data |>
      arrange(col1) |>
      distinct(rho, .keep_all = TRUE)

    ct_rep$data <- ct_rep$data |>
      arrange(col1) |>
      distinct(rho, .keep_all = TRUE)

    p <- plot_grid(
      plotlist = list(ct_exp + scale_fill_gradient2(
        low = bri_color[["low"]],
        mid = bri_color[["mid"]],
        high = bri_color[["high"]],
        na.value = scales::alpha("white", 0)
      ) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ct_rep + theme(
        legend.position = "right",
        plot.margin = unit(c(0, 9.5, 10, 0), "cm")
      )),
      ncol = 1,
      labels = "AUTO",
      rel_heights = c(1, 0.9)
    ) +
      theme(plot.margin = unit(c(0.1, 0.2, -10, 0), "cm"))
  } else {
    p <- plot_grid(
      plotlist = list(ct_exp + scale_fill_gradient2(
        low = bri_color[["low"]],
        mid = bri_color[["mid"]],
        high = bri_color[["high"]],
        na.value = scales::alpha("white", 1)
      ) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ct_rep + theme(legend.position = "right", plot.margin = unit(c(0, 9.5, 10, 0), "cm"))),
      ncol = 1,
      labels = "AUTO",
      rel_heights = c(1, 0.9)
    ) +
      theme(plot.margin = unit(c(0.1, 0.2, -10, 0), "cm"))
  }

  bri_ggsave(paste0(replication_datapath, "/predictors/cortable - predictors - combined spearman.png"), plot = p, width = 10, height = 14)
}

# Plot single replication results (used for individual labs to check)
plot_individual_results <- function(fn, es_label, output_path) {
  EXP_code <- fn |> str_extract("(EPM|ALTMTT|MTT|ALTPCR|PCR)[0-9]+")
  LAB_code <- fn |> str_extract("LAB[0-9]+")

  method <- EXP_code |> str_extract("(EPM|MTT|PCR)")
  altmethod <- EXP_code |> str_extract("(EPM|ALTMTT|MTT|ALTPCR|PCR)")

  if (method == "MTT") {
    meta_table <- meta_data_MTT
  } else if (method == "PCR") {
    meta_table <- meta_data_PCR
  } else if (method == "EPM") {
    meta_table <- meta_data_EPM
  }

  group1_name <- meta_table |>
    filter(EXP == EXP_code |> str_remove("ALT"), LAB == LAB_code) |>
    pull(Group1_name) |>
    unique() |>
    stri_trans_general("Latin-ASCII") |>
    str_replace("---", "\n")

  group2_name <- meta_table |>
    filter(EXP == EXP_code |> str_remove("ALT"), LAB == LAB_code) |>
    pull(Group2_name) |>
    unique() |>
    stri_trans_general("Latin-ASCII") |>
    str_replace("---", "\n")

  if (method == "EPM") {
    outcome_name <- meta_table |>
      filter(EXP == EXP_code, LAB == LAB_code) |>
      pull(outcome_name_short) |>
      unique() |>
      stri_trans_general("Latin-ASCII") |>
      str_replace("---", "\n")
  } else if (altmethod == "ALTPCR") {
    outcome_name <- "Relative Expression (% control)"
  } else {
    outcome_name <- meta_table |>
      filter(EXP == EXP_code |> str_remove("ALT"), LAB == LAB_code) |>
      pull(outcome_name) |>
      unique() |>
      stri_trans_general("Latin-ASCII") |>
      str_replace("---", "\n")

    if (str_length(outcome_name > 40)) {
      outcome_name <- str_replace(outcome_name, "\\(", "\n(")
    }
  }

  if (method == "MTT" | altmethod == "ALTPCR") {
    idf <- read_tsv(fn, show_col_types = F) |>
      select(Group1_Perc, Group2_Perc) |>
      `colnames<-`(c(group1_name, group2_name)) |>
      pivot_longer(cols = everything()) |>
      mutate(name = factor(name, levels = c(group1_name, group2_name))) |>
      filter(!is.na(value))
  } else {
    idf <- read_tsv(fn, show_col_types = F) |>
      select(Group1, Group2) |>
      `colnames<-`(c(group1_name, group2_name)) |>
      pivot_longer(cols = everything()) |>
      mutate(name = factor(name, levels = c(group1_name, group2_name))) |>
      filter(!is.na(value))
  }

  make_plot <- function(idf, suffix = "", note_label = "", scale_idf = NULL, perc_scale) {
    if (is.null(scale_idf)) {
      scale_idf <- idf
    }

    idfs <- scale_idf |>
      group_by(name) |>
      summarise(
        m = mean(value, na.rm = T),
        s = sd(value, na.rm = T),
        maxms = m + s,
        minms = m - s
      )

    p <- ggplot(
      data = idf,
      mapping = aes(x = name, y = value)
    ) +
      geom_point(size = 2.5, alpha = 0.4, position = position_jitter(width = 0.1, height = 0)) +
      geom_errorbar(data = idfs, mapping = aes(x = name, y = m, ymin = m - s, ymax = m + s), color = bri_color[["dark"]], linewidth = 0.68, width = 0.1) +
      geom_segment(data = idfs, mapping = aes(x = as.numeric(name) - 0.085, xend = as.numeric(name) + 0.085, y = m, yend = m), color = bri_color[["dark"]], linewidth = 1.03) +
      labs(
        x = "", y = outcome_name, title = paste0(LAB_code, " ", EXP_code)
        # caption = paste0(es_label, "\n", note_label, "\nPlotted: mean ± SD")
      ) +
      bri_theme +
      theme(
        legend.position = "none",
        plot.caption = element_text(hjust = 0)
      )

    if (perc_scale) {
      y_limits <- c(
        min(-10, min(c(idf$value, idfs$minms) - 10, na.rm = T)),
        max(110, max(c(idf$value, idfs$maxms) + 10, na.rm = T))
      )
    } else {
      y_limits <- c(
        min(c(idf$value, idfs$minms), na.rm = T),
        max(c(idf$value, idfs$maxms), na.rm = T)
      )
    }

    p <- p +
      scale_y_continuous(limits = y_limits, breaks = scales::pretty_breaks(n = 5))

    dir.create(paste0(output_path, "/individual-plots"))
    bri_ggsave(paste0(output_path, "/individual-plots/", LAB_code, " ", EXP_code, suffix, ".png"),
      plot = p,
      width = 6, height = 4
    )
  }

  make_plot(idf, perc_scale = (altmethod != "PCR"))

  lb <- median(idf$value, na.rm = F) - 5 * IQR(idf$value)
  ub <- median(idf$value, na.rm = F) + 5 * IQR(idf$value)
  full_n <- nrow(idf)
  full_idf <- idf
  idf <- idf |> filter(value > lb, value < ub)
  if (nrow(idf) < full_n) {
    make_plot(idf, " no outliers", "* Not showing outliers (outside of median +- 5 * IQR)", scale_idf = full_idf, perc_scale = (method != "PCR"))
  }
}

# Make forest plots
forest_plot <- function(plot_data, EXP_code, plot_label, output_path, scale_to_use) {
  n_reps <- plot_data$LAB |>
    str_detect("LAB") |>
    sum()

  plot_data <- plot_data |> mutate(
    type = as.character(LAB),
    type = ifelse(str_detect(type, "LAB"), "REP", type)
  )

  if (scale_to_use == "ROM") {
    plot_data <- plot_data |>
      mutate(
        yi = 100 * exp(yi),
        lb = 100 * exp(lb),
        ub = 100 * exp(ub)
      )
  } else if (scale_to_use == "MD") {
    plot_data <- rbind(
      plot_data |>
        filter(type == "Original Effect") |>
        mutate(
          yi = 100 * exp(yi),
          lb = 100 * exp(lb),
          ub = 100 * exp(ub)
        ),
      plot_data |>
        filter(type != "Original Effect") |>
        mutate(
          yi = 100 * (yi),
          lb = 100 * (lb),
          ub = 100 * (ub)
        )
    )
  }

  if (scale_to_use %in% c("ROM", "MD")) {
    plot_data$es_label <- paste0(
      round(plot_data$yi, 0), "%",
      " [", round(plot_data$lb, 0), "%, ",
      round(plot_data$ub, 0), "%]"
    )
  } else {
    plot_data$es_label <- paste0(
      round(plot_data$yi, 1),
      " [", round(plot_data$lb, 1), ", ",
      round(plot_data$ub, 1), "]"
    )
  }

  labs <- plot_data$LAB[str_which(plot_data$LAB, "LAB")] |> as.character()
  plot_data$LAB <- factor(
    plot_data$LAB,
    levels = c("Original Effect", "RE Model PI", "RE Model CI", "FE Model CI", labs)
  )

  original_ES <- plot_data |>
    filter(LAB == "Original Effect") |>
    pull(yi)

  x_axis_scale <- max(plot_data$ub) - min(plot_data$lb)

  if (scale_to_use %in% c("ROM", "MD")) {
    label_pos <- max(plot_data$ub) * 1.15
    x_axis_min <- min(-10, min(plot_data$lb))
    x_axis_max <- label_pos + x_axis_scale * 0.5
  } else {
    xsf <- 3 / x_axis_scale
    label_pos <- max(plot_data$ub) + 0.2 * xsf
    x_axis_min <- min(plot_data$lb) - 0.2 * xsf
    x_axis_max <- label_pos + 1 * xsf
  }

  color_vector <- plot_data |>
    arrange(type) |>
    pull(type) |>
    dplyr::recode(
      "Original Effect" = bri_color[["dark"]],
      "RE Model PI" = bri_color[["main"]],
      "RE Model CI" = bri_color[["main"]],
      "FE Model CI" = bri_color[["main"]],
      "REP" = bri_color[["dark"]]
    )

  p <- ggplot(plot_data) +
    aes(
      x = yi, y = LAB, xmin = lb, xmax = ub, color = type,
      shape = type, size = weight, label = es_label
    ) +
    geom_vline(xintercept = original_ES, linetype = "dashed", color = bri_color[["light"]], linewidth = 0.7) +
    geom_errorbarh(linewidth = 0.7, height = 0.2, color = bri_color[["dark"]]) +
    geom_point(shape = 15) +
    geom_textbox(
      aes(x = label_pos),
      fill = NA, box.colour = NA,
      size = 3.2,
      hjust = 0, halign = 0
    ) +
    scale_x_continuous(breaks = pretty_breaks(), limits = c(x_axis_min, x_axis_max)) +
    scale_size(range = c(1, 4)) +
    scale_color_manual(values = color_vector) +
    bri_theme +
    theme(legend.position = "none", axis.text.y = element_text(size = 11), panel.grid = element_blank())

  if (scale_to_use %in% c("ROM", "MD")) {
    p <- p +
      geom_vline(xintercept = 100, linetype = "dashed", color = bri_color[["mid"]], linewidth = 0.7) +
      labs(title = EXP_code, x = "Mean Percentage of Control, from Log Ratio of Means", y = "", caption = plot_label)
  } else {
    p <- p +
      geom_vline(xintercept = 0, linetype = "dashed", color = bri_color[["mid"]], linewidth = 0.7) +
      labs(title = EXP_code, x = "Log Transformed Ratio of Means", y = "", caption = plot_label)
  }

  h <- ifelse(n_reps == 1, 1, n_reps + 1)

  dir.create(paste0(output_path, "/forest-plots"))
  bri_ggsave(paste0(output_path, "/forest-plots/", EXP_code, ".png"), plot = p, width = 6, height = 2.2 + h * 0.5)

  p
}

# Build a table of main numbers to be reported
make_general_numbers_table <- function(inclusion_sets) {
  bri_IS <- inclusion_sets |>
    filter(UNIT == "BRI")

  # Flow of the 180 planned replications (simple, for humans)
  exps_not_done <- bri_IS |>
    filter(done == "No") |>
    count(reason_for_not_done) |>
    ungroup() |>
    mutate(reason_for_not_done = paste0("Not done for: ", reason_for_not_done)) |>
    pivot_wider(names_from = reason_for_not_done, values_from = n)

  exps_done <- bri_IS |>
    filter(done == "Yes") |>
    count(`Validation - Excluded by`) |>
    ungroup() |>
    mutate(`Validation - Excluded by` = paste0("Invalid because: ", `Validation - Excluded by`)) |>
    pivot_wider(names_from = `Validation - Excluded by`, values_from = n) |>
    rename(`Valid experiments` = `Invalid because: NA`)

  IS_all_exps <- bri_IS |> filter(analysis_all_exps_lab_units != "")
  IS_primary <- bri_IS |> filter(analysis_primary != "")

  df <- cbind(
    tibble(Planned = 180),
    exps_not_done,
    exps_done
  ) |> pivot_longer(cols = everything())

  # Flow of the 180 planned replications (with separate columns)
  exps_not_done <- bri_IS |>
    filter(done == "No") |>
    count(reason_for_not_done) |>
    ungroup() |>
    rename(reason = reason_for_not_done) |>
    mutate(status = "Not done", validation = NA)

  exps_done <- bri_IS |>
    filter(done == "Yes") |>
    count(`Validation - Excluded by`) |>
    ungroup() |>
    rename(reason = `Validation - Excluded by`) |>
    mutate(status = "Done")

  exps_done <- rbind(
    exps_done |> filter(!is.na(reason)) |> mutate(validation = "Invalid"),
    exps_done |> filter(is.na(reason)) |> mutate(validation = "Valid")
  )

  dfw <- rbind(
    tibble(status = "Planned", validation = NA, reason = NA, n = 180),
    exps_not_done,
    exps_done
  )

  # Experiments by method and number of copies (all exps)
  copy_count_all <- IS_all_exps |>
    mutate(Method = str_extract(EXP, "(MTT|EPM|PCR)")) |>
    count(Method, EXP) |>
    ungroup() |>
    mutate(Copies = factor(n)) |>
    count(Method, Copies) |>
    ungroup() |>
    pivot_wider(id_cols = Method, names_from = Copies, values_from = n, values_fill = 0) |>
    mutate(Total = `1` + `2` + `3`)

  # Experiments by method and number of copies (primary analysis)
  copy_count_primary <- IS_primary |>
    mutate(Method = str_extract(EXP, "(MTT|EPM|PCR)")) |>
    count(Method, EXP) |>
    ungroup() |>
    mutate(Copies = factor(n)) |>
    count(Method, Copies) |>
    ungroup() |>
    pivot_wider(id_cols = Method, names_from = Copies, values_from = n, values_fill = 0) |>
    mutate(Total = `1` + `2` + `3`)

  paste0("output/", results_path, "/Replication Rate Aggregate.tsv")
  write_tsv(dfw, file = paste0("output/", results_path, "/Experiment Numbers Overview - wide.tsv"))
}

# Plots a single scatterplot of predictor/outcome
plot_scatterplot_predictors <- function(df, corr, c1, c2, c1_cat, c2_cat, fn) {
  dir.create(paste0("./", fn, "/scatterplots"))
  c1_lab <- dplyr::recode(
    c1,
    "Log ES Ratio" = "-logES ratio",
    "Original in replication PI" = "Original in replication's 95% PI",
    "Replication in original CI" = "Replication in original 95% CI",
    "Same sense significance" = "Same-sign significance (p<0.05)",
    "Voting (with ties)" = "≥50% replications significant",
    "Subjective (with ties)" = "≥50% subjectively replicated",
    "Subjective" = "≥50% subjectively replicated"
  )

  c2_lab <- dplyr::recode(
    c2,
    "Log ES Ratio" = "-logES ratio",
    "Original in replication PI" = "Original in replication's 95% PI",
    "Replication in original CI" = "Replication in original 95% CI",
    "Same sense significance" = "Same-sign significance (p<0.05)",
    "Voting (with ties)" = "≥50% replications significant",
    "Subjective (with ties)" = "≥50% subjectively replicated",
    "Subjective" = "≥50% subjectively replicated"
  )

  # Conditional plotting based on whether c1 or c2 are categorical
  if (c1_cat & c2_cat) {
    # If both c1 and c2 are categorical: 2x2 table
    table2x2 <- df |>
      mutate(
        vx = factor(vx),
        vy = factor(vy)
      ) |>
      count(vx, vy) |>
      mutate(
        vx = if_else(vx == 0, "No", "Yes"),
        vy = if_else(vy == 0, "No", "Yes")
      ) |>
      pivot_wider(names_from = vy, values_from = n, values_fill = 0) |>
      rename_with(~c1_lab, vx) |>
      flextable() |>
      add_header_row(
        values = c("", c2_lab),
        colwidths = c(1, 2),
        top = TRUE
      ) |>
      add_footer_lines(
        paste0("Spearman's ρ = ", signif(unname(corr$estimate), 2), " | n = ", nrow(df), " | p = ", signif(corr$p.value, 2))
      ) |>
      set_table_properties(layout = "autofit")

    save_tbl(list(table2x2), paste0("./", fn, "/scatterplots/", c1, " x ", c2, ".docx"), date_hour = F, open_folder = F)
  } else if (c1_cat) {
    # If only c1 is categorical, use a scatterplot 1D for c1 (categorical) vs vy (continuous)
    p <- df |>
      mutate(vx = factor(vx, labels = c("No", "Yes"))) |>
      ggplot() +
      aes(x = vx, y = vy, color = Method) + # x based on c1 (categorical) and y based on vy (continuous)
      geom_point(alpha = 0.6, size = 2, position = position_jitter(width = 0.1, height = 0)) +
      scale_color_manual(values = c(bri_color[["epm"]], bri_color[["mtt"]], bri_color[["pcr"]])) +
      labs(
        x = c1_lab, y = c2_lab,
        caption = paste0("Spearman's ρ = ", signif(unname(corr$estimate), 2), " | n = ", nrow(df), " | p = ", signif(corr$p.value, 2))
      ) +
      bri_theme

    bri_ggsave(paste0("./", fn, "/scatterplots/", c1, " x ", c2, ".png"), plot = p, width = 6, height = 4)
  } else if (c2_cat) {
    # If only c2 is categorical, use a boxplot for c2 (categorical) vs vx (continuous)
    p <- df |>
      mutate(vy = factor(vy, labels = c("No", "Yes"))) |>
      ggplot() +
      aes(x = vx, y = factor(vy), color = Method) + # x based on c2 (categorical) and y based on vx (continuous)
      geom_point(alpha = 0.6, size = 2, position = position_jitter(width = 0, height = 0.1)) +
      scale_color_manual(values = c(bri_color[["epm"]], bri_color[["mtt"]], bri_color[["pcr"]])) +
      labs(
        x = c1_lab, y = c2_lab,
        caption = paste0("Spearman's ρ = ", signif(unname(corr$estimate), 2), " | n = ", nrow(df), " | p = ", signif(corr$p.value, 2))
      ) +
      bri_theme

    bri_ggsave(paste0("./", fn, "/scatterplots/", c1, " x ", c2, ".png"), plot = p, width = 6, height = 4)
  } else {
    # If both are continuous, use a scatter plot
    p <- ggplot(df) +
      aes(x = vx, y = vy, color = Method) +
      scale_color_manual(values = c(bri_color[["epm"]], bri_color[["mtt"]], bri_color[["pcr"]])) +
      labs(
        x = c1_lab, y = c2_lab,
        caption = paste0("Spearman's ρ = ", signif(unname(corr$estimate), 2), " | n = ", nrow(df), " | p = ", signif(corr$p.value, 2))
      ) +
      bri_theme +
      geom_point(size = 2, alpha = 0.6)

    bri_ggsave(paste0("./", fn, "/scatterplots/", c1, " x ", c2, ".png"), plot = p, width = 6, height = 4)
  }
}

# Plots the effect sizes in order to compare variation between experiments
plot_icc_matrix <- function(m, fn, cap_at = NULL, tt = NULL) {
  # Enhanced error handling for required columns
  if (!"EXP" %in% colnames(m)) {
    stop(paste0("[ERROR] plot_icc_matrix: Input matrix is missing required column 'EXP'. Columns present: ", paste(colnames(m), collapse = ", ")))
  }
  
  # Pivot data to long format
  m_long <- m |>
    pivot_longer(cols = -"EXP") |>
    mutate(
      Technique = EXP |> str_extract("(MTT|EPM|PCR)"),
      Type = ifelse(name == "orig", "Original", "Replication"),
      TypeTechnique = ifelse(Type == "Original", Type, Technique)
    )
  
  # Get the desired order based on the 'Original' effect size (reversed for Y-axis display)
  # For Y-axis: arrange desc() to show negative values at top, positive at bottom
  exp_order <- m_long |>
    filter(Type == "Original") |>
    arrange(desc(value)) |>
    pull(EXP)
  
  # Apply this order to the main dataframe by converting EXP to a factor
  m_ordered <- m_long |>
    mutate(EXP = factor(EXP, levels = exp_order))
  
  # Filter out rows not needed for the plot
  m_plot_data <- m_ordered |>
    filter(name != "avg", !is.na(value))
  
  if (!is.null(cap_at)) {
    m_plot_data <- m_plot_data |> filter(value <= cap_at)
  }
  
  m_plot_data$TypeTechnique <- factor(m_plot_data$TypeTechnique, levels = c("Original", "EPM", "MTT", "PCR"))
  
  # --- REVISED GGPLOT CALL ---
  p <- ggplot(m_plot_data, aes(y = EXP, x = value)) + # 'color' aesthetic removed from here
    geom_vline(xintercept = 0, linetype = "dashed", color = bri_color[["dark"]]) +
    
    # Layer 1: Original points, with its own color aesthetic
    geom_point(data = . %>% filter(Type == "Original"), aes(color = TypeTechnique), size = 3.5, alpha = 0.7) +
    scale_color_manual(
      values = c("Original" = bri_color[["dark"]]),
      name = "" # No legend title for the 'Original' points
    ) +
    
    # Reset the color scale for the next layer
    ggnewscale::new_scale_color() +
    
    # Layer 2: Replication points, with its own color aesthetic
    geom_point(data = . %>% filter(Type != "Original"), aes(color = TypeTechnique), size = 3.5, alpha = 0.7) +
    scale_color_manual(
      values = c(
        "EPM" = bri_color[["epm"]],
        "MTT" = bri_color[["mtt"]],
        "PCR" = bri_color[["pcr"]]
      ),
      name = "Replication:" # Legend title for the replication points
    ) +
    
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    facet_wrap(~Technique, nrow = 1, scales = "free") +
    labs(
      y = "",
      x = "Effect size (log ratio of means, normalized, ordered by magnitude)",
      title = tt
    ) +
    bri_theme +
    # Ensure the legend is positioned correctly
    theme(legend.position = "bottom", legend.title = element_text(size = rel(0.8)))
  
  bri_ggsave(fn, plot = p, width = 14, height = 10)
  
  p
}



# Plots a correlation between original and replication effect size
plot_effect_correlation <- function(rep_summary_folder, df, fn, tt) {
  df <- df |>
    mutate(
      Method = str_extract(EXP, "[A-Z]{3}"),
      original_es = log(exp(original_es)),
      replication_es = log(exp(replication_es))
    )

  r1 <- cor.test(x = df$original_es, y = df$replication_es, na.action = na.omit(), method = "pearson")
  r2 <- cor.test(x = df$original_es, y = df$replication_es, na.action = na.omit(), method = "spearman")

  valid_idx <- !is.na(df$original_es) & !is.na(df$replication_es)
  x_clean <- df$original_es[valid_idx]
  y_clean <- df$replication_es[valid_idx]

  spearman_ci <- tryCatch({
    sci <- spearmanCI::spearmanCI(x_clean, y_clean, level = 0.95, plot = FALSE)
    c(sci[1], sci[2])
  }, error = function(e) {
    z <- 0.5 * log((1 + r2$estimate) / (1 - r2$estimate))
    se <- 1 / sqrt(max(length(x_clean) - 3, 1))
    c(tanh(z - 1.96 * se), tanh(z + 1.96 * se))
  })

  r1_ci_str <- paste0("[95% CI: ", round(r1$conf.int[1], 2), ", ", round(r1$conf.int[2], 2), "]")
  r2_ci_str <- paste0("[95% CI: ", round(spearman_ci[1], 2), ", ", round(spearman_ci[2], 2), "]")

  rcor_label <- paste0(
    "Pearson's R = ", round(r1$estimate, 2), " ", r1_ci_str, " (p = ", signif(r1$p.value, 2), ")\n",
    "Spearman's ρ = ", round(r2$estimate, 2), " ", r2_ci_str, " (p = ", signif(r2$p.value, 2), ")"
  )

  p <- ggplot(df) +
    aes(x = original_es, y = replication_es, color = Method) +
    geom_point(size = 3, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.7, color = bri_color[["mid"]]) +
    annotate("text", label = rcor_label, x = Inf, y = -Inf, vjust = -0.5, hjust = 1.05, size = 3) +
    scale_color_manual(values = c(bri_color[["epm"]], bri_color[["mtt"]], bri_color[["pcr"]])) +
    labs(x = "Original Effect Size", y = "Replication Effect Size", title = tt, color = "Method:") +
    bri_theme +
    theme(legend.title = element_text(size = rel(0.8)))

  bri_ggsave(paste0(rep_summary_folder, "/additional-figures/", fn),
    plot = p,
    width = 6, height = 6
  )

  p
}

# Plot the ratios of replications and original effect sizes
plot_log_es_ratios <- function(df, fn, rep_summary_folder) {
  df <- df |>
    mutate(
      Method = str_extract(EXP, "[A-Z]{3}")
    )

  p <- ggplot(df) +
    aes(x = Method, y = log_es_ratio) +
    geom_boxplot(width = 0.5) +
    geom_point(size = 2, alpha = 0.7, position = position_jitter(width = 0.15)) +
    labs(x = "Method", y = "log ES ratio") +
    bri_theme

  bri_ggsave(paste0(rep_summary_folder, "/additional-figures/", fn),
    plot = p,
    width = 6, height = 6
  )

  p
}

# Plot combined figure of effect sizes for presentation
plot_effect_size_figure <- function(p_A, p_B, p_C, rep_summary_folder) {
  top_plots <- plot_grid(plotlist = list(p_B, p_C), labels = c("B", "C"), nrow = 1)

  p <- plot_grid(plotlist = list(p_A, top_plots), labels = c("A", ""), nrow = 2)

  bri_ggsave(paste0(rep_summary_folder, "/additional-figures/", "effect_size_figure_combined.png"), plot = p, width = 14, height = 10)

  p
}

# Plot of all coefficients of variation, for replications and original, flagging the ones that are from badly estimated SDs, in many variations
plot_cvs <- function(rep_df, orig_df, suffix, rep_summary_folder) {
  dir.create(paste0(rep_summary_folder, "/CV plots"))

  # Extract full method label BEFORE collapsing EXP for joins (keeps ALTMTT/ALTPCR distinct)
  rep_df <- rep_df |>
    mutate(
      Method = str_extract(EXP, "(EPM|ALTMTT|MTT|ALTPCR|PCR)"),
      # Map ALT variants to base technique names for display (ALTPCR -> PCR, ALTMTT -> MTT)
      Method = case_when(
        Method == "ALTPCR" ~ "PCR",
        Method == "ALTMTT" ~ "MTT",
        TRUE ~ Method
      ),
      # Collapse ALT from EXP only for joining with original data and x-axis labeling
      EXP = str_remove(EXP, "ALT")
    )
  
  # Set order of techniques to match Fig. 2 (specification curve): EPM, MTT, PCR
  rep_df$Method <- factor(rep_df$Method, levels = c("EPM", "MTT", "PCR"))

  by_rep <- "LAB" %in% colnames(rep_df)

  if (by_rep) {
    df <- rep_df |>
      left_join(orig_df |> select(EXP, `SD Estimate`)) |>
      rename(sd_estimate = `SD Estimate`) |>
      select(LAB, EXP, Method, original_cv, replication_cv, sd_estimate) |>
      pivot_longer(cols = -c(LAB, EXP, Method, sd_estimate)) |>
      rename(cv = value) |>
      mutate(
        type = ifelse(name == "original_cv", "Original", "Replications"),
        # Map ALT variants to base techniques for coloring/legend, but keep Method for faceting
        TypeTechnique = ifelse(
          type == "Original",
          type,
          case_when(
            Method %in% c("ALTPCR", "PCR") ~ "PCR",
            Method %in% c("ALTMTT", "MTT") ~ "MTT",
            Method == "EPM" ~ "EPM",
            TRUE ~ Method
          )
        ),
        sd_estimate_ok = ifelse(sd_estimate != "reported SEM and range for N", "Reported SD", "Estimated SD"),
        flag = ifelse(
          type == "Original",
          paste0(type, "/\n", sd_estimate_ok),
          "Replications"
        )
      ) |>
      select(-name)
  } else {
    df <- rep_df |>
      left_join(orig_df |> select(EXP, `SD Estimate`)) |>
      rename(sd_estimate = `SD Estimate`) |>
      select(EXP, Method, original_cv, replication_cv, sd_estimate) |>
      pivot_longer(cols = -c(EXP, Method, sd_estimate)) |>
      rename(cv = value) |>
      mutate(
        type = ifelse(name == "original_cv", "Original", "Replications"),
        # Map ALT variants to base techniques for coloring/legend, but keep Method for faceting
        TypeTechnique = ifelse(
          type == "Original",
          type,
          case_when(
            Method %in% c("ALTPCR", "PCR") ~ "PCR",
            Method %in% c("ALTMTT", "MTT") ~ "MTT",
            Method == "EPM" ~ "EPM",
            TRUE ~ Method
          )
        ),
        sd_estimate_ok = ifelse(sd_estimate != "reported SEM and range for N", "Reported SD", "Estimated SD"),
        flag = ifelse(
          type == "Original",
          paste0(type, "/\n", sd_estimate_ok),
          "Replications"
        )
      ) |>
      select(-name)
  }

  df1 <- df |>
    filter(type == "Original") |>
    mutate(LAB = "Original") |>
    distinct(EXP, .keep_all = T)
  df2 <- df |> filter(type == "Replications")

  if (!by_rep) {
    df2 <- df2 |>
      mutate(LAB = "Mean of Replications") |>
      distinct(EXP, .keep_all = T)
  }

  df <- rbind(df1, df2)
  
  # Ensure Method factor order is preserved after rbind (EPM, MTT, PCR to match Fig. 2)
  df$Method <- factor(df$Method, levels = c("EPM", "MTT", "PCR"))

  # Criar ordenação baseada no CV original (maior para menor) para que apareça menor para maior após coord_flip()
  exp_order_by_cv <- df |>
    filter(type == "Original") |>
    arrange(desc(cv)) |>
    pull(EXP) |>
    unique()
  
  df$EXP <- factor(df$EXP, levels = exp_order_by_cv)

  df$TypeTechnique <- factor(df$TypeTechnique, levels = c("Original", "EPM", "MTT", "PCR"))

  df <- df |>
    mutate(
      technique_flag = case_when(
        flag == "Original/\nEstimated SD" ~ "Original/\nEstimated SD",
        flag == "Original/\nReported SD" ~ "Original/\nReported SD",
        TypeTechnique == "EPM" ~ "EPM",
        TypeTechnique == "MTT" ~ "MTT",
        TypeTechnique == "PCR" ~ "PCR"
      ),
      technique_flag = factor(technique_flag, levels = c(
        "Original/\nEstimated SD",
        "Original/\nReported SD",
        "EPM",
        "MTT",
        "PCR"
      ))
    )

  p <- ggplot(df) +
    aes(x = EXP, y = cv, color = technique_flag, shape = technique_flag) +
    geom_point(aes(color = technique_flag, shape = technique_flag), size = 2.5, alpha = 0.8, data = df |> filter(str_detect(technique_flag, "Original"))) +
    scale_color_manual(
      values = c(
        "Original/\nEstimated SD" = bri_color[["dark"]],
        "Original/\nReported SD" = bri_color[["dark"]]
      ),
      labels = c(
        "Original/\nEstimated SD" = "Original/\nEstimated SD",
        "Original/\nReported SD" = "Original/\nReported SD"
      ),
      name = "",
      guide = guide_legend(order = 1)
    ) +
    scale_shape_manual(
      values = c(
        "Original/\nEstimated SD" = 17,
        "Original/\nReported SD" = 16
      ),
      labels = c(
        "Original/\nEstimated SD" = "Original/\nEstimated SD",
        "Original/\nReported SD" = "Original/\nReported SD"
      ),
      name = "",
      guide = guide_legend(order = 1)
    ) +
    ggnewscale::new_scale("color") +
    ggnewscale::new_scale("shape") +
    geom_point(aes(color = technique_flag, shape = technique_flag), size = 2.5, alpha = 0.8, data = df |> filter(!str_detect(technique_flag, "Original"))) +
    geom_blank(data = data.frame(EXP = "EPM1", cv = 0, Method = "EPM", technique_flag = "EPM"), aes(x = EXP, y = cv), inherit.aes = FALSE) +
    scale_color_manual(
      values = c(
        "EPM" = bri_color[["epm"]],
        "MTT" = bri_color[["mtt"]],
        "PCR" = bri_color[["pcr"]]
      ),
      labels = c(
        "EPM" = "EPM",
        "MTT" = "MTT",
        "PCR" = "PCR"
      ),
      name = "Replication:",
      guide = guide_legend(order = 2)
    ) +
    scale_shape_manual(
      values = c(
        "EPM" = 16,
        "MTT" = 16,
        "PCR" = 16
      ),
      labels = c(
        "EPM" = "EPM",
        "MTT" = "MTT",
        "PCR" = "PCR"
      ),
      name = "Replication:",
      guide = guide_legend(order = 2)
    ) +
    facet_wrap(~Method, scales = "free") +
    labs(x = "", y = "Coefficient of Variation", color = "", shape = "") +
    coord_flip() +
    theme(legend.title = element_text(size = rel(0.7))) +
    bri_theme

  p1 <- p

  bri_ggsave(paste0(rep_summary_folder, "/CV plots/CV plot - individual", suffix, ".png"), plot = p, width = 9, height = 6)

  write_tsv(df |> select(-Method, -flag), paste0(rep_summary_folder, "/CV plots/CV plot data.tsv"))

  dfd <- df |>
    mutate(
      flag = ifelse(
        type == "Original",
        sd_estimate_ok,
        "Replications"
      )
    )

  p <- ggplot(dfd) +
    aes(x = type, y = cv, color = TypeTechnique, group = type, shape = flag) +
    geom_boxplot(outliers = F) +
    geom_point(
      size = 2, alpha = 0.7,
      position = position_jitter(width = 0.1)
    ) +
    scale_color_manual(values = c(bri_color[["dark"]], bri_color[["epm"]], bri_color[["mtt"]], bri_color[["pcr"]])) +
    scale_shape_manual(values = c(17, 16, 16)) +
    facet_wrap(~Method) +
    labs(x = "", y = "Coefficient of Variation", color = "") +
    bri_theme

  p2 <- p

  bri_ggsave(paste0(rep_summary_folder, "/CV plots/CV plot - aggregate color", suffix, ".png"), plot = p, width = 10, height = 4)


  df <- rbind(
    df,
    df |> mutate(
      sd_estimate_ok = "Everything",
      flag = ifelse(
        type == "Original",
        paste0(type, "/\n", sd_estimate_ok),
        "Replications"
      )
    )
  )

  p <- ggplot(df) +
    aes(x = flag, y = cv) +
    geom_boxplot(outliers = F) +
    geom_point(size = 2.5, alpha = 0.8, position = position_jitter(width = 0.15)) +
    facet_wrap(~Method) +
    labs(x = "", y = "Coefficient of Variation", color = "") +
    bri_theme

  bri_ggsave(paste0(rep_summary_folder, "/CV plots/CV plot - aggregate", suffix, ".png"), plot = p, width = 10, height = 4)

  p <- ggplot(df) +
    aes(x = flag, y = cv) +
    geom_boxplot(outliers = F) +
    geom_point(size = 2.5, alpha = 0.8, position = position_jitter(width = 0.15)) +
    facet_wrap(~Method, nrow = 1) +
    labs(x = "", y = "Coefficient of Variation", color = "") +
    bri_theme

  bri_ggsave(paste0(rep_summary_folder, "/CV plots/CV plot - aggregate, same scale", suffix, ".png"), plot = p, width = 4, height = 10)

  # Combined figure
  p2 <- p2 + theme(legend.position = "none")
  p <- plot_grid(plotlist = list(p1, p2), ncol = 1, labels = "AUTO", rel_heights = c(0.6, 0.4))

  bri_ggsave(paste0(rep_summary_folder, "/CV plots/CV plot - combined", suffix, ".png"), plot = p, width = 8, height = 9)
}

# Original-study SD versus mean regression ---------------------------------

plot_original_sd_mean_regression <- function(
    data_path,
    output_dir) {
  if (!file.exists(data_path)) {
    stop("Original-study summary file not found: ", data_path)
  }

  required_columns <- c(
    "Prediction Markets ID",
    "EXP",
    "Selected Figure/Table",
    "Original Article DOI",
    "Data Processing",
    "Error Type",
    "SD Estimate",
    "Control Mean",
    "Control Error",
    "Control SD",
    "Reported Control Sample Size",
    "Assumed Control Sample Size",
    "Control Sample Size for DFs",
    "Treated Mean",
    "Treated Error",
    "Treated SD",
    "Reported Treated Sample Size",
    "Assumed Treated Sample Size",
    "Sample size (range for both groups)",
    "Treated Sample Size for DFs"
  )

  original_data <- readr::read_tsv(data_path, show_col_types = FALSE)
  missing_columns <- setdiff(required_columns, colnames(original_data))

  if (length(missing_columns) > 0) {
    stop(
      "Original-study summary is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  unexpected_error_types <- setdiff(
    unique(stats::na.omit(original_data[["Error Type"]])),
    c("SD", "SEM")
  )

  if (length(unexpected_error_types) > 0) {
    stop(
      "Unexpected original error-bar types: ",
      paste(unexpected_error_types, collapse = ", ")
    )
  }

  expected_sd_estimate_bases <- c(
    "reported SD",
    "reported SEM and N",
    "reported SEM and range for N"
  )
  unexpected_sd_estimate_bases <- setdiff(
    unique(stats::na.omit(original_data[["SD Estimate"]])),
    expected_sd_estimate_bases
  )

  if (length(unexpected_sd_estimate_bases) > 0) {
    stop(
      "Unexpected SD estimate bases: ",
      paste(unexpected_sd_estimate_bases, collapse = ", ")
    )
  }

  if (nrow(original_data) != 60 || n_distinct(original_data$EXP) != 60) {
    stop("Figure S10 expects exactly 60 unique source experiments")
  }

  expected_mtt_percentage_experiments <- c(
    "MTT4",
    "MTT13",
    "MTT26",
    "MTT33",
    "MTT37",
    "MTT42",
    "MTT45",
    "MTT60",
    "MTT71",
    "MTT83"
  )
  mtt_scale_audit <- original_data |>
    filter(stringr::str_detect(EXP, "^MTT")) |>
    mutate(
      percentage_scale = (
        .data[["Control Mean"]] > 10 &
          .data[["Treated Mean"]] > 10
      ),
      ratio_or_other_scale = (
        .data[["Control Mean"]] <= 10 &
          .data[["Treated Mean"]] <= 10
      )
    )

  if (any(!mtt_scale_audit$percentage_scale &
          !mtt_scale_audit$ratio_or_other_scale)) {
    stop(
      "MTT unit classification is ambiguous because an experiment spans ",
      "the percentage- and ratio-scale ranges"
    )
  }

  mtt_percentage_experiments <- mtt_scale_audit |>
    filter(percentage_scale) |>
    pull(EXP)

  if (!setequal(
    mtt_percentage_experiments,
    expected_mtt_percentage_experiments
  )) {
    stop(
      "The percentage-scale MTT experiments no longer match the audited set"
    )
  }

  expected_normalized_treated_experiments <- c(
    "MTT86",
    "MTT87",
    "MTT96",
    "PCR16",
    "PCR69",
    "PCR147",
    "PCR175",
    "PCR184"
  )
  normalized_zero_control_experiments <- original_data |>
    filter(
      is.finite(.data[["Control Mean"]]),
      dplyr::near(.data[["Control Mean"]], 1),
      is.finite(.data[["Control SD"]]),
      .data[["Control SD"]] == 0
    ) |>
    pull(EXP)

  if (!setequal(
    normalized_zero_control_experiments,
    expected_normalized_treated_experiments
  )) {
    stop(
      "The normalized zero-SD control experiments no longer match the ",
      "audited set"
    )
  }

  sanitize_tsv_text <- function(value) {
    value |>
      stringr::str_replace_all("[\\r\\n\\t]+", " ") |>
      stringr::str_squish()
  }

  make_arm_data <- function(arm) {
    mean_column <- paste(arm, "Mean")
    error_column <- paste(arm, "Error")
    sd_column <- paste(arm, "SD")
    reported_sample_size_column <- paste("Reported", arm, "Sample Size")
    assumed_sample_size_column <- paste("Assumed", arm, "Sample Size")
    sample_size_column <- paste(arm, "Sample Size for DFs")

    original_data |>
      transmute(
        study_id = .data[["Prediction Markets ID"]],
        EXP,
        selected_figure_or_table = sanitize_tsv_text(
          .data[["Selected Figure/Table"]]
        ),
        original_article_doi = .data[["Original Article DOI"]],
        data_processing = sanitize_tsv_text(.data[["Data Processing"]]),
        method = stringr::str_extract(EXP, "(EPM|MTT|PCR)"),
        arm = .env$arm,
        reported_error_type = .data[["Error Type"]],
        original_reported_error = .data[[error_column]],
        reported_sample_size = .data[[reported_sample_size_column]],
        sample_size_range = sanitize_tsv_text(
          .data[["Sample size (range for both groups)"]]
        ),
        assumed_sample_size = .data[[assumed_sample_size_column]],
        sample_size_used = .data[[sample_size_column]],
        sd_estimate_basis = .data[["SD Estimate"]],
        original_mean = .data[[mean_column]],
        original_sd = .data[[sd_column]]
      )
  }

  arm_data <- bind_rows(
    make_arm_data("Control"),
    make_arm_data("Treated")
  ) |>
    mutate(
      method = factor(method, levels = c("EPM", "MTT", "PCR")),
      arm = factor(arm, levels = c("Control", "Treated")),
      study_arm = case_when(
        as.character(arm) == "Control" ~ "Control",
        EXP %in% normalized_zero_control_experiments ~ "Treated (Normalized)",
        TRUE ~ "Treated"
      ),
      study_arm = factor(
        study_arm,
        levels = c("Control", "Treated", "Treated (Normalized)")
      ),
      sd_provenance = case_when(
        reported_error_type == "SD" &
          sd_estimate_basis == "reported SD" ~ "SD",
        reported_error_type == "SEM" &
          sd_estimate_basis == "reported SEM and N" ~
            "SEM (exact)",
        reported_error_type == "SEM" &
          sd_estimate_basis == "reported SEM and range for N" ~
            "SEM (estimated)",
        TRUE ~ NA_character_
      ),
      sd_provenance = factor(
        sd_provenance,
        levels = c(
          "SD",
          "SEM (exact)",
          "SEM (estimated)"
        )
      ),
      mtt_original_scale = case_when(
        as.character(method) != "MTT" ~ "Not applicable",
        EXP %in% mtt_percentage_experiments ~
          "Percentage (converted to ratio)",
        TRUE ~ "Ratio/other (unchanged)"
      ),
      scale_factor = if_else(
        EXP %in% mtt_percentage_experiments,
        0.01,
        1
      ),
      analysis_reported_error = original_reported_error * scale_factor,
      analysis_mean = original_mean * scale_factor,
      analysis_sd = original_sd * scale_factor,
      exclusion_reason = case_when(
        is.na(analysis_mean) | !is.finite(analysis_mean) ~
          "Missing or non-finite mean",
        analysis_mean <= 0 ~ "Nonpositive mean cannot be log-transformed",
        is.na(analysis_sd) | !is.finite(analysis_sd) ~
          "Missing or non-finite SD",
        analysis_sd <= 0 ~ "Nonpositive SD cannot be log-transformed",
        TRUE ~ NA_character_
      )
    )

  if (any(is.na(arm_data$method))) {
    stop("Could not derive EPM, MTT, or PCR method for every experiment")
  }

  if (any(is.na(arm_data$sd_provenance))) {
    stop("Could not derive SD provenance for every source arm")
  }

  if (
    nrow(arm_data) != 120 ||
      n_distinct(paste(arm_data$EXP, arm_data$arm)) != 120
  ) {
    stop("Figure S10 expects exactly 120 unique source arms")
  }

  included_data <- arm_data |>
    filter(is.na(exclusion_reason)) |>
    mutate(
      log_mean = log(analysis_mean),
      log_sd = log(analysis_sd)
    )

  excluded_data <- arm_data |>
    filter(!is.na(exclusion_reason)) |>
    select(
      study_id,
      EXP,
      selected_figure_or_table,
      original_article_doi,
      data_processing,
      method,
      arm,
      study_arm,
      reported_error_type,
      sd_provenance,
      original_reported_error,
      analysis_reported_error,
      reported_sample_size,
      sample_size_range,
      assumed_sample_size,
      sample_size_used,
      sd_estimate_basis,
      mtt_original_scale,
      scale_factor,
      original_mean,
      original_sd,
      analysis_mean,
      analysis_sd,
      exclusion_reason
    )

  normalized_zero_control <- (
    as.character(excluded_data$arm) == "Control" &
      is.finite(excluded_data$original_mean) &
      dplyr::near(excluded_data$original_mean, 1) &
      is.finite(excluded_data$original_sd) &
      excluded_data$original_sd == 0
  )

  if (
    nrow(included_data) != 112 ||
      nrow(excluded_data) != 8 ||
      !all(normalized_zero_control)
  ) {
    stop(
      "Figure S10 expects 112 included arms and eight normalized ",
      "zero-SD control exclusions"
    )
  }

  expected_study_arm_counts <- c(
    "Control" = 52,
    "Treated" = 52,
    "Treated (Normalized)" = 8
  )
  study_arm_counts <- table(included_data$study_arm)

  if (!identical(
    unname(as.integer(study_arm_counts)),
    unname(as.integer(expected_study_arm_counts))
  )) {
    stop("The Figure S10 study-arm counts no longer match the audited data")
  }

  excluded_experiment_ids <- excluded_data |>
    mutate(
      experiment_number = stringr::str_extract(EXP, "[0-9]+") |>
        as.integer()
    ) |>
    arrange(method, experiment_number) |>
    distinct(EXP) |>
    pull(EXP)

  format_experiment_list <- function(experiment_ids) {
    number_of_experiments <- length(experiment_ids)

    if (number_of_experiments == 1) {
      return(experiment_ids)
    }

    if (number_of_experiments == 2) {
      return(paste(experiment_ids, collapse = " and "))
    }

    paste0(
      paste(experiment_ids[-number_of_experiments], collapse = ", "),
      ", and ",
      experiment_ids[[number_of_experiments]]
    )
  }

  figure_caption <- paste0(
    "MTT percentage-scale means and SDs were divided by 100. ",
    "\"Treated (Normalized)\" denotes treated arms whose normalized ",
    "control had mean = 1 and SD = 0.\n",
    "The corresponding controls from ",
    format_experiment_list(excluded_experiment_ids),
    " were omitted because ln(SD) is undefined."
  )

  panel_levels <- c("EPM", "MTT", "PCR")
  panel_data <- included_data |>
    mutate(
      panel = factor(as.character(method), levels = panel_levels)
    )

  fit_panel <- function(panel_name) {
    data <- panel_data |>
      filter(as.character(panel) == panel_name)

    if (nrow(data) < 3) {
      stop("Not enough observations to fit panel: ", panel_name)
    }

    model <- stats::lm(log_sd ~ log_mean, data = data)
    model_summary_object <- summary(model)
    point_prediction <- suppressWarnings(
      stats::predict(
        model,
        newdata = data,
        interval = "prediction",
        level = 0.95
      )
    )
    residual_log_sd <- stats::residuals(model)
    fitted_log_sd <- stats::fitted(model)
    externally_studentized_residual <- stats::rstudent(model)

    diagnostics <- data |>
      mutate(
        fitted_log_sd = unname(fitted_log_sd),
        fitted_sd = exp(fitted_log_sd),
        prediction_lower_log_sd = unname(point_prediction[, "lwr"]),
        prediction_upper_log_sd = unname(point_prediction[, "upr"]),
        prediction_lower_sd = exp(prediction_lower_log_sd),
        prediction_upper_sd = exp(prediction_upper_log_sd),
        residual_log_sd = unname(residual_log_sd),
        externally_studentized_residual = unname(
          externally_studentized_residual
        ),
        abs_externally_studentized_residual = abs(
          externally_studentized_residual
        ),
        residual_rank_within_method = dplyr::min_rank(
          dplyr::desc(abs_externally_studentized_residual)
        ),
        meets_abs_rstudent_2 = abs_externally_studentized_residual >= 2,
        outside_95_prediction_interval = (
          log_sd < prediction_lower_log_sd |
            log_sd > prediction_upper_log_sd
        ),
        diagnostic_direction = case_when(
          residual_log_sd < 0 ~ "Below fitted line",
          residual_log_sd > 0 ~ "Above fitted line",
          TRUE ~ "On fitted line"
        ),
        sd_to_fitted_ratio = analysis_sd / fitted_sd
      )

    prediction_grid <- tibble(
      log_mean = seq(
        min(data$log_mean),
        max(data$log_mean),
        length.out = 200
      )
    )
    prediction_interval <- stats::predict(
      model,
      newdata = prediction_grid,
      interval = "prediction",
      level = 0.95
    )
    prediction_data <- prediction_grid |>
      mutate(
        panel = factor(panel_name, levels = panel_levels),
        fitted_log_sd = unname(prediction_interval[, "fit"]),
        prediction_lower_log_sd = unname(prediction_interval[, "lwr"]),
        prediction_upper_log_sd = unname(prediction_interval[, "upr"])
      )

    model_summary <- tibble(
      panel = factor(panel_name, levels = panel_levels),
      intercept = unname(stats::coef(model)[[1]]),
      slope = unname(stats::coef(model)[[2]]),
      r_squared = unname(model_summary_object$r.squared),
      n_points = nrow(data),
      n_studies = n_distinct(data$EXP)
    )

    list(
      model = model,
      diagnostics = diagnostics,
      prediction_data = prediction_data,
      model_summary = model_summary
    )
  }

  fit_results <- lapply(panel_levels, fit_panel)
  names(fit_results) <- panel_levels

  diagnostics <- purrr::map_dfr(fit_results, function(result) result$diagnostics) |>
    arrange(panel, residual_rank_within_method, EXP, arm) |>
    select(
      panel,
      residual_rank_within_method,
      study_id,
      EXP,
      selected_figure_or_table,
      original_article_doi,
      data_processing,
      method,
      arm,
      study_arm,
      reported_error_type,
      sd_provenance,
      original_reported_error,
      analysis_reported_error,
      reported_sample_size,
      sample_size_range,
      assumed_sample_size,
      sample_size_used,
      sd_estimate_basis,
      mtt_original_scale,
      scale_factor,
      original_mean,
      original_sd,
      analysis_mean,
      analysis_sd,
      log_mean,
      log_sd,
      fitted_log_sd,
      fitted_sd,
      prediction_lower_log_sd,
      prediction_upper_log_sd,
      prediction_lower_sd,
      prediction_upper_sd,
      residual_log_sd,
      externally_studentized_residual,
      abs_externally_studentized_residual,
      meets_abs_rstudent_2,
      outside_95_prediction_interval,
      diagnostic_direction,
      sd_to_fitted_ratio
    )
  prediction_data <- purrr::map_dfr(
    fit_results,
    function(result) result$prediction_data
  )
  model_summary <- purrr::map_dfr(
    fit_results,
    function(result) result$model_summary
  )

  p <- ggplot(diagnostics, aes(x = log_mean, y = log_sd)) +
    geom_ribbon(
      data = prediction_data,
      aes(
        x = log_mean,
        ymin = prediction_lower_log_sd,
        ymax = prediction_upper_log_sd
      ),
      inherit.aes = FALSE,
      fill = bri_color[["light"]],
      alpha = 0.45
    ) +
    geom_line(
      data = prediction_data,
      aes(x = log_mean, y = fitted_log_sd),
      inherit.aes = FALSE,
      color = bri_color[["dark"]],
      linewidth = 0.8
    ) +
    geom_point(
      aes(color = study_arm, shape = sd_provenance),
      size = 2.5,
      alpha = 0.8
    ) +
    facet_wrap(~panel, nrow = 1, scales = "free") +
    scale_color_manual(
      values = c(
        "Control" = bri_color[["none"]],
        "Treated" = bri_color[["main"]],
        "Treated (Normalized)" = bri_color[["q3"]]
      ),
      drop = FALSE,
      name = "Study arm"
    ) +
    scale_shape_manual(
      values = c(
        "SD" = 16,
        "SEM (exact)" = 17,
        "SEM (estimated)" = 2
      ),
      drop = FALSE,
      name = "SD derivation"
    ) +
    guides(
      color = guide_legend(order = 1, nrow = 1),
      shape = guide_legend(order = 2, nrow = 1)
    ) +
    labs(
      x = expression(ln(mean)),
      y = expression(ln(SD))
    ) +
    bri_theme +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical"
    )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  writeLines(figure_caption, file.path(output_dir, "figure-caption.txt"))

  figure_path <- file.path(
    output_dir,
    "original-study-log-sd-vs-log-mean.png"
  )
  bri_ggsave(
    figure_path,
    plot = p,
    width = 12,
    height = 6,
    dpi = 300
  )
  readr::write_tsv(
    model_summary,
    file.path(output_dir, "model-summary.tsv")
  )
  readr::write_tsv(
    diagnostics,
    file.path(output_dir, "point-diagnostics.tsv")
  )
  readr::write_tsv(
    excluded_data,
    file.path(output_dir, "excluded-points.tsv")
  )

  invisible(list(
    plot = p,
    models = lapply(fit_results, function(result) result$model),
    model_summary = model_summary,
    diagnostics = diagnostics,
    excluded_points = excluded_data,
    figure_path = figure_path,
    figure_caption = figure_caption
  ))
}

# Plot histograms of post hoc power
plot_power_histograms <- function(dfi) {
  dir.create(paste0("./output"))
  dir.create(paste0("./output/", results_path, "/Power Histograms"))


  plot_single_histogram <- function(df, varname, xlab, ylab) {
    df$value <- df[[varname]]
    df$Method <- str_extract(df$EXP, "(MTT|PCR|EPM)")
    df <- df |> filter(!is.na(value))
    ggplot(df) +
      aes(x = value) +
      geom_histogram(binwidth = 0.1, center = 0.05, aes(fill = Method), color = bri_color[["none"]]) +
      scale_x_continuous(limits = c(0, 1), n.breaks = 10) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1), add = 0)) +
      facet_wrap(~Method) +
      labs(y = ylab, x = xlab) +
      scale_fill_manual(values = c("EPM" = bri_color[["epm"]],
                        "MTT" = bri_color[["mtt"]],
                        "PCR" = bri_color[["pcr"]]),
                        guide="none") +
      bri_theme +
      theme(
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position="none"
      )
  }

  dfs <- dfi |>
    group_by(EXP) |>
    summarise(
      power_z = mean(`post hoc power Z`, na.rm = T),
      power_t = mean(`post hoc power T`, na.rm = T),
      power_knha = mean(`post hoc power knha`, na.rm = T),
      power_bigexp = mean(`post hoc power bigexp`, na.rm = T)
    )


  p_power_z <- plot_single_histogram(dfs, "power_z", "Post-hoc power (z distribution)", "# of experiments")
  bri_ggsave(paste0("./output/", results_path, "/Power Histograms/post hoc power, by experiment, z.png"), plot = p_power_z, width = 12, height = 4)

  p_power_t <- plot_single_histogram(dfs, "power_t", "Post-hoc power (t distribution - # of units)", "# of experiments")
  bri_ggsave(paste0("./output/", results_path, "/Power Histograms/post hoc power, by experiment, t.png"), plot = p_power_t, width = 12, height = 4)

  p_power_knha <- plot_single_histogram(dfs, "power_knha", "Post-hoc power (t distribution - # of replications)", "# of experiments")
  bri_ggsave(paste0("./output/", results_path, "/Power Histograms/post hoc power, by experiment, knha.png"), plot = p_power_knha, width = 12, height = 4)

  p_power_bigexp <- plot_single_histogram(dfs, "power_bigexp", "Post-hoc power (single mean)", "# of experiments")
  bri_ggsave(paste0("./output/", results_path, "/Power Histograms/post hoc power, by experiment, bigexp.png"), plot = p_power_bigexp, width = 12, height = 4)

  p_power_combined <- ggarrange(
    plotlist = list(
      p_power_t + xlab("Post-hoc power (t distribution - # of units)"),
      p_power_z + xlab("Post-hoc power (z distribution)"), 
      p_power_knha + xlab("Post-hoc power (t distribution - # of replications)"),
      p_power_bigexp + xlab("Post-hoc power (single mean)")), 
    ncol = 1, 
    labels = "AUTO", common.legend = T)
  bri_ggsave(paste0("./output/", results_path, "/Power Histograms/post hoc power, by experiment, combined.png"), plot = p_power_combined, width = 12, height = 12)
}

# Plots specification curves
plot_specification_curve <- function(results_path, include_method, suffix = "") {
  # Data prep
  if (include_method) {
    AGGDATA <- read_tsv(paste0("output/", results_path, "/Replication Rate Aggregate.tsv"))
  } else {
    AGGDATA <- read_tsv(paste0("output/", results_path, "/Replication Rate Aggregate - without Method.tsv"))
  }

  AGGDATA <- AGGDATA |>
    select(-successful) |> 
    filter(
      IsPercentage == T,
      !(Metric %in% c("Significant", "t_score_individual", "OriginalES", "SignalErrorAll", "ExaggerationAll", "mean_abs_diff_reps", "mean_abs_diff_reps_orig", "OriginalCV", "REP_Most_Criteria_WithTies")),
      !(Inclusion_Set %in% c("all_exps", "only_80_single_rep_T", "only_80_single_rep_Z", "only_80_single_rep_KNHA"))
    ) |>
    mutate(
      Metric = MetricShortName,
      index = 1:n()
    ) |>
    rename(
      repro_rate = Value
    ) |>
    mutate(
      MA_Dist = factor(
        case_match(
          MA_Dist,
          "t" ~ "t distribution (# of units)",
          "z" ~ "z distribution",
          "knha" ~ "t distribution (# of replications)",
          "bigexp" ~ "Single mean"
        ),
        levels = c("t distribution (# of units)", "z distribution", "t distribution (# of replications)", "Single mean")
      )
    )

  desired_inclusion_levels <- c(
    "Primary",
    "Lab's choice",
    "All experiments",
    "All experiments (BRI UNIT)",
    "≥2 replications",
    "3 replications",
    "≥80% power (t - # of units)",
    "≥80% power (z)",
    "≥80% power (t - # of replications)",
    "≥80% power (single mean)"
  )

  AGGDATA <- AGGDATA |>
    mutate(
      Inclusion_Set = case_match(
        Inclusion_Set,
        "all_exps" ~ "All experiments (BRI UNIT)",
        "all_exps_lab_units" ~ "All experiments",
        "included_by_lab" ~ "Lab's choice",
        "primary" ~ "Primary",
        "only_3_reps" ~ "3 replications",
        "at_least_2_reps" ~ "≥2 replications",
        "only_80_power_a_posteriori_T" ~ "≥80% power (t - # of units)",
        "only_80_power_a_posteriori_Z" ~ "≥80% power (z)",
        "only_80_power_a_posteriori_KNHA" ~ "≥80% power (t - # of replications)",
        "only_80_power_a_posteriori_bigexp" ~ "≥80% power (single mean)"
      ),
      Inclusion_Set = factor(
        Inclusion_Set,
        levels = desired_inclusion_levels
      )
    )

  # IMPORTANT: We avoid mislabeling "ALL" as a specific choice.
  # For the specification curve with method toggles, we'll filter to the
  # four explicit combinations (ALL_*), and derive two binary toggles
  # instead of using the combined "Method" levels. This yields 50/50 splits.
  AGGDATA <- AGGDATA |>
    mutate(Method = as.character(Method))

  AGGDATA <- AGGDATA |>
    mutate(
      Level = case_match(
        Level,
        "Set of replications" ~ "Aggregate of replications",
        .default = Level
      )
    )

  # Prepare data for the top plot consistent with the bottom ticks
  AGGDATA_USED <- AGGDATA

  # Tick graph
  # Making the specification curve graph
  # If include_method = TRUE: only individual methods (single method band, black).
  # If include_method = FALSE: only ALL_* + toggles (PCR_Scale, MTT_Pairing) to show 50/50 splits.
  if (include_method) {
    # Only individual methods; no toggles
    AGGDATA_TOG <- AGGDATA |>
      filter(Method %in% c("EPM","PCR","ALTPCR","MTT","ALTMTT")) |>
      mutate(
        Method_Individual = Method
      )

    spec_parameters <- c("Inclusion_Set", "Level", "MA_Dist", "Metric", "Method_Individual")

    method_levels <- na.omit(unique(AGGDATA_TOG$Method_Individual))
    # Ensure stable order of method levels
    method_levels <- intersect(c("EPM","PCR","ALTPCR","MTT","ALTMTT"), method_levels)

    available_inclusion_levels <- unique(as.character(AGGDATA_TOG$Inclusion_Set))
    inclusion_levels <- desired_inclusion_levels[desired_inclusion_levels %in% available_inclusion_levels]
    
    ma_dist_levels <- c("t distribution (# of units)", "z distribution", "t distribution (# of replications)", "Single mean")[c("t distribution (# of units)", "z distribution", "t distribution (# of replications)", "Single mean") %in% AGGDATA_TOG$MA_Dist]

    # Fixed order for specification curve y-axis (user-specified)
    metric_order <- c(
      "Original in replication 95% PI",
      "Replication in original 95% CI",
      "Same-sign significance (p < 0.05)", "Same-sign significance (p<0.05)",
      "> 50% subjectively replicated", ">50% subjectively replicated", ">= 50% subjectively replicated", "≥50% subjectively replicated",
      "> 50% significant", ">50% significant", ">= 50% significant", "≥50% significant"
    )
    metric_order <- metric_order[metric_order %in% unique(AGGDATA_TOG$Metric)]
    
    indiv_metric_order <- c(
      "Subjectively replicated (Indiv)",
      "Replication in original CI (Indiv)",
      "Same-sign significance (Indiv)"
    )
    indiv_metric_order <- indiv_metric_order[indiv_metric_order %in% unique(AGGDATA_TOG$Metric)]
    
    # Filter data to remove any metrics not in our valid list (fixes NA issue)
    AGGDATA_TOG <- AGGDATA_TOG |> filter(Metric %in% c(metric_order, indiv_metric_order))
    
    inclusion_order <- c(
      "All experiments",
      "≥2 replications",
      "Lab's choice",
      "3 replications",
      "≥80% power (t - # of replications)",
      "≥80% power (t - # of units)",
      "≥80% power (z)",
      "≥80% power (single mean)",
      "Primary"
    )
    inclusion_order <- inclusion_order[inclusion_order %in% inclusion_levels]
    
    ma_dist_order <- c(
      "t distribution (# of replications)",
      "t distribution (# of units)",
      "z distribution",
      "Single mean"
    )
    ma_dist_order <- ma_dist_order[ma_dist_order %in% ma_dist_levels]
    
    spec_parameters_levels <- rev(c(
      metric_order,
      "Aggregate of replications",
      indiv_metric_order,
      "Individual replication",
      inclusion_order,
      method_levels,
      ma_dist_order
    ))

    # Methods in black, as requested
    method_colors <- rep("black", length(method_levels))
    spec_colors <- rev(c(
      rep(RColorBrewer::brewer.pal(5, "Dark2")[1], length(1:7)),
      rep(RColorBrewer::brewer.pal(7, "Greens")[7], 1),
      rep(RColorBrewer::brewer.pal(5, "Dark2")[1], length(8:10)),
      rep(RColorBrewer::brewer.pal(7, "Greens")[7], 1),
      rep(RColorBrewer::brewer.pal(7, "Blues")[7], length(inclusion_levels)),
      method_colors,
      rep(RColorBrewer::brewer.pal(7, "Oranges")[7], length(ma_dist_levels))
    ))

    # Deduplicate breaks while keeping order and map colors to labels
    spec_breaks <- unique(spec_parameters_levels)
    spec_colors_named <- setNames(spec_colors, spec_parameters_levels)

    # Re-index so both panels align
    AGGDATA_TOG <- AGGDATA_TOG |>
      mutate(index = row_number())

    SPECIFICATION <- AGGDATA_TOG |>
      select(all_of(c("index", "repro_rate", spec_parameters))) |>
      pivot_longer(cols = -c(index, repro_rate)) |>
      filter(!is.na(value)) |>
      arrange(repro_rate) |>
      mutate(value = factor(value, levels = spec_parameters_levels))
    AGGDATA_USED <- AGGDATA_TOG
  } else {
    # No individual methods; keep only ALL_* and add toggles to show 50/50 splits
    AGGDATA_TOG <- AGGDATA |>
      filter(str_detect(Method, "^ALL_")) |>
      mutate(
        PCR_Scale = if_else(str_detect(Method, "ALTPCR"), "PCR (linear)", "PCR (log)"),
        MTT_Pairing = if_else(str_detect(Method, "ALTMTT"), "MTT (original)", "MTT (paired)")
      )

    spec_parameters <- c("Inclusion_Set", "Level", "MA_Dist", "Metric", "PCR_Scale", "MTT_Pairing")

    available_inclusion_levels <- unique(as.character(AGGDATA_TOG$Inclusion_Set))
    inclusion_levels <- desired_inclusion_levels[desired_inclusion_levels %in% available_inclusion_levels]

    ma_dist_levels <- c("t distribution (# of units)", "z distribution", "t distribution (# of replications)", "Single mean")[c("t distribution (# of units)", "z distribution", "t distribution (# of replications)", "Single mean") %in% AGGDATA_TOG$MA_Dist]
    
    # Fixed order for specification curve y-axis (user-specified)
    metric_order <- c(
      "Original in replication 95% PI",
      "Replication in original 95% CI",
      "Same-sign significance (p < 0.05)", "Same-sign significance (p<0.05)",
      "> 50% subjectively replicated", ">50% subjectively replicated", ">= 50% subjectively replicated", "≥50% subjectively replicated",
      "> 50% significant", ">50% significant", ">= 50% significant", "≥50% significant"
    )
    metric_order <- metric_order[metric_order %in% unique(AGGDATA_TOG$Metric)]
    
    indiv_metric_order <- c(
      "Subjectively replicated (Indiv)",
      "Replication in original CI (Indiv)",
      "Same-sign significance (Indiv)"
    )
    indiv_metric_order <- indiv_metric_order[indiv_metric_order %in% unique(AGGDATA_TOG$Metric)]
    
    # Filter data to remove any metrics not in our valid list (fixes NA issue)
    AGGDATA_TOG <- AGGDATA_TOG |> filter(Metric %in% c(metric_order, indiv_metric_order))
    
    inclusion_order <- c(
      "All experiments",
      "≥2 replications",
      "Lab's choice",
      "3 replications",
      "≥80% power (t - # of replications)",
      "≥80% power (t - # of units)",
      "≥80% power (z)",
      "≥80% power (single mean)",
      "Primary"
    )
    inclusion_order <- inclusion_order[inclusion_order %in% inclusion_levels]
    
    pcr_order <- c("PCR (log)", "PCR (linear)")
    mtt_order <- c("MTT (paired)", "MTT (original)")
    
    ma_dist_order <- c(
      "t distribution (# of replications)",
      "t distribution (# of units)",
      "z distribution",
      "Single mean"
    )
    ma_dist_order <- ma_dist_order[ma_dist_order %in% ma_dist_levels]
    
    spec_parameters_levels <- rev(c(
      metric_order,
      "Aggregate of replications",
      indiv_metric_order,
      "Individual replication",
      inclusion_order,
      pcr_order,
      mtt_order,
      ma_dist_order
    ))

    spec_colors <- rev(c(
      rep(RColorBrewer::brewer.pal(5, "Dark2")[1], length(1:7)),
      rep(RColorBrewer::brewer.pal(7, "Greens")[7], 1),
      rep(RColorBrewer::brewer.pal(5, "Dark2")[1], length(8:10)),
      rep(RColorBrewer::brewer.pal(7, "Greens")[7], 1),
      rep(RColorBrewer::brewer.pal(7, "Blues")[7], length(inclusion_levels)),
      rep(RColorBrewer::brewer.pal(7, "Greys")[7], 2),  # PCR_Scale (2 levels)
      rep(RColorBrewer::brewer.pal(7, "Greys")[5], 2),  # MTT_Pairing (2 levels)
      rep(RColorBrewer::brewer.pal(7, "Oranges")[7], length(ma_dist_levels))
    ))

    spec_breaks <- unique(spec_parameters_levels)
    spec_colors_named <- setNames(spec_colors, spec_parameters_levels)

    # Re-index after filtering so both panels align
    AGGDATA_TOG <- AGGDATA_TOG |>
      mutate(index = row_number())

    SPECIFICATION <- AGGDATA_TOG |>
      select(all_of(c("index", "repro_rate", spec_parameters))) |>
      pivot_longer(cols = -c(index, repro_rate)) |>
      filter(!is.na(value)) |>
      arrange(repro_rate) |>
      mutate(value = factor(value, levels = spec_parameters_levels))
    AGGDATA_USED <- AGGDATA_TOG
  }

  # Top/line graph built from AGGDATA_USED to match SPECIFICATION
  # Find the quartiles and median
  quartiles <- AGGDATA_USED |>
    pull(repro_rate) |>
    quantile(na.rm = T, probs = c(0.25, 0.5, 0.75)) |>
    round(2)

  # Making the plot with reproducibility rates
  plot_repro <- ggplot(AGGDATA_USED) +
    aes(x = reorder(index, repro_rate), y = repro_rate) +
    geom_point(size = 1) +
    geom_hline(yintercept = quartiles[1], linetype = "dashed", color = bri_color[["none"]], alpha = 0.6) +
    geom_hline(yintercept = quartiles[2], linetype = "solid", color = bri_color[["low"]]) +
    geom_hline(yintercept = quartiles[3], linetype = "dashed", color = bri_color[["none"]], alpha = 0.6) +
    annotate(
      geom = "text", x = 2, y = quartiles[3] + 0.1,
      label = paste0("Median: ", quartiles[2] * 100, "%\nIQR: [", quartiles[1] * 100, "%, ", quartiles[3] * 100, "%]"), hjust = 0
    ) +
    labs(x = "", y = "Reproducibility Rate") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::label_percent()) + 
    theme_minimal() +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      panel.grid = element_blank(), plot.title = element_blank()
    )

  plot_spec2 <- ggplot(SPECIFICATION) +
    aes(x = reorder(index, repro_rate), y = value, fill = value) +
    geom_tile(height = 0.4) +
    scale_alpha(range = c(0, 1)) +
    scale_fill_manual(breaks = spec_breaks, values = spec_colors_named[spec_breaks]) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.5, color = "#999999"),
      legend.position = "none", panel.grid = element_blank(),
      axis.text.y = element_text(size = 11), plot.title = element_blank()
    )

  # Full graph
  # Plotting the aligned graphs
  p <- plot_grid(plot_repro, plot_spec2,
    align = "v", axis = "l",
    ncol = 1, nrow = 2, rel_heights = c(1, 1.5, 1.5)
  )

  w <- 18
  if (!include_method) w <- 9
  ggsave(paste0("output/", results_path, "/specification-curve", suffix, ".png"),
    plot = p,
    width = w, height = 12, dpi = 250, bg = "white"
  )
}

# Kappa -------------------------------------------------------------------

format_log10_p_value_labels <- function(log10_breaks) {
  labels <- scales::label_math(10^.x)(-log10_breaks)
  labels[log10_breaks == 0] <- expression(1)
  labels
}

plot_kappa_exp <- function(input_path, output_path) {
  # Import the data
  predictor_correlations_exp <- read_tsv(input_path) |>
    select(
      `Original in replication PI`,
      `Replication in original CI`,
      `Same sense significance`,
      `Voting (with ties)`,
      `Subjective (with ties)`
    )

  # Convert variables to integers (0/1)
  data_exp <- predictor_correlations_exp |>
    mutate(across(everything(), ~ as.integer(. == TRUE)))

  # Calculate kappa for variable pairs
  kappa_results_exp <- combn(names(data_exp), 2, simplify = FALSE) |>
    map_dfr(~ {
      # Remove missing values
      paired_data <- data_exp |>
        select(all_of(.x)) |>
        drop_na()
      # Skip if insufficient data for kappa calculation
      if (nrow(paired_data) < 2) {
        return(tibble(
          var1 = .x[1],
          var2 = .x[2],
          kappa = NA_real_,
          p_value = NA_real_
        ))
      }
      # Calculate kappa
      kappa_result <- kappa2(paired_data)
      tibble(
        var1 = .x[1],
        var2 = .x[2],
        kappa = kappa_result$value,
        p_value = kappa_result$p.value
      )
    })

  # Reformat data for visualization
  kappa_matrix_exp <- kappa_results_exp |>
    pivot_longer(c(kappa, p_value), names_to = "metric", values_to = "value") |>
    pivot_wider(
      names_from = metric,
      values_from = value
    ) |>
    mutate(
      var1 = if_else(row_number() == n(), "Subjective (with ties)", var1),
      var2 = if_else(row_number() == n(), "Voting (with ties)", var2),
      var1 = fct_recode(var1,
        "Original in replication's 95% PI" = "Original in replication PI",
        "Replication in original 95% CI" = "Replication in original CI",
        "Same-sign significance (p<0.05)" = "Same sense significance",
        "≥50% replications significant" = "Voting (with ties)",
        "≥50% subjectively replicated" = "Subjective (with ties)",
        "Subjectively replicated" = "Subjective"
      ),
      var2 = fct_recode(var2,
        "Original in replication's 95% PI" = "Original in replication PI",
        "Replication in original 95% CI" = "Replication in original CI",
        "Same-sign significance (p<0.05)" = "Same sense significance",
        "≥50% replications significant" = "Voting (with ties)",
        "≥50% subjectively replicated" = "Subjective (with ties)",
        "Subjectively replicated" = "Subjective"
      )
    )

  kappa_vals_exp <- -log10(kappa_matrix_exp$p_value)
  kappa_max_exp <- suppressWarnings(max(kappa_vals_exp, na.rm = TRUE))
  if (!is.finite(kappa_max_exp) || kappa_max_exp <= 0) {
    kappa_breaks_exp <- c(0, 0.1)
    kappa_limits_exp <- c(0, 0.1)
  } else {
    kappa_breaks_exp <- scales::pretty_breaks(n = 3)(c(0, kappa_max_exp))
    kappa_breaks_exp <- unique(kappa_breaks_exp[kappa_breaks_exp >= 0 & kappa_breaks_exp <= kappa_max_exp])
    if (length(kappa_breaks_exp) == 0) {
      kappa_breaks_exp <- c(0, kappa_max_exp)
    }
    kappa_limits_exp <- c(0, kappa_max_exp)
  }

  # Create visualization with geom_tile
  kappa_plot_exp <- kappa_matrix_exp |>
    ggplot(aes(x = var1, y = fct_rev(var2), fill = -log10(p_value), label = round(kappa, 2))) +
    geom_tile(color = "white") +
    geom_text() +
    scale_fill_gradient(
      low = bri_color[["mid"]],
      high = bri_color[["high"]],
      na.value = bri_color[["none"]],
      breaks = kappa_breaks_exp,
      limits = kappa_limits_exp,
      labels = format_log10_p_value_labels(kappa_breaks_exp)
    ) +
    labs(x = "", y = "", title = "Agreement between criteria (experiment)", fill = "p-value") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    bri_theme +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      title = element_text(size = 9),
      panel.border = element_blank()
    )

  # Save
  bri_ggsave(paste0(output_path, "cortable - experiment - kappa.png"), plot = kappa_plot_exp, width = 8, height = 8)


  # Return the visualization
  return(kappa_plot_exp)
}

plot_kappa_rep <- function(input_path, output_path) {
  # Import and preprocess the data
  predictor_correlations_rep <- read_tsv(input_path) |>
    filter(EXP != "(Protocol Team)") |>
    filter(EXP != "(Data Collection Team)") |>
    select(
      `Replication in original CI`,
      `Same sense significance`,
      `Subjective`
    )

  # Convert variables to integers (0/1)
  data_rep <- predictor_correlations_rep |>
    mutate(across(everything(), ~ as.integer(. == TRUE)))

  # Calculate kappa for variable pairs
  kappa_results_rep <- combn(names(data_rep), 2, simplify = FALSE) |>
    map_dfr(~ {
      # Remove missing values
      paired_data <- data_rep |>
        select(all_of(.x)) |>
        drop_na()
      # Skip if insufficient data for kappa calculation
      if (nrow(paired_data) < 2) {
        return(tibble(
          var1 = .x[1],
          var2 = .x[2],
          kappa = NA_real_,
          p_value = NA_real_
        ))
      }
      # Calculate kappa
      kappa_result <- kappa2(paired_data)
      tibble(
        var1 = .x[1],
        var2 = .x[2],
        kappa = kappa_result$value,
        p_value = kappa_result$p.value
      )
    })

  # Reformat data for visualization
  kappa_matrix_rep <- kappa_results_rep |>
    pivot_longer(c(kappa, p_value), names_to = "metric", values_to = "value") |>
    pivot_wider(
      names_from = metric,
      values_from = value
    ) |>
    mutate(
      var1 = fct_recode(var1,
        "Replication in original 95% CI" = "Replication in original CI",
        "Same-sign significance (p<0.05)" = "Same sense significance",
        "Subjectively replicated" = "Subjective"
      ),
      var2 = fct_recode(var2,
        "Replication in original 95% CI" = "Replication in original CI",
        "Same-sign significance (p<0.05)" = "Same sense significance",
        "Subjectively replicated" = "Subjective"
      )
    )

  kappa_vals_rep <- -log10(kappa_matrix_rep$p_value)
  kappa_max_rep <- suppressWarnings(max(kappa_vals_rep, na.rm = TRUE))
  if (!is.finite(kappa_max_rep) || kappa_max_rep <= 0) {
    kappa_breaks_rep <- c(0, 0.1)
    kappa_limits_rep <- c(0, 0.1)
  } else {
    kappa_breaks_rep <- scales::pretty_breaks(n = 3)(c(0, kappa_max_rep))
    kappa_breaks_rep <- unique(kappa_breaks_rep[kappa_breaks_rep >= 0 & kappa_breaks_rep <= kappa_max_rep])
    if (length(kappa_breaks_rep) == 0) {
      kappa_breaks_rep <- c(0, kappa_max_rep)
    }
    kappa_limits_rep <- c(0, kappa_max_rep)
  }

  # Create visualization with geom_tile
  kappa_plot_rep <- kappa_matrix_rep |>
    ggplot(aes(x = var1, y = fct_rev(var2), fill = -log10(p_value), label = round(kappa, 2))) +
    geom_tile(color = "white") +
    geom_text() +
    scale_fill_gradient(
      low = bri_color[["mid"]],
      high = bri_color[["high"]],
      na.value = bri_color[["none"]],
      breaks = kappa_breaks_rep,
      limits = kappa_limits_rep,
      labels = format_log10_p_value_labels(kappa_breaks_rep)
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(x = "", y = "", title = "Agreement between criteria (individual replication)", fill = "p-value") +
    bri_theme +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      title = element_text(size = 9),
      panel.border = element_blank()
    )

  # Save
  bri_ggsave(paste0(output_path, "cortable - replication - kappa.png"), plot = kappa_plot_rep, width = 8, height = 8)

  # Return the visualization
  return(kappa_plot_rep)
}

# Define the function
combine_kappa_plots <- function(plot_kappa_exp, plot_kappa_rep, output_path) {
  combined_kappa <- plot_grid(
    plotlist = list(
      plot_kappa_exp + theme(
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      ),
      plot_kappa_rep + theme(
        legend.position = "right",
        plot.margin = unit(c(0, 3, 0, 0), "cm")
      )
    ),
    ncol = 1,
    labels = "AUTO",
    rel_heights = c(1, 1)
  ) +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0, 0), "cm"), # 0.5, 0.5, -2, 0
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  bri_ggsave(paste0(output_path, "cortable - combined - kappa.png"), plot = combined_kappa, width = 5.5, height = 7)

  return(combined_kappa)
}

# Função para normalizar e ordenar a matriz ICC (com debug)
normalize_icc_matrix_effect_sizes <- function(m) {
  # Debug: estrutura inicial
  cat("[DEBUG] normalize_icc_matrix_effect_sizes: Initial columns:", paste(colnames(m), collapse = ", "), "\n")
  cat("[DEBUG] normalize_icc_matrix_effect_sizes: Initial dimensions:", nrow(m), "x", ncol(m), "\n")
  print(head(m, 3))
  
  # Verificar se EXP existe
  if (!"EXP" %in% colnames(m)) {
    stop("[ERROR] normalize_icc_matrix_effect_sizes: Coluna EXP não encontrada no data frame. Colunas presentes: ", paste(colnames(m), collapse = ", "))
  }
  
  # Verificar se orig existe
  if (!"orig" %in% colnames(m)) {
    stop("[ERROR] normalize_icc_matrix_effect_sizes: Coluna orig não encontrada no data frame. Colunas presentes: ", paste(colnames(m), collapse = ", "))
  }
  
  # Criar cópia do data frame
  result <- m
  
  # Preservar sinais originais - não aplicar abs() ao orig
  # Ordenar do maior para o menor valor original (mais positivo → mais negativo)
  order_idx <- order(result$orig, decreasing = TRUE)
  result <- result[order_idx, ]
  
  # Debug: estrutura final
  cat("[DEBUG] normalize_icc_matrix_effect_sizes: Final columns:", paste(colnames(result), collapse = ", "), "\n")
  cat("[DEBUG] normalize_icc_matrix_effect_sizes: Final dimensions:", nrow(result), "x", ncol(result), "\n")
  print(head(result, 3))
  
  # Adicionar atributo de normalização
  attr(result, "normalized_and_ordered") <- TRUE
  
  return(result)
}
