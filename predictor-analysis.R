### ### ### ### ### ### ### ### ### ### ### ###

# predictor-analysis.r
# This script contains the functions to run the analyses of predictors of reproducibility success.

### ### ### ### ### ### ### ### ### ### ### ###

# Function to run the experiment level predictor analysis
run_predictor_analysis_exp_level <- function(replication_datapath, make_plots) {
  # Gets replication outcomes
  REP_DATA <- read_tsv(paste0(replication_datapath, "/Replication Success by Experiment.tsv"))

  REP_DATA <- REP_DATA |>
    select(
      EXP,
      REP_Orig_in_REMA_PI,
      REP_REMA_in_Orig_CI,
      REP_FEMA_SSS_Orig,
      REP_Voting_WithTies,
      REP_Subjective_WithTies,
      SignalErrorAll,
      t_score,
      log_es_ratio,
      corrected_sign_replication_es,
      smd_effect_size,
      OriginalCV
    )

  colnames(REP_DATA) <- c(
    "EXP",
    "Original in replication PI",
    "Replication in original CI",
    "Same sense significance",
    "Voting (with ties)",
    "Subjective (with ties)",
    "Signal error",
    "t Value",
    "Log ES Ratio",
    "Replication Effect Size",
    "SMD Effect Size",
    "Original CV"
  )

  # Merge with predictor data
  dir.create(paste0(replication_datapath, "/predictors/"))
  dir.create(paste0(replication_datapath, "/predictors/by experiment"))
  FDATA <- REP_DATA |> left_join(PRED_DATA_EXP_LEVEL, by = "EXP")
  write_tsv(FDATA, paste0(replication_datapath, "/predictors/by experiment/Table for Predictor Correlations.tsv"))

  # List predictor and outcome variables
  predictors <- c(unique(colnames(PRED_DATA_EXP_LEVEL |> select(-EXP, -DOI))))
  predictors <- append(predictors, "Original CV", after = which(predictors == "p-value"))

  categorical_rep_outcome <- c(
    "Original in replication PI",
    "Replication in original CI",
    "Same sense significance",
    "Voting (with ties)",
    "Subjective (with ties)"
  )
  continuous_rep_outcome <- c(
    "t Value",
    "Log ES Ratio",
    "Replication Effect Size",
    "SMD Effect Size"
  )

  # Correlation table and plot between predictors
  print("Building correlation table, predictors x predictors ...")
  vars1 <- predictors
  vars2 <- vars1

  p_pred <- plot_cortable(FDATA, vars1, vars2, "Correlation between predictors", paste0(replication_datapath, "/predictors/by experiment/predictor x predictor"), show_label = T, make_individual_plots = F)


  # Correlation table and plot between replication criteria
  print("Building correlation table, replication outcomes x replication outcomes ...")
  vars1 <- c(continuous_rep_outcome, categorical_rep_outcome)
  vars2 <- vars1

  p_rep <- plot_cortable(FDATA, vars1, vars2, "Correlation between replication criteria", paste0(replication_datapath, "/predictors/by experiment/rep outcome x rep outcome"), show_label = T, make_individual_plots = F)


  # Correlation table and plot between predictors and replication outcomes
  print("Building correlation table, predictors x replication outcomes ...")
  vars1 <- predictors
  vars2 <- c(continuous_rep_outcome, categorical_rep_outcome)

  p_pred_rep <- plot_cortable_cluster(FDATA, vars1, vars2, "Correlation between predictors and replication", paste0(replication_datapath, "/predictors/by experiment/predictor x rep outcome"), show_label = T, make_individual_plots = make_plots)

  return(list(p_pred, p_rep, p_pred_rep))
}

# Function to run the replication level predictor analysis
run_predictor_analysis_rep_level <- function(replication_datapath, make_plots) {
  # Gets replication outcomes
  REP_DATA <- read_tsv(paste0(replication_datapath, "/Replication Success by Replication.tsv"))

  REP_DATA <- REP_DATA |>
    select(
      EXP,
      LAB,
      REP_Individual_Rep_in_Orig_CI,
      REP_Individual_SSS,
      REP_Individual_IndivSubjective,
      SignalErrorAll_Individual,
      t_score_individual,
      log_es_ratio_individual,
      corrected_sign_replication_es,
      smd_effect_size
    )

  colnames(REP_DATA) <- c(
    "EXP",
    "LAB",
    "Replication in original CI",
    "Same sense significance",
    "Subjective",
    "Signal error",
    "t Value",
    "Log ES Ratio",
    "Replication Effect Size",
    "SMD Effect Size"
  )

  # Merge with predictor data
  dir.create(paste0(replication_datapath, "/predictors/"))
  dir.create(paste0(replication_datapath, "/predictors/by replication/"))
  FDATA <- REP_DATA |> left_join(PRED_DATA_REP_LEVEL, by = c("EXP", "LAB"))
  write_tsv(FDATA, paste0(replication_datapath, "/predictors/by replication/Table for Predictor Correlations.tsv"))

  # List predictor and outcome variables
  predictors <- unique(colnames(PRED_DATA_REP_LEVEL |> select(-EXP, -LAB)))

  categorical_rep_outcome <- c(
    "Replication in original CI",
    "Same sense significance",
    "Subjective"
  )
  continuous_rep_outcome <- c(
    "t Value",
    "Log ES Ratio",
    "Replication Effect Size",
    "SMD Effect Size"
  )

  # Correlation table and plot between predictors
  print("Building correlation table, predictors x predictors ...")
  vars1 <- predictors
  vars2 <- vars1

  p_pred <- plot_cortable(FDATA, vars1, vars2, "Correlation between predictors", paste0(replication_datapath, "/predictors/by replication/predictor x predictor"), show_label = T, make_individual_plots = F)


  # Correlation table and plot between replication criteria
  print("Building correlation table, replication outcomes x replication outcomes ...")
  vars1 <- c(continuous_rep_outcome, categorical_rep_outcome)
  vars2 <- vars1

  p_rep <- plot_cortable(FDATA, vars1, vars2, "Correlation between replication criteria", paste0(replication_datapath, "/predictors/by replication/rep outcome x rep outcome"), show_label = T, make_individual_plots = F)


  # Correlation table and plot between predictors and replication outcomes
  print("Building correlation table, predictors x replication outcomes ...")
  vars1 <- predictors
  vars2 <- c(continuous_rep_outcome, categorical_rep_outcome)

  p_pred_rep <- plot_cortable_cluster(FDATA, vars1, vars2, "Correlation between predictors and replication", paste0(replication_datapath, "/predictors/by replication/predictor x rep outcome"), show_label = T, make_individual_plots = make_plots)

  return(list(p_pred, p_rep, p_pred_rep))
}

# Function to run the predictor analysis in both levels (EXP and REP)
run_predictor_analysis_both_level <- function(replication_datapath, make_plots) {
  dir.create(paste0(replication_datapath, "/predictors/"))
  dir.create(paste0(replication_datapath, "/predictors/both"))

  # EXP - Gets predictor data
  FDATA_EXP <- read_tsv(paste0(replication_datapath, "/predictors/by experiment/Table for Predictor Correlations.tsv")) |>
    select(-any_of(c("Original in replication PI", "Replication in original CI", "Same sense significance", "Voting (with ties)", "Subjective (with ties)", "Signal error", "t Value", "Log ES Ratio", "Replication Effect Size", "SMD Effect Size", "DOI"))) |>
    relocate(`Original CV`, .after = `p-value`) |>
    rename(`University Ranking (Experiment)` = `Institution Ranking`)

  # REP - Gets predictor data
  FDATA_REP <- read_tsv(paste0(replication_datapath, "/predictors/by replication/Table for Predictor Correlations.tsv")) |>
    select(-any_of(c("Replication in original CI", "Same sense significance", "Subjective", "Signal error", "t Value", "Log ES Ratio", "Replication Effect Size", "SMD Effect Size"))) |>
    rename(`University Ranking (Replication)` = `Institution Ranking`)

  # Join
  FDATA <- full_join(FDATA_EXP, FDATA_REP, by = "EXP")
  write_tsv(FDATA, paste0(replication_datapath, "/predictors/both/Table for Predictor Correlations.tsv"))

  # List predictor and outcome variables
  global_vars <- setdiff(colnames(FDATA_EXP), c("EXP", "Original CV"))
  global_vars <- append(global_vars, "Original CV", after = which(global_vars == "p-value"))
  individual_vars <- setdiff(colnames(FDATA_REP), c("EXP", "LAB"))

  # Correlation table and plot between predictors
  print("Building correlation table with both levels, predictors x predictors ...")
  vars1 <- c(global_vars, individual_vars)
  vars2 <- c(global_vars, individual_vars)

  # Make the plot
  plot_cortable_alternative(FDATA, FDATA_EXP, FDATA_REP, vars1, vars2, "Correlation between predictors", paste0(replication_datapath, "/predictors/both/predictor x predictor"), show_label = T, make_individual_plots = F)
}

# Function to run predictor x predictor analysis for ALL experiments (Done=Yes, UNIT=BRI)
# This ensures correlations between predictors include all 143 replications, not just those in analysis sets.
# Accepts a vector of paths and saves output to all of them (runs analysis only once)
# Uses global inclusion_sets variable (must be already loaded in main.R)
run_predictor_analysis_all_exps <- function(replication_datapaths) {
  # Read CV data from replication assessment file (contains original_cv column)
  # Use first path to read the data
  cv_data <- fread(paste0(replication_datapaths[1], "/Replication Assessment by Experiment.tsv")) |>
    select(EXP, original_cv) |>
    distinct() |>
    filter(str_detect(EXP, "^(MTT|EPM|PCR)")) |>
    rename(`Original CV` = original_cv)
    
  # Filter to get all 143 replications
  all_exps <- inclusion_sets |>
    filter(UNIT == "BRI") |>
    filter(done == "Yes") |>
    select(EXP, LAB) |>
    distinct()
  
  # Prepare experiment-level predictor data (all 60 experiments)
  FDATA_EXP <- PRED_DATA_EXP_LEVEL |>
    filter(EXP %in% unique(all_exps$EXP)) |>
    select(-DOI) |>
    rename(`University Ranking (Experiment)` = `Institution Ranking`) |>
    left_join(cv_data, by = "EXP")
    
  # Prepare replication-level predictor data (all 143 replications)
  FDATA_REP <- PRED_DATA_REP_LEVEL |>
    inner_join(all_exps, by = c("EXP", "LAB")) |>
    rename(`University Ranking (Replication)` = `Institution Ranking`)
    
  # Combine for correlation analysis
  FDATA <- full_join(FDATA_EXP, FDATA_REP, by = "EXP")
    
  # List predictor variables (position "Original CV" after "p-value" like other predictor analyses)
  global_vars <- setdiff(colnames(FDATA_EXP), c("EXP", "Original CV"))
  global_vars <- append(global_vars, "Original CV", after = which(global_vars == "p-value"))
  individual_vars <- setdiff(colnames(FDATA_REP), c("EXP", "LAB"))
  vars1 <- c(global_vars, individual_vars)
  vars2 <- vars1
  
  # Generate predictor x predictor plot (compute once)
  print("Building correlation table for ALL experiments, predictors x predictors ...")
  
  # Use first path for plot generation (save to 'both' folder alongside other combined predictor analyses)
  first_path <- replication_datapaths[1]
  dir.create(paste0(first_path, "/predictors/both"), showWarnings = FALSE, recursive = TRUE)
  
  p <- plot_cortable_alternative(FDATA, FDATA_EXP, FDATA_REP, vars1, vars2, 
    "Correlation between predictors (all experiments)", 
    paste0(first_path, "/predictors/both/predictor x predictor (all experiments)"), 
    show_label = TRUE, make_individual_plots = FALSE)
  
  # Save TSV to first path
  write_tsv(FDATA, paste0(first_path, "/predictors/both/Table for Predictor Correlations (all experiments).tsv"))
  
  # Copy outputs to all other paths
  if (length(replication_datapaths) > 1) {
    for (path in replication_datapaths[-1]) {
      dir.create(paste0(path, "/predictors/both"), showWarnings = FALSE, recursive = TRUE)
      
      # Copy TSV
      file.copy(
        paste0(first_path, "/predictors/both/Table for Predictor Correlations (all experiments).tsv"),
        paste0(path, "/predictors/both/Table for Predictor Correlations (all experiments).tsv"),
        overwrite = TRUE
      )
      
      # Copy plots
      file.copy(
        paste0(first_path, "/predictors/both/predictor x predictor (all experiments) - spearman.png"),
        paste0(path, "/predictors/both/predictor x predictor (all experiments) - spearman.png"),
        overwrite = TRUE
      )
      file.copy(
        paste0(first_path, "/predictors/both/predictor x predictor (all experiments) - pearson.png"),
        paste0(path, "/predictors/both/predictor x predictor (all experiments) - pearson.png"),
        overwrite = TRUE
      )
      file.copy(
        paste0(first_path, "/predictors/both/predictor x predictor (all experiments) cordata.tsv"),
        paste0(path, "/predictors/both/predictor x predictor (all experiments) cordata.tsv"),
        overwrite = TRUE
      )
    }
    message(paste0("Copied outputs to ", length(replication_datapaths) - 1, " additional folders"))
  }
  
  return(p)
}

# Function to run a correlation between two predictors/outcomes

run_correlation = function (FDATA, c1, c2, fn, make_plot) {
  print(paste(c1, c2, sep = " X "))
  
  col1 = FDATA[[c1]]
  col2 = FDATA[[c2]]
  
  c1_cat = is.character(col1) | is.logical(col1)
  c2_cat = is.character(col2) | is.logical(col2)
  
  if (c1_cat) col1 = as.integer(as.factor(col1)) - 1
  if (c2_cat) col2 = as.integer(as.factor(col2)) - 1
  
  df = tibble(Method = FDATA$EXP |> str_extract("(MTT|EPM|PCR)"), vx = col1, vy = col2) |>
    filter(!is.na(vx), !is.na(vy))
  
  if (length(unique(df$vx)) < 2 | length(unique(df$vy)) < 2 | nrow(df) < 3) {
    return (
      tibble(
        col1 = c1,
        col2 = c2,
        rho = NA,
        p.value = NA,
        pearson_r = NA,
        pearson_p = NA,
        n = nrow(df)
      )
    )
  }
  
  corr = cor.test(col1, col2, method = "spearman", use = "pairwise.complete.obs")
  corr_pearson = cor.test(col1, col2, method = "pearson", use = "pairwise.complete.obs")
  
  if (make_plot) {
    plot_scatterplot_predictors(df, corr, c1, c2, c1_cat, c2_cat, fn)
  }
  
  tibble(
    col1 = c1,
    col2 = c2,
    rho = unname(corr$estimate),
    p.value = corr$p.value,
    pearson_r = unname(corr_pearson$estimate),
    pearson_p = corr_pearson$p.value,
    n = nrow(df)
  )
}

# Function to run a correlation between predictors and outcomes simultaneously (this function is used inside the function run_predictor_analysis_both_level)
run_correlation_alternative = function (FDATA, c1, c2, fn, make_plot, global_vars) {
  print(paste(c1, c2, sep = " X "))
  
  if(c1 %in% global_vars & c2 %in% global_vars){
    FDATA <-  FDATA |>
      group_by(EXP) |>
      summarise(across(all_of(global_vars), ~ mean(.x, na.rm = TRUE))) |>
      ungroup()
    
    col1 = FDATA[[c1]]
    col2 = FDATA[[c2]]
    
    c1_cat = is.character(col1) | is.logical(col1)
    c2_cat = is.character(col2) | is.logical(col2)
    
    if (c1_cat) col1 = as.integer(as.factor(col1)) - 1
    if (c2_cat) col2 = as.integer(as.factor(col2)) - 1
    
    df = tibble(Method = FDATA$EXP |> str_extract("(MTT|EPM|PCR)"), vx = col1, vy = col2) |>
      filter(!is.na(vx), !is.na(vy))

    if (length(unique(df$vx)) < 2 | length(unique(df$vy)) < 2 | nrow(df) < 3) {
      return (
        tibble(
          col1 = c1,
          col2 = c2,
          rho = NA,
          p.value = NA,
          pearson_r = NA,
          pearson_p = NA,
          n = nrow(df)
        )
      )
    }
    
    corr = cor.test(col1, col2, method = "spearman", use = "pairwise.complete.obs")
    corr_pearson = cor.test(col1, col2, method = "pearson", use = "pairwise.complete.obs")
    
    if (make_plot) {
      plot_scatterplot_predictors(df, corr, c1, c2, c1_cat, c2_cat, fn)
    }
    
    tibble(
      col1 = c1,
      col2 = c2,
      rho = unname(corr$estimate),
      p.value = corr$p.value,
      pearson_r = unname(corr_pearson$estimate),
      pearson_p = corr_pearson$p.value,
      n = nrow(df)
    )
  } else{
    col1 = FDATA[[c1]]
    col2 = FDATA[[c2]]
    
    c1_cat = is.character(col1) | is.logical(col1)
    c2_cat = is.character(col2) | is.logical(col2)
    
    if (c1_cat) col1 = as.integer(as.factor(col1)) - 1
    if (c2_cat) col2 = as.integer(as.factor(col2)) - 1
    
    df = tibble(Method = FDATA$EXP |> str_extract("(MTT|EPM|PCR)"), vx = col1, vy = col2) |>
      filter(!is.na(vx), !is.na(vy))
    
    if (length(unique(df$vx)) < 2 | length(unique(df$vy)) < 2 | nrow(df) < 3) {
      return (
        tibble(
          col1 = c1,
          col2 = c2,
          rho = NA,
          p.value = NA,
          pearson_r = NA,
          pearson_p = NA,
          n = nrow(df)
        )
      )
    }
    
    corr = cor.test(col1, col2, method = "spearman", use = "pairwise.complete.obs")
    corr_pearson = cor.test(col1, col2, method = "pearson", use = "pairwise.complete.obs")
    
    if (make_plot) {
      plot_scatterplot_predictors(df, corr, c1, c2, c1_cat, c2_cat, fn)
    }
    
    tibble(
      col1 = c1,
      col2 = c2,
      rho = unname(corr$estimate),
      p.value = corr$p.value,
      pearson_r = unname(corr_pearson$estimate),
      pearson_p = corr_pearson$p.value,
      n = nrow(df)
    )
  }
}
