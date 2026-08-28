### ### ### ### ### ### ### ### ### ### ### ###

# aggregate.r
# This script contains the functions to aggregate the results of the individual analyses into measures of replication success and overall replication rates.

### ### ### ### ### ### ### ### ### ### ### ###

# Given a set of replications, calculates success for each replication
make_summary_subset <- function(
  rep_summaries_in,
  rep_summaries_individual_in,
  filtering = NULL,
  filtering_label = filtering
) {
  # Replication success by experiment (group of replications)
  rep_summaries <- rep_summaries_in
  if (!is.null(filtering)) {
    rep_summaries <- rep_summaries |>
      filter(
        str_detect(EXP, paste0("^", filtering))
      )
  }

  #From the summary of the meta-analyses, generate the large table of replication success by applying each criterion.
  rep_summaries <- rep_summaries |>
    mutate(
      # 1. Original estimate is within the 95% prediction interval of a random-effects meta-analysis of the replications.
      REP_Orig_in_REMA_PI = ifelse(
        is.na(rema_es),
        NA,
        original_es >= rema_pi_lower & original_es <= rema_pi_upper
      ),

      # 2. Estimate of a random-effects meta-analysis of the replications is within the 95% confidence interval of the original result.
      REP_REMA_in_Orig_CI = ifelse(
        is.na(rema_es),
        replication_es > original_ci_lower & replication_es < original_ci_upper,
        rema_es > original_ci_lower & rema_es < original_ci_upper
      ),

      # 3. Point estimate of a fixed-effects meta-analysis of the replications is significant and has the same signal as the original result.
      REP_FEMA_SSS_Orig = ifelse(
        is.na(fema_es),
        significant_reps > n_completed_reps / 2,
        (fema_pvalue < 0.05) & ((fema_es / original_es) > 0)
      ),

      # 4. A majority of the replications are significant and in the same signal as the original result.
      REP_Voting = (significant_reps > n_voting_reps / 2),
      # i.e. 3/3, 2/3, 2/2, 1/1

      REP_Voting_WithTies = (significant_reps >= n_voting_reps / 2),

      REP_Voting_Ties = (significant_reps == n_voting_reps / 2),

      # 5. Subjective assessment of the results, as judged by the teams who conducted each replication.
      REP_Subjective = (subjective_success_reps > n_voting_reps / 2),
      # i.e. 3/3, 2/3, 2/2, 1/1

      REP_Subjective_WithTies = (subjective_success_reps >= n_voting_reps / 2),

      REP_Subjective_Ties = (subjective_success_reps == n_voting_reps / 2),

      REP_Most_Criteria = ifelse(
        is.na(REP_Orig_in_REMA_PI),
        (REP_Voting +
          REP_Subjective +
          REP_FEMA_SSS_Orig +
          REP_REMA_in_Orig_CI) >=
          2,
        (REP_Voting +
          REP_Subjective +
          REP_FEMA_SSS_Orig +
          REP_REMA_in_Orig_CI +
          REP_Orig_in_REMA_PI) >=
          3
      ),

      REP_Most_Criteria_WithTies = ifelse(
        is.na(REP_Orig_in_REMA_PI),
        (REP_Voting_WithTies +
          REP_Subjective_WithTies +
          REP_FEMA_SSS_Orig +
          REP_REMA_in_Orig_CI) >=
          2,
        (REP_Voting_WithTies +
          REP_Subjective_WithTies +
          REP_FEMA_SSS_Orig +
          REP_REMA_in_Orig_CI +
          REP_Orig_in_REMA_PI) >=
          3
      ),

      # Is significant?
      RE_Significant = (rema_pvalue < 0.05),

      # Signal error: whether or not the original effect size has a different signal from that of the meta-analysis.
      SignalError = ifelse(
        RE_Significant,
        (replication_es / original_es) < 0,
        NA
      ),
      SignalErrorAll = ifelse(is.na(SignalError), F, SignalError),

      # Magnitude error: ratio between the original effect and the replication effect size.
      Exaggeration = ifelse(!SignalError, original_es / replication_es, NA),

      # Magnitude error: ratio between the original effect and the replication effect size, for all ES, even if not significant.
      ExaggerationAll = ifelse(
        !SignalErrorAll,
        original_es / replication_es,
        NA
      ),

      OriginalES = corrected_sign_original_es,
      OriginalCV = original_cv
    )

  # _____________________________________________

  # Replication success by individual replication
  rep_summaries_individual <- rep_summaries_individual_in |>
    filter(!ma_only)

  if (!is.null(filtering)) {
    rep_summaries_individual <- rep_summaries_individual |>
      filter(
        str_detect(EXP, paste0("^", filtering))
      )
  }

  # From the summary of the individual replications, generate the large table of replication success by applying each criterion.

  rep_summaries_individual <- rep_summaries_individual |>
    mutate(
      # 1. Replication estimate is within the 95% confidence interval of the original result.
      REP_Individual_Rep_in_Orig_CI = replication_es > original_ci_lower &
        replication_es < original_ci_upper,

      # 4. "Individual voting", i.e. significant and in the same signal as the original result.
      REP_Individual_SSS = (pvalue < 0.05) & (replication_es / original_es > 0),

      # 5. Individual replication subjective assessment of the results, as judged by the teams who conducted each replication.
      REP_Individual_IndivSubjective = subjective_success,

      # Is significant?
      Significant_Individual = (pvalue < 0.05),

      # Signal error: whether or not the original effect size has a different signal from that of the meta-analysis.
      SignalError_Individual = ifelse(
        Significant_Individual,
        (replication_es / original_es) < 0,
        NA
      ),

      SignalErrorAll_Individual = ifelse(
        is.na(SignalError_Individual),
        F,
        SignalError_Individual
      ),

      # Magnitude error: ratio between the original effect and the replication effect size, for significant ones.
      Exaggeration_Individual = ifelse(
        !SignalError_Individual,
        original_es / replication_es,
        NA
      ),

      # Magnitude error: ratio between the original effect and the replication effect size, for all ES, even if not significant.
      ExaggerationAll_Individual = ifelse(
        !SignalErrorAll_Individual,
        original_es / replication_es,
        NA
      ),

      OriginalES_Individual = original_es,
      OriginalCV_Individual = original_cv
    )

  print(paste0(
    "Calculating replication rates ",
    ifelse(is.null(filtering), "ALL_PCR_MTT", filtering),
    " ..."
  ))

  # Summarizing the reproducibility rates, by experiment and by individual replicate
  rep_success_complete <- rep_summaries |>
    select(
      EXP,
      REP_Orig_in_REMA_PI,
      REP_REMA_in_Orig_CI,
      REP_FEMA_SSS_Orig,
      REP_Voting,
      REP_Voting_WithTies,
      REP_Subjective,
      REP_Subjective_WithTies,
      REP_Most_Criteria_WithTies,
      SignalErrorAll,
      es_ratio,
      log_es_ratio,
      t_score_mean,
      t_score_FEMA,
      t_score,
      OriginalES,
      OriginalCV,
      mean_abs_diff_reps,
      mean_abs_diff_reps_orig,
      corrected_sign_replication_es,
      smd_effect_size
    )

  rep_success_individual_complete <- rep_summaries_individual |>
    select(
      LAB,
      EXP,
      REP_Individual_Rep_in_Orig_CI,
      REP_Individual_SSS,
      REP_Individual_IndivSubjective,
      Significant_Individual,
      SignalErrorAll_Individual,
      es_ratio_individual,
      log_es_ratio_individual,
      t_score_individual,
      OriginalES_Individual,
      OriginalCV_Individual,
      corrected_sign_replication_es,
      smd_effect_size
    )

  # Calculate agreement between the rates
  calc_agreement <- function(REP_AGREE) {
    # Agreement between the pairs of rates
    combs <- expand_grid(x = colnames(REP_AGREE), y = colnames(REP_AGREE))
    pair_agreement <- pmap_dfr(combs, function(x, y) {
      amat <- cbind(REP_AGREE[[x]], REP_AGREE[[y]])
      agreement <- irr::kappa2(amat, weight = "unweighted")
      tibble(
        `Outcome 1` = x,
        `Outcome 2` = y,
        Agreement = agreement$value,
        Method = "Cohen's Kappa",
        Call = 'irr::kappa2(DATA_FRAME, weight = "unweighted")'
      )
    })

    # General agreement
    general_agreement <- irr::kappam.fleiss(REP_AGREE, exact = F)

    rbind(
      tibble(
        `Outcome 1` = "General Agreement",
        `Outcome 2` = "General Agreement",
        Agreement = general_agreement$value,
        Method = "Fleiss' Kappa",
        Call = 'irr::kappam.fleiss(DATA_FRAME, exact = F)'
      ),
      pair_agreement
    )
  }

  # Agreement by experiment
  REP_AGREE <- rep_success_complete |>
    select(
      REP_Orig_in_REMA_PI,
      REP_REMA_in_Orig_CI,
      REP_FEMA_SSS_Orig,
      REP_Voting,
      REP_Subjective
    ) |>
    `colnames<-`(c(
      "Orig. in REMA PI",
      "REMA in Orig. CI",
      "FEMA SSS",
      "Voting",
      "Subjective"
    ))

  REP_AGREE <- REP_AGREE |>
    select(where(~ any(!is.na(.))))

  agreements_exp <- calc_agreement(REP_AGREE)

  # Agreement by replication
  REP_AGREE <- rep_success_individual_complete |>
    select(
      REP_Individual_Rep_in_Orig_CI,
      REP_Individual_SSS,
      REP_Individual_IndivSubjective
    ) |>
    `colnames<-`(c("Rep. in orig. CI", "Rep. SSS", "Rep. Subjective"))

  REP_AGREE <- REP_AGREE |>
    select(where(~ any(!is.na(.))))

  agreements_individual <- calc_agreement(REP_AGREE)

  # Binding individual and experiment level measures
  rep_success <- rbind(
    rep_success_complete |>
      select(EXP, where(is.logical)) |>
      pivot_longer(cols = -c(EXP)) |>
      group_by(name) |>
      filter(!is.na(value)) |>
      summarise(
        N = n(),
        successful = sum(value, na.rm = T),
        Value = mean(value, na.rm = T)
      ),
    rep_success_complete |>
      select(EXP, where(is.numeric)) |>
      pivot_longer(cols = -c(EXP)) |>
      group_by(name) |>
      filter(!is.na(value)) |>
      summarise(
        N = n(),
        successful = NA,
        Value = median(value, na.rm = T)
      )
  ) |>
    ungroup() |>
    mutate(Value = round(Value, 2)) |>
    `colnames<-`(c("Metric", "N", "successful", "Value"))

  rep_success_individual <- rbind(
    rep_success_individual_complete |>
      select(EXP, LAB, where(is.logical)) |>
      pivot_longer(cols = -c(EXP, LAB)) |>
      group_by(name) |>
      filter(!is.na(value)) |>
      summarise(
        N = n(),
        successful = sum(value, na.rm = T),
        Value = mean(value, na.rm = T)
      ),
    rep_success_individual_complete |>
      select(EXP, LAB, where(is.numeric)) |>
      pivot_longer(cols = -c(EXP, LAB)) |>
      group_by(name) |>
      filter(!is.na(value)) |>
      summarise(
        N = n(),
        successful = NA,
        Value = median(value, na.rm = T)
      )
  ) |>
    ungroup() |>
    mutate(Value = round(Value, 2)) |>
    `colnames<-`(c("Metric", "N", "successful", "Value")) |>
    mutate(
      Metric = ifelse(
        Metric %in% c("corrected_sign_replication_es", "smd_effect_size"),
        paste0(Metric, "_individual"),
        Metric
      )
    )

  success_summary <- rbind(rep_success, rep_success_individual)
  if (!is.null(filtering)) {
    success_summary <- success_summary |> mutate(Method = filtering_label)
  } else {
    success_summary <- success_summary |> mutate(Method = "ALL_PCR_MTT")
  }

  # Build matrix for calculating agreement
  r1 <- rep_summaries_individual |>
    group_by(EXP) |>
    mutate(R = 1:n()) |>
    pivot_wider(
      id_cols = c(EXP),
      names_from = R,
      values_from = replication_es
    )
  n_max_reps <- ncol(r1) - 1
  r1$avg <- rowMeans(r1[2:(2 + n_max_reps - 1)], na.rm = T)

  r2 <- rep_summaries_individual |>
    select(EXP, original_es) |>
    distinct(.keep_all = T) |>
    rename(orig = original_es)

  # Before left_join for r1 and r2, validate EXP column exists
  if (!"EXP" %in% colnames(r1)) {
    stop("[ERROR] r1 is missing the EXP column before left_join.")
  }
  if (!"EXP" %in% colnames(r2)) {
    stop("[ERROR] r2 is missing the EXP column before left_join.")
  }

  raters_matrix <- r1 |> left_join(r2)

  icc_summary <- map_dfr(c("EPM", "MTT", "PCR", "ALL_PCR_MTT"), function(mtd) {
    mrm <- raters_matrix
    if (mtd == "ALL_PCR_MTT") {
      mrm <- mrm |> filter(str_detect(EXP, "^(PCR|MTT|EPM)"))
    } else {
      mrm <- mrm |> filter(str_detect(EXP, mtd))
    }

    icc_reps_only <- icc(
      ratings = mrm[, 2:(2 + n_max_reps - 1)],
      model = "oneway"
    )
    icc_rep_orig <- icc(
      ratings = mrm[, (2 + n_max_reps):(2 + n_max_reps + 1)],
      model = "oneway"
    )

    rbind(
      tibble(
        ICC = "ICC Between Replications Only",
        method = mtd,
        model = icc_reps_only$model,
        type = icc_reps_only$type,
        unit = icc_reps_only$unit,
        name = icc_reps_only$name,
        value = icc_reps_only$value,
        ci95_lower = icc_reps_only$lbound,
        ci95_upper = icc_reps_only$ubound
      ),
      tibble(
        ICC = "ICC Average of Replications x Original",
        method = mtd,
        model = icc_rep_orig$model,
        type = icc_rep_orig$type,
        unit = icc_rep_orig$unit,
        name = icc_rep_orig$name,
        value = icc_rep_orig$value,
        ci95_lower = icc_rep_orig$lbound,
        ci95_upper = icc_rep_orig$ubound
      )
    )
  })

  # Return output to be combined
  list(
    complete_exp = rep_success_complete,
    complete_individual = rep_success_individual_complete,
    success = success_summary,
    agreements_exp = agreements_exp,
    agreements_individual = agreements_individual,
    icc_matrix = raters_matrix,
    icc_summary = icc_summary
  )
}

# Main function, summarizes rates overall (alternatively using PCR or ALTPCR) and by method
summarize_rates <- function(rep_summary_folder) {
  print("Evaluating replication results ...")

  rep_summaries <- fread(
    file = paste0(
      rep_summary_folder,
      "/Replication Assessment by Experiment.tsv"
    )
  )

  rep_summaries_individual <- fread(
    file = paste0(
      rep_summary_folder,
      "/Replication Assessment by Replication.tsv"
    )
  )

  dir.create(paste0(rep_summary_folder, "/additional-figures"))

  # Calculates aggregate rates summaries for every method subset
  rep_calculated <- list(
    ALL_PCR_MTT = make_summary_subset(
      rep_summaries,
      rep_summaries_individual,
      "(MTT|EPM|PCR)",
      "ALL_PCR_MTT"
    ),
    ALL_ALTPCR_MTT = make_summary_subset(
      rep_summaries,
      rep_summaries_individual,
      "(MTT|EPM|ALTPCR)",
      "ALL_ALTPCR_MTT"
    ),
    ALL_PCR_ALTMTT = make_summary_subset(
      rep_summaries,
      rep_summaries_individual,
      "(ALTMTT|EPM|PCR)",
      "ALL_PCR_ALTMTT"
    ),
    ALL_ALTPCR_ALTMTT = make_summary_subset(
      rep_summaries,
      rep_summaries_individual,
      "(ALTMTT|EPM|ALTPCR)",
      "ALL_ALTPCR_ALTMTT"
    ),
    MTT = make_summary_subset(rep_summaries, rep_summaries_individual, "MTT"),
    ALTMTT = make_summary_subset(
      rep_summaries,
      rep_summaries_individual,
      "ALTMTT"
    ),
    PCR = make_summary_subset(rep_summaries, rep_summaries_individual, "PCR"),
    ALTPCR = make_summary_subset(
      rep_summaries,
      rep_summaries_individual,
      "ALTPCR"
    ),
    EPM = make_summary_subset(rep_summaries, rep_summaries_individual, "EPM")
  )

  rep_success <- rbind(
    rep_calculated$ALL_PCR_MTT$success,
    rep_calculated$ALL_ALTPCR_MTT$success,
    rep_calculated$ALL_PCR_ALTMTT$success,
    rep_calculated$ALL_ALTPCR_ALTMTT$success,
    rep_calculated$MTT$success,
    rep_calculated$ALTMTT$success,
    rep_calculated$PCR$success,
    rep_calculated$ALTPCR$success,
    rep_calculated$EPM$success
  ) |>
    pivot_wider(
      id_cols = Metric,
      names_from = Method,
      values_from = c(Value, N, successful)
    )

  # Gets nicer names for metrics
  dict_criteria <- read_excel("./other-data/replication-criteria-dict.xlsx")
  rep_success <- rep_success |>
    left_join(dict_criteria, by = "Metric") |>
    mutate(
      OldMetric = Metric,
      Metric = Metric |> str_remove("_Individual")
    ) |>
    mutate(Metric = OldMetric) |>
    select(-OldMetric)

  # Plots effect size correlations (with variations on including outliers, and PCR alternatives)
  print("Plotting correlation between effects ...")

  p_B <- plot_effect_correlation(
    rep_summary_folder,
    rep_summaries |> filter(!str_detect(EXP, "ALT")),
    "effect size correlation by exps.png",
    "Effect size correlation, by experiment"
  )

  plot_effect_correlation(
    rep_summary_folder,
    rep_summaries_individual |> filter(!str_detect(EXP, "ALT")),
    "effect size correlation by indiv reps.png",
    "Effect size correlation, by individual replication"
  )

  p_C <- plot_effect_correlation(
    rep_summary_folder,
    rep_summaries |>
      filter(original_es < 5, replication_es < 5) |>
      filter(!str_detect(EXP, "ALT")),
    "effect size correlation by exps small ES.png",
    "Effect size correlation, by experiment (ES < 5)"
  )

  plot_effect_correlation(
    rep_summary_folder,
    rep_summaries_individual |>
      filter(original_es < 5, replication_es < 5) |>
      filter(!str_detect(EXP, "ALT")),
    "effect size correlation by indiv reps small ES.png",
    "Effect size correlation, by individual replication (ES < 5)"
  )

  plot_effect_correlation(
    rep_summary_folder,
    rep_summaries |> filter(!str_detect(EXP, "^PCR")),
    "effect size correlation by exps with ALTPCR.png",
    "Effect size correlation, by experiment"
  )

  plot_effect_correlation(
    rep_summary_folder,
    rep_summaries_individual |> filter(!str_detect(EXP, "^PCR")),
    "effect size correlation by indiv reps with ALTPCR.png",
    "Effect size correlation, by individual replication"
  )

  plot_effect_correlation(
    rep_summary_folder,
    rep_summaries |>
      filter(original_es < 5, replication_es < 5) |>
      filter(!str_detect(EXP, "^PCR")),
    "effect size correlation by exps small ES with ALTPCR.png",
    "Effect size correlation, by experiment (ES < 5)"
  )

  plot_effect_correlation(
    rep_summary_folder,
    rep_summaries_individual |>
      filter(original_es < 5, replication_es < 5) |>
      filter(!str_detect(EXP, "^PCR")),
    "effect size correlation by indiv reps small ES with ALTPCR.png",
    "Effect size correlation, by individual replication (ES < 5)"
  )

  # Plot summary of effect sizes
  print("Plotting log ES ratios ...")

  plot_log_es_ratios(
    rep_summaries |> filter(!str_detect(EXP, "ALT")),
    "Log ES Ratios by Method.png",
    rep_summary_folder
  )
  plot_log_es_ratios(
    rep_summaries |> filter(!str_detect(EXP, "^PCR")),
    "Log ES Ratios by Method with ALTPCR.png",
    rep_summary_folder
  )

  # Plot coefficient of variations
  print("Plotting CVs ...")

  plot_cvs(
    rep_summaries_individual |> filter(!str_detect(EXP, "ALT")),
    original_experiment_data,
    " - by REP",
    rep_summary_folder
  )
  # ALTPCR-only view (não incluir ALTMTT para não misturar com MTT)
  plot_cvs(
    rep_summaries_individual |> filter(!str_detect(EXP, "^(PCR|ALTMTT)")),
    original_experiment_data,
    " - by REP, ALTPCR",
    rep_summary_folder
  )
  # ALTMTT-only view (não incluir MTT base nem ALTPCR)
  plot_cvs(
    rep_summaries_individual |> filter(!str_detect(EXP, "^(MTT|ALTPCR)")),
    original_experiment_data,
    " - by REP, ALTMTT",
    rep_summary_folder
  )

  plot_cvs(
    rep_summaries |> filter(!str_detect(EXP, "ALT")),
    original_experiment_data,
    " - by EXP",
    rep_summary_folder
  )
  # ALTPCR-only view (não incluir ALTMTT para não misturar com MTT)
  plot_cvs(
    rep_summaries |> filter(!str_detect(EXP, "^(PCR|ALTMTT)")),
    original_experiment_data,
    " - EXP, ALTPCR",
    rep_summary_folder
  )
  # ALTMTT-only view (não incluir MTT base nem ALTPCR)
  plot_cvs(
    rep_summaries |> filter(!str_detect(EXP, "^(MTT|ALTPCR)")),
    original_experiment_data,
    " - EXP, ALTMTT",
    rep_summary_folder
  )

  # Plot a matrix of estimates to visually represent intraclass correlations
  print("Plotting ICC ...")

  # Normaliza e ordena a matriz ICC para a figura A (effect sizes originais sempre positivos e ordenados do menor para o maior)
  # Debug: verificar estrutura antes da normalização
  print("Estrutura antes da normalização:")
  print(colnames(rep_calculated$ALL_PCR_MTT$icc_matrix))
  print(head(rep_calculated$ALL_PCR_MTT$icc_matrix))
  # Normaliza e ordena a matriz ICC para a figura A (effect sizes originais sempre positivos e ordenados do menor para o maior)
  normalized_ordered_icc_matrix <- normalize_icc_matrix_effect_sizes(
    rep_calculated$ALL_PCR_MTT$icc_matrix
  )
  # --- START DEBUGGING BLOCK ---
  # Add this code to inspect the object before the error
  cat("\n\n--- DEBUGGING INFO ---\n")
  cat("Object Class:\n")
  print(class(normalized_ordered_icc_matrix))
  cat("\nColumn Names:\n")
  print(colnames(normalized_ordered_icc_matrix))
  cat("\nFirst 3 Rows:\n")
  print(head(normalized_ordered_icc_matrix, 3))
  cat("--- END DEBUGGING ---\n\n")
  # --- END DEBUGGING BLOCK ---

  p_A <- plot_icc_matrix(
    m = normalized_ordered_icc_matrix,
    fn = paste0(
      rep_summary_folder,
      "/additional-figures/All Replication Estimates.png"
    ),
    cap_at = NULL,
    tt = "Effect size comparison, by replication"
  )

  plot_icc_matrix(
    normalized_ordered_icc_matrix,
    cap_at = 5,
    paste0(
      rep_summary_folder,
      "/additional-figures/All Replication Estimates (ES less than 5).png"
    )
  )

  # Plot combined figure for better presentation
  plot_effect_size_figure(p_A, p_B, p_C, rep_summary_folder)

  # Write output
  write_tsv(
    rep_calculated$ALL_PCR_MTT$complete_exp,
    paste0(rep_summary_folder, "/Replication Success by Experiment.tsv")
  )
  write_tsv(
    rep_calculated$ALL_PCR_MTT$complete_individual,
    paste0(rep_summary_folder, "/Replication Success by Replication.tsv")
  )

  write_tsv(
    rep_success,
    paste0(rep_summary_folder, "/Replication Rate Summary.tsv")
  )
}

# Function to run through each analysis subfolder and gather the replication rates, then put it in a long-format - used for specification curves
gather_all_rep_rates <- function(list_of_analyses, pn) {
  # All variables
  all_replication_rates <- map_dfr(list_of_analyses, function(p) {
    rep_rates <- fread(paste0(p, "/Replication Rate Summary.tsv"))

    rep_rates <- rep_rates |>
      select(-MetricLongName) |>
      mutate(
        Level = ifelse(
          str_detect(Metric, "Individual"),
          "Individual replication",
          "Set of replications"
        ),
        Metric = Metric |> str_remove("_Individual")
      ) |>
      pivot_longer(
        cols = -c(Metric, MetricShortName, Level)
      ) |>
      mutate(
        Method = name |> str_remove("^(Value|N|successful)_"),
        What = name |> str_extract("^(Value|N|successful)_") |> str_remove("_")
      ) |>
      select(-name) |>
      pivot_wider(
        id_cols = c(Level, Metric, MetricShortName, Method),
        names_from = What,
        values_from = value
      ) |>
      mutate(
        Inclusion_Set = (basename(p) |>
          str_trim() |>
          str_split(" ") |>
          unlist())[1],
        MA_Dist = (basename(p) |> str_trim() |> str_split(" ") |> unlist())[2],
        IsPercentage = !(Metric %in%
          c(
            "Difference",
            "Exaggeration",
            "ExaggerationAll",
            "ESRatio",
            "SignalError",
            "SignalErrorAll",
            "t_score_mean",
            "t_score_mean_transf",
            "t_score_FEMA",
            "t_score",
            "t_score_transf",
            "t_score_individual",
            "OriginalES",
            "es_ratio",
            "log_es_ratio",
            "es_ratio_individual",
            "log_es_ratio_individual",
            "mean_abs_diff_reps",
            "mean_abs_diff_reps_orig",
            "original_cv",
            "OriginalCV"
          ))
      ) |>
      select(
        Inclusion_Set,
        Level,
        MA_Dist,
        Metric,
        MetricShortName,
        IsPercentage,
        Method,
        Value,
        N,
        successful
      ) |>
      filter(
        IsPercentage,
        Metric != "REP_Most_Criteria_WithTies",
        Metric != "Significant"
      )

    rep_rates
  })

  write_tsv(
    all_replication_rates,
    paste0(pn, "/Replication Rate Aggregate.tsv")
  )

  # Without aggregation by method, only analysis specifications
  all_replication_rates <- map_dfr(list_of_analyses, function(p) {
    rep_rates <- fread(paste0(p, "/Replication Rate Summary.tsv"))
    rep_rates <- rep_rates |>
      select(
        -any_of(c(
          "Value_ALL",
          "Value_MTT",
          "Value_PCR",
          "Value_EPM",
          "Value_ALTPCR",
          "Value_ALTMTT",
          "N_ALL",
          "N_MTT",
          "N_PCR",
          "N_EPM",
          "N_ALTPCR",
          "N_ALTMTT",
          "successful_ALL",
          "successful_MTT",
          "successful_PCR",
          "successful_EPM",
          "successful_ALTPCR",
          "successful_ALTMTT"
        ))
      )

    rep_rates <- rep_rates |>
      select(-MetricLongName) |>
      mutate(
        Level = ifelse(
          str_detect(Metric, "Individual"),
          "Individual replication",
          "Set of replications"
        ),
        Metric = Metric |> str_remove("_Individual")
      ) |>
      pivot_longer(
        cols = -c(Metric, MetricShortName, Level)
      ) |>
      mutate(
        Method = name |> str_remove("^(Value|N|successful)_"),
        What = name |> str_extract("^(Value|N|successful)_") |> str_remove("_")
      ) |>
      select(-name) |>
      pivot_wider(
        id_cols = c(Level, Metric, MetricShortName, Method),
        names_from = What,
        values_from = value
      ) |>
      mutate(
        Inclusion_Set = (basename(p) |>
          str_trim() |>
          str_split(" ") |>
          unlist())[1],
        MA_Dist = (basename(p) |> str_trim() |> str_split(" ") |> unlist())[2],
        IsPercentage = !(Metric %in%
          c(
            "Difference",
            "Exaggeration",
            "ExaggerationAll",
            "ESRatio",
            "SignalError",
            "SignalErrorAll",
            "t_score_mean",
            "t_score_mean_transf",
            "t_score_FEMA",
            "t_score",
            "t_score_transf",
            "t_score_individual",
            "OriginalES",
            "es_ratio",
            "log_es_ratio",
            "es_ratio_individual",
            "log_es_ratio_individual",
            "mean_abs_diff_reps",
            "mean_abs_diff_reps_orig",
            "original_cv",
            "OriginalCV"
          ))
      ) |>
      select(
        Inclusion_Set,
        Level,
        MA_Dist,
        Metric,
        MetricShortName,
        IsPercentage,
        Method,
        Value,
        N,
        successful
      ) |>
      filter(
        IsPercentage,
        Metric != "REP_Most_Criteria_WithTies",
        Metric != "Significant"
      )

    rep_rates
  })

  write_tsv(
    all_replication_rates,
    paste0(pn, "/Replication Rate Aggregate - without Method.tsv")
  )
}
