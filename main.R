# Analysis protocol is available here: https://osf.io/9rnuj

# From "Extracted Data from Original Experiments.xlsx", generate "Original Experiment Statistical Summaries.tsv". This is calculated from the values extracted in the original papers. This is so that we have intervals for the original experiments. Some N are inputted manually in the table because they were not available in the original report.
source("original-precision-calculation.R")

# Data from every Data Collection table was manually gathered in one single table beforehand. In these tables, the experimental units and outcomes are also indicated. These are in the folder /gathered-data.

# For each experimental technique, run the summarizers, which joins each experiment and the original data in the same file, summarizing the data from each in the way that will be used to run the meta-analyses.

dir.create("./replication-results")
dir.create("./replication-results/For labs")

# MTT
dir.create("./replication-results/MTT")
dir.create("./replication-results/For labs/MTT")
source("mtt-summarizer.R")

dir.create("./replication-results/MTT_ALT")
dir.create("./replication-results/For labs/MTT_ALT")
source("mtt-summarizer-alt.R")

# PCR
dir.create("./replication-results/PCR")
dir.create("./replication-results/For labs/PCR")
source("pcr-summarizer.R")

# PCR Alternative
dir.create("./replication-results/PCR_ALT")
dir.create("./replication-results/For labs/PCR_ALT")
source("pcr-summarizer-alt.R")

# EPM
dir.create("./replication-results/EPM")
dir.create("./replication-results/For labs/EPM")
dir.create("./replication-results/EPM/primary_outcomes")
dir.create("./replication-results/For labs/EPM/primary_outcomes")
dir.create("./replication-results/EPM/secondary_outcomes")
dir.create("./replication-results/For labs/EPM/secondary_outcomes")
source("epm-summarizer.R")

##### Run the analysis script to perform the first part, the meta-analysis to estimate replication rates. Plots and datasets with the results will be saved in the folder output, with a subfolder for each subset. ####

source("analysis.R")

results_path <- paste0(today())
dir.create(paste0("./output/", results_path))

plot_original_sd_mean_regression(
  data_path = "./other-data/Original Experiments Statistical Summaries.tsv",
  output_dir = file.path("./output", results_path, "original-error-bar-regression")
)

# Every column that starts with "analysis" will be considered an inclusion set column and must only have either cells marked "INCLUDE" or empty cells (not included).
inclusion_sets <- read_excel("./other-data/inclusion_sets.xlsx", 1)

list_of_sets <- colnames(inclusion_sets)[str_detect(
  colnames(inclusion_sets),
  "^analysis_"
)]

for (setcol in list_of_sets) {
  message("\n\n\n\n")
  message(setcol)

  if (
    !(setcol %in%
      c(
        "analysis_only_80_power_a_posteriori_T",
        "analysis_only_80_power_a_posteriori_KNHA",
        "analysis_only_80_power_a_posteriori_bigexp"
      ))
  ) {
    # Run analysis with z distribution only with inclusion sets also using z
    run_all_meta_analyses(
      inclusion_set_column = setcol,
      save_results_to = results_path,
      params = list(exclude_outliers = F, ma_dist = "z"),
      simulated = F
    )
  }

  if (
    !(setcol %in%
      c(
        "analysis_only_80_power_a_posteriori_Z",
        "analysis_only_80_power_a_posteriori_KNHA",
        "analysis_only_80_power_a_posteriori_bigexp"
      ))
  ) {
    # Run analysis with t distribution only with inclusion sets also using t
    run_all_meta_analyses(
      inclusion_set_column = setcol,
      save_results_to = results_path,
      params = list(exclude_outliers = F, ma_dist = "t"),
      simulated = F
    )
  }

  if (
    !(setcol %in%
      c(
        "analysis_only_80_power_a_posteriori_Z",
        "analysis_only_80_power_a_posteriori_T",
        "analysis_only_80_power_a_posteriori_bigexp"
      ))
  ) {
    # Run analysis with t distribution only with inclusion sets also using t
    run_all_meta_analyses(
      inclusion_set_column = setcol,
      save_results_to = results_path,
      params = list(exclude_outliers = F, ma_dist = "knha"),
      simulated = F
    )
  }

  if (
    !(setcol %in%
      c(
        "analysis_only_80_power_a_posteriori_Z",
        "analysis_only_80_power_a_posteriori_KNHA"
      ))
  ) {
    # Run analysis with bigexp and only with inclusion sets also using t
    run_all_meta_analyses(
      inclusion_set_column = setcol,
      save_results_to = results_path,
      params = list(exclude_outliers = F, ma_dist = "bigexp"),
      simulated = F
    )
  }
}

# For every folder specified below with an analysis in the output folder, calculates the replications rates
list_of_analyses <- list.dirs(
  path = paste0("./output/", results_path),
  recursive = F,
  full.names = T
)
list_of_analyses <- list_of_analyses[file.exists(paste0(list_of_analyses, "/Replication Assessment by Experiment.tsv"))]

for (analysis in list_of_analyses) {
  message("\n\n\n\n")
  message(analysis)
  summarize_rates(analysis)
}

# Gather data for the large replication rate table
gather_all_rep_rates(list_of_analyses, paste0("./output/", results_path))

# Plot specification curve for all analyses
plot_specification_curve(results_path, include_method = T, suffix = "")
plot_specification_curve(
  results_path,
  include_method = F,
  suffix = "-without-method"
)

# Build table of general numbers based on inclusion sets
make_general_numbers_table(inclusion_sets)

##### Run the second part of the analysis, the correlations between predictors and replication success, for the same specified folders in output #####
list_of_analyses <- list.dirs(
  path = paste0("./output/", results_path),
  recursive = F,
  full.names = T
)
list_of_analyses <- list_of_analyses[file.exists(paste0(list_of_analyses, "/Replication Assessment by Experiment.tsv"))]

# Prepares the predictor data, adjusting names and preprocessing predictor calculations
source("predictor-data-prep.R")

make_individual_graphs <- T # Will make graphs of predictors x outcomes only

# Run predictor x predictor analysis for ALL experiments (Done=Yes, UNIT=BRI) - only once
# This ensures predictor correlations include all 143 replications, not just those in analysis sets
# Outputs are copied to all analysis folders
run_predictor_analysis_all_exps(list_of_analyses)

for (analysis in list_of_analyses) {
  message("\n\n\n\n")
  message(analysis)

  pct1 <- run_predictor_analysis_exp_level(analysis, make_individual_graphs)
  pct2 <- run_predictor_analysis_rep_level(analysis, make_individual_graphs)
  run_predictor_analysis_both_level(analysis, make_individual_graphs)

  # Makes figures with table combined for better presentation
  plot_combined_cortable_predictor(pct1[[1]], pct2[[1]], analysis)
  plot_combined_cortable_rep_outcome(pct1[[2]], pct2[[2]], analysis)
  plot_combined_cortable(pct1[[3]], pct2[[3]], analysis)

  # Makes kappa figures
  pkexp <- plot_kappa_exp(
    paste0(
      analysis,
      "/predictors/by experiment/Table for Predictor Correlations.tsv"
    ),
    paste0(analysis, "/predictors/by experiment/")
  )
  pkrep <- plot_kappa_rep(
    paste0(
      analysis,
      "/predictors/by replication/Table for Predictor Correlations.tsv"
    ),
    paste0(analysis, "/predictors/by replication/")
  )
  combine_kappa_plots(pkexp, pkrep, paste0(analysis, "/predictors/"))
}

##### Make higher-level figures #####

# Runs post hoc power analysis and plots results
#source("post-hoc-power-sim.R")
plot_power_histograms(inclusion_sets)

# Run self-assessment analyses
source("self-assessment.R")

# Run survey analyses
source("survey-plot.R")

# Generate tables and numbers that are cited in text
source("tables.R")

# Organize manuscript figures and tables
source("prepare-figures.R")

# Dump session info for reproducibility
p_load(sessioninfo)
writeLines(capture.output(session_info()), "R session info.txt")
