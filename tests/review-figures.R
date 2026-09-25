# Run from the repository root after original-precision-calculation.R.
# Outputs go to an ignored directory unless an output root is supplied.
suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(readxl)
})
source("plots.R")
args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args)) args[[1]] else "output/review-figure-checks"
result <- plot_original_sd_mean_regression(
  "other-data/Original Experiments Statistical Summaries.tsv",
  file.path(output_root, "original-error-bar-regression")
)
diagnostics <- result$diagnostics
stopifnot(
  nrow(diagnostics) == 112L,
  n_distinct(paste(diagnostics$EXP, diagnostics$arm)) == 112L,
  nrow(result$excluded_points) == 8L,
  identical(as.integer(table(diagnostics$panel)), c(40L, 37L, 35L)),
  all(result$excluded_points$arm == "Control"),
  all(result$excluded_points$analysis_sd == 0),
  setequal(result$excluded_points$EXP,
    c("MTT86", "MTT87", "MTT96", "PCR16", "PCR69", "PCR147", "PCR175", "PCR184")),
  all(result$excluded_points$EXP %in% diagnostics$EXP[diagnostics$arm == "Treated"]),
  identical(levels(diagnostics$study_arm), c("Control", "Treated", "Treated (Normalized)")),
  identical(levels(diagnostics$sd_provenance), c("SD", "SEM (exact)", "SEM (estimated)")),
  is.null(result$plot$labels$caption)
)

# Changing both logarithms' base must preserve slopes, R-squared,
# original-scale fitted values and prediction limits, and residual flags.
for (panel in names(result$models)) {
  data <- filter(diagnostics, as.character(.data$panel) == .env$panel)
  old <- lm(log10(analysis_sd) ~ log10(analysis_mean), data = data)
  new <- result$models[[panel]]
  interval <- 10 ^ predict(old, newdata = data, interval = "prediction")
  equal <- function(x, y) isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-9))
  stopifnot(
    equal(coef(old)[[2]], coef(new)[[2]]),
    equal(coef(old)[[1]] * log(10), coef(new)[[1]]),
    equal(summary(old)$r.squared, summary(new)$r.squared),
    equal(interval[, "fit"], data$fitted_sd),
    equal(interval[, "lwr"], data$prediction_lower_sd),
    equal(interval[, "upr"], data$prediction_upper_sd),
    equal(rstudent(old), data$externally_studentized_residual),
    identical(unname(abs(rstudent(old)) >= 2), data$meets_abs_rstudent_2)
  )
}

# S4 labels show probabilities while the color mapping remains -log10(p).
labels <- format_log10_p_value_labels(c(0, 2, 4, 6))
stopifnot(identical(labels[[1]], 1))
for (i in 2:4) {
  stopifnot(isTRUE(all.equal(eval(labels[[i]]), 10 ^ (-2 * (i - 1)))))
}

# S3 must account for every article once and preserve the standalone maps.
articles <- read_excel("other-data/intermediate_steps.xlsx", sheet = "Original Articles")
spec <- jsonlite::fromJSON("maps/figure_s3.vega", simplifyVector = FALSE)
stopifnot(nrow(articles) == 60L, n_distinct(articles$EXP) == 60L,
  n_distinct(articles$DOI) == 60L,
  all(!is.na(articles$selected_correspondence_affiliation)),
  all(nzchar(articles$selected_correspondence_affiliation)),
  all(!is.na(articles$state)))
counts <- count(articles, state)
for (row in spec$datasets$panel_d_counts) {
  expected <- counts$n[counts$state == row$id]
  stopifnot(row$n == if (length(expected)) expected else 0)
}
stopifnot(all(articles$state %in% vapply(spec$datasets$panel_d_counts, `[[`, "", "id")))
standalone <- c("registered_labs", "included_labs", "labs_contributing_experiments")
for (i in seq_along(standalone)) {
  lines <- readLines(paste0("maps/", standalone[[i]], ".vega"), warn = FALSE)
  original <- jsonlite::fromJSON(paste(lines[!grepl("^//", lines)], collapse = "\n"), simplifyVector = FALSE)
  values <- original$transform[[1]]$from$data$values
  stopifnot(identical(values, spec$datasets[[paste0("panel_", letters[[i]], "_counts")]]),
    sum(vapply(values, `[[`, 0, "n")) == c(97, 75, 56)[[i]])
}

mapping <- read_excel("other-data/Manuscript Figure Correspondence.xlsx")
s10 <- filter(mapping, type == "figure", manuscript_name == "S10")
stopifnot(nrow(s10) == 1L,
  identical(s10$code_generated_filename,
    "/original-error-bar-regression/original-study-log-sd-vs-log-mean.png"))
manuscript_dir <- file.path(output_root, "_manuscript figures and tables/manuscript")
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.copy(result$figure_path, file.path(manuscript_dir, "Figure_S10.png"), overwrite = TRUE))
cat("PASS: S3 article counts and maps; S4 labels; S10 counts, exclusions, natural-log equivalence, legends and registration\n")
