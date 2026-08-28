if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(tidyverse, readxl, EnvStats)

# Load data extracted from original articles

D = read_excel("./other-data/Extracted Data from Original Experiments.xlsx", na = c("NA", "NI", ""))

# Calculates columns from extracted data
D = D |>
  mutate(
    `Control SD` = ifelse(`Error Type` == "SD", `Control Error`, `Control Error` * `Control Sample Size for DFs` ^ 0.5),
    `Treated SD` = ifelse(`Error Type` == "SD", `Treated Error`, `Treated Error` * `Treated Sample Size for DFs` ^ 0.5),
    Delta = `Treated Mean` - `Control Mean`,
    `Pooled SD` = (((`Control SD` ^ 2) + (`Treated SD` ^ 2)) / 2) ^ 0.5,
    
    `Treated Mean (Relative to Control)` = `Treated Mean` / `Control Mean`,
    `Delta (Relative to Control)` =  Delta / `Control Mean`,
    `Control SD (Relative to Control)` =  `Control SD` / `Control Mean`, 
    `Treated SD (Relative to Control)` = `Treated SD` / `Control Mean`, 
    `Pooled SD (Relative to Control)` = `Pooled SD` / `Control Mean`
  )

# Calculate sample size for 95% power using log ratio of means
calc_lnorm_n <- function(t_mean, c_mean, t_sd, c_sd) {
  if (is.na(t_mean) || is.na(c_mean) || is.na(t_sd) || is.na(c_sd)) return(NA_real_)
  if (c_mean <= 0 || t_mean <= 0) return(NA_real_)
  ratio <- t_mean / c_mean
  cv1 <- abs(c_sd / c_mean)
  cv2 <- abs(t_sd / t_mean)
  pooled_cv <- sqrt((cv1^2 + cv2^2) / 2)
  if (pooled_cv <= 0 || ratio <= 0) return(NA_real_)
  res <- tryCatch({
    tTestLnormAltN(ratio.of.means = ratio, cv = pooled_cv, alpha = 0.05, power = 0.95, sample.type = "two.sample")
  }, error = function(e) NA_real_)
  as.numeric(res)
}

D$`Estimated Sample Size for 95% power (Log Ratio)` <- pmap_dbl(
  list(D$`Treated Mean`, D$`Control Mean`, D$`Treated SD`, D$`Control SD`),
  calc_lnorm_n
)

# There are multiple cases where we assumed the tests used: one-sample and two-sample (with unequal variances - the most common). Then, next step is to calculate the ts and dfs for each experiment, to pick the correct t distribution to estimate the corresponding p-value.

# One sample:
D_ONE_SAMPLE = D |> filter(
  `Assumed Test for Sample Size Calculation` == "One-sample t test power calculation"
) |>
  mutate(
    degrees_of_freedom = `Treated Sample Size for DFs` - 1
  ) |>
  select(EXP, degrees_of_freedom)

# Two-sample, unequal variances:
D_TWO_SAMPLE_UNEQUAL = D |> filter(
  `Assumed Test for Sample Size Calculation` == "Two-sample t test power calculation with unequal variances"
) |>
  mutate(
    sd1_n1 = (`Control SD (Relative to Control)` ^ 2) / `Control Sample Size for DFs`,
    sd2_n2 = (`Treated SD (Relative to Control)` ^ 2) / `Treated Sample Size for DFs`,
    degrees_of_freedom_numerator = ((sd1_n1 + sd2_n2) ^ 2),
    degrees_of_freedom_denominator = (sd1_n1 ^ 2) / (`Control Sample Size for DFs` - 1) + (sd2_n2 ^ 2) / (`Treated Sample Size for DFs` - 1),
    degrees_of_freedom = degrees_of_freedom_numerator / degrees_of_freedom_denominator
  ) |>
  select(EXP, degrees_of_freedom)


D = full_join(D, rbind(D_ONE_SAMPLE, D_TWO_SAMPLE_UNEQUAL), by = "EXP")


write_tsv(file = "./other-data/Original Experiments Statistical Summaries.tsv", D, quote = "all")

rm(list = ls(all = TRUE))
