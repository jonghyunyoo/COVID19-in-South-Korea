################################################################################
# Replication Code for: Heterogeneous Impacts of Fear and Policy on Building Energy Use during COVID-19 in South Korea
# Journal: iScience
# Authors: Jonghyun Yoo, Daewon Kim, and Minjung Kim
# Date: 2025.11.01
#
# This script replicates the main regression results reported in the manuscript.
# It performs Two-Stage Least Squares (2SLS) instrumental variable regressions
# examining the impact of COVID-19 deaths on electricity consumption across
# different building types.
#
# Required data file: df_cleaned_replication.rds (cleaned dataset)
################################################################################

rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(bit64)
  library(data.table)
  library(dplyr)
  library(lfe)
  library(fixest)
})

# ------------------------------ LOAD DATA -------------------------------------
df <- readRDS("df_cleaned_replication.rds")

cat("Dataset loaded successfully\n")
cat(sprintf("Observations: %d | Variables: %d\n", nrow(df), ncol(df)))

# ----------------------------- SANITY CHECKS ----------------------------------
need <- c("type3","elec_size","covid_death","death","ym",
          "social_distance_lag","unemployment","temp","prcp",
          "identification","dong","month","sido")
miss <- setdiff(need, names(df))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse=", "))

# ---------------------------- TRANSFORMATIONS ---------------------------------
if (!is.integer(df$ym)) df$ym <- as.integer(df$ym)

# --------------------------- SPLIT BY BUILDING TYPE ---------------------------
bt_levels <- c("Commercial_Retail","Community","Civic_Business",
               "Industrial","Cultural_Religious","Residential")
building_types <- setNames(lapply(bt_levels, function(b) df[type3 == b]), bt_levels)

cat("\nSample sizes by building type:\n")
for (nm in names(building_types)) {
  cat(sprintf("  %-22s %10d\n", nm, nrow(building_types[[nm]])))
}

# --------------------------- FIRST STAGE (feols) ------------------------------
run_first_stage <- function(dt, label){
  # First-stage formula: log(covid_death) ~ log(death) + controls | month + sido
  fml_fs <- log(covid_death) ~ log(death) +
    social_distance_lag * ym +
    ym + I(ym^2) +
    unemployment +
    temp + I(temp^2) +
    prcp + I(prcp^2) | month + sido
  
  m_fs <- feols(fml_fs, data = dt, cluster = ~dong)
  # Predicted regressor for second stage (same name as in your original)
  dt$total_case_no_sido_hat <- as.numeric(predict(m_fs, newdata = dt))
  list(data = dt, fs = m_fs)
}

# -------------------------- SECOND STAGE (felm) -------------------------------
run_second_stage <- function(dt, label){
  fml_ss <- as.formula(
    "log(elec_size) ~ total_case_no_sido_hat * ym +
                      social_distance_lag * ym +
                      ym + I(ym^2) +
                      temp + I(temp^2) +
                      prcp + I(prcp^2) +
                      unemployment | month + identification | 0 | dong"
  )
  m_ss <- felm(fml_ss, data = dt)
  list(model = m_ss, summary = summary(m_ss))
}

# ------------------------- ESTIMATE & COLLECT RESULTS -------------------------
main_results <- list()
for (nm in names(building_types)) {
  cat("\n=== First stage:", nm, "===\n")
  fs_out <- run_first_stage(building_types[[nm]], nm)
  
  cat("First-stage F (approx, instrument log(death)):\n")
  sfs <- summary(fs_out$fs)
  if ("log(death)" %in% rownames(sfs$coeftable)) {
    t_iv <- sfs$coeftable["log(death)", "t value"]
    cat(sprintf("  t = %.2f, F ≈ %.2f\n", t_iv, t_iv^2))
  } else {
    cat("  (instrument term not found in coefficient table)\n")
  }
  
  cat("=== Second stage:", nm, "===\n")
  ss_out <- run_second_stage(fs_out$data, nm)
  main_results[[nm]] <- ss_out
  print(ss_out$summary)
}

# --------------------------- FORMAT OUTPUT TABLE -------------------------------
var_names <- c(
  "total_case_no_sido_hat",
  "total_case_no_sido_hat:ym",
  "social_distance_lag",
  "ym:social_distance_lag",
  "unemployment",
  "ym",
  "I(ym^2)",
  "temp",
  "I(temp^2)",
  "prcp",
  "I(prcp^2)"
)

var_labels <- c(
  "COVID-19 Deaths",
  "Deaths × Time",
  "Anti-contagion Policy",
  "Policy × Time",
  "Unemployment Rate",
  "Time",
  "Time²",
  "Temperature",
  "Temperature²",
  "Precipitation",
  "Precipitation²"
)

mk_cell <- function(est, se, p){
  stars <- ifelse(p < 0.01,"***", ifelse(p < 0.05,"**", ifelse(p < 0.10,"*","")))
  sprintf("%.4f (%.4f)%s", est, se, stars)
}

extract_col <- function(x){
  ct <- x$summary$coefficients
  rn <- rownames(ct)
  # order-insensitive matcher for interactions
  find_row <- function(term){
    if (!grepl(":", term, fixed = TRUE)) return(match(term, rn))
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    a <- paste0(parts[1], ":", parts[2])
    b <- paste0(parts[2], ":", parts[1])
    i <- match(a, rn); if (is.na(i)) i <- match(b, rn); i
  }
  vapply(var_names, function(v){
    i <- find_row(v)
    if (is.na(i)) return(NA_character_)
    mk_cell(ct[i,"Estimate"], ct[i,"Cluster s.e."], ct[i,"Pr(>|t|)"])
  }, character(1))
}

results_table <- data.frame(Variable = var_labels, stringsAsFactors = FALSE)
for (nm in names(main_results)) {
  results_table[[nm]] <- extract_col(main_results[[nm]])
}

results_table <- subset(results_table, !(Variable %in% c("Observations","R²")))

Ns  <- sapply(main_results, function(x) x$model$N)
R2  <- sapply(main_results, function(x) summary(x$model)$r.squared)
bt_cols <- c("Commercial_Retail","Community","Civic_Business",
             "Industrial","Cultural_Religious","Residential")

obs_row <- data.frame(Variable = "Observations", check.names = FALSE)
r2_row  <- data.frame(Variable = "R²",          check.names = FALSE)

for (cc in bt_cols) {
  obs_row[[cc]] <- format(Ns[[cc]], big.mark = ",")
  r2_row[[cc]]  <- sprintf("%.2f", R2[[cc]])
}

results_table <- rbind(results_table, obs_row, r2_row, make.row.names = FALSE)

write.csv(results_table, "REGRESSION_RESULTS_TABLE.csv", row.names = FALSE)
cat("\nSaved: REGRESSION_RESULTS_TABLE.csv\n")
