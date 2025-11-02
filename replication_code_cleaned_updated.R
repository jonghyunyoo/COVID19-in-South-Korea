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
# Required data file: df_cleaned_replication.csv (cleaned dataset)
################################################################################

# Clear workspace
rm(list = ls())
gc()

# Load required packages
library(bit64)
library(data.table)
library(dplyr)
library(lfe)        # For high-dimensional fixed effects models
library(fixest)     # Alternative for fixed effects estimation

################################################################################
# 1. DATA LOADING AND PREPARATION
################################################################################

# Load cleaned dataset
df <- fread("df_cleaned_replication.csv")

cat("Dataset loaded successfully\n")
cat(sprintf("Observations: %d\n", nrow(df)))
cat(sprintf("Variables: %d\n", ncol(df)))
cat(sprintf("Building types: %s\n", 
            paste(unique(df$type3), collapse = ", ")))

# Verify no missing values
missing_check <- sapply(df, function(x) sum(is.na(x)))
if (any(missing_check > 0)) {
  warning("Missing values detected in dataset!")
  print(missing_check[missing_check > 0])
}

################################################################################
# 2. CREATE BUILDING TYPE SUBSETS
################################################################################

# Building type categories (now in English)
building_types <- list(
  df_sample_etc    = df[type3 == "Commercial_Retail"],
  df_sample_store  = df[type3 == "Community"],
  df_sample_office = df[type3 == "Civic_Business"],
  df_sample_ind    = df[type3 == "Industrial"],
  df_sample_cult   = df[type3 == "Cultural_Religious"],
  df_sample_res    = df[type3 == "Residential"]
)

# Display sample sizes
cat("\nSample sizes by building type:\n")
for (name in names(building_types)) {
  cat(sprintf("  %s: %d observations\n", 
              name, nrow(building_types[[name]])))
}

# Assign to global environment for easier access
list2env(building_types, envir = .GlobalEnv)

################################################################################
# 3. FIRST STAGE: INSTRUMENTAL VARIABLE REGRESSION
################################################################################

# This section implements the first stage of 2SLS, using total deaths as an
# instrument for COVID-19 deaths, clustered at the district (dong) level.

cat("\n", rep("=", 70), "\n", sep = "")
cat("RUNNING FIRST STAGE REGRESSIONS\n")
cat(rep("=", 70), "\n\n")

# Initialize storage lists for first-stage results
first_stage_summaries <- list()
first_stage_r2 <- list()
weak_iv_tests <- list()

# Function to run first-stage IV regression for each building type
run_first_stage <- function(dt, name) {
  
  cat(sprintf("Processing %s...\n", name))
  
  # First-stage model: Regress log(COVID deaths) on log(total deaths)
  # controlling for social distancing policies, time trends, weather, and unemployment
  first_stage_model <- feols(
    log(covid_death) ~ log(death) + 
      social_distance_lag * ym + 
      I(ym^2) + 
      unemployment + 
      temp + I(temp^2) + 
      prcp + I(prcp^2) | 
      month + sido,                    # Fixed effects for month and region
    data = dt,
    cluster = "dong"                    # Cluster standard errors at district level
  )
  
  # Store model diagnostics
  first_stage_summaries[[name]] <<- summary(first_stage_model)
  first_stage_r2[[name]] <<- fitstat(first_stage_model, type = "ar2")
  weak_iv_tests[[name]] <<- fitstat(first_stage_model, type = "f")
  
  # Display F-statistic (weak instrument test)
  f_stat <- fitstat(first_stage_model, type = "f")$f
  cat(sprintf("  F-statistic: %.2f (weak IV if < 10)\n", f_stat))
  
  # Generate predicted values for second stage
  dt[, total_case_no_sido_hat := predict(first_stage_model, newdata = .SD)]
  
  return(dt)
}

# Apply first-stage regression to all building types
building_types <- Map(
  function(dt, name) run_first_stage(dt, name), 
  building_types, 
  names(building_types)
)

# Update global environment with IV-instrumented datasets
list2env(building_types, envir = .GlobalEnv)

cat("\n✓ First stage completed\n")

################################################################################
# 4. SECOND STAGE: MAIN REGRESSION ANALYSIS
################################################################################

cat("\n", rep("=", 70), "\n", sep = "")
cat("RUNNING SECOND STAGE REGRESSIONS (MAIN RESULTS)\n")
cat(rep("=", 70), "\n\n")

# Function to run second-stage regression for each building type
run_second_stage <- function(data, building_type_name) {
  
  cat(sprintf("Processing %s...\n", building_type_name))
  
  # Main 2SLS specification
  second_stage_model <- felm(
    log(elec_size) ~ 
      total_case_no_sido_hat * ym +        # Instrumented COVID deaths and interaction
      social_distance_lag * ym +            # Social distancing policy and interaction
      I(ym^2) +                             # Quadratic time trend
      temp + I(temp^2) +                    # Temperature controls (linear and quadratic)
      prcp + I(prcp^2) +                    # Precipitation controls (linear and quadratic)
      unemployment |                        # Unemployment rate
      month + identification |              # Fixed effects: month and building ID
      0 |                                   # No IV in this stage (already instrumented)
      dong,                                 # Cluster standard errors at district level
    data = data
  )
  
  return(list(
    model = second_stage_model,
    summary = summary(second_stage_model),
    building_type = building_type_name
  ))
}

# Apply second-stage regression to all building types with descriptive names
main_results <- list(
  Commercial_Retail      = run_second_stage(df_sample_etc, "Commercial_Retail"),
  Community              = run_second_stage(df_sample_store, "Community"),
  Civic_Business         = run_second_stage(df_sample_office, "Civic_Business"),
  Industrial             = run_second_stage(df_sample_ind, "Industrial"),
  Cultural_Religious     = run_second_stage(df_sample_cult, "Cultural_Religious"),
  Residential            = run_second_stage(df_sample_res, "Residential")
)

################################################################################
# 5. DISPLAY AND SAVE RESULTS
################################################################################

cat("\n", rep("=", 70), "\n", sep = "")
cat("REGRESSION RESULTS\n")
cat(rep("=", 70), "\n")

# Extract and print summary statistics for each building type
main_summaries <- lapply(main_results, function(x) {
  cat("\n\n===", x$building_type, "===\n")
  print(x$summary)
  return(x$summary)
})

################################################################################
# 6. CREATE RESULTS TABLE
################################################################################

cat("\n", rep("=", 70), "\n", sep = "")
cat("FORMATTED RESULTS TABLE\n")
cat(rep("=", 70), "\n\n")

# Function to extract coefficient and standard error
extract_coef <- function(model_result, var_name) {
  summ <- model_result$summary
  coef_table <- summ$coefficients
  
  if (var_name %in% rownames(coef_table)) {
    coef <- coef_table[var_name, "Estimate"]
    se <- coef_table[var_name, "Cluster s.e."]
    pval <- coef_table[var_name, "Pr(>|t|)"]
    
    # Add significance stars
    stars <- ifelse(pval < 0.01, "***",
                    ifelse(pval < 0.05, "**",
                           ifelse(pval < 0.10, "*", "")))
    
    return(sprintf("%.4f (%.4f)%s", coef, se, stars))
  } else {
    return(NA)
  }
}

# Variable names mapping
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

# Create results table
results_table <- data.frame(
  Variable = var_labels,
  stringsAsFactors = FALSE
)

building_type_names <- c("Commercial_Retail", "Community", "Civic_Business", 
                         "Industrial", "Cultural_Religious", "Residential")

for (bt in building_type_names) {
  results_table[[bt]] <- sapply(var_names, function(v) {
    extract_coef(main_results[[bt]], v)
  })
}

# Add model statistics
for (bt in building_type_names) {
  model <- main_results[[bt]]$model
  results_table[nrow(results_table) + 1, "Variable"] <- "Observations"
  results_table[nrow(results_table), bt] <- format(model$N, big.mark = ",")
  
  results_table[nrow(results_table) + 1, "Variable"] <- "R²"
  results_table[nrow(results_table), bt] <- sprintf("%.2f", summary(model)$r.squared)
}

# Print table
cat("\nFormatted Results Table:\n")
print(results_table, row.names = FALSE)

# Save table
results_file <- "REGRESSION_RESULTS_TABLE.csv"
write.csv(results_table, results_file, row.names = FALSE)
cat(sprintf("\nResults table saved to: %s\n", results_file))

################################################################################
# 7. DIAGNOSTICS AND VALIDATION
################################################################################

cat("\n", rep("=", 70), "\n", sep = "")
cat("MODEL DIAGNOSTICS\n")
cat(rep("=", 70), "\n\n")

cat("First Stage F-Statistics (Weak Instrument Test):\n")
cat("(F-statistic should be > 10 to rule out weak instruments)\n\n")

for (name in names(weak_iv_tests)) {
  f_stat <- weak_iv_tests[[name]]$f
  cat(sprintf("  %s: %.2f %s\n", 
              name, 
              f_stat,
              ifelse(f_stat > 10, "✓", "⚠")))
}

################################################################################
# END OF SCRIPT
################################################################################

cat("\n", rep("=", 70), "\n", sep = "")
cat("REPLICATION COMPLETED SUCCESSFULLY\n")
cat(rep("=", 70), "\n\n")

cat("Generated files:\n")
cat(sprintf("  1. %s\n", results_file))

cat("\nNote: Results should match the published manuscript table.\n")
cat("Any discrepancies should be investigated carefully.\n")