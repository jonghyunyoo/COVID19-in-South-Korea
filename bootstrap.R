################################################################################
# Cluster bootstrap SEs
# Reads the same df_cleaned_replication.rds and writes REGRESSION_RESULTS_TABLE_bootstrap.csv
################################################################################

rm(list = ls()); gc()
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(lfe); library(fixest) })

# ------------------------------ SETTINGS --------------------------------------
B_BOOT      <- 200
SEED        <- 2025
IN_DATA     <- "df_cleaned_replication.rds"       # same file used for manuscript
IN_MAIN_CSV <- "REGRESSION_RESULTS_TABLE.csv"     # produced by main two-step script
OUT_BOOT    <- "REGRESSION_RESULTS_TABLE_bootstrap.csv"

set.seed(SEED)

# ------------------------------ LOAD & PREP -----------------------------------
df <- readRDS(IN_DATA)

# Use EXACT transforms/names as in manuscript code:
# - log(elec_size) outcome
# - log(covid_death) and log(death)  (no +1)
# - Stage 1 FE: month + sido
# - Stage 2 FE: month + identification
# - Cluster (stage 2): dong
if (!is.integer(df$ym)) df$ym <- as.integer(df$ym)

bt_levels <- c("Commercial_Retail","Community","Civic_Business",
               "Industrial","Cultural_Religious","Residential")

RHS_STR <- paste(c("social_distance_lag*ym",
                   "ym","I(ym^2)",
                   "temp","I(temp^2)",
                   "prcp","I(prcp^2)",
                   "unemployment"), collapse = " + ")

# Coef names to extract (must match the second-stage model’s names)
coef_names <- c("total_case_no_sido_hat",
                "total_case_no_sido_hat:ym",
                "social_distance_lag",
                "ym:social_distance_lag",
                "ym","I(ym^2)",
                "unemployment",
                "temp","I(temp^2)",
                "prcp","I(prcp^2)")

all_dongs <- unique(df$dong)

# ------------------------------ HELPERS ---------------------------------------
# First stage (building panel), FE: month + sido; cluster for robustness not needed for prediction
first_stage_fit <- function(dtb){
  feols(
    log(covid_death) ~ log(death) +
      social_distance_lag*ym + ym + I(ym^2) +
      unemployment + temp + I(temp^2) + prcp + I(prcp^2) |
      month + sido,
    data = dtb
  )
}

# Second stage (building panel), FE: month + identification; cluster by dong
second_stage_fit <- function(dtb){
  felm(as.formula(paste0(
    "log(elec_size) ~ total_case_no_sido_hat*ym + ", RHS_STR,
    " | month + identification | 0 | dong"
  )), data = dtb)
}

# Find row by coefficient name (order-insensitive for interactions)
find_row <- function(term, rn){
  if (!grepl(":", term, fixed = TRUE)) return(match(term, rn))
  parts <- strsplit(term, ":", fixed = TRUE)[[1]]
  a <- paste0(parts[1],":",parts[2]); b <- paste0(parts[2],":",parts[1])
  i <- match(a, rn); if (is.na(i)) i <- match(b, rn); i
}

# One bootstrap replicate for a given building type
fit_two_step_bt <- function(smp_dongs, bt){
  dt <- df[type3 == bt & dong %in% smp_dongs]
  
  # Stage 1 on the building panel (matches manuscript weighting)
  fs <- try(first_stage_fit(dt), silent = TRUE)
  if (inherits(fs, "try-error")) return(setNames(rep(NA_real_, length(coef_names)), coef_names))
  
  # Predicted regional COVID (generated regressor, exact name)
  dt$total_case_no_sido_hat <- as.numeric(predict(fs, newdata = dt))
  
  # Stage 2 with clustering by dong
  ss <- try(second_stage_fit(dt), silent = TRUE)
  if (inherits(ss, "try-error")) return(setNames(rep(NA_real_, length(coef_names)), coef_names))
  
  s  <- summary(ss)
  co <- s$coefficients
  rn <- rownames(co)
  
  out <- vapply(coef_names, function(v){
    i <- find_row(v, rn); if (is.na(i)) NA_real_ else co[i, "Estimate"]
  }, numeric(1))
  setNames(out, coef_names)
}

# ------------------------------ BOOTSTRAP -------------------------------------
boot_se <- setNames(vector("list", length(bt_levels)), bt_levels)

for (bt in bt_levels) {
  cat(sprintf("Bootstrap %s (B=%d)\n", bt, B_BOOT))
  mat <- matrix(NA_real_, nrow = B_BOOT, ncol = length(coef_names))
  colnames(mat) <- coef_names
  
  for (b in seq_len(B_BOOT)) {
    smp <- sample(all_dongs, length(all_dongs), replace = TRUE)  # cluster resample
    mat[b, ] <- fit_two_step_bt(smp, bt)
    if (b %% 20 == 0) cat(".")
  }
  cat(" done\n")
  
  # bootstrap SEs: SD across replicates (ignore NA)
  boot_se[[bt]] <- apply(mat, 2, function(x) sd(x, na.rm = TRUE))
}

# ------------------------ INJECT SEs INTO MAIN CSV ----------------------------
main_tab <- fread(IN_MAIN_CSV)
K <- length(coef_names)

inject_boot <- function(col_vec, se_named) {
  out <- col_vec
  for (k in seq_len(K)) {
    nm <- coef_names[k]; se_k <- se_named[[nm]]
    if (is.finite(se_k)) out[k] <- sub("\\(([-0-9\\.eE]+)\\)", sprintf("(%.4f)", se_k), out[k])
  }
  out
}

res_tab_boot <- copy(main_tab)
for (bt in bt_levels) {
  col <- res_tab_boot[[bt]]
  res_tab_boot[[bt]][seq_len(K)] <- inject_boot(col, boot_se[[bt]])
}

fwrite(res_tab_boot, OUT_BOOT)
cat(sprintf("Saved bootstrap-annotated table: %s\n", OUT_BOOT))
