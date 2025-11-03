# Replication Package

## Heterogeneous Impacts of Fear and Policy on Building Energy Use during COVID-19 in South Korea
**Authors:** Jonghyun Yoo, Daewon Kim, and Minjung Kim  
**Journal:** iScience  
**Year:** 2025

## Overview

This repository contains the replication code and cleaned data for the main regression analysis reported in our iScience paper. The analysis examines the impact of COVID-19 deaths on electricity consumption across different building types in [Location], using a Two-Stage Least Squares (2SLS) instrumental variable approach.

## Requirements

### Software
- R version 4.0 or higher

### Required R Packages
```r
install.packages(c(
  "bit64",
  "data.table",
  "dplyr",
  "lfe",
  "fixest"
))
```

## Repository Structure

```
.
├── README.md
├── replication.R                 # Main replication (writes table)
├── bootstrap.R                  # Optional cluster bootstrap of SEs
├── df_cleaned_replication.rds            # Cleaned dataset (input)
├── REGRESSION_RESULTS_TABLE.csv           # Main results (output)
└── REGRESSION_RESULTS_TABLE_bootstrap.csv # Bootstrap SEs (optional output)

```

## Quick Start
```r
# 1) Main results
source("replication.R")
# -> writes REGRESSION_RESULTS_TABLE.csv

# 2) (Optional) Cluster bootstrap (by dong)
source("bootstrap.R")
# -> writes REGRESSION_RESULTS_TABLE_bootstrap.csv
```

This will replicate all results from the manuscript.

### Building Type Categories

| Building Type | N (observations) |
| ------------------- | -------- |
| Commercial_Retail | ~716,000 |
| Community | ~14,328,000 |
| Civic_Business | ~1,509,000 |
| Industrial | ~3,578,000 |
| Cultural_Religious | ~717,000 |
| Residential | ~2,167,000 |

### Data Dictionary

| Variable | Description | Type | Use |
|----------|-------------|------|-----|
| `elec_size` | Electricity consumption (kWh) | numeric | Dependent variable |
| `type3` | Building type category | character | Subsetting |
| `covid_death` | COVID-19 deaths | numeric | First stage (dep. var.) |
| `death` | Total deaths | numeric | Instrument |
| `social_distance_lag` | Lagged social distancing policy | numeric | Treatment variable |
| `ym` | Time (months from baseline) | numeric | Time trend |
| `unemployment` | Unemployment rate (%) | numeric | Control |
| `temp` | Temperature (°C) | numeric | Control |
| `prcp` | Precipitation (mm) | numeric | Control |
| `month` | Month (1-12) | integer | Fixed effects |
| `identification` | Building ID | character | Fixed effects |
| `sido` | Region identifier | integer | Fixed effects |
| `dong` | District identifier | integer | Clustering |
