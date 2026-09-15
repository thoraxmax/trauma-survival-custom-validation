# Trauma Survival Analysis

This repository contains the R code used for the survival analyses in the associated manuscript, together with a synthetic validation framework for the custom survival-analysis functions.

Participant-level study data are not included.

## Analysis code

The complete analysis code is in the `analysis/` folder:

```text
analysis/
├── custom_functions.R
├── 01_RQ1_incidence.R
├── 02_RQ2_time_to_onset.R
├── 03_RQ3_risk_factors.R
└── 04_sensitivity_analyses.R
```

The scripts cover:

* incidence and cumulative-risk analyses;
* time-to-onset analyses;
* delayed-entry Cox models;
* proportional-hazards diagnostics;
* trauma-type and trauma-count models;
* demographic and age-at-trauma models; and
* prespecified sensitivity analyses.

Run the scripts in this order:

```text
01_RQ1_incidence.R
02_RQ2_time_to_onset.R
03_RQ3_risk_factors.R
04_sensitivity_analyses.R
```

RQ2 should be run before RQ3 and the sensitivity analyses because it generates the sex-handling decision used by downstream Cox models.

## Input data

The analysis scripts expect a prepared dataset named:

```text
cox_survival_dataset_full.csv
```

Alternatively, the input path can be specified using:

```r
Sys.setenv(SURVIVAL_DATA_PATH = "/path/to/cox_survival_dataset_full.csv")
```

The participant-level dataset is not publicly distributed because of data-governance and confidentiality requirements.

## Requirements

The main analyses require R and the `survival` package.

Optional packages used for formatted Word output are:

```text
officer
flextable
```

## Validation

`testing_run.R` and `custom_functions_minimal.R` provide a synthetic-data validation framework for key custom survival-analysis procedures without using study data.

## Reproducibility

The repository provides the complete statistical analysis code required to reproduce the analyses when access to the prepared study dataset is available.
