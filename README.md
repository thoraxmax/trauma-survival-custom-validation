# TRAUMA Survival Analysis — Custom Functions & Validation

This repository contains the custom R survival-analysis functions and standalone
simulation-based validation code accompanying the TRAUMA analyses.

The repository is intended to support reproducibility and transparency of the
custom statistical procedures used in the study. It does not contain the full
study analysis scripts or participant-level study data.

## Contents

### `custom_functions_test_only.R`

A minimal version of the custom survival-analysis functions required for the
validation tests. It contains only the functions needed for:

- delayed-entry Kaplan-Meier estimation and RMST;
- participant-level bootstrap estimation and percentile confidence intervals;
- midpoint-onset reconstruction; and
- post-trauma start-stop interval construction.

### `custom_functions_top3_demo_v7.R`

Reviewer-facing simulation code that generates a synthetic cohort and exercises
the three most custom analysis components:

1. delayed-entry KM/RMST and participant bootstrap;
2. midpoint-onset reconstruction; and
3. post-trauma start-stop Cox analysis.

The script includes numerical and structural validation checks and reports a
clear pass/fail status.

## Requirements

R 4.6.0 or later.

Required package:

```r
install.packages("survival")
