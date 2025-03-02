# Survival Analysis of Chilean Annuitants

## Project Overview
This project analyses the mortality of Chilean Annuitants from 2014 to 2018 with the aim to identify whether the current assumptions used to price annuities are valid. The analysis includes descriptive statistics, survival analysis, life table graduation, and an ethical evaluation of gender-based pricing.

## Data Analysis Stages 

### Exploratory Data Analysis (EDA) 
**Objective:** 
- Within the dataset, several annuitants ages reached up to roughly 130 years old. This unrealistic/bad data has to be cleaned as the oldest person recorded Jeanne Calment was aged 122. Due to the scractiy of data and the absence of verification for these supercenternarions.
- Generate visualisations and summary statistics to uncover patterns and relationships within the dataset.

**Key Findings** 
- Potential predictors for future analysis include:
  - Age of the annuitant
  - Gender of the annuitant

### Survival Analysis 
**Approach:** Application of both semi-parametric and non-parametric techniques, including Cox regression and Kaplan-Meier(KM) estimations to identify the relationships predictors have on mortality rates.
**Findings:**
- Survival probability decreases with age
- Male annuitants and disabled annuitants experience higher mortality rates compared to their female and healthy annuitants respectively
**Challenges:** Since the covariates change with time, the proportional hazard assumption of the Cox regression model, the results provided by the model maybe an inaccurate representation, thus there is a need to rely on other models.

### Graduation of Unisex Life Table 
**Methodology:**
- Calculation of the crude mortality rates.
- Graduating the mortality rates via various parametric and non-parametric techniques, including:
 - Gompertz
 - Makeham
 - Cubic Splines
 - Smoothing Splines

**Recommended Model:**
