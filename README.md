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

**Recommended Model:** By applying a series of statistical tests, it is deduced that the smoothing splines align most closely the data. 

### Key Findings & Implications
The graduated unisex mortality rates are consistently underestimating the mortality rates for males and overestimating the mortality rates for females. This observation suggests that graduating unisex mortality rates may not accurately capture the mortality experience.

## Ethical Considerations 
**Focus:** Analysis of the ethical implications of using gender as a rating factor for annuity pricing. 
**Pros:**
- Insurers are able to properly price their insurance products, which ensures sustainability as by disregarding gender insurers would suffer huge losses if unisex life tables were used to price their annuities products as the life tables would consistently overestimate the mortality rates for females which means the insurer would have to pay more benefits.
- Estimation of mortality rates may be more inaccurate due to the inability to consider gender-based diseases.
**Cons:**
- Under a utilitarian lens, in Chile females on average earn 22% less than males. The pay disparity is unfair as fewer females may have the financial means to afford annuity products thereby hindering their ability to receive adequate pension benefits.
- Under a deontological perspective, it is morally wrong to penalise females for a factor out their control.
**Recommendation:** Gender should be used as a rating factor but to provide reassurance to policyholders, it is best to implement legislation similar to "The Anti-Discrimination Act 1977" in New South
Wales.

## Technical details 
**Tools Used**
- R


