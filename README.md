# Survival Analysis of Chilean Annuitants  

## Project Overview  
This project analyses the mortality of Chilean annuitants from 2014 to 2018 to assess whether the current assumptions used for annuity pricing are valid. The analysis includes **descriptive statistics, survival analysis, life table graduation, and an ethical evaluation** of gender-based pricing.  

## Data Analysis Stages  

### Exploratory Data Analysis (EDA)  
**Objective:**  
- Identify and clean unrealistic data points, such as annuitants with reported ages exceeding **130 years**. Given that the oldest recorded person, **Jeanne Calment**, lived to **122 years**, these extreme values likely indicate data entry errors.  
- Generate visualisations and summary statistics to uncover patterns and relationships within the dataset.  

**Key Findings:**  
- The primary predictors for future analysis include:  
  - **Age of the annuitant**  
  - **Gender of the annuitant**  

### Survival Analysis  
**Approach:**  
- Apply **semi-parametric and non-parametric survival analysis techniques**, including:  
  - **Cox Regression**  
  - **Kaplan-Meier (KM) Estimations**  

**Findings:**  
- **Survival probability decreases with age.**  
- **Male annuitants** and **disabled annuitants** experience **higher mortality rates** compared to their **female** and **healthy** counterparts.  

**Challenges:**  
- Since **covariates change over time**, the **proportional hazard assumption** of the Cox regression model does not hold, potentially leading to inaccurate results. As a result, alternative models may be needed.  

### Graduation of Unisex Life Table  
**Methodology:**  
- Calculate **crude mortality rates**.  
- Graduate mortality rates using various **parametric and non-parametric techniques**, including:  
  - **Gompertz Model**  
  - **Makeham Model**  
  - **Cubic Splines**  
  - **Smoothing Splines**  

**Recommended Model:**  
- After applying a series of statistical tests, the **smoothing splines method** was found to best align with the data.  

### Key Findings & Implications  
- The **graduated unisex mortality rates underestimate male mortality rates** and **overestimate female mortality rates**.  
- This suggests that **unisex mortality tables may not accurately represent real mortality experiences**, potentially leading to mispricing in annuity products.  

## Ethical Considerations  

**Focus:**  
- Evaluate the ethical implications of using **gender** as a rating factor for annuity pricing.  

**Pros of Gender-Based Pricing:**  
- Allows insurers to **accurately price annuity products**, ensuring long-term financial sustainability.  
- Helps account for **gender-based health differences**, leading to more precise mortality estimations.  

**Cons of Gender-Based Pricing:**  
- **Utilitarian Perspective:** In Chile, **women earn 22% less than men on average**. This wage disparity means fewer women may be able to afford annuities, limiting their financial security in retirement.  
- **Deontological Perspective:** Penalising women for a factor **beyond their control** (i.e., longer life expectancy) is **morally questionable**.  

**Recommendation:**  
- **Gender should be used as a rating factor**, but insurers should provide transparency and reassurance to policyholders.  
- Implementation of **legislation similar to "The Anti-Discrimination Act 1977" (NSW, Australia)** may help balance fairness and financial viability.  

## Technical Details  
**Tools Used:**  
- **R**  

---

