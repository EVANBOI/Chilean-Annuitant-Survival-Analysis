# Packages Used
library(dplyr)
library(ggplot2)
library(survival)
library(KMsurv)
library(readxl)
library(survminer)
library(ldatools)
library(lifecontingencies)
library(demography)
library(splines)

# Loading the Data 
annuitant_data <- read.csv("C:/Users/Evan/Downloads/ChileanMortality.csv")

# Obtaining Elementary Data
str(annuitant_data)
summary(annuitant_data)


# Exploratory Data Analysis

# Annuitants Distributed by Gender
gender_prop <- annuitant_data %>% group_by(SEX) %>% 
  summarise(proportions = n()/nrow(annuitant_data))

gender_prop

# Annuitants Distributed by Health

disability_prop <- annuitant_data %>% group_by(HEALTH) %>%
  summarise(proportions = n()/nrow(annuitant_data))

disability_prop


# Annuitants Distributed by Death

death_prop <- annuitant_data %>% group_by(DEATH) %>%
  summarise(proportions = n()/nrow(annuitant_data))

death_prop


# Data Cleaning 
annuitant_data$DEATH <- as.numeric(annuitant_data$DEATH)
annuitant_data$SEX <- as.factor(annuitant_data$SEX)
annuitant_data$HEALTH <- as.factor(annuitant_data$HEALTH)
annuitant_data$PERSON_TYPE <- as.factor(annuitant_data$PERSON_TYPE)


# Add Age of Death, Start Age and time of exposure
annuitant_data$BIRTHDATE <- as.Date(annuitant_data$BIRTHDATE)
annuitant_data$DATE_START <- as.Date(annuitant_data$DATE_START)
annuitant_data$DATE_END <- as.Date(annuitant_data$DATE_END)

annuitant_data <- annuitant_data %>% 
  mutate(START_AGE = difftime(annuitant_data$DATE_START, 
                              annuitant_data$BIRTHDATE, units = "days")/365.25,
         END_AGE = difftime(annuitant_data$DATE_END, 
                             annuitant_data$BIRTHDATE, units = "days")/365.25,
         EXPOSURE_TIME = annuitant_data$DATE_END - annuitant_data$DATE_START)

annuitant_data$START_AGE <- as.numeric(annuitant_data$START_AGE)
annuitant_data$END_AGE <- as.numeric(annuitant_data$END_AGE)
annuitant_data$EXPOSURE_TIME <- as.numeric(annuitant_data$EXPOSURE_TIME)


# Don't need to round the ages up; the floor function can handle it 
# annuitant_data$START_AGE <- floor(annuitant_data$START_AGE)
# annuitant_data$END_AGE <- floor(annuitant_data$END_AGE)

# Ensure individuals start at age 60
# To start at age 60

aged_less_than <- annuitant_data %>% filter(START_AGE < 60)

aged_less_than$DATE_START <- aged_less_than$BIRTHDATE + 60*365.25
aged_less_than$START_AGE <- aged_less_than$DATE_START - aged_less_than$BIRTHDATE
aged_less_than$START_AGE <- as.numeric(aged_less_than$START_AGE/365.25)

annuitant_data <- annuitant_data %>% filter(START_AGE >= 60)

annuitant_data <- rbind(annuitant_data, aged_less_than)

# Omit if the end age is 60 

annuitant_data <- annuitant_data %>% filter(START_AGE != END_AGE)

# Mostly unverified
annuitant_data <- annuitant_data %>% filter(END_AGE <= 110)



# Correlation Between Variables

# Proportions of Healthy Annuitants by Gender
health_proportions <- annuitant_data %>% group_by(SEX) %>% 
  summarise(total_count = n(),
            health_count = sum(HEALTH == "Healthy"),
            health_proportion = health_count/total_count)

# Proportion of Beneficiaries by Gender

# Calculate the count of male and female beneficiaries separately
beneficiary_counts <- annuitant_data %>% 
  filter(PERSON_TYPE == "Beneficiary") %>%
  group_by(SEX) %>%
  summarise(total_count = n())

# Calculate the total count of beneficiaries
total_beneficiaries <- sum(beneficiary_counts$total_count)

# Calculate the proportion of male and female beneficiaries
beneficiary_proportions <- beneficiary_counts %>%
  mutate(beneficiary_proportion = total_count / total_beneficiaries)


# Proportions of Deaths by Gender

death_by_gender <- annuitant_data %>% group_by(SEX, DEATH)

death_proportions <- annuitant_data %>% 
  group_by(SEX) %>%
  summarise(total_count = n(),
            death_count = sum(DEATH == 1),
            death_proportion = death_count / total_count)

ggplot(data = death_by_gender, mapping = aes(x = SEX)) + 
  geom_bar(aes(fill = factor(DEATH)), position = "fill") +
  geom_text(data = death_proportions, 
            aes(x = SEX, y = death_proportion, 
                label = scales::percent(death_proportion)),
            color = "black", size = 4, vjust = 1.2) +
  xlab("Sex") +
  ylab("Proportion of Deaths") +
  ggtitle("Proportion of Deaths by Gender") +
  scale_fill_manual(values = c("0" = "skyblue", "1" = "red"), 
                    labels = c("0" = "Censored", "1" = "Death")) +
  labs(fill = "Outcome") +
  theme_minimal()




# Distribution of the Age of Death

# Histogram (Do it by gender as well)
death_data <- annuitant_data %>% filter(DEATH == 1)

hist(death_data$END_AGE, 
     col = 'lightblue', 
     ylab = "Number of Deaths", 
     xlab = "Age", 
     main = "Distribution of the Age of Death")

ggplot(death_data, aes(x = END_AGE, fill = SEX)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  labs(
    x = "Age",
    y = "Number of Deaths",
    title = "Distribution of the Age of Death by Gender"
  ) +
  theme_minimal()

ggplot(death_data, aes(x = END_AGE, fill = PERSON_TYPE)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  labs(
    x = "Age",
    y = "Number of Deaths",
    title = "Distribution of the Age of Death by Person Type"
  ) +
  theme_minimal()


# Task 2 

# Creating Survival Objects
cens_annuitant <- Surv(annuitant_data$START_AGE - 60, 
                       annuitant_data$END_AGE - 60, annuitant_data$DEATH)

summary(cens_annuitant)

cens_male <- Surv(annuitant_data$START_AGE[annuitant_data$SEX == "M"] - 60,
                  annuitant_data$END_AGE[annuitant_data$SEX == "M"] - 60,
                  annuitant_data$DEATH[annuitant_data$SEX == "M"])

cens_female <- Surv(annuitant_data$START_AGE[annuitant_data$SEX == "F"] - 60,
                    annuitant_data$END_AGE[annuitant_data$SEX == "F"] - 60,
                    annuitant_data$DEATH[annuitant_data$SEX == "F"])

cens_disabled <- Surv(annuitant_data$START_AGE[annuitant_data$HEALTH == "Disabled"] - 60,
                      annuitant_data$END_AGE[annuitant_data$HEALTH == "Disabled"] - 60,
                      annuitant_data$DEATH[annuitant_data$HEALTH == "Disabled"])

cens_healthy <- Surv(annuitant_data$START_AGE[annuitant_data$HEALTH == "Healthy"] - 60,
                     annuitant_data$END_AGE[annuitant_data$HEALTH == "Healthy"] - 60,
                     annuitant_data$DEATH[annuitant_data$HEALTH == "Healthy"])

cens_main <- Surv(annuitant_data$START_AGE[annuitant_data$PERSON_TYPE == "Main Annuitant"] - 60,
                         annuitant_data$END_AGE[annuitant_data$PERSON_TYPE == "Main Annuitant"] - 60,
                         annuitant_data$DEATH[annuitant_data$PERSON_TYPE == "Main Annuitant"])

cens_beneficiary <- Surv(annuitant_data$START_AGE[annuitant_data$PERSON_TYPE == "Beneficiary"] - 60,
                  annuitant_data$END_AGE[annuitant_data$PERSON_TYPE == "Beneficiary"] - 60,
                  annuitant_data$DEATH[annuitant_data$PERSON_TYPE == "Beneficiary"])

km_survival_prob <- survfit(cens_annuitant ~ 1, conf.int = 0.95, 
                             conf.type = "log")
# Non-Parametric Approaches

# K-M Estimation
km_male_sp <- survfit(cens_male ~ 1, conf.int = 0.95, conf.type = "log")
km_female_sp <- survfit(cens_female ~ 1, conf.int = 0.95, conf.type = "log")
km_disabled_sp <- survfit(cens_disabled ~ 1, conf.int = 0.95, conf.type = "log")
km_healthy_sp <- survfit(cens_healthy ~ 1, conf.int = 0.95, conf.type = "log")
km_main_sp <- survfit(cens_main ~ 1, conf.int = 0.95, conf.type = "log")
km_beneficiary_sp <- survfit(cens_beneficiary ~ 1, conf.int = 0.95, conf.type = "log")


# Total KM
plot(km_survival_prob, main = "KM Estimate with 95% Confidence Intervals", 
     xlab = "Time", ylab = "S(X)")

# Survival probability of males 
plot(km_female_sp, main = "KM Estimate with 95% Confidence Intervals", 
     xlab = "Time", ylab = "S(X)", col = "skyblue")


# Survival probability of females
lines(km_male_sp, col = "red")
legend("bottomleft", legend = c("Female", "Male"), col = c("skyblue", "red"), 
       lty = 1, cex = 0.6, x.intersp = 0.5, y.intersp = 0.5)

# Log-rank Test 
sex_cox <- coxph(cens_annuitant ~ SEX, data = annuitant_data)
summary(sex_cox)

# Survival Proability of disabled 
plot(km_disabled_sp, main = "KM Estimate with 95% Confidence Intervals", 
     xlab = "Time", ylab = "S(X)", col = "skyblue")

#Survival probability of healthy
lines(km_healthy_sp, col = "red")
legend("bottomleft", legend = c("Disabled", "Healthy"), col = c("skyblue", "red"), 
       lty = 1, cex = 0.6, x.intersp = 0.5, y.intersp = 0.5)

# Log-rank Test 
health_cox <- coxph(cens_annuitant ~ HEALTH, data = annuitant_data)
summary(health_cox)

#Survival probability of main annuitant 
plot(km_main_sp, main = "KM Estimate with 95% Confidence Intervals", 
     xlab = "Time", ylab = "S(X)", col = "skyblue")

#Survival probability of beneficiaries
lines(km_beneficiary_sp, col = "red")
legend("bottomleft", legend = c("Main Annuitant", "Beneficiary"), col = c("skyblue", "red"), 
       lty = 1, cex = 0.6, x.intersp = 0.5, y.intersp = 0.5)

# Log-rank Test 
person_type_cox <- coxph(cens_annuitant ~ PERSON_TYPE, data = annuitant_data)
summary(person_type_cox)

# Need for a parametric approach
# The cumulative hazard against t
plot(x = km_survival_prob$time, y = -log(km_survival_prob$surv), type = "l" ,
     col = "skyblue", main = "Cumulative Hazard over Time", xlab = "Time", 
     ylab = "Cumulative Hazard")


# Semi-Parametric Approach
# If cox regression estimate has a 1 within CI
# Cox Regression 
cox_reg <- coxph(cens_annuitant ~ HEALTH + SEX + PERSON_TYPE, 
                 data = annuitant_data, method = "breslow")

# Regression Results
summary(cox_reg)

# Baseline Hazard
baseline_hazard <- basehaz(cox_reg, centered = FALSE) 
plot(baseline_hazard$time, baseline_hazard$hazard, xlab = "Time", 
     ylab = "Baseline Hazard Rate", main = "Baseline Hazard Rate against Time", 
     type = "s", col = "skyblue")

# Anova Table
anova(cox_reg)

# Testing Proportionality Assumption 
test_ph <- cox.zph(cox_reg)
test_ph

ggcoxzph(test_ph)
ggcoxzph(test_ph[1]) + ylab("Health")
ggcoxzph(test_ph[2]) + ylab("Sex")
ggcoxzph(test_ph[3]) + ylab("Person Type")

# Cox-Snell Residuals
gg_coxsnell(cox_reg, type = c("cumu_hazard", "cdf")) +
  ylim(0,2.2) + xlim(0,2.2) +
  geom_abline(intercept = 0, slope = 1, col = 2) + 
  theme_minimal()

# Life Tables

chilean_lifetables <- read_excel("C:/Users/Evan/Downloads/ChileanLifeTables.xlsx")

FLT_1 <- probs2lifetable(probs = chilean_lifetables$`qx 2020 (RV-M-2020)`, 
                         type = "qx", radix = 100000, 
                         name = "Female Main Annuitants not Disabled LT")

FLT_2 <- probs2lifetable(probs = chilean_lifetables$`qx 2020 (B-M-2020)`, 
                         type = "qx", radix = 100000, 
                         name = "Female Beneficiary not Disabled LT")

FLT_3 <- probs2lifetable(probs = chilean_lifetables$`qx 2020 (MI-M-2020)`, 
                         type = "qx", radix = 100000, 
                         name = "Female Annuitants who are Disabled LT")

MLT_1 <- probs2lifetable(probs = chilean_lifetables$`qx 2020 (CB-H-2020)`, 
                       type = "qx", radix = 100000, 
                       name = "Male not Disabled LT")

MLT_2 <- probs2lifetable(probs = chilean_lifetables$`qx 2020 (MI-H-2020)`, 
                         type = "qx", radix = 100000, 
                         name = "Males who are Disabled LT")


# Plot of Lx vs x
plot(FLT_1@x, FLT_1@lx, type = "l", col = "red", xlab = "x (age)", 
     ylab = "lx (initial population)", main = "lx vs x")

lines(FLT_2@lx, type = "l", col = "green")
lines(FLT_3@lx, type = "l", col = "skyblue")
lines(MLT_1@lx, type = "l", col = "black")
lines(MLT_2@lx, type = "l", col = "purple")

legend("bottomleft", legend = c("Female Main Annuitants who are not diasbled", 
                                "Female Beneficiaries who are not Disabled", 
                                "Female Annuitants who are Disabled", 
                                "Males Annuitants who are not Disabled", 
                                "Males Annuitants who are Disabled" ),
       col = c("red", "green", "skyblue", "black", "purple"),
       lty = 1, cex = 0.45, x.intersp = 1.2, y.intersp = 0.5)

# Graduation 

# Further data cleaning
annuitant_data <- annuitant_data %>% 
  mutate(age_last_birthday = floor(annuitant_data$START_AGE),
         End_age_last_birthday = floor(annuitant_data$END_AGE))

healthy_annuitant <- annuitant_data %>% filter(HEALTH == "Healthy")


# Account for death and censored observations

# Finding Ex
Ex <- rep(0, 40)
p = 1

for (x in 60:99) {
  i <- 0
  for (row in 1:nrow(healthy_annuitant)) {
    if (healthy_annuitant$age_last_birthday[row] == x && 
        healthy_annuitant$DEATH[row] == 1 &&
        healthy_annuitant$age_last_birthday[row] == floor(healthy_annuitant$END_AGE[row])) {
      i = i + (x + 1 - healthy_annuitant$START_AGE[row]) 
    } else if (healthy_annuitant$age_last_birthday[row] == x && 
               healthy_annuitant$DEATH[row] != 1 &&
               healthy_annuitant$age_last_birthday[row] == floor(healthy_annuitant$END_AGE[row])) {
      i = i + (healthy_annuitant$END_AGE[row] - healthy_annuitant$START_AGE[row])
    } else if (healthy_annuitant$age_last_birthday[row] == x && 
               healthy_annuitant$age_last_birthday[row] != floor(healthy_annuitant$END_AGE[row])) {
      i = i + (x + 1 - healthy_annuitant$START_AGE[row])
    } else if (healthy_annuitant$age_last_birthday[row] < x && 
               floor(healthy_annuitant$END_AGE[row]) > x) {
      i = i + 1
    } else if (healthy_annuitant$age_last_birthday[row] < x &&
               floor(healthy_annuitant$END_AGE[row]) == x 
               && healthy_annuitant$DEATH[row] == 1) {
      i = i + 1
    } else if (healthy_annuitant$age_last_birthday[row] < x &&
               floor(healthy_annuitant$END_AGE[row]) == x 
               && healthy_annuitant$DEATH[row] != 1) {
      i = i + (healthy_annuitant$END_AGE[row] - x)
    }
  }
  Ex[p] <- i
  p = p + 1
}

print(Ex)

# Finding Dx
Dx <- rep(0, 40)
k = 1
for (x in 60:99) {
  i <- 0
  for (row in 1:nrow(healthy_annuitant)) {
    if (healthy_annuitant$End_age_last_birthday[row] == x 
        && healthy_annuitant$DEATH[row] == 1) {
      i = i + 1
    }
  }
  Dx[k] <- i
  k = k + 1
}

print(Dx)

# Finding the crude rates qx
qx <- Dx/Ex


# Finding the mu_x and central exposed risk using UDD assumption
x <- 60:99
mx <- qx2mx(qx, ax = 0.5)
Exc <- Ex - 1/2*Dx

gompertz <- nls(mx ~ exp(b0 + b1*x), start = list(b0 = 1, b1 = 0), 
                weights = Exc/mx)
gompertz

mx_gompertz <- fitted(gompertz)

plot(x, log(mx), pch = 20, xlab = "age", 
     ylab = "Central mortality rate (log scale)", 
     main = "Australian Men: Gompertz law", ylim = range(-11, 0))
lines(x, log(mx_gompertz), col ='blue')

# qx plot
plot(x, 1 - exp(-mx), pch = 20, xlab = "age", ylab = "qx")
lines(x, 1- exp(-mx_gompertz), col = "red")


# Makeham
# Starting point are the parameter estimates for the gompertz model
makeham <- nls(mx ~ A + exp(b0 + b1*x), 
               start = list(A = 0, b0 = coef(gompertz)[1], 
                            b1 = coef(gompertz)[2]), weights = Exc/mx)
makeham

mx_makeham <- fitted(makeham)
plot(x, log(mx), pch = 20, xlab = "age", 
     ylab = "Central mortality rate (log scale)",
     main = "Australian Men: Makeham law",
     ylim = range(-11, 0))
lines(x,log(mx_makeham), col = 'red')

# qx plot
plot(x, 1 - exp(-mx), pch = 20, xlab = "age", ylab = "qx")
lines(x, 1- exp(-mx_makeham), col = "blue")


# Analysis of parametric laws
plot(x, log(mx), pch = 20, xlab = "age", 
     ylab = "Central mortality rate (log scale)", 
     main = "Annuitant Data", ylim = range(-11, 0))
lines(x, log(mx_gompertz), col ='blue')
lines(x,log(mx_makeham), col = 'red')
legend("bottomleft", legend = c("Gompertz", "Makeham"), 
       col = c("blue", "red"), lty = 1, cex = 0.6, 
       x.intersp = 2.5, y.intersp = 0.5)

# Smoothing Splines 
# The smoothing parameter has been obtained via trial and error
smSpline <- smooth.spline(x, qx, spar = 0.6)
smSpline

qx_smSpline <- fitted(smSpline)
plot(x, qx, pch = 20, xlab = "age", 
     ylab = "Central mortality rate (log scale)",
     main = "Australian Men: Smoothing spline")
lines(x, qx_smSpline, col = 'green')


# Using the knots chosen for Australian Life Tables ; weighted average 
# To handle the imbalance in the data, gender prop
knots <- c(67.866, 77.433, 88.716, 94.711)
cubic_basis <- ns(x, knots = knots)
cubSpline <- lm(qx ~ cubic_basis, weights = Ex/qx)
cubSpline

qx_cubSpline <- fitted(cubSpline)
plot(x, qx, pch = 20, xlab = "age", 
     ylab = "Central mortality rate (log scale)",
     main = "Australian Men: Natural cubic spline")
lines(x, qx_cubSpline, col = 'purple')


# Graduated Plots 
plot(x, qx, pch = 20, xlab = "ages", ylab = "qx", 
     main = "Mortality Rates of Annuitants")
lines(x, 1- exp(-mx_gompertz), col = "red", lwd = 2)
lines(x, 1- exp(-mx_makeham), col = "blue", lwd = 2)
lines(x, qx_smSpline, col = 'green', lwd = 2) 
lines(x, qx_cubSpline, col = 'purple', lwd = 2)
legend("topleft", legend = c("Gompertz", "Makeham","Smoothing Spline",
                             "Natural Cubic Spline", "Crude Rates"), 
       col = c("red", "blue", "green", "purple", "black"), lty = 1, lwd = 2,
       cex = 0.6, x.intersp = 2.5, y.intersp = 0.5)

# Statistical Analysis 
# qx = Dx/Ex

zx_makeham <- (Dx - Ex * (1- exp(-mx_makeham)))/ sqrt(Ex * (1- exp(-mx_makeham)))
zx_gompertz <- (Dx -  Ex * (1- exp(-mx_gompertz))) / sqrt(Ex * (1- exp(-mx_gompertz)))
zx_cubSpline <- (Dx - Ex * qx_cubSpline) / sqrt(Ex * qx_cubSpline)
zx_smSpline <- (Dx - Ex * qx_smSpline) / sqrt(Ex * qx_smSpline)

# Chi Squared Goodness of Fit Test 
chi2Test <- function(O, E, npar, alpha = 0.05) {
  chi2 <- sum((O - E)^2 / E) # Test statistic
  df <- length(O) - npar
  chi2_alpha <- qchisq(1 - alpha, df) # Critical value
  p.value <- 1 - pchisq(chi2, df) # p value
  list(statistic = chi2, c.value = chi2_alpha, df = df, p.value = p.value)
}

# Low p value indicates the mortality experience does not conform to gompertz law
chi2Test(Dx, Ex * (1- exp(-mx_gompertz)), length(coef(gompertz)))

# Low p value --> does not conform to makeham law
chi2Test(Dx, Ex * (1- exp(-mx_makeham)), length(coef(makeham)))

# Low p value --> does not conform 
chi2Test(Dx, Ex * qx_cubSpline, cubSpline$rank)

# Passes the chi squared goodness of fit
chi2Test(Dx, Ex * qx_smSpline, smSpline$df)

# Standardised Deviations Test 
stdTest <- function(zx, breaks = c(-Inf, -1, 0, 1, Inf)) {
  observed <- table(cut(zx, breaks)) # count observation in each interval
  expected.p <- diff(pnorm(breaks)) # expected probabilities for standard normal
  chisq.test(observed, p = expected.p) # apply chi squared test
}


stdTest_gompertz <- stdTest(zx_gompertz)
stdTest_makeham <- stdTest(zx_makeham)
stdTest_cubSpline <- stdTest(zx_cubSpline)
stdTest_smSpline <- stdTest(zx_smSpline)

# Low p values indicates the deviations don't conform to a standard normal distribution
stdTest_gompertz

# High p value indicates the deviations conforms to a standard normal distribution
# Could be a sign of over graduation
stdTest_makeham

# High p value indicates the deviations conforms to a standard normal distribution
stdTest_cubSpline

# High p value indicates the deviations conforms to a standard normal distribution
stdTest_smSpline



# Also use qq plots
# Horrible
qqnorm(zx_gompertz)
qqline(zx_gompertz, col = "steelblue", lwd = 2)

# Alright
qqnorm(zx_makeham)
qqline(zx_makeham, col = "steelblue", lwd = 2)

# Alright
qqnorm(zx_cubSpline)
qqline(zx_cubSpline, col = "steelblue", lwd = 2)

# Alright
qqnorm(zx_smSpline)
qqline(zx_smSpline, col = "steelblue", lwd = 2)

# Complimentary tests to the chi-squared tests

# Signs Test 
nages <- length(x)
signTest_gompertz <- binom.test(sum(zx_gompertz > 0), nages)
signTest_makeham <- binom.test(sum(zx_makeham > 0), nages)
signTest_cubSpline <- binom.test(sum(zx_cubSpline > 0), nages)
signTest_smSpline <- binom.test(sum(zx_smSpline > 0), nages)



# Don't reject --> +ve signs for certain periods can be cancelled out by negative signs
signTest_gompertz

# Don't reject --> +ve signs for certain periods can be cancelled out by negative signs
signTest_makeham

# Don't reject --> +ve signs for certain periods can be cancelled out by negative signs
signTest_cubSpline

# Don't reject --> +ve signs for certain periods can be cancelled out by negative signs
# The number of +ve and -ve deviations balance each other out
signTest_smSpline

# Cumulative Deviations Test  --> goodness of fit test 
cumDevTest <- function(A, E, alpha = 0.05) {
  cumDev <- sum(A - E) / sqrt(sum(E)) # Test statistic
  z_alpha <- qnorm(1 - alpha/2) # Critical value
  p.value <- 2 *(1 - pnorm(cumDev)) # p value (Note it is two-tailed)
list(statistic = cumDev, c.value = z_alpha, p.value = p.value)
}

cumDevTest_gompertz <- cumDevTest(Dx, Ex * (1- exp(-mx_gompertz)))
cumDevTest_makeham <- cumDevTest(Dx, Ex * (1- exp(-mx_makeham)))
cumDevTest_cubSpline <- cumDevTest(Dx, Ex * qx_cubSpline)
cumDevTest_smSpline <- cumDevTest(Dx, Ex * qx_smSpline)

# The mortality experience conforms to the data
cumDevTest_gompertz

# The mortality experience conforms to the data
cumDevTest_makeham

# The mortality experience conforms to the data
cumDevTest_cubSpline

# The mortality experience conforms to the data
cumDevTest_smSpline

# Indicates the data is not biased. 

# Grouping of Signs Test 

groupSignTest <- function(zx, alpha = 0.05){
  
  #Count +'s and -'s
  signs <- sign(zx)
  n1 <- sum(signs == 1)
  n2 <- sum(signs == -1)
  
  #Count runs
  y <- c(-1, sign(zx))
  G <- sum((y[-1] != y[-(n1 + n2 + 1)]) & y[-1]!= -1) # No Runs
  
  #Normal approximation
  mu <- n1 * (n2 + 1) / (n1 + n2)
  s2 <- (n1 * n2)^2 / (n1 + n2)^3
  G_alpha <- qnorm(alpha, mean = mu, sd = sqrt(s2)) #Critical value
  p.value <- (pnorm(G + 0.5, mean = mu, sd = sqrt(s2))) #p value (one sided)
  list(statistic = G, c.value = G_alpha, p.value = p.value)
}

groupSignTest_gompertz <- groupSignTest(zx_gompertz)
groupSignTest_makeham <- groupSignTest(zx_makeham)
groupSignTest_cubSpline <- groupSignTest(zx_cubSpline)
groupSignTest_smSpline <- groupSignTest(zx_smSpline)

# Normal Approximation
# 4 positive runs --> not enough reject hypothesis --> overgraduation
groupSignTest_gompertz

# 8 positive runs --> not enough evidence to reject the null
groupSignTest_makeham

# 8 positive runs --> not enough evidence to reject the null
groupSignTest_cubSpline

# 8 positive runs --> not enough evidence to reject the null
groupSignTest_smSpline

# WHy do we need a large number of positive runs


# Serial Correlation Test
# Shows significant correlation --> the curve tends to be biased for relatively
# long periods --> we also don't want it to be constantly oscillating; they 
# will be negatively correlated
par(mfrow = c(1, 1))
# Consider lag 1 --> overgraduated
acf(zx_gompertz) 

# Shows significant correlation --> overgraduated
acf(zx_makeham)

# Doesn't show significant correlation --> overgraduated
acf(zx_cubSpline)
# Doesn't show significant correlation --> not overgraduated
acf(zx_smSpline)

# Choose the right one compare 
# Graduation Plot 

chilean_lifetables_clean <- chilean_lifetables %>% filter(between(Age, 60, 99))

plot(chilean_lifetables_clean$Age, 
     chilean_lifetables_clean$`qx 2020 (CB-H-2020)`, type = "l", 
     xlab = "x (age)", ylab = "qx (initial rate of mortality)",
     main = "qx vs x")
lines(chilean_lifetables_clean$Age, 
      chilean_lifetables_clean$`qx 2020 (RV-M-2020)`, col = "red")
lines(chilean_lifetables_clean$Age, 
      chilean_lifetables_clean$`qx 2020 (B-M-2020)`, col = "skyblue")
lines(x, qx_cubSpline, col = 'purple')
legend("topleft", legend = c("Males Annuitants who are not Disabled",
                             "Female Main Annuitants who are not diasbled",
                             "Female Beneficiaries who are not Disabled",
                             "Graduated Rates"),
       col = c("black", "red", "skyblue", "purple"),
       lty = 1, cex = 0.45, x.intersp = 0.4, y.intersp = 0.5)




























