# Libraries
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(corrplot)
library(tidyr)
library(patchwork)

# ============================================================
# DATA PRE PROCESSING
# ============================================================

# set wd and file path
setwd("/Users/baste/Documents/GitHub/math2269-bayesian-workforce-analysis")
raw_file <- "data/ANZSCO Occupation data - February 2026.xlsx"

# load separate , relevant sheets
# remove duplicate Occupation columns
overview <- read_excel(raw_file, sheet = "Table_1", skip = 6, na = c("", "N/A")) 
earnings_hours <- read_excel(raw_file, sheet = "Table_4", skip = 6, na = c("", "N/A")) %>%
  select(-Occupation)
states <- read_excel(raw_file, sheet = "Table_6", skip = 6, na = c("", "N/A")) %>%
  select(-Occupation)
age_profile <- read_excel(raw_file, sheet = "Table_7", skip = 6, na = c("", "N/A")) %>%
  select(-Occupation)
education <- read_excel(raw_file, sheet = "Table_8", skip = 6, na = c("", "N/A")) %>%
  select(-Occupation)

# merge sheets
anzsco_clean <- overview %>%
  left_join(earnings_hours, by = "ANZSCO Code") %>%
  left_join(age_profile, by = "ANZSCO Code") %>%
  left_join(education, by = "ANZSCO Code")

#filter to 4 digit anzsco codes
# correct variable types and rename columns
anzsco_clean <- anzsco_clean %>%
  filter(`ANZSCO Code` >= 1000 & `ANZSCO Code` <= 9999) %>%
  mutate(Employed = as.numeric(Employed)) %>%
  rename(
    anzsco_code = `ANZSCO Code`,
    occupation = Occupation,
    employed = Employed,
    part_time_share = `Part-time share (%)`,
    female_share = `Female share (%)`,
    median_weekly_earnings = `Median weekly earnings ($)`,
    median_age = `Median age`,
    annual_employment_growth = `Annual employment growth`,
    full_time_share = `Share of workers who work full-time hours (%)`,
    average_full_time_hours = `Average full-time hours worked per week`,
    median_full_time_weekly_earnings = `Median full-time earnings per week ($)`,
    median_full_time_hourly_earnings = `Median full-time hourly earnings ($)`,
    age_15_19 = `15 - 19\r\n(%)`,
    age_20_24 = `20 - 24\r\n(%)`,
    age_25_34 = `25 - 34\r\n(%)`,
    age_35_44 = `35 - 44\r\n(%)`,
    age_45_54 = `45 - 54\r\n(%)`,
    age_55_59 = `55 - 59\r\n(%)`,
    age_60_64 = `60 - 64\r\n(%)`,
    age_65_plus = `65 +\r\n(%)`,
    postgraduate_share = `Post Graduate/ Graduate Diploma or Graduate Certificate \r\n(%)`,
    bachelor_share = `Bachelor degree\r\n(%)`,
    advanced_diploma_share = `Advanced Diploma/ Diploma \r\n(%)`,
    certificate_iii_iv_share = `Certificate III/ IV \r\n(%)`,
    year_12_share = `Year 12\r\n (%)`,
    year_11_share = `Year 11 \r\n(%)`,
    year_10_below_share = `Year 10 and below \r\n(%)`)

# derive analysis variables
anzsco_clean <- anzsco_clean %>%
  mutate(
    previous_year_employment = employed - annual_employment_growth,
    employment_growth_rate = annual_employment_growth / previous_year_employment * 100,
    young_worker_share = age_15_19 + age_20_24 + age_25_34,
    older_worker_share = age_55_59 + age_60_64 + age_65_plus,
    higher_education_share = postgraduate_share + bachelor_share,
    tertiary_education_share = postgraduate_share + bachelor_share + advanced_diploma_share)

# ============================================================
# INITIAL DESCRIPTIVE ANALYSIS
# ============================================================

# summaries and inspection of data structure
glimpse(anzsco_clean)
summary(anzsco_clean$employment_growth_rate)
colSums(is.na(anzsco_clean))

# examine outcome variable
# occupations with largest declines
anzsco_clean %>%
  select(anzsco_code, occupation, employed, annual_employment_growth, employment_growth_rate) %>%
  arrange(employment_growth_rate) %>%
  head(10)

# occupations with largest growth
anzsco_clean %>%
  select(anzsco_code, occupation, employed, annual_employment_growth, employment_growth_rate) %>%
  arrange(desc(employment_growth_rate)) %>%
  head(10)

# distribution of employment growth rate
ggplot(anzsco_clean, aes(x = employment_growth_rate)) +
  geom_histogram(
    bins = 25,
    fill = "#4C78A8",
    colour = "white"
  ) +
  labs(
    title = "Distribution of Annual Employment Growth Rates",
    x = "Annual employment growth rate (%)",
    y = "Number of occupations"
  ) +
  theme_minimal()

# examine employment growth rate by occupation size
ggplot(anzsco_clean, aes(x = employed, y = employment_growth_rate) +
  geom_point() +
  scale_x_log10() +
  labs(
    title = "Employment Growth Rate by Occupation Size",
    x = "Number employed (log scale)",
    y = "Annual employment growth rate (%)") +
  theme_minimal()

# examine relationship size of occupations and their growth rates
cor(
  log(anzsco_clean$employed),
  abs(anzsco_clean$employment_growth_rate),
  use = "complete.obs")

# visualise
ggplot(anzsco_clean, aes(x = employed, y = abs(employment_growth_rate))) +
  geom_point() +
  scale_x_log10() +
  labs(
    title = "Absolute Employment Growth Rate by Occupation Size",
    x = "Number employed (log scale)",
    y = "Absolute annual employment growth rate (%)") +
  theme_minimal()

# log transform occupation size since employment is right skewed
# for descriptive analysis purposes
anzsco_clean <- anzsco_clean %>%
  mutate(log_employed = log(employed))

# summarise outcome and predictor variables
descriptive_data <- anzsco_clean %>%
  select(
    employment_growth_rate,
    employed,
    log_employed,
    median_age,
    young_worker_share,
    older_worker_share,
    higher_education_share,
    tertiary_education_share,
    female_share,
    part_time_share,
    average_full_time_hours
  )

summary(descriptive_data)

# summary statistics
summary_stats <- descriptive_data %>%
  summarise(
    across(
      everything(),
      list(
        N = ~ sum(!is.na(.)),
        Mean = ~ mean(., na.rm = TRUE),
        SD = ~ sd(., na.rm = TRUE),
        Median = ~ median(., na.rm = TRUE),
        Q1 = ~ quantile(., 0.25, na.rm = TRUE),
        Q3 = ~ quantile(., 0.75, na.rm = TRUE),
        Min = ~ min(., na.rm = TRUE),
        Max = ~ max(., na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "(.*)_(N|Mean|SD|Median|Q1|Q3|Min|Max)"
  ) %>%
  mutate(
    across(
      c(Mean, SD, Median, Q1, Q3, Min, Max),
      ~ round(., 2)
    )
  )

summary_stats

# report summary stats
report_descriptive_data <- anzsco_clean %>%
  select(
    employment_growth_rate,
    log_employed,
    median_age,
    higher_education_share,
    female_share,
    part_time_share,
    average_full_time_hours
  )

# examine distributions of final predictors

# log employed
p1 <- ggplot(anzsco_clean, aes(x = log_employed)) +
  geom_histogram(
    bins = 20,
    fill = "#4C78A8",
    colour = "white"
  ) +
  labs(
    title = "Log Employed",
    x = "Log number employed",
    y = "Number of occupations"
  ) +
  theme_minimal()

# median age
p2 <- ggplot(anzsco_clean, aes(x = median_age)) +
  geom_histogram(
    bins = 20,
    fill = "#4C78A8",
    colour = "white"
  ) +
  labs(
    title = "Median Age",
    x = "Median age",
    y = "Number of occupations"
  ) +
  theme_minimal()

# higher education share
p3 <- ggplot(anzsco_clean, aes(x = higher_education_share)) +
  geom_histogram(
    bins = 20,
    fill = "#4C78A8",
    colour = "white"
  ) +
  labs(
    title = "Higher Education Share",
    x = "Workers with higher education (%)",
    y = "Number of occupations"
  ) +
  theme_minimal()

# female share
p4 <- ggplot(anzsco_clean, aes(x = female_share)) +
  geom_histogram(
    bins = 20,
    fill = "#4C78A8",
    colour = "white"
  ) +
  labs(
    title = "Female Share",
    x = "Female workers (%)",
    y = "Number of occupations"
  ) +
  theme_minimal()

# part-time share
p5 <- ggplot(anzsco_clean, aes(x = part_time_share)) +
  geom_histogram(
    bins = 20,
    fill = "#4C78A8",
    colour = "white"
  ) +
  labs(
    title = "Part-time Share",
    x = "Part-time workers (%)",
    y = "Number of occupations"
  ) +
  theme_minimal()

# average full-time hours
p6 <- ggplot(anzsco_clean, aes(x = average_full_time_hours)) +
  geom_histogram(
    bins = 20,
    fill = "#4C78A8",
    colour = "white"
  ) +
  labs(
    title = "Average Full-time Hours",
    x = "Average full-time hours per week",
    y = "Number of occupations"
  ) +
  theme_minimal()

(p1 + p2 + p3 + p4 + p5 + p6) +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Distribution of Final Predictor Variables"
  )

# examine correlations among candidate predictors
candidate_vars <- anzsco_clean %>%
  select(
    employment_growth_rate,
    employed,
    median_age,
    young_worker_share,
    older_worker_share,
    higher_education_share,
    tertiary_education_share,
    female_share,
    part_time_share,
    average_full_time_hours
  )

cor_matrix <- cor(
  candidate_vars,
  use = "pairwise.complete.obs"
)

round(cor_matrix, 2)

# visualise the correlation matrix
corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  tl.cex = 0.8,
  tl.col = "black",
  tl.srt = 45,
  diag = FALSE
)

# Correlation matrix for the response and final candidate predictors
final_cor_vars <- anzsco_clean %>%
  select(
    employment_growth_rate,
    log_employed,
    median_age,
    higher_education_share,
    female_share,
    part_time_share,
    average_full_time_hours
  )

# Calculate Pearson correlations using available complete pairs
final_cor_matrix <- cor(
  final_cor_vars,
  use = "pairwise.complete.obs"
)

round(final_cor_matrix, 2)

# examine relationships between candidate predictors and employment growth

# occupation size
p1 <- ggplot(
  anzsco_clean,
  aes(x = log_employed, y = employment_growth_rate)
) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Log Employed",
    x = "Log number employed",
    y = "Annual employment growth rate (%)"
  ) +
  theme_minimal()

# median age
p2 <- ggplot(
  anzsco_clean,
  aes(x = median_age, y = employment_growth_rate)
) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Median Age",
    x = "Median age",
    y = "Annual employment growth rate (%)"
  ) +
  theme_minimal()

# higher education share
p3 <- ggplot(
  anzsco_clean,
  aes(x = higher_education_share, y = employment_growth_rate)
) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Higher Education Share",
    x = "Workers with higher education (%)",
    y = "Annual employment growth rate (%)"
  ) +
  theme_minimal()

# female share
p4 <- ggplot(
  anzsco_clean,
  aes(x = female_share, y = employment_growth_rate)
) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Female Share",
    x = "Female workers (%)",
    y = "Annual employment growth rate (%)"
  ) +
  theme_minimal()

# part-time share
p5 <- ggplot(
  anzsco_clean,
  aes(x = part_time_share, y = employment_growth_rate)
) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Part-time Share",
    x = "Part-time workers (%)",
    y = "Annual employment growth rate (%)"
  ) +
  theme_minimal()

# average full-time hours
p6 <- ggplot(
  anzsco_clean,
  aes(x = average_full_time_hours, y = employment_growth_rate)
) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Average Full-time Hours",
    x = "Average full-time hours per week",
    y = "Annual employment growth rate (%)"
  ) +
  theme_minimal()

((p1 | p2 | p3) /
    (p4 | p5 | p6)) +
  plot_annotation(
    title = "Relationships Between Candidate Predictors and Employment Growth"
  )

# ------------------------------------------------------------
# EXPLORATORY CHECKS
# Done to help with choosing the Bayesian model to do
# ------------------------------------------------------------

# create anzsco major groups
anzsco_clean <- anzsco_clean %>%
  mutate(
    major_group_code = floor(anzsco_code / 1000),
    major_group = case_when(
      major_group_code == 1 ~ "Managers",
      major_group_code == 2 ~ "Professionals",
      major_group_code == 3 ~ "Technicians and Trades Workers",
      major_group_code == 4 ~ "Community and Personal Service Workers",
      major_group_code == 5 ~ "Clerical and Administrative Workers",
      major_group_code == 6 ~ "Sales Workers",
      major_group_code == 7 ~ "Machinery Operators and Drivers",
      major_group_code == 8 ~ "Labourers"))

# number of occupations in each group
anzsco_clean %>%
  count(major_group)

# growth summary
group_summary <- anzsco_clean %>%
  group_by(major_group) %>%
  summarise(
    n = sum(!is.na(employment_growth_rate)),
    mean_growth = mean(employment_growth_rate, na.rm = TRUE),
    sd_growth = sd(employment_growth_rate, na.rm = TRUE)
  )

group_summary

# compare growth distribution across major groups
ggplot(
  anzsco_clean,
  aes(x = reorder(major_group, employment_growth_rate, FUN = median), y = employment_growth_rate)) +
  geom_boxplot() +
  coord_flip() +
  labs(
    title = "Employment Growth Rate by ANZSCO Major Group",
    x = "ANZSCO major group",
    y = "Annual employment growth rate (%)") +
  theme_minimal()

# variation between group means
sd(group_summary$mean_growth)

# Prepare candidate datasets for model comparison
# Select the outcome and non-redundant candidate predictors
# identified from the descriptive analysis.
# Occupation code and name are retained for reference only.

model_data <- anzsco_clean %>%
  select(
    anzsco_code,
    occupation,
    employment_growth_rate,
    log_employed,
    median_age,
    higher_education_share,
    female_share,
    part_time_share,
    average_full_time_hours
  )

# check missing values
colSums(is.na(model_data))

# check occupations that would be excluded bc of missing info
model_data %>%
  filter(
    is.na(employment_growth_rate) |
      is.na(female_share) |
      is.na(part_time_share)
  ) %>%
  select(
    occupation,
    employment_growth_rate,
    female_share,
    part_time_share
  )

# check complete case sample size
model_data %>%
  tidyr::drop_na() %>%
  nrow()

# ============================================================
# BAYESIAN MODEL PREPARATION
# ============================================================

# Prepare complete-case modelling dataset
# Keeps ANZSCO identifiers for reference and later prediction
model_data <- anzsco_clean %>%
  select(
    anzsco_code,
    occupation,
    employment_growth_rate,
    log_employed,
    median_age,
    higher_education_share,
    female_share,
    part_time_share,
    average_full_time_hours
  ) %>%
  tidyr::drop_na()

# Check number of complete observations used for modelling
nrow(model_data)


# ============================================================
# EXPORT COMPLETE MODEL DATA
# ============================================================

# Save dataset with ANZSCO code and occupation name
write.csv(
  model_data,
  "data/model_data.csv",
  row.names = FALSE
)


# ============================================================
# EXPORT DATA FOR JAGS
# ============================================================

# Keep only the response and predictors required by JAGS
jags_data <- model_data %>%
  select(
    employment_growth_rate,
    log_employed,
    median_age,
    higher_education_share,
    female_share,
    part_time_share,
    average_full_time_hours
  )

# Check JAGS dataset
nrow(jags_data)
names(jags_data)

# Save JAGS modelling dataset
write.csv(
  jags_data,
  "data/jags_model_data.csv",
  row.names = FALSE
)

# ============================================================
# RESIDUAL DIAGNOSTICS
# ============================================================

# Fit a preliminary multiple linear regression
# This is used only to examine the residuals and help choose
# between a Normal and robust Student-t likelihood.
preliminary_model <- lm(
  employment_growth_rate ~
    log_employed +
    median_age +
    higher_education_share +
    female_share +
    part_time_share +
    average_full_time_hours,
  data = model_data
)

# Create a separate dataset for residual diagnostics
diagnostic_data <- model_data %>%
  mutate(
    fitted_value = fitted(preliminary_model),
    residual = resid(preliminary_model)
  )


# Residual distribution
p1 <- ggplot(diagnostic_data, aes(x = residual)) +
  geom_histogram(
    bins = 25,
    fill = "#4C78A8",
    colour = "white"
  ) +
  labs(
    title = "Residual Distribution",
    x = "Residual",
    y = "Number of occupations"
  ) +
  theme_minimal()


# Normal Q-Q plot
p2 <- ggplot(diagnostic_data, aes(sample = residual)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Normal Q-Q Plot",
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  ) +
  theme_minimal()


# Residuals against fitted values
p3 <- ggplot(
  diagnostic_data,
  aes(
    x = fitted_value,
    y = residual
  )
) +
  geom_point() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Residuals vs Fitted Values",
    x = "Fitted employment growth rate (%)",
    y = "Residual"
  ) +
  theme_minimal()


patchwork::wrap_plots(
  p1, p2, p3,
  ncol = 2
) +
  plot_annotation(
    title = "Preliminary Residual Diagnostics for Likelihood Selection"
  )
