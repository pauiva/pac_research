# Title: Perceived Accessibility Analysis
#Related publication: Perceived Accessibility Scale adapted to cycling – What insights can it provide in the context of Stockholm? 
  # Link: 
# Description: This script computes descriptive statistics, CFA, non-parametric tests, 
#              odds ratios, and logistic regression for perceived accessibility (PAC) data.
# Author: Ivana Paulusova
# Date: 2025-08-02

# Required packages
library(readr)
library(questionr)
library(dunn.test)
library(lavaan)
library(psych)
library(epitools)
library(pscl)
library(dplyr)
library(lavaan)


# Data 
df_access <- read_csv("individuals_with_poi_densities.csv")


# PAC descriptives and CFA ------------------------------------------------
pac <- df_access[ , c("q10_easy", 
                      "q11_live_life_as_want",
                      "q12_able_to_do_activities", 
                      "q13_satisfactory_access")]

corr.test(pac)
psych::describe(pac)
psych::alpha(pac)


### CFA
cfa_model_pac <- '
  PAC_factor =~ q10_easy + q11_live_life_as_want + q12_able_to_do_activities + q13_satisfactory_access
'
cfa_pac_fit <- cfa(cfa_model_pac, data = df_access)
summary(cfa_pac_fit, fit.measures = TRUE, standardized = TRUE)

### Compute factor scores
factor_scores_pac <- lavPredict(cfa_pac_fit)
factor_scores_df_pac <- as.data.frame(factor_scores_pac)

df_access <- cbind(df_access, factor_scores_df_pac)

### Obtain and categorize Above and Below PAC and spatial values 
mean_PAC_factor <- mean(df_access$PAC_factor, na.rm = TRUE)

df_access$PAC_cat <- ifelse(df_access$PAC_factor <= mean_PAC_factor, 
                                                "Below Average PB", 
                                                "Above Average PB")


mean_bike_infra <- mean(df_access$high_low_ratio, na.rm = TRUE)

df_access$bike_infra_cat <- ifelse(
  !is.na(df_access$high_low_ratio) & df_access$high_low_ratio > mean_bike_infra,
  "Above average",
  "Below average"
)



dens_vars <- c("dens_shops", "dens_amenity", "dens_leisure", "dens_office")

dens_means <- sapply(df_access[dens_vars], function(x) mean(x, na.rm = TRUE))

for (v in dens_vars) {
  thr <- dens_means[[v]]
  newcol <- paste0(v, "_cat")
  df_access[[newcol]] <- ifelse(
    !is.na(df_access[[v]]) & df_access[[v]] > thr, "Above average", "Below average"
  )
}


# Desriptives - means -----------------------------------------------------
group_vars <- c("gender", "age_group", "education", "living_children", 
                "car_household_bin", "car_access", "driving_license", 
                "sl_card", "cycle_confidence", "commute_lenght_bin", 
                "peer_commute", "peer_commute_bin")

for (var in group_vars) {
  cat("===", var, "===\n")
  print(tapply(df_access$PAC_factor, df_access[[var]], mean, na.rm = TRUE))
  cat("\n")
}

# Non-parametric tests ----------------------------------------------------
#Mann-Whitney U test for 2 categories
binary_vars <- c("gender", "education", "living_children", "car_household_bin", 
                 "driving_license", "sl_card", "commute_lenght_bin", "peer_commute_bin")

for (var in binary_vars) {
  cat("\nVariable:", var, "\n")
  print(wilcox.test(df_access$PAC_factor ~ df_access[[var]]))
}


#Kruskal test if more than 2 categories
multi_vars <- c("age_group", "car_access", "cycle_confidence", "peer_commute")

for (var in multi_vars) {
  cat("\nKruskal-Wallis test for:", var, "\n")
  print(kruskal.test(df_access$PAC_factor ~ df_access[[var]]))
}

##to test differences between the groups
posthoc_vars <- c("cycle_confidence", "peer_commute")

for (var in posthoc_vars) {
  cat("\nDunn test for:", var, "\n")
  print(dunn.test(df_access$PAC_factor, df_access[[var]], method = "bonferroni"))
}

### Spearman correlation for spatial variables 
cor.test(df_access$PAC_factor, df_access$high_low_ratio,
         method = "spearman", use = "complete.obs")


cor.test(df_access$PAC_factor, df_access$dens_shops,
         method = "spearman", use = "complete.obs")

cor.test(df_access$PAC_factor, df_access$dens_amenity,
         method = "spearman", use = "complete.obs")

cor.test(df_access$PAC_factor, df_access$dens_office,
         method = "spearman", use = "complete.obs")

cor.test(df_access$PAC_factor, df_access$dens_leisure,
         method = "spearman", use = "complete.obs")




# Conditional correlations ------------------------------------------------
# AGE 
cat("\n=== Age Group ===\n")
contingency_table <- table(df_access$age_group, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

fisher.test(contingency_table[,,1])                #low bike frequency
chisq.test(contingency_table[,,2], correct=FALSE)  #high bike frequency

# EDUCATION
cat("\n=== Education ===\n")
contingency_table <- table(df_access$education, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

fisher.test(contingency_table[,,1])
chisq.test(contingency_table[,,2], correct = FALSE)

# GENDER + ODDS RATIO
cat("\n=== Gender ===\n")
contingency_table <- table(df_access$gender, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

chisq.test(contingency_table[,,1], correct=FALSE)
chisq.test(contingency_table[,,2], correct=FALSE)

#Odds ratios - for high frequency
high_frequency_table <- contingency_table[,,2]

a <- high_frequency_table[1,1]  # Female (0) and Above Average PB
b <- high_frequency_table[1,2]  # Female (0) and Below Average PB
c <- high_frequency_table[2,1]  # Male (1) and Above Average PB
d <- high_frequency_table[2,2]  # Male (1) and Below Average PB

oddsratio(matrix(c(c, d, a, b), ncol = 2))

cat("\n=== Living with Children ===\n")
contingency_table <- table(df_access$living_children, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

chisq.test(contingency_table[,,1], correct=FALSE)
chisq.test(contingency_table[,,2], correct=FALSE)


# CAR HOUSEHOLD + ODDS RATIO
cat("\n=== Car Household ===\n")
contingency_table <- table(df_access$car_household_bin, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

chisq.test(contingency_table[,,1], correct=FALSE) 
chisq.test(contingency_table[,,2], correct=FALSE) 

#Odds ratios
high_frequency_table <- contingency_table[,,2]

a <- high_frequency_table[1,1]  # No car (0) and Above average
b <- high_frequency_table[1,2]  # No car(0) and Below average
c <- high_frequency_table[2,1]  # Car (1) and Above average
d <- high_frequency_table[2,2]  # Car (1) and Below average

oddsratio(matrix(c(c, d, a, b), ncol = 2))


# CAR ACCESS + ODDS RATIO
cat("\n=== Car Access ===\n")
contingency_table <- table(df_access$car_access, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

fisher.test(contingency_table[,,1])
chisq.test(contingency_table[,,2], correct=FALSE)

#Odds ratios
high_frequency_table <- contingency_table[,,2]

a <- high_frequency_table[1,1]  # No access and Above Average PB 
b <- high_frequency_table[1,2]  # No access and Below Average PB 
c <- high_frequency_table[3,1]  # Own access and Above Average PB
d <- high_frequency_table[3,2]  # Own access and Below Average PB

oddsratio(matrix(c(c, d, a, b), ncol = 2))

# DRIVING LICENSE
cat("\n=== Driving License ===\n")
contingency_table <- table(df_access$driving_license, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

fisher.test(contingency_table[,,1])
chisq.test(contingency_table[,,2], correct=FALSE)


# SL public transit CARD + ODDS RATIO
cat("\n=== SL Card ===\n")
contingency_table <- table(df_access$sl_card, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

chisq.test(contingency_table[,,1], correct=FALSE)
chisq.test(contingency_table[,,2], correct = FALSE)

#Odds ratios
high_frequency_table <- contingency_table[,,2]

a <- high_frequency_table[1,1]  # No SL (0) and Above Average PB
b <- high_frequency_table[1,2]  # No SL (0) and Below Average PB
c <- high_frequency_table[2,1]  # SL (1) and Above Average PB
d <- high_frequency_table[2,2]  # SL (1) and Below Average PB

oddsratio(matrix(c(c, d, a, b), ncol = 2))

# COMMUTE LENGTH + ODDS RATIO
cat("\n=== Commute Length ===\n")
contingency_table <- table(df_access$commute_lenght_bin, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

chisq.test(contingency_table[,,1], correct=FALSE)
chisq.test(contingency_table[,,2], correct=FALSE)

#Odds ratios
high_frequency_table <- contingency_table[,,2]

a <- high_frequency_table[1,1]  # Short (0) and Above Average PB 
b <- high_frequency_table[1,2]  # Short (0) and Below Average PB
c <- high_frequency_table[2,1]  # Long (1) and Above Average PB
d <- high_frequency_table[2,2]  # Long (1) and Below Average PB

oddsratio(matrix(c(c, d, a, b), ncol = 2))

# PEER COMMUTE (BIN)
cat("\n=== Peer Commute (Binary) ===\n")
contingency_table <- table(df_access$peer_commute_bin, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

fisher.test(contingency_table[,,1])
chisq.test(contingency_table[,,2], correct=FALSE)


# PEER COMMUTE (NUMBER) + ODDS RATIO
cat("\n=== Peer Commute (Number) ===\n")
contingency_table <- table(df_access$peer_commute, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

fisher.test(contingency_table[,,1])
chisq.test(contingency_table[,,2], correct=FALSE)

#Odds ratios
#Low Frequency group
low_frequency_table <- contingency_table[,,1]

a <- low_frequency_table[2,1]  # 1 peer and Above Average PB
b <- low_frequency_table[2,2]  # 1 peer and Below Average PB 
c <- low_frequency_table[4,1]  # 3 peers and Above Average PB 
d <- low_frequency_table[4,2]  # 3 peers and Below Average PB 

oddsratio(matrix(c(c, d, a, b), ncol = 2))

# CYCLE CONFIDENCE
cat("\n=== Cycle Confidence ===\n")
contingency_table <- table(df_access$cycle_confidence, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

fisher.test(contingency_table[,,1])
chisq.test(contingency_table[,,2], correct=FALSE)



# BIKE INFRASTRUCTURE
cat("\n=== Bike infrastructure ===\n")
contingency_table <- table(df_access$bike_infra_cat, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

fisher.test(contingency_table[,,1])
chisq.test(contingency_table[,,2], correct=FALSE)

# DESTINATION DENSITIES

# Example for dens_shops
cat("\n=== Density: shops ===\n")
contingency_table <- table(df_access$dens_shops_cat, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

fisher.test(contingency_table[,,1])
chisq.test(contingency_table[,,2], correct = FALSE)

# Example for dens_amenities
cat("\n=== Density: amenities ===\n")
contingency_table <- table(df_access$dens_amenity_cat, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

fisher.test(contingency_table[,,1])
chisq.test(contingency_table[,,2], correct = FALSE)

# Example for dens_office
cat("\n=== Density: office ===\n")
contingency_table <- table(df_access$dens_office_cat, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

fisher.test(contingency_table[,,1])
chisq.test(contingency_table[,,2], correct = FALSE)

# Example for dens_leisure
cat("\n=== Density: leisure ===\n")
contingency_table <- table(df_access$dens_leisure_cat, 
                           df_access$PAC_cat, 
                           df_access$bike_frequency_bin)

fisher.test(contingency_table[,,1])
chisq.test(contingency_table[,,2], correct = FALSE)




# Binomial logistic regression -----------------------------------------------------
model <- glm(bike_frequency_bin ~ PAC_factor, data = df_access, family = binomial)
summary(model)

# Odds ratios with 95% CI
exp_coef <- exp(cbind(OR = coef(model), confint(model)))
print(round(exp_coef, 3))

# McFadden R²
pR2(model)
