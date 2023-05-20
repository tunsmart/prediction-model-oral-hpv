################################################################################
################################################################################
## Predicition Model for Prognostic enrichment of oral HPV prevention trials ###
################################################################################
################################################################################

# Author: Babatunde Y. Alli
# Date: May 19, 2022

################################
#### Load required packages ####
################################

library(haven)
library(naniar)
library(recipes)
library(tidyverse)
library(rms)
library(gridExtra)
library(tableone)
library(CalibrationCurves)
library(rmda)
library(flextable)
source('perf.R') # Loads in custom performance functions, adapted from van den Goorbergh et al, doi.org/10.1093/jamia/ocac093



# Data wrangling for raw NHANES data obtained from https://wwwn.cdc.gov/nchs/nhanes/Default.aspx
################################################################################

cycles = c('F', 'G', 'H', 'I')
cycle_years = c('2010', '2012', '2014', '2016')
hpv_data = data.frame()
demo_data = data.frame()
sex_behaviour_data = data.frame()
smoking_data = data.frame()

hr_hpv = c("ORXH16", "ORXH18", "ORXH26", "ORXH31", "ORXH33", "ORXH35", "ORXH39", 
           "ORXH45", "ORXH51", "ORXH52", "ORXH53", "ORXH56", "ORXH58", "ORXH59", 
           "ORXH66", "ORXH68", "ORXH73", "ORXH82")

y=0

for (cycle in cycles) {
  y = y + 1
  return = read_xpt(paste("ORHPV_", cycle, ".XPT",sep = "")) %>%
    select(SEQN,HPV=ORXHPV, all_of(hr_hpv))
  return$year = cycle_years[y]
  hpv_data = rbind(hpv_data, return)
  
  return = read_xpt(paste("DEMO_", cycle, ".XPT",sep = "")) %>%
    select(SEQN, gender=RIAGENDR, age=RIDAGEYR, race=RIDRETH1, 
           marital_status=DMDMARTL, cluster=SDMVPSU, stratum=SDMVSTRA, weight=WTINT2YR)
  return$year = cycle_years[y]
  demo_data = rbind(demo_data, return)
  
  return = read_xpt(paste("SXQ_", cycle, ".XPT",sep = "")) %>%
    select(SEQN, ever_sex=SXD021, age_first_oral=SXD621, 
           age_first_sex=SXD031, 
           num_males_performed_oral_sex_on_life=SXQ624, 
           num_females_performed_oral_sex_on_life=SXQ636,
           male_num_sex_with_men_life=SXQ836,
           male_num_sex_with_women_life=SXD171,
           female_num_sex_with_men_life=SXD450,
           female_num_sex_with_women_life=SXQ130
    )
  return$year = cycle_years[y]
  sex_behaviour_data = rbind(sex_behaviour_data, return)
  
  return = read_xpt(paste("SMQ_", cycle, ".XPT",sep = "")) %>%
    select(SEQN, ever_smoker=SMQ020,
           start_age=SMD030,current_smoker=SMQ040, stop_how_long=SMQ050Q,
           unit_of_stop=SMQ050U, past_smoker_dly_cig=SMD057,
           current_smoker_dly_cig=SMD650)
  return$year = cycle_years[y]
  smoking_data = rbind(smoking_data, return)
}

sex_behaviour_data = sex_behaviour_data %>%
  replace_with_na_all(condition = ~.x %in% c(7, 77, 777, 77777, 9, 99, 999, 99999, 2000))

sex_behaviour_data$lifetime_num_oral_sex_partner =  select(sex_behaviour_data, num_males_performed_oral_sex_on_life,
                                                           num_females_performed_oral_sex_on_life) %>%
  rowSums(na.rm = T)
sex_behaviour_data$lifetime_num_oral_sex_partner = ifelse(sex_behaviour_data$lifetime_num_oral_sex_partner < 1 & 
                                                            !sex_behaviour_data$ever_sex %in% c(1,2), 
                                                          NA, sex_behaviour_data$lifetime_num_oral_sex_partner) #make all zero NA

sex_behaviour_data$lifetime_num_sex_partner =  select(sex_behaviour_data, male_num_sex_with_men_life,
                                                      male_num_sex_with_women_life,
                                                      female_num_sex_with_men_life,
                                                      female_num_sex_with_women_life) %>%
  rowSums(na.rm = T)
sex_behaviour_data$lifetime_num_sex_partner = ifelse(sex_behaviour_data$lifetime_num_sex_partner < 1 & 
                                                       !sex_behaviour_data$ever_sex %in% c(1,2), 
                                                     NA, sex_behaviour_data$lifetime_num_sex_partner)



sex_behaviour_data = select(sex_behaviour_data, 
                            SEQN, ever_sex, age_first_oral, age_first_sex,
                            lifetime_num_oral_sex_partner, lifetime_num_sex_partner,
                            year)

sex_behaviour_train = sex_behaviour_data %>%
  filter(year != '2016') %>%
  select(-year)

sex_behaviour_val = sex_behaviour_data %>%
  filter(year == '2016') %>%
  select(-year)



smoking_data = smoking_data %>%
  replace_with_na_all(condition = ~.x %in% c(7, 77, 777, 7777, 77777, 9, 99, 999, 9999, 99999))

smoking_data$stop_how_long = ifelse(smoking_data$stop_how_long == 66666, 50, smoking_data$stop_how_long) 
smoking_data$stop_how_long_years = with(smoking_data, ifelse(unit_of_stop == 4, stop_how_long, 
                                                             ifelse(unit_of_stop == 3, stop_how_long/12, 
                                                                    ifelse(unit_of_stop == 2, stop_how_long/52,
                                                                           ifelse(unit_of_stop == 1, stop_how_long/365, NA)))))
smoking_train = smoking_data %>%
  filter(year != '2016') %>%
  select(-year)

smoking_val = smoking_data %>%
  filter(year == '2016') %>%
  select(-year)

demo_data = demo_data %>%
  replace_with_na_all(condition = ~.x %in% c(7, 77, 777, 7777, 77777, 9, 99, 999, 9999, 99999))

demo_train = demo_data %>%
  filter(year != '2016') %>%
  select(-year)

demo_val = demo_data %>%
  filter(year == '2016') %>%
  select(-year)

hpv_train = hpv_data %>%
  filter(year != '2016') %>%
  select(-year)

hpv_val = hpv_data %>%
  filter(year == '2016') %>%
  select(-year)

data_train = inner_join(hpv_train, sex_behaviour_train, by='SEQN')
data_train = inner_join(data_train, demo_train, by='SEQN')
train_data = inner_join(data_train, smoking_train, by='SEQN')

data_val = inner_join(hpv_val, sex_behaviour_val, by='SEQN')
data_val = inner_join(data_val, demo_val, by='SEQN')
val_data = inner_join(data_val, smoking_val, by='SEQN')

train_data$smoke_status = with(train_data, ifelse(ever_smoker == 2, "Never smoker", 
                                                  ifelse(current_smoker %in% c(1,2), "Current smoker", 
                                                         ifelse(current_smoker == 3, "Former smoker", NA))))
train_data$smoke_years = with(train_data, ifelse(smoke_status == "Current smoker", age - start_age, 
                                                 ifelse(smoke_status == "Former smoker", age - start_age - stop_how_long_years, 
                                                        ifelse(smoke_status == "Never smoker", 0, NA))))
train_data$smk_pk_yrs = with(train_data, ifelse(smoke_status == "Current smoker", (current_smoker_dly_cig * smoke_years)/20, 
                                                ifelse(smoke_status == "Former smoker", (past_smoker_dly_cig * smoke_years)/20, 
                                                       ifelse(smoke_status == "Never smoker", 0, NA))))
train_data_final = train_data %>% 
  select(HPV, starts_with("ORX"), age_first_oral:marital_status, smk_pk_yrs)
train_data_final$HPV = ifelse(train_data_final$HPV == 1, 1, 0)
train_data_final$gender = ifelse(train_data_final$gender == 1, 1, 0)
train_data_final$marital_status = train_data_final$marital_status %>% 
  as.character() %>%
  recode('1'='Married/Living with partner', '2'='Widowed/Divorced/Separated', 
         '3'='Widowed/Divorced/Separated', '4'='Widowed/Divorced/Separated', 
         '5'='Never married', '6'='Married/Living with partner') %>%
  as.factor()
train_data_final$race = train_data_final$race %>% 
  as.character() %>%
  recode('1'='Mexican American', '2'='Other Hispanic', '3'='Non-Hispanic White',
         '4'='Non-Hispanic Black', 
         '5'='Other Race')%>%
  as.factor()



val_data$smoke_status = with(val_data, ifelse(ever_smoker == 2, "Never smoker", 
                                              ifelse(current_smoker %in% c(1,2), "Current smoker", 
                                                     ifelse(current_smoker == 3, "Former smoker", NA))))
val_data$smoke_years = with(val_data, ifelse(smoke_status == "Current smoker", age - start_age, 
                                             ifelse(smoke_status == "Former smoker", age - start_age - stop_how_long_years, 
                                                    ifelse(smoke_status == "Never smoker", 0, NA))))
val_data$smk_pk_yrs = with(val_data, ifelse(smoke_status == "Current smoker", (current_smoker_dly_cig * smoke_years)/20, 
                                            ifelse(smoke_status == "Former smoker", (past_smoker_dly_cig * smoke_years)/20, 
                                                   ifelse(smoke_status == "Never smoker", 0, NA))))
val_data_final = val_data %>% 
  select(HPV, starts_with("ORX"), age_first_oral:marital_status, smk_pk_yrs)
val_data_final$HPV = ifelse(val_data_final$HPV == 1, 1, 0)
val_data_final$gender = ifelse(val_data_final$gender == 1, 1,0)
val_data_final$marital_status = val_data_final$marital_status %>% 
  as.character() %>%
  recode('1'='Married/Living with partner', '2'='Widowed/Divorced/Separated', 
         '3'='Widowed/Divorced/Separated', '4'='Widowed/Divorced/Separated', 
         '5'='Never married', '6'='Married/Living with partner') %>%
  as.factor()
val_data_final$race = val_data_final$race %>% 
  as.character() %>%
  recode('1'='Mexican American', '2'='Other Hispanic', '3'='Non-Hispanic White', 
         '4'='Non-Hispanic Black', 
         '5'='Other Race')%>%
  as.factor()

train_data_final2 = train_data_final %>% select(-HPV)
train_data_final =  train_data_final %>% select(!starts_with("ORX"))

val_data_final2 = val_data_final %>% select(-HPV)
val_data_final =  val_data_final %>% select(!starts_with("ORX"))

############################## End of data wranglng ############################


# Evaluation of missing values
##############################

df = train_data_final
missing.values <- df %>%
  gather(key = "key", value = "val") %>%
  mutate(is.missing = is.na(val)) %>%
  group_by(key, is.missing) %>%
  summarise(num.missing = n()) %>%
  filter(is.missing==T) %>%
  select(-is.missing) %>%
  arrange(desc(num.missing)) 

missing.values %>%
  ggplot() +
  geom_bar(aes(x=key, y=num.missing), stat = 'identity') +
  labs(x='variable', y="number of missing values", title='Number of missing values') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

missing.values <- df %>%
  gather(key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, isna) %>%
  summarise(num.isna = n()) %>%
  mutate(pct = num.isna / total * 100)

levels <- (missing.values  %>% filter(isna == T) %>% arrange(desc(pct)))$key

percentage.plot <- missing.values %>%
  ggplot() +
  geom_bar(aes(x = reorder(key, desc(pct)), 
               y = pct, fill=isna), 
           stat = 'identity', alpha=0.8) +
  scale_x_discrete(limits = levels) +
  scale_fill_manual(name = "", 
                    values = c('steelblue', 'tomato3'), labels = c("Present", "Missing")) +
  coord_flip() +
  labs(title = "Percentage of missing values", x =
         'Variable', y = "% of missing values")

percentage.plot #Plot percentage missing values


df_train <- train_data_final %>% 
  select(-age_first_oral) #remove age_first_oral from further considerations based
                          #on the very high percentage of missing values

df_train2 <- train_data_final2 %>% 
  select(-age_first_oral)


# Table 1 (combines tab1 and tab2 below)
#########

df_table <- df_train %>%
  mutate_at(.vars="gender", factor)
tab1 = CreateTableOne(data=df_table, 
                      vars=c("age", "age_first_sex", "smk_pk_yrs", 
                             "lifetime_num_sex_partner", "lifetime_num_oral_sex_partner",
                             "gender", "race", "marital_status"), strata = "HPV",
                      includeNA=TRUE, test=FALSE)
print(tab1,showAllLevels=TRUE, quote = T, nonnormal=c("age", "age_first_sex", 
                                                      "smk_pk_yrs", 
                                                      "lifetime_num_sex_partner", 
                                                      "lifetime_num_oral_sex_partner"))

df_test <- val_data_final %>% 
  select(-age_first_oral) #remove age_first_oral from further considerations
df_test2 <- val_data_final2 %>% 
  select(-age_first_oral)

df_table2 <- df_test %>%
  mutate_at(.vars="gender", factor)
tab2 = CreateTableOne(data=df_table2, 
                      vars=c("age", "age_first_sex", "smk_pk_yrs", "lifetime_num_sex_partner", 
                             "lifetime_num_oral_sex_partner", "gender", "race", "marital_status"), 
                      strata = "HPV",includeNA=TRUE, test=FALSE)
print(tab2,showAllLevels=TRUE, quote = T, nonnormal=c("age", "age_first_sex", "smk_pk_yrs", 
                                                      "lifetime_num_sex_partner", 
                                                      "lifetime_num_oral_sex_partner"))


# Prep and impute data. The blueprint will also be used for the external temporal validation data later
################################################################################

set.seed(12345)
nhanes_recipe <- recipe(HPV ~ ., data = df_train)

blueprint = nhanes_recipe %>%
  step_impute_knn(all_predictors()) %>%
  step_normalize(age, age_first_sex, 
                 lifetime_num_oral_sex_partner, lifetime_num_sex_partner, smk_pk_yrs)

prepare <- prep(blueprint, training = df_train)

baked_train <- bake(prepare, new_data = df_train)



#######################################################
#######################################################
#### Model with any oral HPV status as the outcome ####
#######################################################
#######################################################

## Fit standard logistic regression for the Any-HPV Model
#########################################################

set.seed(234)

standard_model <- lrm(HPV ~ rcs(age, 4) + rcs(age_first_sex, 4) + 
                        rcs(lifetime_num_oral_sex_partner, 4) + 
                        rcs(lifetime_num_sex_partner, 4) + rcs(smk_pk_yrs, 4) + 
                        gender + race + marital_status,
                      data=baked_train,x=T,y=T, linear.predictors = T)

## Fit and internally validate penalized logistic model
#######################################################

# determine performance over range of penalties
p_trace <- pentrace(standard_model, c(0,2,4,6,8,10,12,14,16,18,20,22, 24, 28, 32,40))
# fit penalized model
penalized_model <- update(standard_model, penalty=p_trace$penalty)
# Internal validation by bootstrapping of penalized model
val_penalized <- validate(penalized_model, B=250)

# Look at the results
penalized_model
val_penalized

## Temporal validation
######################

set.seed(12345)
baked_test <- bake(prepare, new_data = df_test) #Use existing blueprint for data prep

#Get outcome imbalance
n_minor <- train_data_final %>% 
  filter(HPV == 1) %>% 
  nrow()
n_major <- train_data_final %>% 
  filter(HPV == 0) %>% 
  nrow()
imb = n_minor/n_major


baked_test = as.data.frame(baked_test)
test_outcome = baked_test$HPV

ruleLP1	<- predict(penalized_model, baked_test, type="lp")
probLP1  	<- plogis(ruleLP1) # from linear predictor to probabilities

## Calibration plots and model performance
##########################################

pred1 <- rep(0, length(probLP1)) # empty vector for class prediction
pred1[probLP1 > imb] = 1 #predictions based on the class imbalance threshold and not the arbitrary 0.5

# Creating table with all estimated probabilities
probs_table <- cbind(probLP1)

# Creating table with all predicted classes
pred_table <- cbind(pred1)

# Creating empty vectors for performance measures
accuracy_vector <- rep(NA, ncol(pred_table))
sensitivity_vector <- rep(NA, ncol(pred_table))
specificity_vector <- rep(NA, ncol(pred_table))

# Loop over all models and imbalance approaches to get performance measures
for (i in 1:ncol(pred_table)){
  accuracy_vector[i] <- accuracy(pred_table[,i], test_outcome)
  sensitivity_vector[i] <- sensitivity(pred_table[,i], test_outcome)
  specificity_vector[i] <- specificity(pred_table[,i], test_outcome)
}

# Create matrix to store calibration measures (intercept & slope + CI)
calibration_matrix <- matrix(ncol = 6, nrow = ncol(probs_table))

# Calculate calibration measures for all models
for (i in 1:ncol(probs_table)) {
  calibration_matrix[i,] <- calibration_ci(probs = probs_table[,i], outcome = test_outcome)
}

# Create matrix to store c statistics + CI
cstat_matrix <- matrix(nrow = ncol(probs_table), ncol = 3)

# Calculate c-statistic + CI
for (i in 1:ncol(probs_table)){
  mat <- cbind(probs_table[,i], test_outcome)
  x <- data.frame(mat) %>% 
    filter(test_outcome == 1)
  y <- data.frame(mat) %>% 
    filter(test_outcome == 0)
  cstat_matrix[i,] <- auRoc::auc.nonpara.mw(x = x[,1], 
                                            y = y[,1])
}

# Bind all results together
results <- cbind(accuracy_vector, sensitivity_vector, specificity_vector,
                 cstat_matrix, calibration_matrix)

# Round results to 2 digits
results <- format(round(results, digits = 2), nsmall = 2)

# Name rows and columns of the results object
colnames(results) <- c("Accuracy", "Sensitivity", "Specificity", "C-statistic",
                       "C-statistic lower", "C-statistic upper", "CIL", 
                       "Lower CIL", "Upper CIL", "Calibration slope", 
                       "Lower slope", "Upper slope")
rownames(results) <- c("Penalized")


results <- data.frame(results)


# Loop over results object to get all CI's in parentheses
for (i in 1:nrow(results)){
  results$C.statistic[i] <- str_c(results$C.statistic[i],
                                  " (", 
                                  results$C.statistic.lower[i],
                                  " to ",
                                  results$C.statistic.upper[i],
                                  ")")
  results$CIL[i] <- str_c(results$CIL[i],
                          " (", 
                          results$Lower.CIL[i],
                          " to ",
                          results$Upper.CIL[i],
                          ")")
  results$Calibration.slope[i] <- str_c(results$Calibration.slope[i],
                                        " (", 
                                        results$Lower.slope[i],
                                        " to ",
                                        results$Upper.slope[i],
                                        ")")
}

# Remove old CI measure columns
results <- results %>% 
  select(!c(C.statistic.lower, C.statistic.upper, 
            Upper.CIL, Lower.CIL, 
            Upper.slope, Lower.slope))


val.prob.ci.2(p = probs_table[,1], y = baked_test$HPV, smooth="loess", lwd.smooth=2, 
              lwd.ideal=2, lty.ideal=2, xlim = c(0, 0.5),ylim = c(0, 0.5),
              legendloc = c(0.35, 0.1), statloc = c(0, 0.4))
mtext("Calibration Curve for ANY-HPV Model", side=3, line=1.3)

## Plot Decision Curve
######################

dcdata = as.data.frame(cbind(baked_test$HPV,probs_table)) 
colnames(dcdata) = c("Y", "Penalized")

pen <- decision_curve(Y~Penalized, data = dcdata, fitted.risk = TRUE, 
                      thresholds = seq(0, 0.4, by = .01), population.prevalence = 0.069, 
                      study.design = 'case-control', bootstraps = 250) 

plot_decision_curve( list(pen), 
                     curve.names = c("Recruiting with model"),
                     col = c("red"), 
                     ylim = c(-0.05, 0.1), #set ylim
                     #xlim = c(0,0.5),
                     lty = c(1), lwd = c(2), confidence.intervals = FALSE,
                     standardize = F, #plot Net benefit instead of standardized net benefit
                     legend.position = "topright",xlab="Probability Threshold (Pt)",
                     ylab = "Net Benefit", cost.benefit.axis=F) 
mtext("Decision Curve for ANY-HPV Model")



#######################################################
#######################################################
#### Model with high-risk ora HPV status as the outcome
#######################################################
#######################################################

# bit of data wrangling for that
hr_vars = df_train2 %>% 
  select(starts_with("ORX")) %>%
  mutate(across(everything(), ~ recode(.x, `1`=1L, `2`=0L, `3`=0L))) %>%
  as.matrix()
df_train2$HPV = ifelse(rowSums(hr_vars) > 0, 1, 0)

df_train_hr = df_train2 %>% 
  select(HPV, !starts_with("ORX"))


hr_vars_val = df_test2 %>% 
  select(starts_with("ORX")) %>%
  mutate(across(everything(), ~ recode(.x, `1`=1L, `2`=0L, `3`=0L))) %>%
  as.matrix()
df_test2$HPV = ifelse(rowSums(hr_vars_val) > 0, 1, 0)

df_test_hr = df_test2 %>% 
  select(HPV, !starts_with("ORX"))

## Prep and impute data.
#######################

set.seed(123456)
nhanes_recipe2 <- recipe(HPV ~ ., data = df_train_hr)

blueprint2 = nhanes_recipe2 %>%
  step_impute_knn(all_predictors()) %>%
  step_normalize(age, age_first_sex, 
                 lifetime_num_oral_sex_partner, lifetime_num_sex_partner, smk_pk_yrs)

prepare2 <- prep(blueprint2, training = df_train_hr)

baked_train2 <- bake(prepare2, new_data = df_train_hr)

## Fit standard logistic regression for hr HPV
#############################################

set.seed(2345)

standard_model2 <- lrm(HPV ~ rcs(age, 4) + rcs(age_first_sex, 4) + 
                         rcs(lifetime_num_oral_sex_partner, 4)
                       + rcs(lifetime_num_sex_partner, 4) + rcs(smk_pk_yrs, 4) + 
                         gender + race + marital_status,
                       data=baked_train2,x=T,y=T, linear.predictors = T)


## Fit and internally validate penalized logistic model for hr HPV
##################################################################

set.seed(2345)
# determine performance over range of penalties
p_trace2 <- pentrace(standard_model2, c(0,2,4,6,8,10,12,14,16,18,20,22, 24, 28, 32,40))
# fit penalized model
penalized_model2 <- update(standard_model2, penalty=p_trace2$penalty)
# Internal validation by bootstrapping of penalized model
val_penalized2 <- validate(penalized_model2, B=250)

# Look at the results
penalized_model2
val_penalized2

## Temporal validation for hr HPV
#################################

set.seed(123456)
baked_test2 <- bake(prepare2, new_data = df_test_hr) #Use existing blueprint for data prep

#Get outcome imbalance
n_minor2 <- df_train2 %>% 
  filter(HPV == 1) %>% 
  nrow()
n_major2 <- df_train2 %>% 
  filter(HPV == 0) %>% 
  nrow()
imb2 = n_minor2/n_major2


## Temporal Validation
######################

baked_test2 = as.data.frame(baked_test2)
test_outcome2 = baked_test2$HPV
ruleLP_hr	<- predict(penalized_model2, baked_test2, type="lp")
probLP_hr  	<- plogis(ruleLP_hr)


## Calibration plots and model performance
##########################################

pred_hr <- rep(0, length(probLP_hr)) # empty vector for class prediction
pred_hr[probLP_hr > imb2] = 1 #predictions based on the class imbalance threshold and not 0.5

# Creating table with all estimated probabilities
probs_table2 <- cbind(probLP_hr)

# Creating table with all predicted classes
pred_table2 <- cbind(pred_hr)

# Creating empty vectors for performance measures
accuracy_vector2 <- rep(NA, ncol(pred_table2))
sensitivity_vector2 <- rep(NA, ncol(pred_table2))
specificity_vector2 <- rep(NA, ncol(pred_table2))

# Loop over all models
for (i in 1:ncol(pred_table2)){
  accuracy_vector2[i] <- accuracy(pred_table2[,i], test_outcome2)
  sensitivity_vector2[i] <- sensitivity(pred_table2[,i], test_outcome2)
  specificity_vector2[i] <- specificity(pred_table2[,i], test_outcome2)
}

# Create matrix to store calibration measures (intercept & slope + CI)
calibration_matrix2 <- matrix(ncol = 6, nrow = ncol(probs_table2))

# Calculate calibration measures for all models
for (i in 1:ncol(probs_table2)) {
  calibration_matrix2[i,] <- calibration_ci(probs = probs_table2[,i], outcome = test_outcome2)
}

# Create matrix to store c statistics + CI
cstat_matrix2 <- matrix(nrow = ncol(probs_table2), ncol = 3)

# Calculate c-statistic + CI
for (i in 1:ncol(probs_table2)){
  mat2 <- cbind(probs_table2[,i], test_outcome2)
  x2 <- data.frame(mat2) %>% 
    filter(test_outcome2 == 1)
  y2 <- data.frame(mat2) %>% 
    filter(test_outcome2 == 0)
  cstat_matrix2[i,] <- auRoc::auc.nonpara.mw(x = x2[,1], 
                                             y = y2[,1])
}

# Bind all results together
results2 <- cbind(accuracy_vector2, sensitivity_vector2, specificity_vector2,
                  cstat_matrix2, calibration_matrix2)

# Round results to 2 digits
results2 <- format(round(results2, digits = 2), nsmall = 2)

# Name rows and columns of the results object
colnames(results2) <- c("Accuracy", "Sensitivity", "Specificity", "C-statistic",
                        "C-statistic lower", "C-statistic upper", "CIL", 
                        "Lower CIL", "Upper CIL", "Calibration slope", 
                        "Lower slope", "Upper slope")
rownames(results2) <- c("Penalized")


results2 <- data.frame(results2)


# Loop over results object to get all CI's in parentheses
for (i in 1:nrow(results2)){
  results2$C.statistic[i] <- str_c(results2$C.statistic[i],
                                   " (", 
                                   results2$C.statistic.lower[i],
                                   " to ",
                                   results2$C.statistic.upper[i],
                                   ")")
  results2$CIL[i] <- str_c(results2$CIL[i],
                           " (", 
                           results2$Lower.CIL[i],
                           " to ",
                           results2$Upper.CIL[i],
                           ")")
  results2$Calibration.slope[i] <- str_c(results2$Calibration.slope[i],
                                         " (", 
                                         results2$Lower.slope[i],
                                         " to ",
                                         results2$Upper.slope[i],
                                         ")")
}

# Remove old CI measure columns
results2 <- results2 %>% 
  select(!c(C.statistic.lower, C.statistic.upper, 
            Upper.CIL, Lower.CIL, 
            Upper.slope, Lower.slope))


val.prob.ci.2(p = probs_table2[,1], y = baked_test2$HPV, smooth="loess", lwd.smooth=2,                          
              lwd.ideal=2, lty.ideal=2,  xlim = c(0, 0.5), ylim = c(0, 0.5), 
              legendloc = c(0.35, 0.1), statloc = c(0, 0.4))
mtext("Calibration Curve for High-Risk HPV Model", side=3, line=1.3)

## Decision curve for HR-HPV
#############################

dcdata2 = as.data.frame(cbind(baked_test2$HPV,probs_table2)) 
colnames(dcdata2) = c("Y", "Penalized")

pen2 <- decision_curve(Y~Penalized, data = dcdata2, fitted.risk = TRUE, 
                       thresholds = seq(0, 0.4, by = .01), population.prevalence = 0.042, 
                       study.design = 'case-control', bootstraps = 250) 

plot_decision_curve( list(pen2), 
                     curve.names = c("Recruiting with high-risk HPV model"),
                     col = c("red"), 
                     ylim = c(-0.05, 0.1), #set ylim
                     #xlim = c(0,0.5),
                     lty = c(1), lwd = c(2), confidence.intervals = FALSE,
                     standardize = F,
                     legend.position = "topright",xlab="Probability Threshold (Pt)",ylab = "Net Benefit", cost.benefit.axis=F) 
mtext("Decision Curve for High-Risk Oral HPV Model")

