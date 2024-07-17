library(readr)
library(naniar)
library(dplyr)
 
library(glmnet)
library(ggfortify)
library(ISLR)
library(pls)
library(car)
library(boot)
 
# ______________________________________________________________________________
clinical_drug <- read.csv("Data/clinical_drug.csv")
clinical_data <- read.csv("Data/clinical_data.csv")
protein_data <- read.csv("Data/protein_data.csv")
 
length(unique(clinical_drug$bcr_patient_barcode)) # 285 VS 695
length(unique(clinical_data$bcr_patient_barcode))
length(unique(protein_data$bcr_patient_barcode))
 
# Combine the day, month, and year into a single date column for sorting
clinical_drug <- clinical_drug %>% mutate(date = as.Date(paste(year_of_form_completion,
                                                               month_of_form_completion,
                                                               day_of_form_completion, sep = "-"), "%Y-%m-%d"))
 
clinical_drug = clinical_drug %>% relocate(date, .after = bcr_patient_barcode)

df <- clinical_drug %>%
  group_by(bcr_patient_barcode) %>%
  arrange(is.na(days_to_drug_therapy_start), days_to_drug_therapy_start, date, .by_group = TRUE) %>%
  mutate(observation_count = row_number()) # Add a column for the count of observations


# Print the sorted and counted dataframe
print(df)

df <- df[df$observation_count == 1, ]

merged = merge(df, clinical_data, by = 'bcr_patient_barcode')
merged_pr = merge(merged, protein_data, by = 'bcr_patient_barcode')
 
table(merged_pr$measure_of_response)


#merged_pr$outcome = ifelse(merged_pr$measure_of_response == 'Clinical Progressive Disease' |
                             merged_pr$measure_of_response == 'Stable Disease', 0, 1)
 
#merged_pr$outcome = ifelse(merged_pr$measure_of_response == 'Clinical Progressive Disease', 0, 1)

merged_pr$outcome = ifelse(merged_pr$neoplasm_histologic_grade == 'G3', 0, 1)
 
merged_pr = merged_pr[!is.na(merged_pr$outcome), ]
merged_pr$outcome
 
colnames(merged_pr)
 
df = merged_pr %>% select(c('outcome', names(protein_data)[-1]))
df = df[ , colSums(is.na(df))==0]

x <- model.matrix(outcome ~ ., data = df)[, -1]

glmmod = glmnet(x, df$outcome, alpha=1, family="binomial")
plot(glmmod, xvar="lambda")
 
# Come scegliere lambda? Cross validazione!
set.seed(4)
lambdas = c(10^seq(10,-2,length.out = 100), 0)
cv.model_ridge = cv.glmnet(x, df$outcome, lambda = lambdas, alpha = 1, nfolds = 10)
autoplot(cv.model_ridge)
 
str(cv.model_ridge)
lambda <- cv.model_ridge$lambda.1se
 
mod <- glmnet(x, df$outcome, alpha=1, family="binomial", lambda = lambda)
summary(mod)
mod$beta@i # variables indexes
mod$beta@x # coefficients related

bestlam = cv.model_ridge$lambda.min # prima linea tratteggiata
# (lambda in corrispondenza del minor MSE in CV)
bestlam
# selezioniamo il modello corrispondente a quel lambda
indbest = which(cv.model_ridge$lambda == bestlam)
indbest
 
reglam = cv.model_ridge$lambda.1se  # seconda linea tratteggiata
# (il modello piu' regolarizzato,
# che abbia MSE CV entro una dev standard dal minimo MSE)
 
# selezioniamo il modello corrispondente a quel lambda
indreg = which(cv.model_ridge$lambda == reglam)
indreg
 
 
 
 
 
# Note alpha=1 for lasso only and can blend with ridge penalty down to
# alpha=0 ridge only.
glmmod <- glmnet(x, y=as.factor(asthma), alpha=1, family="binomial")
 
# Plot variable coefficients vs. shrinkage parameter lambda.
plot(glmmod, xvar="lambda")
 
# ______________________________________________________________________________
 
 
df <- clinical_drug %>%
  group_by(bcr_patient_barcode) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(observation_count = row_number()) # Add a column for the count of observations
 
df = df %>% relocate(observation_count, .after = bcr_patient_barcode)
 
