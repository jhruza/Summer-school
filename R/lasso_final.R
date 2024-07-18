library(glmnet)
library(ggplot2)
library(dplyr)
library(ISLR)
library(ggfortify)
library(fastDummies)
library(cobalt)
set.seed(7)


load('Data/df_final.Rda')
protein_data <- read.csv("Data/protein_data.csv")
confounders <- c("postoperative_rx_tx", "gender", "white",
                "age_at_initial_pathologic_diagnosis",
                "year_of_initial_pathologic_diagnosis",
                "seizure", "neoplasm_histologic_grade", "icd10",
                "tumor_location", "histological_type")


protein_names <- names(protein_data)[-1]
columns_to_remove <- colnames(df_final)[sapply(df_final, function(col) all(is.na(col)))]
protein_names <- setdiff(protein_names, columns_to_remove)


df_final_cleaned <- df_final %>%
      filter(!if_any(all_of(c(confounders, "outcome")), is.na))



df <- df_final %>% select(c('outcome', names(protein_data)[-1]))
df <- df[!is.na(df$outcome), ]
df <- df[ , colSums(is.na(df)) == 0]


#create covariate matrix
x <- model.matrix(outcome ~ ., data = df)[ ,-1]

glmmod = glmnet(x, df$outcome, alpha=1, family="binomial")
#plot(glmmod, xvar="lambda")

# Cross validazione!

lambdas = c(10^seq(10, -2, length.out = 100), 0)
cv.model_lasso = cv.glmnet(x, df$outcome, lambda = lambdas, alpha = 1, nfolds = 10)
#autoplot(cv.model_lasso)

str(cv.model_lasso)
lambda <- cv.model_lasso$lambda.min

mod <- glmnet(x, df$outcome, alpha=1, family="binomial", lambda = lambda)
summary(mod)
mod$beta@i # variables indexes
mod$beta@x # coefficients values

####################
# addaptive lasso
# Ridge weights with gamma = 1
g = 1
modelr <- cv.glmnet(x, df$outcome, alpha = 0)
coefr <- as.matrix(coef(modelr, s = modelr$lambda.min))
w.r <- 1/(abs(coefr[-1,]))^g

#adaptive elastic net
set.seed(24)
cv.alasso <- cv.glmnet(x, df$outcome, alpha=0.5, family="binomial", penalty.factor = w.r)
# cbind(coef(mod, s="lambda.min"),
#            coef(cv.alasso, s="lambda.min"))

alasso <- glmnet(x, df$outcome, alpha=0.5, family="binomial", lambda = cv.alasso$lambda.min, penalty.factor = w.r)
alasso$beta@i # variables indexes
alasso$beta@x # coefficients values
alasso$beta@Dimnames[[1]][alasso$beta@i] #names



# Initialize a list to store the variable indices for each iteration
variable_indices_list <- vector("list", 1000)

# Run the loop 100 times
for (i in 1:1000) {
  # Perform cross-validation with glmnet
  cv.alasso <- cv.glmnet(x, df$outcome, alpha = 0.5, family = "binomial", penalty.factor = w.r)
  
  # Fit the model using the lambda value from cross-validation
  alasso <- glmnet(x, df$outcome, alpha = 0.5, family = "binomial", lambda = cv.alasso$lambda.min, penalty.factor = w.r)
  
  # Store the variable indices
  variable_indices_list[[i]] <- alasso$beta@i
}

# Initialize a vector to store the containment count for each list
containment_count <- numeric(1000)

# Compare each set of indices with all other sets
for (i in 1:1000) {
  for (j in 1:1000) {
    if (i != j) {
      if (all(variable_indices_list[[i]] %in% variable_indices_list[[j]])) {
        containment_count[i] <- containment_count[i] + 1
      }
    }
  }
}

# Create a data frame for plotting
containment_data <- data.frame(
  List = 1:1000,
  ContainmentCount = containment_count
)
# Plot the results
ggplot(containment_data, aes(x = List, y = ContainmentCount)) +
  geom_bar(stat = "identity") +
  labs(title = "Containment of Variable Indices in Other Lists",
       x = "List Index",
       y = "Containment Count") +
  theme_minimal()

#threshold is set choose that appears half of the time
index <- unique(unlist(variable_indices_list[containment_count>500]))
protein_selected <- alasso$beta@Dimnames[[1]][index]
protein_selected <- c("BCLXL", "CIAP", "CLAUDIN7", "Histone.H3", "RICTOR", "SHC_pY317", "TRIM25","XPA")


###protein prognostic score PPS
formula_string <- paste("outcome ~", paste(protein_selected, collapse = " + "))
formula <- as.formula(formula_string)
PSS.model  <- glm(formula_string, data = df, family = "binomial")
PPS <- predict(PSS.model, df_final_cleaned[ ,protein_selected], type= "response")
hist(PPS, breaks=30)




subset_df = df_final_cleaned[, confounders]

subset_df$tumor_location[subset_df$tumor_location == "Supratentorial, Occipital Lobe"] <- "Supratentorial_Occipital_Lobe"
subset_df$tumor_location[subset_df$tumor_location == "Supratentorial, Parietal Lobe"] <- "Supratentorial_Parietal_Lobe"
subset_df$tumor_location[subset_df$tumor_location == "Supratentorial, Temporal Lobe"] <- "Supratentorial_Temporal_Lobe"
subset_df$tumor_location[subset_df$tumor_location == "Supratentorial, Frontal Lobe"] <- "Supratentorial_Frontal_Lobe"

str(subset_df)

to_be_transformed = c("tumor_location", "histological_type", "icd10")

# Apply model.matrix to convert factors to dummy variables
dummy_matrix <- model.matrix(~ ., data = subset_df[, to_be_transformed])[,-1]

# Convert the result to a dataframe for easier handling if needed
dummy_df <- as.data.frame(dummy_matrix)

results <- cbind(df_final_cleaned[, setdiff(confounders, to_be_transformed)], dummy_df)
colnames(results)

set.cobalt.options(binary = "std")


formula_str <- paste("postoperative_rx_tx ~", paste(colnames(results)[-1], collapse = " + "))
formula_obj <- as.formula(formula_str)


# Load the love.plot package if not already loaded
if (!requireNamespace("love.plot", quietly = TRUE)) {
  install.packages("love.plot")
  library(love.plot)
}

W = weightit(formula_obj, data = results, 
            method = "glm", estimand = "ATE", stabilize = TRUE)

# Use love.plot with the formula object
love.plot(W, data = results, abs = TRUE, thresholds = 0.1)



# CODE FOR PPS, CATE and so on

library(WeightIt)
library(marginaleffects)

#Logstic regression with postoperative radiation as outcome (PROPENSITY SCORE)

df_final_cleaned$PPS = PPS

fit1 <- glm(postoperative_rx_tx ~ gender + white + age_at_initial_pathologic_diagnosis + year_of_initial_pathologic_diagnosis + 
              seizure + neoplasm_histologic_grade + icd10  + 
               tumor_location + histological_type, family=binomial(), 
               data = df_final_cleaned)
summary(fit1)

df_final_cleaned$PS <- predict(fit1, type = "response")


W = weightit(postoperative_rx_tx ~ gender + white + age_at_initial_pathologic_diagnosis + year_of_initial_pathologic_diagnosis + 
              seizure + neoplasm_histologic_grade + icd10  + 
               tumor_location + histological_type, data = df_final_cleaned, 
            method = "glm", estimand = "ATE", stabilize = TRUE)

fit = glm_weightit(outcome ~ postoperative_rx_tx * (gender + white + age_at_initial_pathologic_diagnosis + year_of_initial_pathologic_diagnosis + 
              seizure + neoplasm_histologic_grade + icd10  + 
               tumor_location + histological_type) + 
               postoperative_rx_tx * PPS, 
                  data = df_final_cleaned, weighit = W, family = binomial)
summary(fit)


# _____________

# OTHER_DX
table(df_final_cleaned$other_dx)

ggplot(df_final_cleaned, aes(x = other_dx, y = PPS, fill = other_dx)) +
  geom_boxplot() +
  labs(title = "Boxplot of Values Stratified by Group",
       x = "Group",
       y = "Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Set3")



# _____________


# tissue_source_site
table(df_final_cleaned$tissue_source_site)

ggplot(df_final_cleaned, aes(x = tissue_source_site, y = PPS, fill = tissue_source_site)) +
  geom_boxplot() +
  labs(title = "Boxplot of Values Stratified by Group",
       x = "Group",
       y = "Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),  # Title font size and style
    axis.title.x = element_text(size = 50),  # X-axis label font size
    axis.title.y = element_text(size = 14),  # Y-axis label font size
    axis.text.x = element_text(size = 50, angle = 45, hjust = 1),  # X-axis tick label font size and rotation
    axis.text.y = element_text(size = 12),  # Y-axis tick label font size
    legend.title = element_text(size = 14, face = "bold"),  # Legend title font size and style
    legend.text = element_text(size = 12),  # Legend text font size,
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Set3")


# Compute the analysis of variance
res.aov <- aov(PPS ~ tissue_source_site, data = df_final_cleaned)
# Summary of the analysis
summary(res.aov)


# _____________

table(df_final_cleaned$tissue_prospective_collection_indicator)

ggplot(df_final_cleaned, aes(x = tissue_prospective_collection_indicator, y = PPS, fill = tissue_prospective_collection_indicator)) +
  geom_boxplot() +
  labs(title = "Boxplot of Values Stratified by Group",
       x = "Group",
       y = "Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),  # Title font size and style
    axis.title.x = element_text(size = 50),  # X-axis label font size
    axis.title.y = element_text(size = 14),  # Y-axis label font size
    axis.text.x = element_text(size = 50, angle = 45, hjust = 1),  # X-axis tick label font size and rotation
    axis.text.y = element_text(size = 12),  # Y-axis tick label font size
    legend.title = element_text(size = 14, face = "bold"),  # Legend title font size and style
    legend.text = element_text(size = 12),  # Legend text font size,
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Set3")


# Compute the analysis of variance
res.aov <- aov(PPS ~ tissue_prospective_collection_indicator, data = df_final_cleaned)
# Summary of the analysis
summary(res.aov)
