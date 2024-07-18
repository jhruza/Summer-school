library(glmnet)
library(ggplot2)
library(dplyr)
library(ISLR)
library(ggfortify)
set.seed(7)


load('Data/df_final.Rda')
protein_data <- read.csv("Data/protein_data.csv")

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
variable_indices_list <- vector("list", 100)

# Run the loop 100 times
for (i in 1:100) {
  # Perform cross-validation with glmnet
  cv.alasso <- cv.glmnet(x, df$outcome, alpha = 0.5, family = "binomial", penalty.factor = w.r)
  
  # Fit the model using the lambda value from cross-validation
  alasso <- glmnet(x, df$outcome, alpha = 0.5, family = "binomial", lambda = cv.alasso$lambda.min, penalty.factor = w.r)
  
  # Store the variable indices
  variable_indices_list[[i]] <- alasso$beta@i
}

# Initialize a vector to store the containment count for each list
containment_count <- numeric(100)

# Compare each set of indices with all other sets
for (i in 1:100) {
  for (j in 1:100) {
    if (i != j) {
      if (all(variable_indices_list[[i]] %in% variable_indices_list[[j]])) {
        containment_count[i] <- containment_count[i] + 1
      }
    }
  }
}

# Create a data frame for plotting
containment_data <- data.frame(
  List = 1:100,
  ContainmentCount = containment_count
)
containment_count>50
# Plot the results
ggplot(containment_data, aes(x = List, y = ContainmentCount)) +
  geom_bar(stat = "identity") +
  labs(title = "Containment of Variable Indices in Other Lists",
       x = "List Index",
       y = "Containment Count") +
  theme_minimal()

#threshold is set choose that appear half of the time
protein_selected <- alasso$beta@Dimnames[[1]][index]


###protein prognostic score PPS
formula_string <- paste("outcome ~", paste(protein_selected, collapse = " + "))
formula <- as.formula(formula_string)
PSS.model  <- glm(formula_string, data = df, family = "binomial")
PPS <- predict(PSS.model, df[ ,protein_selected], type= "response" )
hist(PPS, breaks=30)








# CODE FOR PPS, CATE and so on

library(WeightIt)
library(marginaleffects)

# we assume 
# tr = post_operative something
# confounders = X1 + X2 + .... to be replaced with correct ones
# data = data
# outcome = outcome

data$PPS = PPS

W = weightit(tr ~ confounders, data = data, method = "glm", estimand = "ATE")

fit = glm_weightit(outcome ~ A * (confounders) + A * PPS, data = data, weighit = W, family = binomial)
