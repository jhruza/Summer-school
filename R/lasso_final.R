library(glmnet)
library(dplyr)
library(ISLR)
library(ggfortify)
set.seed(4)


load('Data/df_final.Rda')
protein_data <- read.csv("Data/protein_data.csv")

df <- df_final %>% select(c('outcome', names(protein_data)[-1]))
df <- df[!is.na(df$outcome), ]
df <- df[ , colSums(is.na(df)) == 0]

#create covariate matrix
x <- model.matrix(outcome ~ ., data = df)[ ,-1]

glmmod = glmnet(x, df$outcome, alpha=1, family="binomial")
plot(glmmod, xvar="lambda")

# Cross validazione!

lambdas = c(10^seq(10, -2, length.out = 100), 0)
cv.model_lasso = cv.glmnet(x, df$outcome, lambda = lambdas, alpha = 1, nfolds = 10)
autoplot(cv.model_lasso)

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
cv.alasso <- cv.glmnet(x, df$outcome, alpha=0.5, family="binomial", penalty.factor = w.r)
cbind(coef(mod, s="lambda.min"),
           coef(cv.alasso, s="lambda.min"))

alasso <- glmnet(x, df$outcome, alpha=0.5, family="binomial", lambda = cv.alasso$lambda.min, penalty.factor = w.r)
alasso$beta@i # variables indexes
alasso$beta@x # coefficients values
alasso$beta@Dimnames[[1]][alasso$beta@i] #names


###protein prognostic score PPS
formula_string <- paste("outcome ~", paste(alasso$beta@Dimnames[[1]][alasso$beta@i], collapse = " + "))
formula <- as.formula(formula_string)
PSS.model  <- glm(formula_string, data = df, family = "binomial")
PPS <- predict(PSS.model, df[ ,alasso$beta@Dimnames[[1]][alasso$beta@i]], type= "response" )
hist(PPS, breaks=30)
