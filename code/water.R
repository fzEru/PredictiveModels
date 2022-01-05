library(tidyverse)
library(mice)
library(naniar)
library(leaps)
library(glmnet)
library(MASS)
library(tree)
library(rpart)

waterData <- read.csv("Drinking_water.csv")

### tidying waterData
## cleaning column names
waterData <- waterData %>%
  dplyr::rename(
    hardness = Hardness,
    solids = Solids,
    chloramines = Chloramines,
    sulfate = Sulfate,
    conductivity = Conductivity,
    organic.carbon = Organic_carbon,
    trihalomethanes = Trihalomethanes,
    turbidity = Turbidity,
    potability = Potability) %>%
  mutate(
    potability = as.factor(potability)) %>%
  dplyr::select(c(-X, -Carcinogenics, -medical_waste))

## handling NAs
sapply(waterData, function(x) sum(is.na(x)))

## MICE imputation
water.pred.mat <- quickpred(waterData)
water.impute <- mice(waterData, m = 5, meth = "norm", seed = 5, predictorMatrix = water.pred.mat)
waterData <- complete(water.impute, 1)

## train/test definitions
attach(waterData)
train <- sample(1:nrow(waterData), 2457)
train.water <- waterData[train, ]
test.water <- waterData[-train, ]

## classification tree
set.seed(25)
tree.water <- rpart(potability ~ ., train.water, control = rpart.control('minsplit' = 1))
summary(tree.water)
tree.water
plot(tree.water)
text(tree.water, pretty = 0)
pred.water <- predict(tree.water, test.water, type = 'class')
table.water <- table(test.water$potability, pred.water)
table.water
mean(!test.water$potability == pred.water)

## best subset selection
regfit.fwd <- regsubsets(potability ~ ., data = waterData, nvmax = 9, method = "forward")
summary(regfit.fwd)

## logistic regression
dim(test.water)
test.potable <- potability[-train]
glm.fit <- glm(potability ~ solids + organic.carbon + chloramines, 
               data = waterData, family = binomial, subset = train)
glm.probs <- predict(glm.fit, test.water, type = "response")
glm.pred <- rep("0", 819)
glm.pred[glm.probs > 0.4] = "1"
table(glm.pred, test.potable)
mean(glm.pred == test.potable)
plot(glm.pred)


## LDA
lda.fit <- lda(potability ~ solids + organic.carbon + chloramines, 
               data = waterData, subset = train)
lda.fit
lda.pred <- predict(lda.fit, test.water)
names(lda.pred)
lda.class <- lda.pred$class
table(lda.class,test.potable)
mean(lda.class == test.potable)
plot(lda.class)

sum(lda.pred$posterior[,1] >= 0.6)
sum(lda.pred$posterior[,1] < 0.6)
lda.pred$posterior[1:20,1]
lda.class[1:50]

# QDA
qda.fit <- qda(potability ~ solids + organic.carbon + chloramines, 
               data = waterData, subset = train)
qda.fit
qda.class <- predict(qda.fit, test.water)$class
table(qda.class, test.potable)
mean(qda.class == test.potable)
plot(qda.class)

# KNN
library(class)
train.X <- cbind(solids ,organic.carbon, chloramines)[train ,]
test.X <- cbind(solids ,organic.carbon, chloramines)[-train ,]
train.potable <- potability[train]

set.seed(25)
knn.pred <- knn(train.X, test.X, train.potable, k = 1)
table(knn.pred, test.potable)
mean(knn.pred == test.potable)
plot(knn.pred)

