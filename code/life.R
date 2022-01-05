library(devtools)
library(tidyverse)
library(mice)
library(naniar)
library(leaps)
library(glmnet)
library(MASS)
library(tree)
library(rpart)
library(scales)
install_github("vqv/ggbiplot")
library(ggbiplot)

lifeData <- read.csv("Life Expectancy Data.csv")
lifeData$Status <- recode(lifeData$Status, Developing = 0, Developed = 1)
### tidying lifeData
## cleaning column names
lifeData <- lifeData %>%
  dplyr::select(-Country) %>%
  dplyr::rename(
    year = Year,
    status = Status,
    life.exp = Life.expectancy,
    adult.mortality = Adult.Mortality,
    alcohol = Alcohol,
    hep.B = Hepatitis.B,
    measles = Measles,
    bmi = BMI,
    polio = Polio,
    total.expenditure = Total.expenditure,
    diphtheria = Diphtheria,
    hiv.aids = HIV.AIDS,
    gdp = GDP,
    population = Population,
    thinness.10to19 = thinness..1.19.years,
    thinness.5to9 = thinness.5.9.years,
    income.composition.resources = Income.composition.of.resources,
    schooling = Schooling) %>%
  mutate(
    status = as.factor(status),
    year = as.numeric(year),
    infant.deaths = as.numeric(infant.deaths),
    measles = as.numeric(measles),
    under.five.deaths = as.numeric(under.five.deaths)
  )
## checking NAs
sapply(lifeData, function(x) sum(is.na(x)))

## MICE imputation
set.seed(25)
life.pred.mat <- quickpred(lifeData)
life.impute <- mice(lifeData, m = 5, meth = "cart", seed = 5, predictorMatrix = life.pred.mat)
lifeData2 <- complete(life.impute, 1)
lifeData <- lifeData2
sapply(lifeData, function(x) sum(is.na(x)))

## data diligence
pairs(lifeData)
summary(lifeData)
sapply(lifeData, class)

## regression tree
set.seed(25)
train.tree <- sample(1:nrow(lifeData), nrow(lifeData)/2)
tree.life <- tree(life.exp ~ ., lifeData, subset = train.tree)
summary(tree.life)
plot(tree.life)
text(tree.life, pretty = 0)
cv.life <- cv.tree(tree.life)
plot(cv.life$size, cv.life$dev, type = 'b')
prune.life <- prune.tree(tree.life, best = 6)
plot(prune.life)
text(prune.life ,pretty = 0)
yhat <- predict(tree.life, newdata = lifeData[-train.tree,])
life.test <- lifeData[-train.tree,"life.exp"]
plot(yhat, life.test)
abline(0,1)
mean((yhat - life.test)^2)

## simple linear regression
# best subset selection for simple linear regression
regfit.fwd <- regsubsets(life.exp ~ ., data = lifeData, nvmax = 19, method = "forward")
summary(regfit.fwd)
# According to the forward selection algorithm, the variable 'schooling' is the 
# best predictor to use in a simple linear regression 
lm.fit <- lm(life.exp ~ schooling, data = lifeData)
lm.fit
attach(lifeData) 
plot(schooling, life.exp)
abline(lm.fit, lwd = 3, col = 'red')
par(mfrow = c(2, 2))
plot(lm.fit)
# check non linearity 
plot(predict(lm.fit), residuals(lm.fit))
plot(predict(lm.fit), rstudent(lm.fit))
# leverage stats 
plot(hatvalues(lm.fit))
which.max(hatvalues(lm.fit))
lm.pred <- predict(lm.fit, newdata = lifeData[-train.tree,])
mean((lm.pred - life.test)^2)

## multiple linear regression
# best subset selection for multiple linear regression
regfit.full <- regsubsets(life.exp ~ ., lifeData, nvmax = 19)
summary(regfit.full)
reg.summary <- summary(regfit.full)
names(reg.summary)
reg.summary$rsq
par(mfrow=c(2,2))
plot(reg.summary$rss, xlab = "Number of Variables", ylab="RSS", type="l")
which.min(reg.summary$rss)
plot(reg.summary$adjr2, xlab = "Number of Variables", ylab="Adjusted RSq", type="l")
which.max(reg.summary$adjr2)
points(15, reg.summary$adjr2[15], col="red", cex = 2, pch = 20)
plot(reg.summary$cp, xlab = "Number of Variables", ylab = "Cp", type = "l")
which.min(reg.summary$cp)
points(14, reg.summary$cp[14], col = "red", cex = 2, pch = 20)
which.min(reg.summary$bic )
plot(reg.summary$bic, xlab="Number of Variables", ylab="BIC", type = "l")
points(12, reg.summary$bic[12], col = "red", cex = 2, pch = 20)
plot(regfit.full,scale="r2")
plot(regfit.full,scale="adjr2") 
plot(regfit.full,scale="Cp")
plot(regfit.full,scale="bic")
coef(regfit.full, 15)
# best model is with 15 predictors

# Choosing Among Models Using the Validation Set Approach and Cross-Validation
set.seed(25)
train <- sample(c(TRUE,FALSE), nrow(lifeData), rep = TRUE)
test <- (!train)
regfit.best <- regsubsets(life.exp ~ ., data = lifeData[train,], nvmax = 19)
test.mat <- model.matrix(life.exp ~ ., data = lifeData[test,])

k = 10
set.seed(25)
folds <- sample(1:k, nrow(lifeData), replace = TRUE)
cv.errors <- matrix(NA, k, 19, dimnames = list(NULL, paste(1:19)))

## ridge regression
x <- model.matrix(life.exp ~ ., lifeData)[,-1] 
y <- lifeData$life.exp
grid <- 10^seq(10, -2, length = 100)
ridge.mod <- glmnet(x, y, alpha = 0, lambda = grid)
dim(coef(ridge.mod))
ridge.mod$lambda[50]
coef(ridge.mod)[,50]
sqrt(sum(coef(ridge.mod)[-1,50]^2))
ridge.mod$lambda[60]
coef(ridge.mod)[,60]
sqrt(sum(coef(ridge.mod)[-1,60]^2))
predict(ridge.mod, s = 50, type = "coefficients")[1:20,]

set.seed(25)
train <- sample(1:nrow(x), nrow(x)/2)
test <- (-train)
y.test <- y[test]
ridge.mod <- glmnet(x[train,], y[train], alpha = 0, lambda = grid)
ridge.pred <- predict(ridge.mod, s = 4, newx = x[test,])
mean((ridge.pred-y.test)^2)
mean((mean(y[train])-y.test)^2)
ridge.pred <-predict(ridge.mod,s=1e10,newx=x[test,])
mean((ridge.pred-y.test)^2)
ridge.pred <- predict(ridge.mod,s=0,newx=x[test,],exact=T,x=x[train,],y=y[train])
mean((ridge.pred-y.test)^2)
lm(y ~ x, subset = train)
predict(ridge.mod, s = 0, exact = T, type = "coefficients", x = x[train,], y = y[train])[1:20,]
set.seed(25)
cv.out <- cv.glmnet(x[train,], y[train], alpha = 0)
plot(cv.out)
bestlam <- cv.out$lambda.min
bestlam
ridge.pred <- predict(ridge.mod,s = bestlam,newx = x[test,])
mean((ridge.pred-y.test)^2)
out <- glmnet(x, y, alpha = 0)
predict(out,type = "coefficients",s = bestlam)[1:20,]

## lasso 
lasso.mod <- glmnet(x[train,],y[train],alpha = 1,lambda = grid)
plot(lasso.mod)
set.seed(25)
cv.out <- cv.glmnet(x[train,], y[train], alpha = 1)
plot(cv.out)
bestlam <- cv.out$lambda.min
lasso.pred <- predict(lasso.mod, s = bestlam, newx = x[test,])
mean((lasso.pred-y.test)^2)
out <- glmnet(x, y, alpha = 1, lambda = grid)
lasso.coef <- predict(out,type = "coefficients",s = bestlam)[1:20,]
lasso.coef
lasso.coef[lasso.coef!=0]

## PCR 
library(pls)
set.seed(24)
pcr.fit <- pcr(life.exp ~., data = lifeData , scale = TRUE, validation = "CV")
summary(pcr.fit)
validationplot(pcr.fit, val.type = "MSEP")
set.seed(25)
pcr.fit <- pcr(life.exp ~., data= lifeData, subset = train, scale = TRUE, validation = "CV")
validationplot(pcr.fit, val.type = "MSEP")
pcr.pred <- predict(pcr.fit, x[test,], ncomp = 7)
mean((pcr.pred-y.test)^2)
pcr.fit <- pcr(y~x, scale = TRUE, ncomp = 7)
summary(pcr.fit)
life.pca <- prcomp(lifeData[, -2], center = TRUE, scale. = TRUE)
summary(life.pca)
ggbiplot(life.pca, alpha = 0, scale = 0)
ggbiplot(life.pca, scale = 0, labels = life.exp)

## PLS 
set.seed(25)
pls.fit <- plsr(life.exp ~., data = lifeData, subset = train, scale = TRUE, validation = "CV")
summary(pls.fit)
validationplot(pls.fit, val.type = "MSEP")
pls.pred <- predict(pls.fit, x[test,], ncomp = 2)
mean((pls.pred-y.test)^2)
pls.fit <- plsr(life.exp ~., data = lifeData, scale = TRUE, ncomp = 2)
summary(pls.fit)

## bagging
library(randomForest)
set.seed(25)
bag.life <- randomForest(life.exp ~ .,data = lifeData,subset = train, mtry = 19, importance = TRUE)
bag.life
yhat.bag <- predict(bag.life, newdata = lifeData[-train,])
plot(yhat.bag, life.test)
abline(0,1)
mean((yhat.bag - life.test)^2)
bag.life <- randomForest(life.exp ~ .,data = lifeData, subset = train, mtry = 19, ntree = 25)
yhat.bag <- predict(bag.life, newdata = lifeData[-train,])
mean((yhat.bag - life.test)^2)
set.seed(25)
rf.life <- randomForest(life.exp ~. ,data = lifeData, subset = train ,mtry = 6, importance = TRUE)
yhat.rf <- predict(rf.life, newdata = lifeData[-train,])
mean((yhat.rf - life.test)^2)
importance(rf.life)
varImpPlot(rf.life)

## boosting 
library(gbm)
set.seed(25)
boost.life <- gbm(life.exp ~., data = lifeData[train,],distribution = "gaussian", n.trees = 5000, interaction.depth = 4)
summary(boost.life)
par(mfrow = c(1,2))
plot(boost.life, i = "hiv.aids")
plot(boost.life, i = "income.composition.resources")
yhat.boost <- predict(boost.life, newdata=  lifeData[-train,], n.trees = 5000)
mean((yhat.boost - life.test)^2)
boost.life <- gbm(life.exp ~., data = lifeData[train,], distribution="gaussian", n.trees=5000, interaction.depth=4, shrinkage=0.2, verbose=F)
yhat.boost <-predict(boost.life,newdata = lifeData[-train,], n.trees = 5000)
mean((yhat.boost - life.test)^2)
