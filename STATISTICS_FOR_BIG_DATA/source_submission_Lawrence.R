
## Linear regression####

library(leaps)
library(ISLR)
library(corrplot)
library(plotmo)

# Explore Dataset
#fix(College)

names(College)
dim(College)
sum(is.na(College))


C=cor(College[-1])
corrplot(C, method = "square")

# a)  ----------------------  split data

set.seed(1234)
x = model.matrix(Apps~.,College)[-1,]
y = College$Apps
dims(x)
colnames(x)

train = sample(1:nrow(x), nrow(x)/2)
test = (-train)
y.test = y[test]


# Models

# b) Linear Model

regfit.full = regsubsets(Apps~.,data=College, nvmax=17)
reg.summary <- summary(regfit.full)
reg.summary$rsq
reg.summary$adjr2

par(mfrow=c(2,2))

plot(reg.summary$rss, xlab="Number of Variables", ylab="RSS", type="l")
rss.min <- which.min(reg.summary$rss)
rss.min
points(rss.min, reg.summary$rss[rss.min], col="red", cex=2, pch=20)

plot(reg.summary$adjr2, xlab="Number of Variables", ylab="Adjusted", type="l")
adjr2.max <- which.max(reg.summary$adjr2)
points(adjr2.max, reg.summary$adjr2[adjr2.max], col="red", cex=2, pch=20)
adjr2.max

plot(reg.summary$cp, xlab="Number of Variables", ylab="Cp", type="l")
min.cp.vars <- which.min(reg.summary$cp)
points(min.cp.vars, reg.summary$cp[min.cp.vars], col='red', cex=2, pch=20)
min.cp.vars

plot(reg.summary$bic, xlab="Number Of Variables", ylab='bic', type='l')
min.bic.vars = which.min(reg.summary$bic)
points(min.bic.vars, reg.summary$bic[min.bic.vars], col='red', cex=2, pch=20)
min.bic.vars

plot(regfit.full, scale="r2")
plot(regfit.full, scale="adjr2")
plot(regfit.full, scale="Cp")
plot(regfit.full, scale="bic")


## Forward and Backward Stepwise Selection
regfit.fwd = regsubsets(Apps~.,data=College,nvmax=17, method='forward')
summary(regfit.fwd)
regfit.bwd = regsubsets(Apps~.,data=College,nvmax=17, method='backward')
summary(regfit.bwd)


##
regfit.best = regsubsets(Apps~.,data=College[train,], nvmax=17)
test.mat = model.matrix(Apps~.,data=College[test,])
val.errors = rep(NA, 17)

for (i in 1:17){
coefi = coef(regfit.best, id=i)
pred = test.mat[,names(coefi)]%*%coefi
val.errors[i] = mean((College$Apps[test]-pred)^2)
}
#######

predict.regsubsets = function(object, newdata, id,...){
form = as.formula(object$call[[2]])
mat = model.matrix(form, newdata)
coefi=coef(object, id=id)
xvars = names(coefi)
mat[,xvars]%*%coefi
}

coef(regfit.best,13)
coef(regfit.best,5)

# CROSS VALIDATION
k = 10
set.seed(10)
folds = sample(1:k, nrow(College), replace=T)
cv.errors = matrix(NA, k, 17, dimnames=list(NULL, paste(1:17)))
for (j in 1:k){
best.fit=regsubsets(Apps~.,data=College[folds!=j,], nvmax=17)
for(i in 1:17){
pred = predict(best.fit, College[folds==j,], id=i)
cv.errors[j,i] = mean((College$Apps[folds==j]-pred)^2)
}
}
cv.errors
mean.cv.errors = apply(cv.errors, 2, mean)
mean.cv.errors
par(mfrow = c(1,1))
plot(mean.cv.errors, type='b')

# Fit on full date to inspect coefficients
reg.best =regsubsets(Apps~.,College, nvmax=17)
coef(reg.best, 11)




####
#### Ridge Regression
library(glmnet)

grid = 10^seq(10, -2, length=100)
ridge.mod<-glmnet(x[train,], y[train], alpha=0, lambda=grid, thresh=1e-12, standardize = T)
summary(ridge.mod)
plot(ridge.mod)
set.seed(100)
cv.out = cv.glmnet(x[train,], y[train], alpha=0)
plot(cv.out)
bestlam = cv.out$lambda.min
bestlam
ridge.pred  = predict(ridge.mod, s=bestlam, newx=x[test,])
mean((ridge.pred-y.test)^2)
out =glmnet(x,y,alpha=0)
ridge.ceof=predict(out, type="coefficients", s=bestlam)[1:18,]
ridge.ceof


################ Lasso #####################
lasso.mod <- glmnet(x[train,], y[train], alpha=1, lambda=grid, standardize=T)
plot(lasso.mod)
set.seed(100)
cv.out = cv.glmnet(x[train,], y[train], alpha=1)
plot(cv.out)
bestlam = cv.out$lambda.min
lasso.pred = predict(lasso.mod, s=bestlam, newx=x[test,])
mean((lasso.pred - y.test)^2)
out =glmnet(x, y, alpha=1, lambda=grid)
lasso.coef = predict(out, type="coefficients", s=bestlam)[1:18]
lasso.coef
lasso.coef[lasso.coef != 0]



