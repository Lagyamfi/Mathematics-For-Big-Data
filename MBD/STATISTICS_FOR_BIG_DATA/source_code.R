## MATHEMATICS FOR BIG DATA - FUNCTIONAL DATA ANALYSIS 
## ASSIGNMENT
## AUTHOR: LAWRENCE ADU-GYAMFI
## DATE: 10/06/2019


require(fda)
load('medfly.Rdata') #if working directory is set to source file???s.Otherwise put  path.
lifetime = as.numeric(medfly$lifetime)
medbasis = create.bspline.basis(c(0,25),norder=4,breaks=0:25)
lambdas = exp(-10:10)
gcvs = rep(0,length(lambdas))
for(i in 1:length(lambdas)){
  mPar = fdPar(medbasis,int2Lfd(2),lambdas[i])
  gcvs[i] = mean(smooth.basis(0:25,medfly$eggcount,mPar)$gcv)
}
best = which.min(gcvs)
lambda = lambdas[best]
mPar = fdPar(medbasis,int2Lfd(2),lambda)
medfd = smooth.basis(0:25,medfly$eggcount,mPar)

cat("The best lambda is: ", lambda)


# Question 1
plot(medfd )


# Question 2

# calculate the Principal components
med.pca = pca.fd(medfd$fd,nharm=10)
med.var = med.pca$varprop

# view variability explained by each of the PCs
data.table(PC=seq(length(med.var)), Variability_Explained=round(med.var, digits=3), CUMSUM=round(cumsum(med.var), digits=3) )

# Plot the variability explained with the PCs
plot(med.var,type='b')

# view the Principal components
plot(med.pca$harmonics[1:3], legend())

# plot the inidividual PCs
par(mfrow=c(3,1))
plot(med.pca,harm=1:3)

# Question 3

xfdlist = list(rep(1,50),medfd$fd)
xfdlist

# create basic functional parameter object
cPar = fdPar(fdobj=create.constant.basis(),
             Lfdobj=0,
             lambda=0)

# Perform leave one out cross validation to get best smoothing parameter
cross.vals = rep(0, length(lambdas))
for (i in 1:length(lambdas)){
  aux.par=fdPar(fdobj=medbasis, int2Lfd(m=2), lambdas[i])
  cross.vals[i] = fRegress(lifetime,xfdlist,list(cPar,aux.par))$OCV
}

best_lambda = lambdas[which.min(cross.vals)]

# return final functional parameter object for the regression
aux.par = fdPar(fdobj=medbasis, Lfdobj=int2Lfd(m=2),
                best_lambda)

med.reg = fRegress(y=lifetime, xfdlist=xfdlist, betalist=list(cPar, aux.par))

# view predictions against actual values
plot(med.reg$yhatfdobj,lifetime,xlab='predictions',ylab='actuals', main="Predictions vs. Actuals", pty="s", pch=5)
abline(c(0,1), col="red", lwd=2)


# view coefficient function
med.var = sum( (lifetime-med.reg$yhatfdobj)^2 )/(50-med.reg$df)

med.std = fRegress.stderr(med.reg,NULL,med.var*diag(rep(1,50)))

med.coeffun = med.reg$betaestlist[[2]]

plot(med.coeffun$fd,ylim=c(-5,5))
lines(med.reg$betaestlist[[2]]$fd+2*sqrt(med.std$betastderrlist[[2]]),lty=2, col="red")
lines(med.reg$betaestlist[[2]]$fd-2*sqrt(med.std$betastderrlist[[2]]),lty=2, col="red")


# Question 4

# Verify significance of model by permutation testing
med.permtest = Fperm.fd(yfdPar=lifetime, xfdlist=xfdlist,
                        betalist=list(cPar, aux.par), q=0.95)

summary(med.permtest)
med.permtest$pval


# Calculate R_squared of Model
med.rsq = 1 - med.var / var(lifetime)

cat("R_squared of regression model is: ", med.rsq)

# Question 5

# get scores for all the principal components calculated previously
pca.scores <- med.pca$scores

#create linear model
med.pca.mod <- lm(lifetime ~ pca.scores)

summary(med.pca.mod)

#compare results with functional linear regression model
par(mfrow=c(2,1))
plot(med.pca.mod$fitted, lifetime)
abline(0,1)

plot(med.reg$yhatfdobj,lifetime,xlab='predictions',ylab='actuals', main="Predictions vs. Actuals", pty="s", pch=5)
abline(c(0,1), col="red", lwd=2)


# reconstruct coefficient function

med.coef = med.pca.mod$coefficients
med.coefvar = summary(med.pca.mod)$coefficients[,2]

med.pcabeta = med.coef[2]*med.pca$harmonics[1]
med.pcavar = med.coefvar[2]^2 * med.pca$harmonics[1]^2

# iterate through coefficients and add up all harmonics
for (i in 2:dim(pca.scores)[2]){
  med.pcabeta = med.pcabeta + med.coef[i+1]*med.pca$harmonics[i]
  med.pcavar = med.pcavar + med.coefvar[i+1]^2 * med.pca$harmonics[i]^2
}

# compare plot of pca reconstructed coeffficient function with that of functional regression 
par(mfrow=c(2,1))
plot(med.pcabeta, main="PCA coefficient Function", ylim=c(-20,20))
lines(med.pcabeta + 2*sqrt(med.pcavar), lty=2, col="red")
lines(med.pcabeta - 2*sqrt(med.pcavar), lty=2, col="red")

plot(med.coeffun$fd,ylim=c(-5,5), main="Functional Linear Regression Coefficient Function")
lines(med.reg$betaestlist[[2]]$fd+2*sqrt(med.std$betastderrlist[[2]]),lty=2, col="red")
lines(med.reg$betaestlist[[2]]$fd-2*sqrt(med.std$betastderrlist[[2]]),lty=2, col="red")
