library(MASS)
data = read.table('DataTP.txt', header = TRUE)
lm.out = lm(FFo~., data)
# Residuals analysis

X11()
par(mfrow=c(2, 2))
plot(fitted(lm.out), residuals(lm.out), 
     main = "Hypothese d'homoscedasticite", xlab="fitted values", 
     ylab = "residuals") 
qqnorm(residuals(lm.out)) # we can see that the majority of the theoritical quantiles and the sample quantiles 
# match up. this test helps to verify the hypothesis that the residuals follow the normal distribution.
hist(residuals(lm.out)) # same as before we plot the distribution of the residuals to compare it to a normal distribution.
acf(residuals(lm.out)) # this test is to verify that the residuals are independant which is the case because we see a fall off
# in the Acf after the third lag.
summary(lm.out) # we can keep HU P v heure FFp as predictors since they have a p-value below 0.05 so their impact on FFo
# cant be ignored.

lm.outint=lm(FFo~.*.,data)

X11()
par(mfrow=c(2,2))
plot(fitted(lm.outint),residuals(lm.outint),main="Hypothese d'homoscedasticite",xlab="Valeurs ajustees (Y*)",ylab="Residus")
hist(residuals(lm.outint))
qqnorm(residuals(lm.outint))
acf(residuals(lm.outint))

summary(lm.outint)
Bias.out = mean(data$FFo - fitted(lm.out))
RMSE.out = sqrt(mean((data$FFo - fitted(lm.out))**2))
Bias.outint = mean(data$FFo - fitted(lm.outint))
RMSE.outint = sqrt(mean((data$FFo - fitted(lm.outint))**2))
cat("Bias of lm.out :", Bias.out, "\n")
cat("Bias of lm.outint  :", Bias.outint, "\n")
cat("RMSE of lm.out :", RMSE.out, "\n")
cat("RMSE of lm.outint  :", RMSE.outint, "\n")
# the model lm.outint fits the data better in fact the slight U shape in the residuals plot
# suggest that the relation between FFo and the predictors is not really linear.
lm.outAIC    = stepAIC(lm.out)
lm.outintAIC = stepAIC(lm.outint)
lm.outBIC    = stepAIC(lm.out, k=log(nrow(data)))
lm.outintBIC = stepAIC(lm.outint, k=log(nrow(data)))
# The AIC penalizes the number of parameters less strongly than the BIC
eval=function(obs, prev) { 
rmse  = sqrt(mean((prev - obs)**2))
biais = mean(prev - obs)
return(c(biais,rmse))
}
eval(data$FFo, fitted(lm.outBIC))       # RMSE : 2.69
eval(data$FFo, fitted(lm.outAIC))       # RMSE : 2.68
eval(data$FFo, fitted(lm.outintBIC))    # RMSE : 2.61
eval(data$FFo, fitted(lm.outintAIC))    # RMSE : 2.58
# we cant really conclude anything we need to cross validate
nappr = ceiling(0.8*nrow(data))
ii = sample(1:nrow(data),nappr)
jj = setdiff(1:nrow(data),ii)
datatest = data[jj,]
datapp   = data[ii,]
lm.outAIC    = lm(formula(lm.outAIC), datapp)
lm.outBIC    = lm(formula(lm.outBIC), datapp)
lm.outintAIC = lm(formula(lm.outintAIC), datapp)
lm.outintBIC = lm(formula(lm.outintBIC), datapp)
eval(datapp$FFo,fitted(lm.outAIC))              # RMSE : 2.70
eval(datapp$FFo,fitted(lm.outBIC))              # RMSE : 2.71   
eval(datapp$FFo,fitted(lm.outintAIC))           # RMSE : 2.60
eval(datapp$FFo,fitted(lm.outintBIC))           # RMSE : 2.63
eval(datatest$FFo, predict(lm.outAIC, datatest))              # RMSE : 2.61
eval(datatest$FFo, predict(lm.outBIC, datatest))              # RMSE : 2.63   
eval(datatest$FFo, predict(lm.outintAIC, datatest))           # RMSE : 2.54
eval(datatest$FFo, predict(lm.outintBIC, datatest))           # RMSE : 2.54
# we can see that in our simple test the int models outperform the simple models on the train data and the test data.
# K-fold cross validation
set.seed(10)
K = 10
rows = sample(nrow(data))
datashuffle = data[rows,]
ntest = ceiling((1/K) * nrow(data))
endtest = 0
a.aic = c()
a.bic = c()
a.aint = c()
a.bint = c()
a.INT = c()
t.aic = c()
t.bic = c()
t.aint = c()
t.bint = c()
t.INT = c()
for (i in 1:K) {
    starttest = endtest + 1
    endtest   = i * ntest
    if (endtest > nrow(data)) {
       endtest = nrow(data)
    }
    datatest = datashuffle[starttest:endtest,]
    datapp   = datashuffle[setdiff(1:nrow(data), starttest:endtest),]

    lm.outAIC = lm(formula(lm.outAIC), datapp)
    lm.outBIC = lm(formula(lm.outBIC), datapp)
    lm.outintAIC = lm(formula(lm.outintAIC), datapp)
    lm.outintBIC = lm(formula(lm.outintBIC), datapp)
    lm.INT = lm(FFo~.*.*. ,datapp)

    t.aic = c(t.aic, eval(datatest$FFo, predict(lm.outAIC, datatest))[2])
    t.bic = c(t.bic, eval(datatest$FFo, predict(lm.outBIC, datatest))[2])
    t.aint = c(t.aint, eval(datatest$FFo, predict(lm.outintAIC, datatest))[2])
    t.bint = c(t.bint, eval(datatest$FFo, predict(lm.outintBIC, datatest))[2])
    t.INT  = c(t.INT, eval(datatest$FFo, predict(lm.INT, datatest))[2])

    a.aic = c(a.aic, eval(datapp$FFo, fitted(lm.outAIC))[2])
    a.bic = c(a.bic, eval(datapp$FFo, fitted(lm.outBIC))[2])
    a.aint = c(a.aint, eval(datapp$FFo, fitted(lm.outintAIC))[2])
    a.bint = c(a.bint, eval(datapp$FFo, fitted(lm.outintBIC))[2])
    a.INT  = c(a.INT, eval(datapp$FFo, fitted(lm.INT))[2])
}

X11()
boxplot(a.bic, a.aic, a.bint, a.aint, a.INT, t.bic, t.aic, t.bint, t.aint, t.INT, 
        col   = c("blue", "blue", "blue", "blue", "blue", "red", "red", "red", "red", "red"), 
        names = c("a.bic", "a.aic", "a.bint", "a.aint", "a.INT", "t.bic", "t.aic", "t.bint", "t.aint", "t.INT"), 
        main  = "Modele lineaire gaussien - Score RMSE")

# the best model to choose is lm.outintBIC since it's simpler than lm.outintAIC, gives comparable resultes to lm.outintAIC on the test data
# and it's relatively robuste (not as robust as the lm.outAIC and the lm.outBIC but still gives better results)

X11()
plot(data$FFo[300:400], type ="l", main="Force de vent",xlab="Date",ylab="[FFo]")
points(data$FFp[300:400], col="blue", pch="+")
points(fitted(lm.outintBIC)[300:400], col="red", pch="+")
legend(0, 25, lty=1, col=c("black"), legend=c("FFo"), bty="n")
legend(0, 23, pch="+", col="blue", legend="FFp", bty="n") 
legend(0, 21, pch="+", col="red", legend="régression intbic", bty="n")
