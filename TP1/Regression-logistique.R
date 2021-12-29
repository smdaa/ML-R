library(MASS)
library(verification)
data = read.table('DataTP.txt', header = TRUE)
data$OCC = as.factor(as.numeric(data$FFo>13))
data$OCCp = as.factor(as.numeric(data$FFp>13))
source("scores.R")
glm.out    = glm(OCC~., data[,-11], family=binomial)    # we remove FFo from the dataset because duh 
glm.outint = glm(OCC~.*.,data[,-11],family=binomial)
glm.outAIC = stepAIC(glm.out)
glm.outBIC = stepAIC(glm.out,k=log(nrow(data)))
glm.outintAIC = stepAIC(glm.outint)
glm.outintBIC = stepAIC(glm.outint,k=log(nrow(data)))
formula(glm.outAIC)
formula(glm.outBIC)
formula(glm.outintAIC)
formula(glm.outintBIC)

X11()
roc.plot(as.numeric(data$OCC)-1,fitted(glm.outAIC))
# We choose the probability threshold 0.2 which allows to reach a hit rate of 0.9 and a false rate of 0.1

scores(fitted(glm.outAIC)>0.2, data$OCC)                # PSS = 0.69
brier(as.numeric(data$OCC)-1, fitted(glm.outAIC))$bs    #BS = 0.083 equivalent for RMSE
roc.area(as.numeric(data$OCC)-1, fitted(glm.outAIC))$A  # AUC of 1 is a perfect score
nappr=ceiling(0.8*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]
glm.outAIC = glm(formula(glm.outAIC), datapp, family=binomial)
glm.outBIC = glm(formula(glm.outBIC), datapp, family=binomial)
glm.outintAIC = glm(formula(glm.outintAIC), datapp, family=binomial)
glm.outintBIC =glm(formula(glm.outintBIC), datapp, family=binomial)
scores(fitted(glm.outAIC)>0.2, datapp$OCC)       # PSS = 0.70
scores(fitted(glm.outBIC)>0.2, datapp$OCC)       # PSS = 0.70
scores(fitted(glm.outintAIC)>0.2, datapp$OCC)    # PSS = 0.72
scores(fitted(glm.outintBIC)>0.2, datapp$OCC)    # PSS = 0.72
brier(as.numeric(datapp$OCC)-1, fitted(glm.outAIC))$bs          # bs = 0.079
brier(as.numeric(datapp$OCC)-1, fitted(glm.outBIC))$bs          # bs = 0.080
brier(as.numeric(datapp$OCC)-1, fitted(glm.outintAIC))$bs       # bs = 0.072
brier(as.numeric(datapp$OCC)-1, fitted(glm.outintBIC))$bs       # bs = 0.074
scores(predict(glm.outAIC, datatest, type="response")>0.2, datatest$OCC)       # PSS = 0.64
scores(predict(glm.outBIC, datatest, type="response")>0.2, datatest$OCC)       # PSS = 0.64
scores(predict(glm.outintAIC, datatest, type="response")>0.2, datatest$OCC)    # PSS = 0.68
scores(predict(glm.outintBIC, datatest, type="response")>0.2, datatest$OCC)    # PSS = 0.67
brier(as.numeric(datatest$OCC)-1, predict(glm.outAIC, datatest, type="response")>0.2)$bs             #bs = 0.16
brier(as.numeric(datatest$OCC)-1, predict(glm.outBIC, datatest, type="response")>0.2)$bs             #bs = 0.16
brier(as.numeric(datatest$OCC)-1, predict(glm.outintAIC, datatest, type="response")>0.2)$bs           #bs = 0.14
brier(as.numeric(datatest$OCC)-1, predict(glm.outintBIC, datatest, type="response")>0.2)$bs           #bs = 0.15
# cross validation
set.seed(0)
K = 10
rows = sample(nrow(data))
datashuffle = data[rows,]
ntest = ceiling((1/K) * nrow(data))
endtest = 0
a.aic = matrix(NA, nrow=K,ncol=3)
a.bic = matrix(NA, nrow=K,ncol=3)
a.aint = matrix(NA, nrow=K,ncol=3)
a.bint = matrix(NA, nrow=K,ncol=3)
t.aic = matrix(NA, nrow=K,ncol=3)
t.bic = matrix(NA, nrow=K,ncol=3)
t.aint = matrix(NA, nrow=K,ncol=3)
t.bint = matrix(NA, nrow=K,ncol=3)
for (i in 1:K) {
    starttest = endtest + 1
    endtest   = min(i * ntest, nrow(data)) 
    datatest = datashuffle[starttest:endtest,]
    datapp   = datashuffle[setdiff(1:nrow(data), starttest:endtest),]
    glm.outAIC    = glm(formula(glm.outAIC), datapp, family=binomial)
    glm.outBIC    = glm(formula(glm.outBIC), datapp, family=binomial)
    glm.outintAIC = glm(formula(glm.outintAIC), datapp, family=binomial)
    glm.outintBIC = glm(formula(glm.outintBIC), datapp, family=binomial)
    t.aic[i, 1]  = scores(predict(glm.outAIC, datatest, type="response")>0.2, datatest$OCC, quiet=TRUE)[4]
    t.aic[i, 2]  = brier(as.numeric(datatest$OCC)-1, predict(glm.outAIC, datatest, type="response")>0.2)$bs
    t.aic[i, 3]  = roc.area(as.numeric(datatest$OCC)-1, predict(glm.outAIC, datatest, type="response"))$A
    t.bic[i, 1]  = scores(predict(glm.outBIC, datatest, type="response")>0.2, datatest$OCC, quiet=TRUE)[4]
    t.bic[i, 2]  = brier(as.numeric(datatest$OCC)-1, predict(glm.outBIC, datatest, type="response")>0.2)$bs
    t.bic[i, 3]  = roc.area(as.numeric(datatest$OCC)-1, predict(glm.outBIC, datatest, type="response"))$A
    t.aint[i, 1] = scores(predict(glm.outintAIC, datatest, type="response")>0.2, datatest$OCC, quiet=TRUE)[4]
    t.aint[i, 2] = brier(as.numeric(datatest$OCC)-1, predict(glm.outintAIC, datatest, type="response")>0.2)$bs
    t.aint[i, 3] = roc.area(as.numeric(datatest$OCC)-1, predict(glm.outintAIC, datatest, type="response"))$A
    t.bint[i, 1] = scores(predict(glm.outintBIC, datatest, type="response")>0.2, datatest$OCC, quiet=TRUE)[4]
    t.bint[i, 2] = brier(as.numeric(datatest$OCC)-1, predict(glm.outintBIC, datatest, type="response")>0.2)$bs
    t.bint[i, 3] = roc.area(as.numeric(datatest$OCC)-1, predict(glm.outintBIC, datatest, type="response"))$A
    a.aic[i, 1]  = scores(fitted(glm.outAIC)>0.2, datapp$OCC, quiet=TRUE)[4]
    a.aic[i, 2]  = brier(as.numeric(datapp$OCC)-1, fitted(glm.outAIC))$bs 
    a.aic[i, 3]  = roc.area(as.numeric(datapp$OCC)-1, fitted(glm.outAIC))$A 
    a.bic[i, 1]  = scores(fitted(glm.outBIC)>0.2, datapp$OCC, quiet=TRUE)[4]
    a.bic[i, 2]  = brier(as.numeric(datapp$OCC)-1, fitted(glm.outBIC))$bs 
    a.bic[i, 3]  = roc.area(as.numeric(datapp$OCC)-1, fitted(glm.outBIC))$A 
    a.aint[i, 1] = scores(fitted(glm.outintAIC)>0.2, datapp$OCC, quiet=TRUE)[4]
    a.aint[i, 2] = brier(as.numeric(datapp$OCC)-1, fitted(glm.outintAIC))$bs 
    a.aint[i, 3] = roc.area(as.numeric(datapp$OCC)-1, fitted(glm.outintAIC))$A 
    a.bint[i, 1] = scores(fitted(glm.outintBIC)>0.2, datapp$OCC, quiet=TRUE)[4]
    a.bint[i, 2] = brier(as.numeric(datapp$OCC)-1, fitted(glm.outintBIC))$bs 
    a.bint[i, 3] = roc.area(as.numeric(datapp$OCC)-1, fitted(glm.outintBIC))$A 
}
# score PSS higher is better
X11()
boxplot(a.bic[, 1], a.aic[, 1], a.bint[, 1], a.aint[, 1], t.bic[, 1], t.aic[, 1], t.bint[, 1], t.aint[, 1], 
        col   = c("blue", "blue", "blue", "blue", "red", "red", "red", "red"), 
        names = c("a.bic", "a.aic", "a.bint", "a.aint", "t.bic", "t.aic", "t.bint", "t.aint"), 
        main  = "Regression logistique  Score PSS")

# score bs lower is better
X11()
boxplot(a.bic[, 2], a.aic[, 2], a.bint[, 2], a.aint[, 2], t.bic[, 2], t.aic[, 2], t.bint[, 2], t.aint[, 2], 
        col   = c("blue", "blue", "blue", "blue", "red", "red", "red", "red"), 
        names = c("a.bic", "a.aic", "a.bint", "a.aint", "t.bic", "t.aic", "t.bint", "t.aint"), 
        main  = "Regression logistique  Brier Score")

# score auc higher is better
X11()
boxplot(a.bic[, 3], a.aic[, 3], a.bint[, 3], a.aint[, 3], t.bic[, 3], t.aic[, 3], t.bint[, 3], t.aint[, 3], 
        col   = c("blue", "blue", "blue", "blue", "red", "red", "red", "red"), 
        names = c("a.bic", "a.aic", "a.bint", "a.aint", "t.bic", "t.aic", "t.bint", "t.aint"), 
        main  = "Regression logistique  ROC AREA")
lm.out    = lm(FFo~., data[, 1:11])
lm.outBIC = stepAIC(lm.out, k=log(nrow(data)))
K = 10
set.seed(100)
rows = sample(nrow(data))
datashuffle = data[rows,]
ntest = ceiling((1/K) * nrow(data))
endtest = 0
a.MLG = c()
a.RL  = c()
t.MLG = c()
t.RL  = c()
for (i in 1:K) {
  starttest = endtest + 1
  endtest   = min(i * ntest, nrow(data)) 
  datatest = datashuffle[starttest:endtest,]
  datapp   = datashuffle[setdiff(1:nrow(data), starttest:endtest),]
  glm.outAIC = glm(formula(glm.outAIC), datapp, family=binomial)
  lm.outBIC  = lm(formula(lm.outBIC), datapp)
  t.RL = c(t.RL, scores(predict(glm.outAIC, datatest, type="response")>0.2, datatest$OCC, quiet=TRUE)[4])
  a.RL = c(a.RL, scores(fitted(glm.outAIC)>0.2, datapp$OCC, quiet=TRUE)[4])
  t.MLG = c(t.MLG, scores(predict(lm.outBIC, datatest)>13, datatest$OCC, quiet=TRUE)[4])
  a.MLG = c(a.MLG, scores(fitted(lm.outBIC)>13, datapp$OCC, quiet=TRUE)[4])
}
# score PSS higher is better
X11()
boxplot(a.MLG, a.RL, t.MLG, t.RL , 
        col   = c("blue", "blue", "red", "red"), 
        names = c("a.MLG", "a.RL", "t.MLG", "t.RL"),
        main  = "Regression logistique vs Modele Lineaire Gaussien  Score PSS") 
