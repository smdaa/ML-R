library(MASS)
source("scores.R")
data = read.table('DataTP.txt', header = TRUE)
data$OCC = as.factor(as.numeric(data$FFo>13))
data$OCCp = as.factor(as.numeric(data$FFp>13))
iOCC = which(data$OCC==1)
iNOCC = which(data$OCC==0)

X11()
par(mfrow=c(2,6))
hist(data$FFp[iOCC])
hist(data$HU[iOCC])
hist(data$P[iOCC])
hist(data$u[iOCC])
hist(data$v[iOCC])
hist(data$hel[iOCC])
hist(data$FFp[iNOCC])
hist(data$HU[iNOCC])
hist(data$P[iNOCC])
hist(data$u[iNOCC])
hist(data$v[iNOCC])
hist(data$hel[iNOCC])

X11()
par(mfrow=c(1,6))
boxplot(data$FFp[iOCC],data$FFp[iNOCC],names=c("FFp_OCC","FFp_NOCC"))
boxplot(data$HU[iOCC],data$HU[iNOCC],names=c("HU_OCC","HU_NOCC"))
boxplot(data$P[iOCC],data$P[iNOCC],names=c("P_OCC","P_NOCC"))
boxplot(data$u[iOCC],data$u[iNOCC],names=c("u_OCC","u_NOCC"))
boxplot(data$v[iOCC],data$v[iNOCC],names=c("v_OCC","v_NOCC"))
boxplot(data$hel[iOCC],data$hel[iNOCC],names=c("hel_OCC","hel_NOCC"))

t.test(data$FFp[iOCC],data$FFp[iNOCC])            # not equal
t.test(data$HU[iOCC],data$HU[iNOCC])              # not equal
t.test(data$P[iOCC],data$P[iNOCC])                # not equal
t.test(data$u[iOCC],data$u[iNOCC])                # equal
t.test(data$v[iOCC],data$v[iNOCC])                # equal
t.test(data$hel[iOCC],data$hel[iNOCC])            # not equal

var.test(data$FFp[iOCC],data$FFp[iNOCC])          # not equal
var.test(data$HU[iOCC],data$HU[iNOCC])            # not equal
var.test(data$P[iOCC],data$P[iNOCC])              # not equal
var.test(data$u[iOCC],data$u[iNOCC])              # not equal
var.test(data$v[iOCC],data$v[iNOCC])              # not equal
var.test(data$hel[iOCC],data$hel[iNOCC])          # not equal

var(data[iOCC,-c(11,12,13)]) - var(data[iNOCC,-c(11,12,13)])
cor(data[iOCC,-c(11,12,13)]) - cor(data[iNOCC,-c(11,12,13)])
# in order to choose between lda or qda we need to assess the difference between the covariance matrices

X11()
pairs(data[iNOCC, c(1,3,4,5,6,10)])

X11()
pairs(data[iOCC, c(1,3,4,5,6,10)])

lda.out = lda(OCC~., data[,-11]) 

X11()
roc.plot(as.numeric(data$OCC)-1, predict(lda.out)$posterior[,2])

scores(predict(lda.out)$posterior[,2]>0.2,data$OCC) 
