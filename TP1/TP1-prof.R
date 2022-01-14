# Script TP1 HPC / Big Data 2021 (sans procedure de validation a coder vous-meme)


# 1 - Chargement librairies et donnees

library(MASS)
data=read.table(file="DataTP.txt",header=TRUE)


# 2 - Etude preliminaire 

summary(data)
dim(data)


sd(data$FFo)
sd(data$FFp)

par(mfrow=c(1,2))
hist(data$FFo,breaks=seq(0,40,2),ylim=c(0,500))
hist(data$FFp,breaks=seq(0,40,2),ylim=c(0,500))

X11()
boxplot(data$FFo,data$FFp,names=c("FFo","FFp"))

var.test(data$FFo,data$FFp)
# -> les variances sont jugees significativement differentes
t.test(data$FFo,data$FFp)
# -> les moyennes sont jugees significativement differentes

plot(data$FFp,data$FFo)
abline(lm(FFo~FFp,data),col="red")
cor(data$FFo,data$FFp)
cor.test(data$FFo,data$FFp)

pairs(data[,-c(2,8,9)])

x11()
par(mfrow=c(3,2))
hist(data$FFp)
hist(data$P)
hist(data$u)
hist(data$v)
hist(data$hel)
hist(data$DD)

plot(data$FFo[1:300],type ="l",lwd=2,main="Concentration d'ozone",xlab="Date",ylab="[O3]")
points(data$FFp[1:300],col="blue",pch="+")
points(fitted(lm(FFo~FFp,data))[1:300],col="red",pch="+")
legend(0,25,lty=1,col=c("black"),legend=c("FFo"),bty="n")
legend(0,23,pch="+",col="blue",legend="                    FFp",bty="n") 
legend(0,21,pch="+",col="red",legend="                    AS",bty="n") 


# 3 - Modele lineaire gaussien

lm.out=lm(FFo~.,data)
summary(lm.out)
model.matrix(lm.out)[300:310,]

par(mfrow=c(2,2))
plot(fitted(lm.out),residuals(lm.out),main="Hypothese d'homoscedasticite",xlab="Valeurs ajustees (Y*)",ylab="Residus")
hist(residuals(lm.out))
qqnorm(residuals(lm.out))
acf(residuals(lm.out))

lm.outint=lm(FFo~.*.,data)
summary(lm.outint)

par(mfrow=c(2,2))
plot(fitted(lm.outint),residuals(lm.outint),main="Hypothese d'homoscedasticite",xlab="Valeurs ajustees (Y*)",ylab="Residus")
hist(residuals(lm.outint))
qqnorm(residuals(lm.outint))
acf(residuals(lm.outint))

lm.outAIC=stepAIC(lm.out)
lm.outintAIC=stepAIC(lm.outint)
lm.outBIC=stepAIC(lm.out,k=log(nrow(data)))
lm.outintBIC=stepAIC(lm.outint,k=log(nrow(data)))

formula(lm.outAIC)
formula(lm.outBIC)
formula(lm.outintAIC)
formula(lm.outintBIC)

summary(lm.outAIC)
summary(lm.outBIC)
summary(lm.outintAIC)
summary(lm.outintBIC)


# 4 - Evaluation des modeles

eval=function(obs,prev) {
rmse=sqrt(mean((prev-obs)**2))
biais=mean(prev-obs)
print("Biais         RMSE") 
return(c(biais,rmse))
}

eval(data$FFo,data$FFp)
eval(data$FFo,fitted(lm.outBIC))
eval(data$FFo,fitted(lm.outAIC))
eval(data$FFo,fitted(lm.outintBIC))
eval(data$FFo,fitted(lm.outintAIC))


nappr=ceiling(0.8*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]

lm.outAIC=lm(formula(lm.outAIC),datapp)
lm.outBIC=lm(formula(lm.outBIC),datapp)
lm.outintAIC=lm(formula(lm.outintAIC),datapp)
lm.outintBIC=lm(formula(lm.outintBIC),datapp)

eval(datapp$FFo,fitted(lm.outBIC))
eval(datapp$FFo,fitted(lm.outAIC))
eval(datapp$FFo,fitted(lm.outintBIC))
eval(datapp$FFo,fitted(lm.outintAIC))

eval(datatest$FFo,predict(lm.outBIC,datatest))
eval(datatest$FFo,predict(lm.outAIC,datatest))
eval(datatest$FFo,predict(lm.outintBIC,datatest))
eval(datatest$FFo,predict(lm.outintAIC,datatest))

plot(datatest$FFo,type ="l",main="Concentration d'ozone",xlab="Indice",ylab="[O3]")
points(datatest$FFp,col="blue",pch="+")
points(predict(lm.outBIC,datatest),col="red",pch="+")
legend(1,30,lty=1,col=c("black"),legend=c("observe"),bty="n")
legend(1,29,pch="+",col=c("blue","red"),legend=c("       FFp","       ASBIC"),bty="n") 


# ---> Coder procedure de validation croisee
# l'objectif est de confronter les modeles sur apprentissage et test, 
# en analysant leur sensibilite a l'echantillonnage

# ---> Confronter les 4 modeles lineaires



# 5 - Prevision du depassement de seuil

library(verification)
library(MASS)

data=read.table(file="DataTP.txt",header=TRUE)
data$OCC=as.factor(as.numeric(data$FFo>13))
data$OCCp=as.factor(as.numeric(data$FFp>13))
summary(data)

source("scores.R")


#####################################
# Regression logistique
#####################################


glm.out=glm(OCC~.,data[,-11],family=binomial)
summary(glm.out)
glm.outint=glm(OCC~.*.,data[,-11],family=binomial)
summary(glm.outint)


glm.outAIC=stepAIC(glm.out)
glm.outBIC=stepAIC(glm.out,k=log(nrow(data)))
glm.outintAIC=stepAIC(glm.outint)
glm.outintBIC=stepAIC(glm.outint,k=log(nrow(data)))

formula(glm.outAIC)
formula(glm.outBIC)
formula(glm.outintAIC)
formula(glm.outintBIC)

scores(data$OCCp,data$OCC)
scores(fitted(glm.outAIC)>0.5,data$OCC)

roc.plot(as.numeric(data$OCC)-1,fitted(glm.outAIC))
scores(fitted(glm.outAIC)>0.2,data$OCC)

brier(as.numeric(data$OCC)-1,fitted(glm.outAIC))$bs
roc.area(as.numeric(data$OCC)-1,fitted(glm.outAIC))$A

nappr=ceiling(0.8*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]

glm.outAIC=glm(formula(glm.outAIC),datapp,family=binomial)
glm.outBIC=glm(formula(glm.outBIC),datapp,family=binomial)
glm.outintAIC=glm(formula(glm.outintAIC),datapp,family=binomial)
glm.outintBIC=glm(formula(glm.outintBIC),datapp,family=binomial)

scores(fitted(glm.outAIC)>0.2,datapp$OCC)
scores(fitted(glm.outBIC)>0.2,datapp$OCC)
scores(fitted(glm.outintBIC)>0.2,datapp$OCC)
scores(fitted(glm.outintAIC)>0.2,datapp$OCC)

scores(predict(glm.outAIC,datatest,type="response")>0.2,datatest$OCC)
scores(predict(glm.outBIC,datatest,type="response")>0.2,datatest$OCC)
scores(predict(glm.outintBIC,datatest,type="response")>0.2,datatest$OCC)
scores(predict(glm.outintAIC,datatest,type="response")>0.2,datatest$OCC)

#---> Confronter les 4 regressions logistiques avec votre procedure de validation en terme de PSS, BS et ROC AREA
#---> Confronter le meilleur modele logistique avec la strategie suivante :
# prevoir la force du vent avec un modele lineraire gaussien et en deduire une prevision de depassement du seuil de 13 m/s.



#####################################
# Analyse discriminante
#####################################


library(MASS)

iOCC=which(data$OCC==1)
iNOCC=which(data$OCC==0)
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

par(mfrow=c(1,6))
boxplot(data$FFp[iOCC],data$FFp[iNOCC],names=c("FFp_OCC","FFp_NOCC"))
boxplot(data$HU[iOCC],data$HU[iNOCC],names=c("HU_OCC","HU_NOCC"))
boxplot(data$P[iOCC],data$P[iNOCC],names=c("P_OCC","P_NOCC"))
boxplot(data$u[iOCC],data$u[iNOCC],names=c("u_OCC","u_NOCC"))
boxplot(data$v[iOCC],data$v[iNOCC],names=c("v_OCC","v_NOCC"))
boxplot(data$hel[iOCC],data$hel[iNOCC],names=c("hel_OCC","hel_NOCC"))

t.test(data$FFp[iOCC],data$FFp[iNOCC])
t.test(data$HU[iOCC],data$HU[iNOCC])
t.test(data$P[iOCC],data$P[iNOCC])
t.test(data$u[iOCC],data$u[iNOCC])
t.test(data$v[iOCC],data$v[iNOCC])
t.test(data$hel[iOCC],data$hel[iNOCC])

var.test(data$FFp[iOCC],data$FFp[iNOCC])
var.test(data$HU[iOCC],data$HU[iNOCC])
var.test(data$P[iOCC],data$P[iNOCC])
var.test(data$u[iOCC],data$u[iNOCC])
var.test(data$v[iOCC],data$v[iNOCC])
var.test(data$hel[iOCC],data$hel[iNOCC])

var(data[iOCC,-c(11,12,13)])-var(data[iNOCC,-c(11,12,13)])
cor(data[iOCC,-c(11,12,13)])-cor(data[iNOCC,-c(11,12,13)])

pairs(data[iNOCC,c(1,3,4,5,6,10)])
X11()
pairs(data[iOCC,c(1,3,4,5,6,10)])


lda.out=lda(OCC~.,data[,-11]) 
roc.plot(as.numeric(data$OCC)-1,predict(lda.out)$posterior[,2])
scores(predict(lda.out)$posterior[,2]>0.2,data$OCC)

lda.out=lda(OCC~FFp+v,data[,-11]) 
roc.plot(as.numeric(data$OCC)-1,predict(lda.out)$posterior[,2])
scores(predict(lda.out)$posterior[,2]>0.2,data$OCC)


qda.out=qda(OCC~.,data[,-11]) 
roc.plot(as.numeric(data$OCC)-1,predict(qda.out)$posterior[,2])
scores(predict(qda.out)$posterior[,2]>0.1,data$OCC)


nappr=ceiling(0.8*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datapp=data[ii,]
datatest=data[jj,]


lda.out=lda(OCC~.,datapp[,-11]) 
roc.plot(as.numeric(datapp$OCC)-1,predict(lda.out)$posterior[,2])
scores(predict(lda.out)$posterior[,2]>0.2,datapp$OCC)
scores(predict(lda.out,datatest)$posterior[,2]>0.2,datatest$OCC)
brier(as.numeric(datapp$OCC)-1,predict(lda.out)$posterior[,2])$bs 
brier(as.numeric(datatest$OCC)-1,predict(lda.out,datatest)$posterior[,2])$bs 

qda.out=qda(OCC~.,datapp[,-11]) 
roc.plot(as.numeric(datapp$OCC)-1,predict(qda.out)$posterior[,2])
scores(predict(qda.out)$posterior[,2]>0.1,datapp$OCC)
scores(predict(qda.out,datatest)$posterior[,2]>0.1,datatest$OCC)
brier(as.numeric(datapp$OCC)-1,predict(qda.out)$posterior[,2])$bs 
brier(as.numeric(datatest$OCC)-1,predict(qda.out)$posterior[,2])$bs 

# --> Confronter analyses discriminantes lineaire et quadratique, puis avec regression logistique


