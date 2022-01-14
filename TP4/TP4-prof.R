# ScriptTP4 HPC / Big Data 2021


# Chargement des packages
library(verification)
library(rpart)

source("scores.R")

eval=function(obs,prev) {
rmse=sqrt(mean((prev-obs)**2))
biais=mean(prev-obs)
print("Biais       RMSE") 
return(c(biais,rmse))
}

RMSE=function(obs,pr){
return(sqrt(mean((pr-obs)^2)))}


################
# Bagging
################

library(ipred)

data=read.table(file="DataTP.txt",header=TRUE)

#Regression

bagging(FFo~.,data,coob=T)
bagging(FFo~.,data,coob=T,control=rpart.control(cp=0))
bagging(FFo~.,data,nbagg=50,coob=T)
bagging(FFo~.,data,nbagg=500,coob=T)
bagging(FFo~.,data,nbagg=10,coob=T)
bagging(FFo~.,data,coob=T,control=rpart.control(cp=0.1))
bagreg.out=bagging(FFo~.,data,coob=T)
eval(data$FFo,predict(bagreg.out,data))


# influence du nb d'arbres nbag :

erreur=NULL
x=seq(10,300,10)
for (i in x) {erreur=c(erreur,bagging(FFo~.,data,nbagg=i,coob=T)$err)}
plot(x,erreur,type="l",ylim=c(2.7,2.8),xlab="nbagg",ylab="Erreur OOB",main="Bagging : influence du nb d'arbres")

# --> Stabilisation de l'erreur OOB a partir de nbag=100


nappr=ceiling(0.5*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]

# influence de l'elagage :

erreur=NULL
erreur2=NULL
x=seq(0,0.02,0.001)
for (i in x) {
bag.out=bagging(FFo~.,datapp,nbag=100,control=rpart.control(cp=i))
erreur=c(erreur,RMSE(datapp$FFo,predict(bag.out,datapp)))
erreur2=c(erreur2,RMSE(datatest$FFo,predict(bag.out,datatest)))
}
plot(x,erreur,type="l",ylim=c(1.6,3),col="blue",main="Bagging : influence de l'elagage",xlab="cp",ylab="RMSE (bleu -> appr , rouge -> test)")
lines(x,erreur2,type="l",col="red")

#--> cp=0.01 (valeur par defaut) convient ici, surapprentissage severe si arbres peu elages



#Discrimination

data$OCC=as.factor(as.numeric(data$FFo>13))
data$OCCp=as.factor(as.numeric(data$FFp>13))

bagging(OCC~.,data[,-11],nbag=50,coob=T)
bagging(OCC~.,data[,-11],nbag=25,coob=T)
bagdis.out=bagging(OCC~.,data[,-11],nbag=50,coob=T)
roc.plot(as.numeric(data$OCC)-1,predict(bagdis.out,type="prob")[,2])
scores(predict(bagdis.out,type="prob")[,2]>0.2,data$OCC)

# influence du nb d'arbres nbagg :

erreur=NULL
x=seq(10,300,10)
for (i in x) {erreur=c(erreur,bagging(OCC~.,data[,-11],nbagg=i,coob=T)$err)}
plot(x,erreur,type="l",xlab="nbag",ylab="Erreur OOB",main="Bagging : influence du nb d'arbres")

# --> Stabilisation de l'erreur a partir de nbag=100 (25 par defaut n'est pas adapte ici)

nappr=ceiling(0.5*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]

# influence de l'elagage :

erreur=NULL
erreur2=NULL
x=seq(0,0.02,0.001)
for (i in x) {
bag.out=bagging(OCC~.,datapp[,-11],nbag=100,control=rpart.control(cp=i))
erreur=c(erreur,brier(as.numeric(datapp$OCC)-1,predict(bag.out,datapp,type="prob")[,2])$bs)
erreur2=c(erreur2,brier(as.numeric(datatest$OCC)-1,predict(bag.out,datatest,type="prob")[,2])$bs)
}
plot(x,erreur,type="l",col="blue",ylim=c(0,0.2),main="Bagging : influence de l'elagage",xlab="cp",ylab="RMSE (bleu -> appr , rouge -> test)")
lines(x,erreur2,type="l",col="red")

#--> cp=0.01 (valeur par defaut) convient ici, surapprentissage severe si arbres peu elages



################
# Random Forest
################

library(randomForest)
data=read.table(file="DataTP.txt",header=TRUE)

#Regression

rfreg.out=randomForest(FFo~.,data,importance=TRUE)
sort(importance(rfreg.out)[,1],dec=T)
eval(data$FFo,predict(rfreg.out,data))

# influence du nb d'arbres ntree :

RMSE=function(obs,pr){
return(sqrt(mean((pr-obs)^2)))}

erreur=NULL
x=seq(10,600,30)
for (i in x) {
rf.out=randomForest(FFo~.,data,ntree=i)
erreur=c(erreur,RMSE(data$FFo,predict(rf.out))) # par defaut predict renvoie les previsions OOB
}
plot(x,erreur,type="l",ylim=c(2.5,2.7),xlab="ntree",main="RandomForest : influence du nb d'arbres",ylab="RMSE")

# --> Stabilisation de l'erreur OOB a partir de ntree=200 ou 300


# influence du nb de predicteurs tires au hasard mtry :

erreur=NULL
x=1:10
for (i in x) {
rf.out=randomForest(FFo~.,data,ntree=200,mtry=i)
erreur=c(erreur,RMSE(data$FFo,predict(rf.out)))
}
plot(x,erreur,type="l",main="RandomForest - influence de mtry",xlab="mtry",ylab="RMSE")

# --> minimum vers 2-3, le choix par defaut E(q/3)=3 est bon ici


nappr=ceiling(0.5*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]

# influence de la complexite des arbres maxnodes (=nb max de feuilles):

erreur=NULL
erreur2=NULL
x=seq(10,200,10)
for (i in x) {
rf.out=randomForest(FFo~.,datapp,ntree=200,maxnodes=i)
erreur=c(erreur,RMSE(datapp$FFo,predict(rf.out,datapp)))
erreur2=c(erreur2,RMSE(datatest$FFo,predict(rf.out,datatest)))
}
plot(x,erreur,type="l",ylim=c(1,3),col="blue",main="RandomForest - influence du nb de feuilles",xlab="maxnodes",ylab="RMSE (bleu -> appr , rouge -> test)")
lines(x,erreur2,type="l",col="red")

#--> qq dizaines suffisent sinon surapprentissage severe


#Discrimination

data$OCC=as.factor(as.numeric(data$FFo>13))
data$OCCp=as.factor(as.numeric(data$FFp>13))

rfdis.out=randomForest(OCC~.,data=data[,-11],importance=TRUE)
sort(importance(rfdis.out)[,1],dec=T)


roc.plot(as.numeric(data$OCC)-1,predict(rfdis.out,type="prob")[,2])
scores(predict(rfdis.out,type="prob")[,2]>0.25,data$OCC)
brier(as.numeric(data$OCC)-1,predict(rfdis.out,type="prob")[,2])$bs

#--> par rapport au Bagging, sur nos donnees, on observe un gain significatif en regression mais pas en discrimination



################
# Boosting
################

library(gbm)
data=read.table(file="DataTP.txt",header=TRUE)

#Regression

boostreg.out=gbm(FFo~.,data,distribution="gaussian",n.trees=5000,cv.folds=10)
boostreg.out
summary(boostreg.out)
plot(boostreg.out$cv.error)
gbm.perf(boostreg.out)
boostreg.out=gbm(FFo~.,data,distribution="gaussian",n.trees=250)
eval(data$FFo,predict(boostreg.out))


nappr=ceiling(0.5*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]


RMSE=function(obs,pr){
return(sqrt(mean((pr-obs)^2)))}

# influence du shrinkage :

erreur=NULL
erreur2=NULL
x=seq(0.01,2,0.1)
for (i in x) {
boostreg.out=gbm(FFo~.,datapp,distribution="gaussian",n.trees=6000,cv.folds=10,shrinkage=i)
erreur=c(erreur,RMSE(datapp$FFo,predict(boostreg.out)))
erreur2=c(erreur2,RMSE(datatest$FFo,predict(boostreg.out,datatest)))
}
plot(x,erreur,type="l",col="blue",ylim=c(2.4,3.5),main="Boosting - influence du shrinkage",xlab="shrinkage",ylab="RMSE (bleu -> appr , rouge -> test)")
lines(x,erreur2,type="l",col="red")

# --> 0.1, valeur par defaut convient ici
# le nb d'arbres optimal depend de la valeur du shrinkage, autour de 200 pour shrinkage=0.1


# influence de l'elagage - interaction.depth :

erreur=NULL
erreur2=NULL
x=seq(1,30,1)
for (i in x) {
boostreg.out=gbm(FFo~.,datapp,distribution="gaussian",interaction.depth=i)
erreur=c(erreur,RMSE(datapp$FFo,predict(boostreg.out,datapp)))
erreur2=c(erreur2,RMSE(datatest$FFo,predict(boostreg.out,datatest)))
}
plot(x,erreur,type="l",ylim=c(1,3),col="blue",main="Boosting - influence elagage",xlab="interaction.depth",ylab="RMSE (bleu -> appr , rouge -> test)")
lines(x,erreur2,type="l",col="red")

# --> les arbres doivent etre tres simples en boosting : interaction.depth=1, valeur par defaut convient


#Discrimination

data$OCC=as.factor(as.numeric(data$FFo>13))
data$OCCp=as.factor(as.numeric(data$FFp>13))


boostdis.out=gbm(as.numeric(data$OCC)-1~.,data[,-11],distribution="bernoulli",n.trees=500,cv.folds=10)
plot(boostdis.out$cv.error)
gbm.perf(boostdis.out)

roc.plot(as.numeric(data$OCC)-1,predict(boostdis.out,type="response"))
scores(predict(boostdis.out,type="response")>0.12,data$OCC)

nappr=ceiling(0.5*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]


# influence de l'elagage - interaction.depth :

erreur=NULL
erreur2=NULL
x=seq(1,30,1)
for (i in x) {
boostdis.out=gbm(as.numeric(datapp$OCC)-1~.,datapp[,-11],distribution="bernoulli",n.trees=200,interaction.depth=i)
erreur=c(erreur,brier(as.numeric(datapp$OCC)-1,predict(boostdis.out,type="response"))$bs)
erreur2=c(erreur2,brier(as.numeric(datatest$OCC)-1,predict(boostdis.out,datatest,type="response"))$bs)
}
plot(x,erreur,type="l",ylim=c(0,0.1),col="blue",main="Boosting - influence elagage",xlab="interaction.depth",ylab="BS (bleu -> appr , rouge -> test)")
lines(x,erreur2,type="l",col="red")


# influence du shrinkage :

erreur=NULL
erreur2=NULL
x=seq(0.01,2,0.1)
for (i in x) {
boostdis.out=gbm(as.numeric(datapp$OCC)-1~.,datapp[,-11],distribution="bernoulli",cv.folds=10,n.trees=6000,shrinkage=i)
erreur=c(erreur,brier(as.numeric(datapp$OCC)-1,predict(boostdis.out,type="response"))$bs)
erreur2=c(erreur2,brier(as.numeric(datatest$OCC)-1,predict(boostdis.out,datatest,type="response"))$bs)
}
plot(x,erreur,type="l",col="blue",main="Boosting - influence du shrinkage",xlab="shrinkage",ylab="BS (bleu -> appr , rouge -> test)")
lines(x,erreur2,type="l",col="red")



########################################################################################################
# POUR LA RENTREE :
# --> Confronter les techniques d'agregation aux methodes des autres TP en regression puis discrimination
########################################################################################################


