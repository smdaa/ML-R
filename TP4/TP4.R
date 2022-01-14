library(verification)
library(rpart)
library(ipred)
library(randomForest)
library(gbm)

source("scores.R")

eval=function(obs,prev) {
  rmse=sqrt(mean((prev-obs)**2))
  biais=mean(prev-obs)
  print("Biais       RMSE") 
  return(c(biais,rmse))
}

RMSE=function(obs,pr){
  return(sqrt(mean((pr-obs)^2)))
}

data=read.table("DataTP.txt",header=TRUE)

# bagging
bagreg.out = bagging(FFo ~ . , data, coob=T)
eval(data$FFo,predict(bagreg.out,data))

# Analysis l’influence du nombre d’arbres
erreur = c()
x=seq(10,300,10)
for (i in x) {
  erreur = c(erreur, bagging(FFo~.,data,nbagg=i,coob=T)$err)  
}

X11()
plot(x, erreur, type="l", ylim=c(2.7,2.8), xlab="nbagg", 
     ylab="Erreur OOB",main="Bagging : influence du nb d'arbres")

# Stabilisation de l'erreur OOB a partir de nbag=100

nappr=ceiling(0.5*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]

erreur  = c()
erreur2 = c()
x=seq(0, 0.02, 0.001)

for (i in x) {
  bag.out = bagging(FFo~., datapp, nbag=100, control=rpart.control(cp=i))
  erreur = c(erreur, RMSE(datapp$FFo, predict(bag.out, datapp)))
  erreur2 = c(erreur2, RMSE(datatest$FFo, predict(bag.out, datatest)))
}

X11()
plot(x,erreur,type="l",ylim=c(1.6,3),col="blue",main="Bagging : influence de l'elagage",xlab="cp",ylab="RMSE (bleu -> appr , rouge -> test)")
lines(x,erreur2,type="l",col="red")

# cp=0.01 (valeur par defaut) convient ici, surapprentissage severe si arbres peu elages
  
# Random Forest

rfreg.out = randomForest(FFo ~ . , data, importance=T)
# Le modèle obtenu est ininterprétable mais des coefficients estiment les contributions
# des différents prédicteurs
sort(importance(rfreg.out)[,1],dec=T)
eval(data$FFo,predict(rfreg.out,data))    # RMSE = 1.14 good score

erreur = c()
x=seq(10,600,30)
for (i in x) {
  rf.out = randomForest(FFo~.,data,ntree=i)
  erreur = c(erreur, RMSE(data$FFo, predict(rf.out))) # par defaut predict renvoie les previsions OOB
}

X11()
plot(x,erreur,type="l",ylim=c(2.5,2.7),xlab="ntree",main="RandomForest : influence du nb d'arbres",ylab="RMSE")

# Stabilisation de l'erreur OOB a partir de ntree=200 ou 300

erreur = c()
x=1:10
for (i in x) {
  rf.out=randomForest(FFo~.,data,ntree=200,mtry=i)
  erreur=c(erreur,RMSE(data$FFo,predict(rf.out)))
}

X11()
plot(x,erreur,type="l",main="RandomForest - influence de mtry",xlab="mtry",ylab="RMSE")
# minimum vers 2-3, le choix par defaut E(q/3)=3 est bon ici


nappr=ceiling(0.5*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]

erreur  = c()
erreur2 = c()
x=seq(10, 200, 10)
for (i in x) {
  rf.out=randomForest(FFo~.,datapp,ntree=200,maxnodes=i)
  erreur=c(erreur,RMSE(datapp$FFo,predict(rf.out,datapp)))
  erreur2=c(erreur2,RMSE(datatest$FFo,predict(rf.out,datatest)))
}

X11()
plot(x,erreur,type="l",ylim=c(1,3),col="blue",main="RandomForest - influence du nb de feuilles",xlab="maxnodes",ylab="RMSE (bleu -> appr , rouge -> test)")
lines(x,erreur2,type="l",col="red")

# qq dizaines suffisent sinon surapprentissage severe


# Boosting
boostreg.out = gbm(FFo ~ . , data, distribution = "gaussian", n.trees = 5000, cv.folds = 10)

X11()
plot(boostreg.out$cv.error)

X11()
gbm.perf(boostreg.out)

boostreg.out = gbm(FFo~.,data,distribution="gaussian",n.trees = 250)
eval(data$FFo, predict(boostreg.out))

nappr=ceiling(0.5*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]


# influence du shrinkage :

erreur=NULL
erreur2=NULL
x=seq(0.01,2,0.1)
for (i in x) {
  boostreg.out=gbm(FFo~.,datapp,distribution="gaussian",n.trees=6000,cv.folds=10,shrinkage=i)
  erreur=c(erreur,RMSE(datapp$FFo,predict(boostreg.out)))
  erreur2=c(erreur2,RMSE(datatest$FFo,predict(boostreg.out,datatest)))
}

X11()
plot(x,erreur,type="l",col="blue",ylim=c(2.4,3.5),main="Boosting - influence du shrinkage",xlab="shrinkage",ylab="RMSE (bleu -> appr , rouge -> test)")
lines(x,erreur2,type="l",col="red")

