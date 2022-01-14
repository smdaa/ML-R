# ScriptTP3 HPC / Big Data 2021

# Chargement des packages
library(verification)
library(rpart)
library(partykit)

source("scores.R")

eval=function(obs,prev) {
rmse=sqrt(mean((prev-obs)**2))
biais=mean(prev-obs)
print("Biais       RMSE") 
return(c(biais,rmse))
}



###############################
# Arbres binaires
###############################

#Regression

data=read.table(file="DataTP.txt",header=TRUE)

rpartreg.out=rpart(FFo~.,data,cp=0)# Arbre max
rpartreg.out2=rpart(FFo~.,data)# Arbre elague aveccp=0.01 par defaut
summary(rpartreg.out)
print(rpartreg.out)
plot(rpartreg.out)
text(rpartreg.out)
eval(data$FFo,predict(rpartreg.out))

plot(rpartreg.out2)
text(rpartreg.out2)
plot(as.party(rpartreg.out2))
eval(data$O3o,predict(rpartreg.out2))


#Elagage
printcp(rpartreg.out)
plotcp(rpartreg.out)
rpartreg.out3=rpart(FFo~.,data,cp=0.003)
plot(as.party(rpartreg.out3))
eval(data$O3o,predict(rpartreg.out3))

# --> comparer performances arbre maximal vs arbre elague sur apprentissage puis test
# --> confronter l'arbre elague au modele lineaire GAM, meilleur modele de regression obtenu jusqu'a present


#################################################


# Discrimination

data=read.table(file="DataTP.txt",header=TRUE)
data$OCC=as.factor(as.numeric(data$FFo>13))
data$OCCp=as.factor(as.numeric(data$FFp>13))
summary(data)

source("scores.R")


rpartdis.out=rpart(OCC~.,data[,-11],cp=0)
plot(rpartdis.out)
text(rpartdis.out)
plot(as.party(rpartdis.out))
roc.plot(as.numeric(data$OCC)-1,predict(rpartdis.out)[,2])
scores(predict(rpartdis.out,data)[,2]>0.2,data$OCC)
brier(as.numeric(data$OCC)-1,predict(rpartdis.out)[,2])$bs

#Elagage   
plotcp(rpartdis.out)
rpartdis.out2=rpart(OCC~.,data[,-11],cp=0.01)
plot(rpartdis.out2)
text(rpartdis.out2)
plot(as.party(rpartdis.out2))
roc.plot(as.numeric(data$OCC)-1,predict(rpartdis.out2)[,2])
scores(predict(rpartdis.out2)[,2]>0.1,data$OCC)
brier(as.numeric(data$OCC)-1,predict(rpartdis.out2)[,2])$bs

# --> comparer les performances arbre maximal vs arbre elague sur apprentissage puis test
# --> confronter l'arbre elague a la regression logistique et l'analyse discriminante 





