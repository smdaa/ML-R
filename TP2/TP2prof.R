# Script TP2 HPC / Big Data 2021

# Chargement des packages
library(verification)
library(gam)
library(akima)

source("scores.R")

eval=function(obs,prev) {
rmse=sqrt(mean((prev-obs)**2))
biais=mean(prev-obs)
print("Biais       RMSE") 
return(c(biais,rmse))
}


###############################
# Modeles GAM
###############################


#Regression

data=read.table(file="DataTP.txt",header=TRUE)

gamreg.out=gam(FFo~lo(heure)+lo(v)+lo(u)+lo(FFp)+lo(HU)+lo(N)+lo(P)+lo(hel)+lo(DD)+lo(mois),family=gaussian,data)
par(mfrow=c(2,5))
plot(gamreg.out,se=T,col="red")

summary(gamreg.out)
#Anova Parametric effects : lien lineaire significatif si pvalue faible
#Anova NonParametric effects : lien non lineaire a privilegier si pvalue faible

gamreg.out=gam(FFo~lo(heure)+lo(v)+u+lo(FFp)+P,family=gaussian,data)
summary(gamreg.out)
eval(data$FFo,predict(gamreg.out))

gamreg.out=gam(FFo~lo(heure)+lo(v)+lo(FFp),family=gaussian,data)
summary(gamreg.out)
eval(data$FFo,predict(gamreg.out))

#Introduction d'interactions 
gamreg.out=gam(FFo~lo(heure)+lo(v)+lo(FFp)+lo(heure,v)+lo(heure,FFp)+lo(v,FFp),family=gaussian,data)
par(mfrow=c(2,3))
plot(gamreg.out)
eval(data$FFo,predict(gamreg.out))

#Essayer avec les splines s() au lieu des regressions locales lo()
gamreg.out=gam(FFo~s(heure)+s(v)+s(FFp),family=gaussian,data)
summary(gamreg.out)
eval(data$FFo,predict(gamreg.out))

nappr=ceiling(0.8*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]

gamreg.out=gam(FFo~lo(heure)+lo(v)+lo(FFp),family=gaussian,datapp)
eval(datatest$FFo,predict(gamreg.out,datatest))

#--> Confronter regression GAM et modele lineaire gaussien (MLG) avec votre procedure d'evaluation



############################################


#Discrimination

data=read.table(file="DataTP.txt",header=TRUE)
data$OCC=as.factor(as.numeric(data$FFo>13))
data$OCCp=as.factor(as.numeric(data$FFp>13))
summary(data)

source("scores.R")

gamdis.out=gam(OCC~lo(FFp)+lo(heure)+lo(v)+lo(u)+OCCp+lo(HU)+lo(N)+lo(P)+lo(hel)+lo(DD)+lo(mois),data,family=binomial)
par(mfrow=c(2,5))
plot(gamdis.out,se=T,col="red")
summary(gamdis.out)
gamdis.out=gam(OCC~lo(v)+u+lo(FFp),data,family=binomial)

roc.plot(as.numeric(data$OCC)-1,predict(gamdis.out,type="response"))
scores(fitted(gamdis.out)>0.2,data$OCC)
brier(as.numeric(data$OCC)-1,fitted(gamdis.out))$bs
roc.area(as.numeric(data$OCC)-1,fitted(gamdis.out))$A

nappr=ceiling(0.8*nrow(data))
ii=sample(1:nrow(data),nappr)
jj=setdiff(1:nrow(data),ii)
datatest=data[jj,]
datapp=data[ii,]

gamdis.out=gam(OCC~lo(v)+lo(FFp),datapp,family=binomial)

scores(predict(gamdis.out,datatest,type="response")>0.2,datatest$OCC)
brier(as.numeric(datatest$OCC)-1,predict(gamdis.out,datatest,type="response"))$bs
roc.area(as.numeric(datatest$OCC)-1,predict(gamdis.out,datatest,type="response"))$A


#--> Comparer les modeles GAM avec regression logistique et analyse discriminante


