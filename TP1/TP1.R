#--------------------------------------------------------------#
#           Chargement des librairies et des données           #
#--------------------------------------------------------------#
library(MASS)
data=read.table('DataTP.txt', header = TRUE)

#--------------------------------------------------------------#
#                    Etude préliminaire                        #
#--------------------------------------------------------------#
summary(data)
dim(data)

sd(data$FFo)
sd(data$FFp)

X11()
par(mfrow=c(1,2))
hist(data$FFo,breaks=seq(0,40,2),ylim=c(0,500))
hist(data$FFp,breaks=seq(0,40,2),ylim=c(0,500))

X11()
boxplot(data$FFo,data$FFp,names=c("FFo","FFp"))

var.test(data$FFo,data$FFp)
t.test(data$FFo,data$FFp)
# -> les variances et les moyennes sont jugees significativement differentes

X11()
plot(data$FFo, type ="l", main="Force de vent",xlab="Date",ylab="[FFo]")
points(data$FFp, col="blue", pch="+")
points(fitted(lm(FFo~FFp, data)), col="red", pch="+")

Biais = mean(fitted(lm(FFo~FFp, data))) - mean(data$FFo)
Biais

RMSE = sqrt(mean((fitted(lm(FFo~FFp, data))) - data$FFo)**2)
RMSE

#--------------------------------------------------------------#
#           Régression - Modèle linéaire gaussien              #
#--------------------------------------------------------------#
model1 = lm(FFo~., data)
summary(model1)
X = model.matrix(model1)
X[300:310,]

X11()
par(mfrow=c(2,2))
plot(fitted(model1), residuals(model1), main="Hypothese d'homoscedasticite",xlab="Valeurs ajustees (Y*)",ylab="Residus")
hist(residuals(model1))
qqnorm(residuals(model1))
acf(residuals(model1))

# -> On conserve les prédicateurs HU P v heure FFp

model2 = lm(FFo~.*., data)
summary(model2)

X11()
par(mfrow=c(2,2))
plot(fitted(model2), residuals(model2), main="Hypothese d'homoscedasticite",xlab="Valeurs ajustees (Y*)",ylab="Residus")
hist(residuals(model2))
qqnorm(residuals(model2))
acf(residuals(model2))

model1AIC = stepAIC(model1)
model2AIC = stepAIC(model2)
model1BIC = stepAIC(model1, k=log(nrow(data)))
model2BIC = stepAIC(model2, k=log(nrow(data)))


formula(model1)
formula(model1AIC)
formula(model1BIC)

formula(model2)
formula(model2AIC)
formula(model2BIC)

summary(model1)
summary(model1AIC)
summary(model1BIC)

summary(model2)
summary(model2AIC)
summary(model2BIC)

#--------------------------------------------------------------#
#                   Evaluation des modèles                     #
#--------------------------------------------------------------#
eval=function(obs, prev) {
rmse  = sqrt(mean((prev - obs)**2))
biais = mean(prev - obs)
print("Biais         RMSE") 
return(c(biais,rmse))
}

eval(data$FFo, data$FFp)
eval(data$FFo, fitted(model1BIC))
eval(data$FFo, fitted(model1AIC))
eval(data$FFo, fitted(model2BIC))
eval(data$FFo, fitted(model2AIC))

nappr = ceiling(0.8*nrow(data))
ii = sample(1:nrow(data), nappr)
jj = setdiff(1:nrow(data), ii)
datatest = data[jj,]
datapp = data[ii,]

model1AIC = lm(formula(model1AIC), datapp)
model1BIC = lm(formula(model1BIC), datapp)
model2AIC = lm(formula(model2AIC),datapp)
model2BIC = lm(formula(model2BIC),datapp)

eval(datapp$FFo, fitted(model1BIC))
eval(datapp$FFo, fitted(model1AIC))
eval(datapp$FFo, fitted(model2BIC))
eval(datapp$FFo, fitted(model2AIC))

eval(datapp$FFo, predict(model1BIC, datatest))
eval(datapp$FFo, predict(model1AIC, datatest))
eval(datapp$FFo, predict(model2BIC, datatest))
eval(datapp$FFo, predict(model2AIC, datatest))

# Validation croisée
# TODO

#--------------------------------------------------------------------#
#   Discrimination - Régression logistique et Analyse discriminante  #
#--------------------------------------------------------------------#

