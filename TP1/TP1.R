#--------------------------------------------------------------------------------------------#
#           Modèle linéaire gaussien - Régression logistique - Analyse discriminante         #
#--------------------------------------------------------------------------------------------#

# 1- Chargement des librairies et des données
library(MASS)
data = read.table('DataTP.txt', header = TRUE)

# 2- Etude préliminaire
#  2.1- analyse des données
summary(data)
dim(data)
sd(data$FFo)
sd(data$FFp)

X11()
par(mfrow=c(1,2))
hist(data$FFo, breaks=seq(0,40,2), ylim=c(0,500))
hist(data$FFp, breaks=seq(0,40,2), ylim=c(0,500))

X11()
boxplot(data$FFo, data$FFp, names=c("FFo","FFp"))

var.test(data$FFo,data$FFp)
t.test(data$FFo,data$FFp)
# -> les variances et les moyennes sont jugees significativement differentes

#  2.2- la série de force de vent observée FFO et les prévisions FFp et la régression de FFo par FFp
X11()
plot(data$FFo[300:400], type ="l", main="Force de vent",xlab="Date",ylab="[FFo]")
points(data$FFp[300:400], col="blue", pch="+")
points(fitted(lm(FFo~FFp, data))[300:400], col="red", pch="+")
legend(0, 25, lty=1, col=c("black"), legend=c("FFo"), bty="n")
legend(0, 23, pch="+", col="blue", legend="FFp", bty="n") 
legend(0, 21, pch="+", col="red", legend="régression", bty="n")
# Les prévisions issues du modèle de régression correspondent mieux à la FFO

#  2.2- le Biais et RMSE
eval=function(obs, prev) {
# Cet indice fournit une indication par rapport à la dispersion ou la variabilité 
# de la qualité de la prédiction    
rmse  = sqrt(mean((prev - obs)**2))
# Le biais permet d’évaluer si les prédictions sont précises ou non et 
# si le modèle a tendance à sur- ou sous-estimer les valeurs de la variable d’intérêt
# Cet indicateur ne prend pas en compte la variabilité des prédictions
biais = mean(prev - obs)
return(c(biais,rmse))
}
eval(data$FFo, fitted(lm(FFo~FFp, data)))

# 3- Régression - Modèle linéaire gaussien
#  3.1- Estimer le modèle linéaire gaussien exploitant l'ensemble des prédicteurs potentiels
lm.out = lm(FFo~., data)
#  3.1- afficher le bilan de l’estimation
summary(lm.out)
#  3.1- les graphiques des diagnostics
X11()
par(mfrow=c(2,2))
plot(fitted(lm.out), residuals(lm.out), main="Hypothese d'homoscedasticite", xlab="Valeurs ajustees (Y*)", ylab="Residus")
hist(residuals(lm.out))
qqnorm(residuals(lm.out))
acf(residuals(lm.out))
# les prédicteurs qu'on conserve : HU P v heure FFp

#  3.2- Estimer le modèle exploitant l'ensemble des prédicteurs mais également toutes les interactions d’ordre 2 possibles
lm.outint=lm(FFo~.*.,data)
#  3.2- afficher le bilan de l’estimation
summary(lm.outint)
#  3.2- les graphiques des diagnostics
X11()
par(mfrow=c(2,2))
plot(fitted(lm.outint),residuals(lm.outint),main="Hypothese d'homoscedasticite",xlab="Valeurs ajustees (Y*)",ylab="Residus")
hist(residuals(lm.outint))
qqnorm(residuals(lm.outint))
acf(residuals(lm.outint))
#  3.2- Comparez au précédent modèle
eval(data$FFo, fitted(lm.out))
eval(data$FFo, fitted(lm.outint))
# le model ouint a un facteur R^2 supérieur a out (en effet R^2 augmente toujours avec la complexité du model)
# on obtient un RMSE inférieur avec le model outint
#  3.3- Sélection automatique des prédicteurs
lm.outAIC    = stepAIC(lm.out)
lm.outintAIC = stepAIC(lm.outint)
lm.outBIC    = stepAIC(lm.out, k=log(nrow(data)))
lm.outintBIC = stepAIC(lm.outint, k=log(nrow(data)))
# L'AIC pénalise le nombre de paramètres moins fortement que le BIC

# 4- Evaluation des modèles
#  4.1- Evaluation des modèles
eval(data$FFo, data$FFp)                 # RMSE : 4.44
eval(data$FFo, fitted(lm.outBIC))        # RMSE : 2.69
eval(data$FFo, fitted(lm.outAIC))        # RMSE : 2.68
eval(data$FFo, fitted(lm.outintBIC))     # RMSE : 2.61
eval(data$FFo, fitted(lm.outintAIC))     # RMSE : 2.58
# outintAIC > outintBIC > outAIC > outBIC mais on peut pas conclure sur la robustesse des models car on
# doit tester sur des données etrangers au model
#  4.2- créer un fichier d’apprentissage datapp et un fichier de test datatest avec les données restantes
nappr = ceiling(0.8*nrow(data))
ii = sample(1:nrow(data),nappr)
jj = setdiff(1:nrow(data),ii)
datatest = data[jj,]
datapp   = data[ii,]
#  4.3- Réestimer les 4 modèles précédents
lm.outAIC    = lm(formula(lm.outAIC), datapp)
lm.outBIC    = lm(formula(lm.outBIC), datapp)
lm.outintAIC = lm(formula(lm.outintAIC), datapp)
lm.outintBIC = lm(formula(lm.outintBIC), datapp)
eval(datapp$FFo,fitted(lm.outBIC))              # RMSE : 2.69   
eval(datapp$FFo,fitted(lm.outAIC))              # RMSE : 2.68
eval(datapp$FFo,fitted(lm.outintBIC))           # RMSE : 2.60
eval(datapp$FFo,fitted(lm.outintAIC))           # RMSE : 2.58
# outintAIC > outintBIC > outAIC > outBIC
#  4.4- Validation croisée
# We split Data into K chunks
K = 4
set.seed(42)
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
        main  = "Modele lineaire gaussien − Score RMSE")
#  4.5- le meilleur modèle
# Les modèles sans interactions aic et bic sont très robustent mais les 
# performances des modèles avec interactions sont significativement meilleures.
# L’écart constaté sur test entre les 2 modèles avec interactions n’est pas jugé
# significatif on choisit donc le modele le plus simple BIC avec interaction
#  4.6- Illustrer le phénomène de sur-apprentissage.
# voir lm.INT
#  4.7- la série de force de vent observée FFO et les prévisions FFp et la régression de FFo par FFp
X11()
plot(data$FFo[300:400], type ="l", main="Force de vent",xlab="Date",ylab="[FFo]")
points(data$FFp[300:400], col="blue", pch="+")
points(fitted(lm.outintBIC)[300:400], col="red", pch="+")
legend(0, 25, lty=1, col=c("black"), legend=c("FFo"), bty="n")
legend(0, 23, pch="+", col="blue", legend="FFp", bty="n") 
legend(0, 21, pch="+", col="red", legend="régression intbic", bty="n")

# 5- Discrimination - Régression logistique et Analyse discriminante
library(MASS)
library(verification)

#  5-1 Ajouter à la data.frame data deux nouvelles variables, OCC et OCCp
data = read.table(file="DataTP.txt",header=TRUE)
data$OCC = as.factor(as.numeric(data$FFo>13))
data$OCCp = as.factor(as.numeric(data$FFp>13))
#  5-2 Exécuter le script scores.R
source("scores.R")
#  5-3 Régression logistique
glm.out    = glm(OCC~., data[,-11], family=binomial)
glm.outint = glm(OCC~.*.,data[,-11],family=binomial)
glm.outAIC = stepAIC(glm.out)
glm.outBIC = stepAIC(glm.out,k=log(nrow(data)))
glm.outintAIC = stepAIC(glm.outint)
glm.outintBIC = stepAIC(glm.outint,k=log(nrow(data)))

formula(glm.outAIC)
formula(glm.outBIC)
formula(glm.outintAIC)
formula(glm.outintBIC)
scores(data$OCCp, data$OCC)
scores(fitted(glm.outAIC)>0.5, data$OCC)

X11()
roc.plot(as.numeric(data$OCC)-1,fitted(glm.outAIC))
