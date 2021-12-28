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
scores(data$OCCp, data$OCC)                 # PSS = 0.15
scores(fitted(glm.outAIC)>0.5, data$OCC)    # PSS = 0.57

X11()
roc.plot(as.numeric(data$OCC)-1,fitted(glm.outAIC))
# On choisit le seuil de probailité 0.2 qui permet d'atteindre un hit rate de 0.9 et un false rate de 0.1
scores(fitted(glm.outAIC)>0.2, data$OCC)    # PSS = 0.69
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
#  5-4 validation croisée
# We split Data into K chunks
K = 4
set.seed(10)
rows = sample(nrow(data))
datashuffle = data[rows,]
ntest = ceiling((1/K) * nrow(data))
endtest = 0

a.aic = matrix(,nrow=K,ncol=3)
a.bic = matrix(,nrow=K,ncol=3)
a.aint = matrix(,nrow=K,ncol=3)
a.bint = matrix(,nrow=K,ncol=3)

t.aic = matrix(,nrow=K,ncol=3)
t.bic = matrix(,nrow=K,ncol=3)
t.aint = matrix(,nrow=K,ncol=3)
t.bint = matrix(,nrow=K,ncol=3)

for (i in 1:K) {
    starttest = endtest + 1
    endtest   = min(i * ntest, nrow(data)) 

    datatest = datashuffle[starttest:endtest,]
    datapp   = datashuffle[setdiff(1:nrow(data), starttest:endtest),]

    glm.outAIC = glm(formula(glm.outAIC), datapp, family=binomial)
    glm.outBIC = glm(formula(glm.outBIC), datapp, family=binomial)
    glm.outintAIC = glm(formula(glm.outintAIC), datapp, family=binomial)
    glm.outintBIC = glm(formula(glm.outintBIC), datapp, family=binomial)

    t.aic[i, 1] = scores(predict(glm.outAIC, datatest, type="response")>0.2, datatest$OCC, quiet=TRUE)[4]
    t.aic[i, 2] = brier(as.numeric(datatest$OCC)-1, predict(glm.outAIC, datatest, type="response")>0.2)$bs
    t.aic[i, 3] = roc.area(as.numeric(datatest$OCC)-1, predict(glm.outAIC, datatest, type="response"))$A

    t.bic[i, 1] = scores(predict(glm.outBIC, datatest, type="response")>0.2, datatest$OCC, quiet=TRUE)[4]
    t.bic[i, 2] = brier(as.numeric(datatest$OCC)-1, predict(glm.outBIC, datatest, type="response")>0.2)$bs
    t.bic[i, 3] = roc.area(as.numeric(datatest$OCC)-1, predict(glm.outBIC, datatest, type="response"))$A

    t.aint[i, 1] = scores(predict(glm.outintAIC, datatest, type="response")>0.2, datatest$OCC, quiet=TRUE)[4]
    t.aint[i, 2] = brier(as.numeric(datatest$OCC)-1, predict(glm.outintAIC, datatest, type="response")>0.2)$bs
    t.aint[i, 3] = roc.area(as.numeric(datatest$OCC)-1, predict(glm.outintAIC, datatest, type="response"))$A

    t.bint[i, 1] = scores(predict(glm.outintBIC, datatest, type="response")>0.2, datatest$OCC, quiet=TRUE)[4]
    t.bint[i, 2] = brier(as.numeric(datatest$OCC)-1, predict(glm.outintBIC, datatest, type="response")>0.2)$bs
    t.bint[i, 3] = roc.area(as.numeric(datatest$OCC)-1, predict(glm.outintBIC, datatest, type="response"))$A

    a.aic[i, 1] = scores(fitted(glm.outAIC)>0.2, datapp$OCC, quiet=TRUE)[4]
    a.aic[i, 2] = brier(as.numeric(datapp$OCC)-1, fitted(glm.outAIC))$bs 
    a.aic[i, 3] = roc.area(as.numeric(datapp$OCC)-1, fitted(glm.outAIC))$A 

    a.bic[i, 1] = scores(fitted(glm.outBIC)>0.2, datapp$OCC, quiet=TRUE)[4]
    a.bic[i, 2] = brier(as.numeric(datapp$OCC)-1, fitted(glm.outBIC))$bs 
    a.bic[i, 3] = roc.area(as.numeric(datapp$OCC)-1, fitted(glm.outBIC))$A 

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
        main  = "Regression logistique  Score PSSE")

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

# Conclusion
# aic is better
# 6- Analyse discriminante linéaire et quadratique
library(MASS)
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

X11()
pairs(data[iNOCC,c(1,3,4,5,6,10)])

X11()
pairs(data[iOCC,c(1,3,4,5,6,10)])
#  6.1- Estimer plusieurs modèles d’analyse discriminante, linéaire puis quadratique
lda.out = lda(OCC~.,data[,-11]) 

X11()
roc.plot(as.numeric(data$OCC)-1, predict(lda.out)$posterior[,2])

scores(predict(lda.out)$posterior[,2]>0.2,data$OCC)

lda.out=lda(OCC~FFp+v,data[,-11]) 

X11()
roc.plot(as.numeric(data$OCC)-1,predict(lda.out)$posterior[,2])

scores(predict(lda.out)$posterior[,2]>0.2,data$OCC)
