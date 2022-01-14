# Syntaxe par methode pour exploiter les previsions sur apprentissage puis test



#################
# En regression :
#################


# Modele lineaire gaussien

predict(lm.out) ou fitted(lm.out)
predict(lm.out,datatest)


# GAM (modele gaussien)

predict(gam.out) ou fitted(gam.out)
predict(gam.out,datatest)


# Arbres binaires 

predict(rpart.out)
predict(rpart.out,datatest)


# Bagging

predict(bag.out,datapp) #Attention pour bagging, si on ne precise pas le tableau datapp, predict par defaut retourne les previsions OOB.
predict(bag.out,datatest)


# Random forest

predict(rf.out,datapp) #Attention pour randomForest, si on ne precise pas le tableau datapp, predict par defaut retourne les previsions OOB.
predict(rf.out,datatest)


# Boosting

predict(gbm.out)
predict(gbm.out,datatest)


# SVM

predict(svm.out)
predict(svm.out,datatest)




#-----------------------------------------------------------------------




###########################################################
# En discrimination (probabilite du depassement de seuil) :
###########################################################


# Regression logistique

fitted(glm.out) ou predict(glm.out,type="response")
predict(glm.out,datatest,type="response")


# Analyses discriminantes lineaire et quadratique

predict(lda.out)$posterior[,2]
predict(lda.out,datatest)$posterior[,2]

predict(qda.out)$posterior[,2]
predict(qda.out,datatest)$posterior[,2]


# GAM (modele logistique)

fitted(gam.out) ou predict(gam.out,type="response")
predict(gam.out,datatest,type="response")


# Arbres binaires

predict(rpart.out)[,2]
predict(rpart.out,datatest)[,2]


# Bagging

predict(bag.out,datapp,type="prob")[,2] #Attention pour bagging, si on ne precise pas le tableau datapp, predict par defaut retourne les previsions OOB.
predict(bag.out,datatest,type="prob")[,2]


# Random forest

predict(rf.out,datapp,type="prob")[,2] #Attention pour randomForest, si on ne precise pas le tableau datapp, predict par defaut retourne les previsions OOB.
predict(rf.out,datatest,type="prob")[,2]


# Boosting

predict(gbm.out,type="response")
predict(gbm.out,datatest,type="response")


# SVM (apres avoir entraine le modele avec prob=T)

attributes(predict(svm.out,datapp,prob=T))$probabilities[,2]
attributes(predict(svm.out,datatest,prob=T))$probabilities[,2]


