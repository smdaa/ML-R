library("gam")
library("akima")
library("verification")
data=read.table("DataTP.txt",header=TRUE)
source("scores.R")
eval=function(obs,prev) {
  rmse=sqrt(mean((prev-obs)**2))
  biais=mean(prev-obs)
  return(c(biais,rmse))
}

gamreg.out = gam(FFo~ lo(heure)+lo(v)+lo(u)+lo(FFp)+lo(HU)+lo(N)+lo(P)+lo(hel)+lo(DD)+lo(mois), data, family=gaussian)
summary(gamreg.out)
eval(data$FFo,predict(gamreg.out))        # RMSE : 2.55

X11()
par(mfrow=c(2,5))
plot(gamreg.out,se=T,col="red")

#Anova Parametric effects    : lien lineaire significatif si pvalue faible
#Anova NonParametric effects : lien non lineaire a privilegier si pvalue faible

gamreg.out = gam(FFo~ lo(heure)+lo(v)+u+lo(FFp)+P, family=gaussian, data)
summary(gamreg.out)
eval(data$FFo,predict(gamreg.out))        # RMSE : 2.59

gamreg.out = gam(FFo~ lo(heure)+lo(v)+lo(FFp),family=gaussian,data)
summary(gamreg.out)
eval(data$FFo,predict(gamreg.out))

# interactions
gamreg.out=gam(FFo ~ lo(heure)+lo(v)+lo(FFp)+lo(heure,v)+lo(heure,FFp)+lo(v,FFp), family=gaussian, data)

X11()
par(mfrow=c(2,3))
plot(gamreg.out)

eval(data$FFo,predict(gamreg.out))    # RMSE : 2.58

gamreg.out = gam(FFo~ s(heure)+s(v)+s(FFp), family=gaussian, data)
summary(gamreg.out)
eval(data$FFo,predict(gamreg.out))    # RMSE : 2.60

# Cross validation GAM vs MLG
lm.outBIC = lm(FFo~., data)
lm.outintBIC = lm(FFo~.*., data)
lm.outBIC = stepAIC(lm.outBIC, k=log(nrow(data)))
lm.outBIC = stepAIC(lm.outBIC, k=log(nrow(data)))

K = 10
endtest = 0
ntest = ceiling((1/K) * nrow(data))

a.gam = c()
a.bic = c()
a.bint = c()
t.gam = c()
t.bic = c()
t.bint = c()

set.seed(0)
datasuffle = data[sample(nrow(data)),]

for (i in 1:K) {
  starttest = endtest + 1
  endtest   = min(nrow(data), i*ntest)
  
  datatest = datasuffle[starttest:endtest,]
  datapp   = datasuffle[setdiff(1:nrow(data), starttest:endtest),]
  
  lm.outBIC = lm(formula(lm.outBIC), datapp)
  lm.outintBIC = lm(formula(lm.outintBIC), datapp)
  gamreg.out = gam(formula(gamreg.out), datapp, family="gaussian")
  
  a.bic = c(a.bic, eval(datapp$FFo, fitted(lm.outBIC))[2])
  a.bint = c(a.bint, eval(datapp$FFo, fitted(lm.outintBIC))[2])
  a.gam = c(a.gam, eval(datapp$FFo, fitted(gamreg.out))[2])
  
  t.bic = c(t.bic, eval(datatest$FFo, predict(lm.outBIC, datatest))[2])
  t.bint = c(t.bint, eval(datatest$FFo, predict(lm.outintBIC, datatest))[2])
  t.gam = c(t.gam, eval(datatest$FFo, predict(gamreg.out, datatest))[2])
}

X11()
boxplot(a.gam, a.bic, a.bint, t.gam, t.bic, t.bint, 
        col   = c("blue", "blue", "blue", "red", "red", "red"), 
        names = c("a.gam", "a.bic", "a.bint", "t.gam", "t.bic", "t.bint"),
        main  = "Modele lineaire gaussien - Score RMSE")





