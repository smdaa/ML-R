library("rpart")
library("partykit")
library("gam")
library("akima")
library("verification")
source("scores.R")

eval=function(obs,prev) {
  rmse=sqrt(mean((prev-obs)**2))
  biais=mean(prev-obs)
  #print("Biais       RMSE") 
  return(c(biais,rmse))
}

data = read.table("DataTP.txt",header=TRUE)
rpartreg.out  = rpart(FFo~., data, cp=0)   # Arbre max
rpartreg.out2 = rpart(FFo~., data)         # Arbre elague aveccp=0.01 par defaut

summary(rpartreg.out)
print(rpartreg.out)

X11()
plot(rpartreg.out)
text(rpartreg.out)

eval(data$FFo, predict(rpartreg.out))       # RMSE = 2.06

summary(rpartreg.out2)
print(rpartreg.out2)

X11()
plot(rpartreg.out2)
text(rpartreg.out2)

X11()
plot(as.party(rpartreg.out2))

eval(data$FFo, predict(rpartreg.out2))      # RMSE = 2.76

# Pruning of the over trained model.
printcp(rpartreg.out)

X11()
plotcp(rpartreg.out)

rpartreg.out3=rpart(FFo~.,data, cp=0.003)

X11()
plot(rpartreg.out3)
text(rpartreg.out3)

eval(data$FFo, predict(rpartreg.out3))      # RMSE = 2.61

# arbre maximal vs prunned vs gam
# K fold cross validation
gamreg.out = gam(FFo~ s(heure)+s(v)+s(FFp), family=gaussian, data)



set.seed(666)
K = 10
ntest = ceiling(nrow(data) / K)
endtest = 0
data_shuffle = data[sample(nrow(data)),]
a.tree    = c()
a.treemax = c()
a.gam     = c()
t.tree    = c()
t.treemax = c()
t.gam     = c()

for (i in 1:K) {
  starttest = endtest + 1
  endtest   = min(i * ntest, nrow(data))

  datatest  = data_shuffle[starttest:endtest,]
  datapp    = data_shuffle[setdiff(1:nrow(data), starttest:endtest),]
  
  rpartreg.out  = rpart(FFo~., datapp, cp=0)
  rpartreg.out3 = rpart(FFo~., datapp, cp=0.003)
  gamreg.out    = gam(formula(gamreg.out), datapp, family="gaussian")
  
  a.tree    = c(a.tree, eval(datapp$FFo, predict(rpartreg.out3))[2])
  a.treemax = c(a.treemax, eval(datapp$FFo, predict(rpartreg.out))[2])
  a.gam     = c(a.gam, eval(datapp$FFo, predict(gamreg.out))[2])
  
  t.tree    = c(t.tree, eval(datatest$FFo, predict(rpartreg.out3, datatest))[2])
  t.treemax = c(t.treemax, eval(datatest$FFo, predict(rpartreg.out, datatest))[2])
  t.gam     = c(t.gam, eval(datatest$FFo, predict(gamreg.out, datatest))[2])
  
}

X11()
boxplot(a.gam, a.tree, a.treemax, t.gam, t.tree, t.treemax, 
        col   = c("blue", "blue", "blue", "red", "red", "red"),
        names = c("a.gam", "a.tree", "a.treemax", "t.gam", "t.tree", "t.treemax"),
        main  = "Regression : Tree vs GAM Score RMSE")

# gam est tjrs plus robuste

data$OCC = as.factor(as.numeric(data$FFo>13))
data$OCCp = as.factor(as.numeric(data$FFp>13))

rpartdis.out=rpart(OCC~.,data[,-11],cp=0)

X11()
plot(rpartdis.out)
text(rpartdis.out)

X11()
roc.plot(as.numeric(data$OCC)-1, predict(rpartdis.out)[,2])

scores(predict(rpartdis.out,data)[,2]>0.2, data$OCC)        # PSS = 0.76
brier(as.numeric(data$OCC)-1,predict(rpartdis.out)[,2])$bs  # bs  = 0.06 

# pruning 

X11()
plotcp(rpartdis.out)

rpartdis.out2=rpart(OCC~.,data[,-11],cp=0.01)

X11()
plot(rpartdis.out2)
text(rpartdis.out2)

X11()
roc.plot(as.numeric(data$OCC)-1,predict(rpartdis.out2)[,2])

scores(predict(rpartdis.out2)[,2]>0.1,data$OCC)             # PSS = 0.65
brier(as.numeric(data$OCC)-1,predict(rpartdis.out2)[,2])$bs # bs  = 0.085

# arbre maximal vs prunned vs logistic regression
# K fold cross validation
glm.out    = glm(OCC~., data[,-11], family=binomial)
glm.outAIC = stepAIC(glm.out)

set.seed(666)
K = 10
ntest = ceiling(nrow(data) / K)
endtest = 0
data_shuffle = data[sample(nrow(data)),]
a.reglog  = c()
a.tree    = c()
a.treemax = c()
t.reglog  = c()
t.tree    = c()
t.treemax = c()

for (i in 1:K) {
  starttest = endtest + 1
  endtest   = min(i * ntest, nrow(data))
  
  datatest  = data_shuffle[starttest:endtest,]
  datapp    = data_shuffle[setdiff(1:nrow(data), starttest:endtest),]
  
  rpartdis.out  = rpart(OCC~., datapp[,-11], cp=0)
  rpartdis.out2 = rpart(OCC~., datapp[,-11], cp=0.01)
  glm.outAIC    = glm(formula(glm.outAIC), datapp[,-11], family=binomial)
  
  a.tree    = c(a.tree,    brier(as.numeric(datapp$OCC)-1, predict(rpartdis.out2)[,2])$bs)
  a.treemax = c(a.treemax, brier(as.numeric(datapp$OCC)-1, predict(rpartdis.out)[,2])$bs)
  a.reglog  = c(a.reglog,  brier(as.numeric(datapp$OCC)-1, fitted(glm.outAIC))$bs)
  
  t.tree    = c(t.tree,    brier(as.numeric(datatest$OCC)-1, predict(rpartdis.out2, datatest)[,2]>0.2)$bs)
  t.treemax = c(t.treemax, brier(as.numeric(datatest$OCC)-1, predict(rpartdis.out, datatest)[,2]>0.2)$bs)
  t.reglog  = c(t.reglog,  brier(as.numeric(datatest$OCC)-1, predict(glm.outAIC, datatest, type="response")>0.2)$bs)
  
}

X11()
boxplot(a.reglog, a.tree, a.treemax, t.reglog, t.tree, t.treemax, 
        col   = c("blue", "blue", "blue", "red", "red", "red"),
        names = c("a.reglog", "a.tree", "a.treemax", "t.reglog", "t.tree", "t.treemax"),
        main  = "Arbres binaires vs Regression Logistique BS")




