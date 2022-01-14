library(MASS)

eval=function(obs, prev) { 
  rmse  = sqrt(mean((prev - obs)**2))
  biais = mean(prev - obs)
  return(c(biais,rmse))
}

data = read.table('DataTP.txt', header = TRUE)
lm.out = lm(FFo~., data)
lm.outAIC    = stepAIC(lm.out)
lm.outBIC    = stepAIC(lm.out, k=log(nrow(data)))

a.aic = NULL
a.bic = NULL

t.aic = NULL
t.bic = NULL
# Randomly shuffle the data
set.seed(0)
data = data[sample(nrow(data)),]

# Create 10 equally size folds
K = 10
folds = cut(seq(1,nrow(data)), breaks=K, labels=FALSE)

# Perform 10 fold cross validation
for(i in 1:K){
	testind  = which(folds==i, arr.ind=TRUE)
	datatest = data[testind,]
	datapp  = data[-testind, ]
	
	lm.outAIC = lm(formula(lm.outAIC), datapp)
	lm.outBIC = lm(formula(lm.outBIC), datapp)
	
	a.aic = c(a.aic, eval(datapp$FFo, fitted(lm.outAIC))[2])
	a.bic = c(a.bic, eval(datapp$FFo, fitted(lm.outBIC))[2])
	
	t.aic = c(t.aic, eval(datatest$FFo, predict(lm.outAIC, datatest))[2])
	t.bic = c(t.bic, eval(datatest$FFo, predict(lm.outBIC, datatest))[2])
	
}

X11()
boxplot(a.bic, a.aic, t.bic, t.aic, 
        col   = c("blue", "blue", "red", "red"), 
        names = c("a.bic", "a.aic", "t.bic", "t.aic"), 
        main  = "Modele lineaire gaussien - Score RMSE")
