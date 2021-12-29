# Chargement des librairies et des données
data = read.table('DataTP.txt', header = TRUE)
summary(data)
dim(data)       # 11 attributes 2180 observations
sd(data$FFo)
sd(data$FFp)

# Histograms of FFo and FFp
X11()
par(mfrow=c(1,2))
hist(data$FFo, breaks=seq(0,40,2), ylim=c(0,500))
hist(data$FFp, breaks=seq(0,40,2), ylim=c(0,500))

# Boxplots of FFo and FFp
X11()
boxplot(data$FFo, data$FFp, names=c("FFo","FFp"))   # generally FFp underestimates FFo

var.test(data$FFo,data$FFp)     # the variance of FFp and FFo are different
t.test(data$FFo,data$FFp)       # the means of FFp and FFo are different

X11()
plot(data$FFo[300:400], type ="l", main="Force de vent",xlab="Date",ylab="[FFo]")
points(data$FFp[300:400], col="blue", pch="+")
points(fitted(lm(FFo~FFp, data))[300:400], col="red", pch="+")
legend(0, 25, lty=1, col=c("black"), legend=c("FFo"), bty="n")
legend(0, 23, pch="+", col="blue", legend="FFp", bty="n") 
legend(0, 21, pch="+", col="red", legend="régression", bty="n")
# we can see that FFp underestimates FFo and the predication that are comming from the regression model  fits the data better

Biaslm = mean(data$FFo - fitted(lm(FFo~FFp, data)))
RMSElm = sqrt(mean((data$FFo - fitted(lm(FFo~FFp, data)))**2))
BiasFFp = mean(data$FFo - data$FFp)
RMSEFFP = sqrt(mean((data$FFo - data$FFp)**2))

cat("Bias of FFp :", BiasFFp, "\n")
cat("Bias of lm model using FFp  :", Biaslm, "\n")
cat("RMSE of FFp :", RMSEFFP, "\n")
cat("RMSE of lm model using FFp  :", RMSElm, "\n")
# we can see that we get a lower BIAS and RMSE using a regression model

