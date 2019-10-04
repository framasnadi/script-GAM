#####GAM###########
##################
#modelling part: vedere quale modello è più appropriato tramite AIC (più basso possibile)
library(mgcv)
str(dati)
summary(dati)
dati$E=as.factor(dati$E)

M <- gam(Y~+s(var1)+s(var2)+s(var3) + s(var4) + E,data=dati)
summary(M)
AIC(M)
step(M)  ###trova il miglior modello basato su AIC
coef(M)
plot(M) 

##or
M2 <- gam(Y~+s(var1)+s(var2)+s(var3) + s(var4) + E,method = "REML",  data=dati)
##or GAM con tweedie dist, when lots of 0
M3 <- gam(Y~+s(var1)+s(var2)+s(var3) + s(var4) +s(E), family=Tweedie(1.5, power(0)), data=dati)
#gam Binomial, link logit per dati 0-1, 0%-100% 
M <- gam(Y~+s(var1)+s(var2)+s(var3) + s(var4), data=dati, family = binomial(link = "logit"),method = "REML")
#GAM con INTERAZIONE####
M <- gam(Y~+s(var1,var2)+E,method="REML", data=dati)
 
# cOOK'S DISTANCE
par(mfrow=c(1,1))
plot(cooks.distance(M), ylab="Cook's distance", type = "h", ylim=c(0,1))
abline(h=1, col=1,lwd=2)
###visualizare il modello 3D
vis.gam(M, view = (c("var1","var2")))
vis.gam(M, view = c("var1","var2"),plot.type = "contour")
#Analisi dei Residui
E.el <- resid(M)
Fit.el <- fitted(M)
plot(x = Fit.el, y = E.el, xlab = "Fitted values",
     ylab = "Residuals", main = "Dati")
abline(0,0)
tmp <- loess(E.el ~Fit.el ,span=0.75)
tmp2 <- predict(tmp,se=T)
I1 <- order(Fit.el)
lines(Fit.el[I1], tmp2$fit[I1], lty=1)
lines(Fit.el[I1], tmp2$fit[I1] + 2*tmp2$se.fit[I1], lty = 2)
lines(Fit.el[I1], tmp2$fit[I1] - 2*tmp2$se.fit[I1], lty = 2)
#residui pearson
E.m4<- resid(M, type = "pearson")
Fit.m4 <- fitted(M)
plot(x = Fit.m4, y = E.m4, xlab = "Fitted values",
     ylab = "Residuals")
plot(x = Fit.el, y = E.m4, xlab = "Fitted values",
     ylab = "Residuals", main = "Y")
abline(0,0)
tmp <- loess(E.m4 ~Fit.m4 ,span=0.75)
tmp2 <- predict(tmp,se=T)
I1 <- order(Fit.m4)
lines(Fit.m4[I1], tmp2$fit[I1], lty=1)
lines(Fit.m4[I1], tmp2$fit[I1] + 2*tmp2$se.fit[I1], lty = 2)
lines(Fit.m4[I1], tmp2$fit[I1] - 2*tmp2$se.fit[I1], lty = 2)

#graficare i residui (validation model)
E1 <- resid(M, type = "pearson")
var1 <- dati$var
plot(x = var1, y = E1, xlab = "var1",
     ylab = "Residuals",main = "var")
abline(0,0)
tmp <- loess(E1 ~var1 ,span=0.75)
tmp2 <- predict(tmp,se=T)
I1 <- order(var1)
lines(var1[I1], tmp2$fit[I1], lty=1)
lines(var1[I1], tmp2$fit[I1] + 2*tmp2$se.fit[I1], lty = 2)
lines(var1[I1], tmp2$fit[I1] - 2*tmp2$se.fit[I1], lty = 2)
####e cosi per tutte le altre variabili significarive del modello......
##Normality
par(mar=c(5,5,2,2))
hist(E.el, ylab = "Frequency", xlab = "Residuals", las=1,breaks=18, cex.lab=1.1, cex.axis=1.1)
#Or qq-plot
par(mar=c(5,5,2,2))
qqnorm(E.el, lwd=1.5,cex.lab=1.1, las=1,cex.axis=1, bty="l", xlab="Theoretical quantiles", ylab="Sample quantiles")