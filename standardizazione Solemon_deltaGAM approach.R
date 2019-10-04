# DELTA GAM approach      #
###########################
library(mgcv)
library(readr)
library(shapefiles)
library(mapplots)

dataset <- read_delim("Solea stnd (20-6-2019).csv", ";", 
                      escape_double = FALSE, trim_ws = TRUE)
View(dataset)
data<-dataset
summary(data)
# computation of mean haul depth
data$depth <- (data$shooting_depth+data$hauling_depth)/2
# hour extraction
data$hour<-as.numeric(as.character(ifelse(nchar(data$shooting_time)==3, substr(data$shooting_time,1,1), substr(data$shooting_time,1,2))))
# transformation of dependent variable
data$n_sqrt <- log(data$n_km2+1)
dotchart(data$n_sqrt,
         ylab = "Order of observations",
         xlab = "num", main = "Cleveland dotplot")


summary(data$n_sqrt==0)

# creo il database con dati presenza/assenza
data$pres <- rep(0, dim(data)[1])
for (i in 1:(dim(data)[1])) {
  if(data$n_sqrt[i] > 1) 
  {data$pres[i] <- 1} else
  {data$pres[i] <- 0}
} 
summary(factor(data$pres))
# creo il database con le cale solo positivi
abb <- data[data$n_sqrt>0,]
summary(abb)
hist(abb$n_km2)
hist(abb$n_sqrt)

# Stepwise forward inclusion for GAM modeling
#mod1 <- gam(n_sqrt ~ s(X), family= gaussian (link="identity"),data=data, select=T) 
#summary(mod1) 
#mod2 <- gam(n_sqrt ~ s(Y), family= gaussian (link="identity"),data=data, select=T)
#summary(mod2) 
#mod3 <- gam(n_sqrt ~ s(depth), family= gaussian (link="identity"),data=data, select=T)
#summary(mod3) 
#mod4 <- gam(n_sqrt ~ s(hour), family= gaussian (link="identity"),data=data, select=T)
#summary(mod4) 
#mod5 <- gam(n_sqrt ~ s(year), family= gaussian (link="identity"),data=data, select=T)
#summary(mod5) 
#mod6 <- gam(n_sqrt ~ s(hauling_duration), family= gaussian (link="identity"),data=data, select=T)
#summary(mod6) 
#mod7 <- gam(n_sqrt ~ factor(vessel), family= gaussian (link="identity"),data=data, select=T)
#summary(mod7) 
#boxplot(data$n_sqrt~data$vessel)

# once that the best basic model was detected, the variables should be included one at the time in the same way
#.....


#####################
##binomial part
####################
modfin0<- gam(pres ~   s(depth) +s(X,Y)+s(year, k=13) +s(hour)+te(Y,year)+te(X,year), 
              family= binomial(),data=data, select=T)

summary(modfin0)
AIC(modfin0)
par(mfrow=c(1,1))
plot(modfin0)
termplot(modfin0)
par(mfrow=c(2,2))
gam.check(modfin0)
res <- residuals(modfin0)
plot(res, ylim=c(-5.5,5.5))
abline(0,0, col="red")
########################
## Gaussuian part
########################

#modfin<- gam(n_sqrt~  + s(depth) +s(X,Y)+s(year, k=13) +s(hour)+te(Y,year)+te(X,year), 
 #            family= gaussian (link="identity"),data=abb, select=T)
modfin<- gam(n_sqrt~  + s(depth) +s(X,Y)+s(year, k=13) +s(hour)+te(Y,year)+te(X,year), 
             family= gaussian (link="identity"),data=abb, select=T)

summary(modfin)
AIC(modfin)
par(mfrow=c(1,1))
plot(modfin)
termplot(modfin)
par(mfrow=c(2,2))
gam.check(modfin)
vis.gam(modfin, view = c("X","Y"),plot.type = "contour",type = "response", asp=1, n.grid = 100)
#vis.gam(modfin, view = c("year","vessel"),plot.type = "contour", cond = list(data$vessel=="DLP"))
#GSA17<-read.shapefile("Med_Poly - Copia/Med_Poly")
draw.shape(GSA17, col = "white")
points(data$X, data$Y , cex=data$n_sqrt/4, col="black")
# inspection of residuals and diagnostic plots
res <- residuals(modfin)
plot(res, ylim=c(-5.5,5.5))
abline(0,0, col="red")


analysis_stratum1 <- T
analysis_stratum2 <- T
analysis_stratum3 <- T

Adriatic_Grid <- read_csv("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Mappe distribuzione Solemon/Adriatic_Grid2.csv")
View(Adriatic_Grid)
#plot(Adriatic_Grid$lon,Adriatic_Grid$lat, cex=0.1)
library(dplyr)
grid<-Adriatic_Grid%>%  filter(depth<100 & depth>1 & GSA==17 & lon<16&lon>12,lat>41.5)
summary(grid)
stratum_1 <- grid[grid$depth >= 0 & grid$depth < 30, ]
stratum_2 <- grid[grid$depth >= 30 & grid$depth < 50, ]
stratum_3 <- grid[grid$depth >= 50 , ]
stratum_1$X<-stratum_1$lon;stratum_1$Y<-stratum_1$lat
stratum_2$X<-stratum_2$lon;stratum_2$Y<-stratum_2$lat
stratum_3$X<-stratum_3$lon;stratum_3$Y<-stratum_3$lat
summary(stratum_1);summary(stratum_2);summary(stratum_3);

area_s1<- 11512  #11361
area_s2<- 8410
area_s3<- 22466

### standardizzo con ora=12 e haul_dur=30
res_table <- data.frame(seq(2005,2018, 1))
res_table$index_1 <- NA
res_table$index_2 <- NA
res_table$index_3 <- NA


i=1
for(i in 1:length(res_table[,1])){
  if (analysis_stratum1 == T){
    stratum_1$year <- res_table[i,1]
    # stratum_1$hauling_duration <- 30
    stratum_1$hour <- 12
    stratum_1$pred0 <- predict(modfin0, newdata = data.frame(stratum_1), type="response")
    stratum_1$predgau <-  exp(predict(modfin, newdata = data.frame(stratum_1)))+1
    stratum_1$pred <-stratum_1$pred0*stratum_1$predgau 
    res_table[i,2] <- mean(stratum_1$pred) *area_s1} else {
      area_s1 = 0
    }
  
  if (analysis_stratum2 == T){
    stratum_2$year <- res_table[i,1]
    #  stratum_2$hauling_duration <- 30
    stratum_2$hour <- 12
    stratum_2$pred0 <- predict(modfin0, newdata = data.frame(stratum_2), type="response")
    stratum_2$predgau <-  exp(predict(modfin, newdata = data.frame(stratum_2)))+1
    stratum_2$pred <-stratum_2$pred0*stratum_2$predgau 
    res_table[i,3] <- mean(stratum_2$pred) *area_s2} else {
      area_s2 = 0
    }
  
  if (analysis_stratum3 == T){
    stratum_3$year <- res_table[i,1]
    #   stratum_3$hauling_duration <- 30
    stratum_3$hour <- 12
    stratum_3$pred0 <- predict(modfin0, newdata = data.frame(stratum_3), type="response")
    stratum_3$predgau <-  exp(predict(modfin, newdata = data.frame(stratum_3)))+1
    stratum_3$pred <-stratum_3$pred0*stratum_3$predgau 
    res_table[i,4] <- mean(stratum_3$pred) *area_s3} else {
      area_s3 = 0
    }
  
  sum_res <- c(res_table[i,2],res_table[i,3],res_table[i,4])
  res_table[i, 5]<- sum(sum_res[!is.na(sum_res)])/sum(area_s1,area_s2,area_s3)
}

colnames(res_table) <- c("year", "stratum 1","stratum 2", "stratum 3", "Indices")
par(mfrow=c(1,1))
plot(res_table[,1], res_table[,5], type="b", xlab="year", ylab= "index")



# solemon confronto
time_series <- c(279,318,377,227,248,269,368,426,711,826,608,608,521,760)
par(mfrow=c(1,1))
plot(res_table[,1], time_series, col="red", xlab="year", ylab= "N/km^2", main = "Solemon Abb Index" , pch=16, ylim=c(0,1300))
lines(res_table[,1], time_series, col="red")
points(res_table[,1], res_table[,5], col="black", pch=16)
lines(res_table[,1], res_table[,5], col="black")
legend("topleft", c("Time series","Std Prediction (deltagam)"), col=c("red", "black"), lwd=1, pch=16)

Index=res_table2[,c(1,5)]
colnames(Index)=c("Year","Index")
write.table(Index,"Index Solemon standardizzato deltagam.csv",sep=";",row.names=F)
