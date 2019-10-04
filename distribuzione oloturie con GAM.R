#### Mappe distribuzione Solemon con GAM
library(fitdistrplus)
library(mgcv)
library(maps)
library(mapdata)
library(mapplots)
library(shapefiles)
library(caret)
setwd("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/progetto Oloturie gianna")

dataset <- read_csv("Olofor_with_grain.csv")
#View(dataset)

data<-dataset
summary(data)
colSums(is.na(data))


### esplorazione dei dati
par(mfrow=c(1,2))
hist(data$TOT);plot(data$TOT~data$year)
hist(data$HOLOFOR);boxplot(data$HOLOFOR~data$year)
hist(data$HOLOTND);boxplot(data$HOLOTND~data$year)
hist(data$HOLOTUB);boxplot(data$HOLOTUB~data$year)
hist(data$OCNUPLA);boxplot(data$OCNUPLA~data$year)
hist(data$STICREG);boxplot(data$STICREG~data$year)
hist(data$TRACELO);plot(data$TRACELO~data$year)
hist(data$TRACTER);plot(data$TRACTER~data$year)
par(mfrow=c(1,1))
# trasformazione log

data$logTOT <-  log(data$TOT+1)
data$logHOLOFOR <-   log(data$HOLOFOR+1)
data$logHOLOTND <-    log(data$HOLOTND+1)
   data$logHOLOTUB <-   log(data$HOLOTUB+1)
   data$logOCNUPLA <-    log(data$OCNUPLA+1)
   data$logSTICREG <-    log(data$STICREG+1)
   data$logTRACELO <-    log(data$TRACELO+1)
   data$logTRACTER <-    log(data$TRACTER+1)
   

###collinearità
Z2<-cbind( data$logTOT, data$logHOLOFOR,data$logHOLOTND,data$logHOLOTUB,data$logOCNUPLA,data$logSTICREG,data$logTRACELO,data$logTRACTER, data$X, data$Y, data$depth, data$year)
## creo i nomi da visualizzare nel grafico
colnames(Z2) <- c("TOT","HOLOFOR","HOLOTND","HOLOTUB","OCNUPLA","STICREG","TRACELO","TRACTER","lon","lat","depth", "year")
pairs(Z2, lower.panel = panel.smooth2,
      upper.panel = panel.cor, diag.panel = panel.hist)

## carico lo shape della GSA
# MAPS 
GSAtable <<- read.table("gsa_coordinates2.csv", sep=";", header=T)
leg_pos <<- as.character(GSAtable[GSAtable$GSA==19,6])
coord_lon <<- as.numeric( GSAtable[GSAtable$GSA==17,c(2,3)] )
coord_lat <<- as.numeric( GSAtable[GSAtable$GSA==17,c(4,5)] )
GSA17<-read.shapefile("Med_Poly - Copia/Med_Poly")
#carico il grid adriatico da Emodnet
library(readr)
Adriatic_Grid <- read_csv("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Mappe distribuzione Solemon/Adriatic_Grid2.csv")
#View(Adriatic_Grid)
#plot(Adriatic_Grid$lon,Adriatic_Grid$lat, cex=0.1)
library(dplyr)
adr<-Adriatic_Grid%>%  filter(depth<100 & depth>1 & GSA==17 & lon<16&lon>12,lat>41.5)
summary(adr)
length(adr$lon) #54950

###############################################
#####            MODELLO oloturie TOT         #
###############################################
data$logTOT<-data$logHOLOFOR  ## non ho voglia di cambuare tutto lo script!!!!!!!!!!!!
data$TOT<-data$HOLOFOR        ## non ho voglia di cambuare tutto lo script!!!!!!!!!!!!

hist(data$logTOT);plot(data$logTOT~data$year)
descdist(data$TOT,boot = 100)
descdist(data$logTOT,boot = 100)
summary(data$logTOT==0)
#188/661
#473/661
descdist(data$TOT,boot = 100)


####################################################
#   Zero inladcted Model                           #
####################################################
library(pscl)
summary(data$TOT)
data$rTOT<-round(data$TOT);data$rTOT
data$sqrt<-round(sqrt(data$TOT))
data$log<-round(log(data$TOT+1))
zi<-zeroinfl(rTOT ~ X*Y + depth, data = data, dist = "negbin", link = "logit")
summary(zi);AIC(zi)

Fit <- fitted(zi)
plot(Fit)
res <- residuals(zi)
plot(res)
abline(0,0, col="red")
##Normality
par(mar=c(5,5,2,2))
hist(res, ylab = "Frequency", xlab = "Residuals", las=1,breaks=18, cex.lab=1.1, cex.axis=1.1)
#Or qq-plot
par(mar=c(5,5,2,2))
qqnorm(res, lwd=1.5,cex.lab=1.1, las=1,cex.axis=1, bty="l", xlab="Theoretical quantiles", ylab="Sample quantiles")

# predict
pzii<-predict(zi, data,type = "response");summary(pzii)
data$pzii<-round(pzii);summary(data$pzii)
data$pzi <- rep(0, dim(data)[1])
for (i in 1:(dim(data)[1])) {
  if(data$pzii[i] > 0) 
  {data$pzi[i] <- 1} else
  {data$pzi[i] <- 0}
}
table(data$pzi , data$pres)
confusionMatrix(table(data$pzi , data$pres))


## DELTA GAM approch
# creo il database con dati presenza/assenza
data$pres <- rep(0, dim(data)[1])
for (i in 1:(dim(data)[1])) {
  if(data$logTOT[i] > 1) 
  {data$pres[i] <- 1} else
  {data$pres[i] <- 0}
} 
summary(factor(data$pres))
descdist(data$pres,boot = 100)
# creo il database con le cale solo positivi
abb<-data %>% dplyr::select("year", "X", "Y","depth", "logTOT", "TOT", "grain") %>%  filter(logTOT > 0)
summary(abb)
hist(abb$logTOT)
plot(abb$logTOT~abb$year)
descdist(abb$logTOT,boot = 100)
descdist(abb$TOT,boot = 100)
#####################
##binomial part TOT
####################
# Model selection  was done through a forward stepwise selection approach 
modd <- gam(pres ~ s(grain), family= binomial(),data=data, select=T)
summary(modd); AIC(modd);plot(modd)
mod <- gam(pres ~ s(Y), family= binomial(),data=data, select=T)
summary(mod); AIC(mod);plot(mod)
mod1 <- gam(pres ~ s(X), family= binomial(),data=data, select=T)
summary(mod1); AIC(mod1);termplot(mod1)
mod2 <- gam(pres ~ s(depth), family= binomial(),data=data, select=T)
summary(mod2); AIC(mod2);plot(mod2)
mod3 <- gam(pres ~ factor(year), family= binomial(),data=data, select=T)
summary(mod3); AIC(mod3);termplot(mod3)

mod4 <- gam(pres ~ te(X,Y), family= binomial(),data=data, select=T)
summary(mod4); AIC(mod4);plot(mod4)
mod5 <- gam(pres ~ s(X,Y)+s(depth), family= binomial(),data=data, select=T)
summary(mod5); AIC(mod5);plot(mod5)
mod6 <- gam(pres ~ s(X,Y)+s(depth)+s(year), family= binomial(),data=data, select=T)
summary(mod6); AIC(mod6);plot(mod6)

mod7 <- gam(pres ~ s(X,Y)+s(year)+s(depth)+s(grain), family= binomial(),data=data, select=T)
summary(mod7); AIC(mod7)
gam.check(mod7)
#plot(mod7);termplot(mod7)
vis.gam(mod7, view = c("X","Y"),plot.type = "contour",too.far = 1,type = "response", asp=1, n.grid = 100, main="binomial; dev=82%")
points(data$X, data$Y , cex=data$logTOT/4, col="black")

mod8<- gam(pres ~   s(depth, k=40) +s(X,Y,k=120) +s(year, k=10)+te(Y,year)+te(X,year),
           family= binomial(),data=data, select=T)
mod8<- gam(pres ~   s(depth, k=40) +s(X,Y)+te(Y,year)+te(X,year),
           family= binomial(),data=data, select=T)
summary(mod8); AIC(mod8)
gam.check(mod8)
# Model selection was done through a backward stepwise selection approach based on statistical significance (Wood, 2006). 
mod9<- gam(pres ~   s(depth) +s(X,Y) +s(year)+te(Y,year)+te(X,year)+te(depth,year)+s(Y)+s(X)+s(grain),
              family= binomial(),data=data, select=T)
summary(mod9)
mod10<- gam(pres ~   s(depth, k=40) +s(X,Y,k=120) +factor(year)+te(X,year)+te(Y,year),
           family= binomial(),data=data, select=T)
summary(mod10); AIC(mod10)

################### 
# mod8 e mod10 verranno testati per vedere quale è il migliore 
###################
par(mfrow=c(3,2))
plot(mod8);plot(mod10)

par(mfrow=c(2,2))
gam.check(mod8);gam.check(mod10)

res <- residuals(mod8)
plot(res, ylim=c(-5.5,5.5))
abline(0,0, col="red")
res <- residuals(mod10)
plot(res, ylim=c(-5.5,5.5))
abline(0,0, col="red")

# cOOK'S DISTANCE
par(mfrow=c(1,1))
plot(cooks.distance(mod8), ylab="Cook's distance", type = "h", ylim=c(0,1))
abline(h=1, col=1,lwd=2)

#visualizazione grafica
par(mfrow=c(1,2))
vis.gam(mod8, view = c("X","Y"),plot.type = "contour",too.far = 1,type = "response", asp=1, n.grid = 100, cond = list(depth=median(data$depth)))
#draw.shape(GSA17, col = "white")
points(data$X, data$Y , cex=data$logTOT/4, col="black")
vis.gam(mod10, view = c("X","Y"),plot.type = "contour",too.far = 1,type = "response", asp=1, n.grid = 100)
#draw.shape(GSA17, col = "white")
points(data$X, data$Y , cex=data$logTOT, col="black")

# predict con mod8
p8<-predict(mod8, data,type = "response");summary(p8)
p8<-round(p8,0);summary(p8)
table(p8 , data$pres)
confusionMatrix(table(p8, data$pres))
# predict con mod10
p10<-predict(mod10, data,type = "response");summary(p10)
p10<-round(p10,0);summary(p10)
table(p10 , data$pres)
confusionMatrix(table(p10, data$pres))
# predict con mod7
p7<-predict(mod7, data,type = "response");summary(p7)
p7<-round(p7,0);summary(p7)
confusionMatrix(table(p7, data$pres))

##### in base alle varie diagnostiche scelgo il modello mod7 come modello di presenza assenza!!
binomial_mod<-mod7
########################

########################
## Gaussuian part
########################

# Model selection was done through a backward stepwise selection approach based on statistical significance (Wood, 2006). 

gau<- gam(logTOT ~   s(depth) +s(X,Y, k=120) +s(year, k=10)+te(X,year)+te(depth,year),
           family= gaussian (link="identity"),data=abb, select=T)
summary(gau)

gau1<- gam(logTOT ~  + s(depth, k=40) +s(X,Y,k=120)+factor(year) +te(Y,year), 
             family= gaussian (link="identity"),data=abb, select=T)
gau1<- gam(logTOT ~  + s(depth) +s(X,Y)+s(year)+s(grain), 
           family= gaussian (link="identity"),data=abb, select=T)
summary(gau1);AIC(gau1)
gam.check(gau1)
plot(gau1)
vis.gam(gau1, view = c("X","Y"),plot.type = "contour",type = "response", too.far=0.05,asp=1, n.grid = 100, main="gauss")
#draw.shape(GSA17, col = "white")
points(abb$X, abb$Y , cex=abb$logTOT/4, col="black")

# predict con gau1
pp1<-predict(gau1, abb,type = "response");summary(exp(pp1))
summary(abb$TOT)


##############################################################
#      vis.gam finali con i modelli selezionati              #
##############################################################
par(mfrow=c(1,2))
vis.gam(mod7, view = c("X","Y"),plot.type = "contour",too.far = 1,type = "response", asp=1, n.grid = 100, main="binomial; dev=82%")
#draw.shape(GSA17, col = "white")
points(data$X, data$Y , cex=data$logTOT/4, col="black")
vis.gam(gau1, view = c("X","Y"),plot.type = "contour",too.far = 1,type = "response", asp=1, n.grid = 100, main="gauss; dev=94.3% ")
#draw.shape(GSA17, col = "white")
points(abb$X, abb$Y , cex=abb$logTOT/4, col="black")
#vis.gam(tw7.5, view = c("X","Y"),plot.type = "contour",too.far =1,type = "response", asp=1, n.grid = 100, main ="Tweedie; dev=68.2%")
#draw.shape(GSA17, col = "white")
#points(data$X, data$Y , cex=data$logTOT/4, col="black")


########################################################
#   PARTIRE DA QUI SE SI HANNO I MODELLI GIA FATTI   #
########################################################
#  MAPPE DISTRIBUZIONE con ggplot
##########################################################

##pred sul grid adriactico con grainsize(indipendentemente dagli anni)
adrgrid$grain<-adrgrid$grain_mask
summary(adrgrid)
# setto l'anno della mappa!!
adrgrid$year<-rep(2017, 39754)
adr<-adrgrid %>% dplyr::select("X", "Y", "depth", "year", "grain")
summary(adr)
#adr1<- adr[1:55008,]
#adr2<- adr[55009:110016,]
### predizione su griglia Tweedie approch
##################
#pred = predict(tw7.5, adr, type = "response")
#predd = predict(modfin0, adr, type = "response")


#adr$pred = pred^2;summary(adr$pred)


#par(mfrow=c(1,2))
#hist(adr$pred);plot(adr$pred)

#adr$pp<-(adr$pred2/68);summary(adr$pp)

# rimetto insieme i due database per fare la mappa finale
#adr<-rbind(adr1,adr2)
#summary(adr$pred)

# si trasformano i nem_km2 in classi di "presenza"
# adr$abbol <- rep(0, dim(adr)[1])
# for (i in 1:(dim(adr)[1])) {
#    if(adr$pred[i] > 0 &  adr$pred[i] <= 100) 
#     {adr$abbol[i] <- "(L) 0-100"} else
#       if(adr$pred[i] > 0 &  adr$pred[i] <= 500) 
#       {adr$abbol[i] <- "(I) 100-500"} else
#         if(adr$pred[i] > 500 &  adr$pred[i] <= 1500) 
#         {adr$abbol[i] <-"(H) 500-1500"} else
#           if(adr$pred[i] > 1500 &  adr$pred[i] <= 3000) 
#           {adr$abbol[i] <- "(G) 1500-3000"} else
#             if(adr$pred[i] > 3000 &  adr$pred[i] <= 6000) 
#             {adr$abbol[i] <- "(F) 3000-6000"} else
#               if(adr$pred[i] > 6000 &  adr$pred[i] <= 10000) 
#               {adr$abbol[i] <- "(E) 6000-10000"} else
#                 if(adr$pred[i] > 10000 &  adr$pred[i] <= 30000) 
#                 {adr$abbol[i] <- "(D) 10000-30000"} else
#                   if(adr$pred[i] > 30000 &  adr$pred[i] <= 50000) 
#                   {adr$abbol[i] <- "(C) 30000-50000"} else
#                     if(adr$pred[i] > 50000 &  adr$pred[i] <= 100000) 
#                     {adr$abbol[i] <- "(B) 50000-100000"} else
#                       if(adr$pred[i] > 100000 ) 
#                       {adr$abbol[i] <- "(A) >100000"} 
# }
# str(adr)
# summary(adr)
# unique(adr$abbol)
# boxplot(adr$pred~ adr$abbol)
# library(ggplot2)
# #library(rgdal)
# library(readr)
 map1<-map_data("worldHires") %>% filter(long > 11 & long < 18 & lat > 40 & lat < 48)
 bt100<-Adriatic_Grid%>%  filter(depth<2000 & depth>100 & GSA==17 & lon<18&lon>12,lat>41.5& lat<44)
# 
# ggplot(adr,aes(X,Y,fill=adr$abbol))+theme_light()+geom_raster(interpolate=T,show.legend = T)+scale_fill_grey()+coord_quickmap(xlim=c(12,15.25),ylim=c(41.9,45.8)) +theme(legend.title = element_text( size = 10))+
#   labs(x=NULL,y=NULL,fill="Abbondanza (Num/Km2)")+geom_polygon(data=bt100,aes(lon,lat),fill=NA,color= NA )+geom_polygon(data=map1,aes(long,lat,group=group),fill="lightgoldenrodyellow",color="black")+ggtitle("Tweedie 2017")


# DELTA Approach moltiplico le abbondanze predette per la probabilità di presenza del modello binomiale!!
pred0 = predict(binomial_mod, adr, type = "response")
predgau = predict(gau1, adr, type = "response")
#predd = predict(modfin0, adr, type = "response")
adr$pred0 = pred0;summary(adr$pred0)
adr$predgau = exp(predgau);summary(adr$predgau)

adr$predfin<-adr$pred0*adr$predgau
summary(adr$predfin); summary(data$TOT)
# si trasformano i nem_km2 in classi di "presenza"
adr$abbb <- rep(0, dim(adr)[1])
for (i in 1:(dim(adr)[1])) {
 if(adr$predfin[i] > 0 &  adr$predfin[i] <= 1) 
 {adr$abbb[i] <- "(M) <1"} else
  if(adr$predfin[i] > 1 &  adr$predfin[i] <= 100) 
  {adr$abbb[i] <- "(L) 1-100"} else
    if(adr$predfin[i] > 100 &  adr$predfin[i] <= 500) 
    {adr$abbb[i] <- "(I) 100-500"} else
      if(adr$predfin[i] > 500 &  adr$predfin[i] <= 1500) 
      {adr$abbb[i] <-"(H) 500-1500"} else
        if(adr$predfin[i] > 1500 &  adr$predfin[i] <= 3000) 
        {adr$abbb[i] <- "(G) 1500-3000"} else
          if(adr$predfin[i] > 3000 &  adr$predfin[i] <= 6000) 
          {adr$abbb[i] <- "(F) 3000-6000"} else
            if(adr$predfin[i] > 6000 &  adr$predfin[i] <= 10000) 
            {adr$abbb[i] <- "(E) 6000-10000"} else
              if(adr$predfin[i] > 10000 &  adr$predfin[i] <= 30000) 
              {adr$abbb[i] <- "(D) 10000-30000"} else
                if(adr$predfin[i] > 30000 &  adr$predfin[i] <= 50000) 
                {adr$abbb[i] <- "(C) 30000-50000"} else
                  if(adr$predfin[i] > 50000 &  adr$predfin[i] <= 100000) 
                  {adr$abbb[i] <- "(B) 50000-100000"} else
                    if(adr$predfin[i] > 100000 ) 
                    {adr$abbb[i] <- "(A) >100000"} 
}
str(adr)
summary(adr)
unique(adr$abbb)
boxplot(adr$predfin~ adr$abbb)

summary(adr$abbb=="(M) <1")
summary(adr$abbb=="(L) 1-100")
summary(adr$abbb=="(I) 100-500")
summary(adr$abbb=="(H) 500-1500")
summary(adr$abbb=="(G) 1500-3000")
summary(adr$abbb=="(F) 3000-6000")
summary(adr$abbb=="(E) 6000-10000")
summary(adr$abbb=="(D) 10000-30000")
summary(adr$abbb=="(C) 30000-50000")
summary(adr$abbb=="(B) 50000-100000")
summary(adr$abbb=="(A) >100000")

# si trasformano i nem_km2 in classi di "presenza"
adr$pres <- rep(0, dim(adr)[1])
for (i in 1:(dim(adr)[1])) {
  if(adr$pred0[i] > 0 &  adr$pred0[i] <= 0.2) 
   {adr$pres[i] <- "(M) <20%"} else
    if(adr$pred0[i] > 0.2 &  adr$pred0[i] <= 0.4) 
    {adr$pres[i] <- "(L) 20%-40%"} else
      if(adr$pred0[i] > 0.4 &  adr$pred0[i] <= 0.6) 
      {adr$pres[i] <- "(I) 40%-60%"} else
        if(adr$pred0[i] > 0.6 &  adr$pred0[i] <= 0.8) 
        {adr$pres[i] <-"(H) 40%-80%"} else
          if(adr$pred0[i] > 0.8 ) 
          {adr$pres[i] <- "(A) >80%"}
}
  boxplot(adr$pred0~ adr$pres)

## zeroinf 
  adr$predzi<-predict(zi, adr, type = "response")
  summary(adr$predzi)
  adr$zi <- rep(0, dim(adr)[1])
  for (i in 1:(dim(adr)[1])) {
    if(adr$predzi[i] > 0 &  adr$predzi[i] <= 1) 
    {adr$zi[i] <- "(M) <1"} else
      if(adr$predzi[i] > 1 &  adr$predzi[i] <= 100) 
      {adr$zi[i] <- "(L) 1-100"} else
        if(adr$predzi[i] > 100 &  adr$predzi[i] <= 500) 
        {adr$zi[i] <- "(I) 100-500"} else
          if(adr$predzi[i] > 500 &  adr$predzi[i] <= 1500) 
          {adr$zi[i] <-"(H) 500-1500"} else
            if(adr$predzi[i] > 1500 &  adr$predzi[i] <= 3000) 
            {adr$zi[i] <- "(G) 1500-3000"} else
              if(adr$predzi[i] > 3000 &  adr$predzi[i] <= 6000) 
              {adr$zi[i] <- "(F) 3000-6000"} else
                if(adr$predzi[i] > 6000 &  adr$predzi[i] <= 10000) 
                {adr$zi[i] <- "(E) 6000-10000"} else
                  if(adr$predzi[i] > 10000 &  adr$predzi[i] <= 30000) 
                  {adr$zi[i] <- "(D) 10000-30000"} else
                    if(adr$predzi[i] > 30000 &  adr$predzi[i] <= 50000) 
                    {adr$zi[i] <- "(C) 30000-50000"} else
                      if(adr$predzi[i] > 50000 &  adr$predzi[i] <= 100000) 
                      {adr$zi[i] <- "(B) 50000-100000"} else
                        if(adr$predzi[i] > 100000 ) 
                        {adr$zi[i] <- "(A) >100000"} 
  }
  
  # setto l'anno
  adr$year<-rep(2017, 54950)
# mappa deltagam
ggplot(adr,aes(X,Y,fill=adr$abbb))+theme_light()+geom_raster(interpolate=T,show.legend = T)+scale_fill_grey()+coord_quickmap(xlim=c(12,15.25),ylim=c(41.9,45.8)) +theme(legend.title = element_text( size = 10))+
  labs(x=NULL,y=NULL,fill="Abbondanza (Num/Km2)")+geom_polygon(data=bt100,aes(lon,lat),fill=NA,color= NA )+geom_polygon(data=map1,aes(long,lat,group=group),fill="lightgoldenrodyellow",color="black")+ggtitle("HOLOFOR DeltaGam 2017")

# mappa pres
ggplot(adr,aes(X,Y,fill=adr$pres))+theme_light()+geom_raster(interpolate=T,show.legend = T)+scale_fill_grey()+coord_quickmap(xlim=c(12,15.25),ylim=c(41.9,45.8)) +theme(legend.title = element_text( size = 10))+
  labs(x=NULL,y=NULL,fill="Abbondanza (Num/Km2)")+geom_polygon(data=bt100,aes(lon,lat),fill=NA,color= NA )+geom_polygon(data=map1,aes(long,lat,group=group),fill="lightgoldenrodyellow",color="black")+ggtitle("HOLOFOR Presenza 2017")

# mappa zeroinf
ggplot(adr,aes(X,Y,fill=adr$zi))+theme_light()+geom_raster(interpolate=T,show.legend = T)+scale_fill_grey()+coord_quickmap(xlim=c(12,15.25),ylim=c(41.9,45.8)) +theme(legend.title = element_text( size = 10))+
  labs(x=NULL,y=NULL,fill="Abbondanza (Num/Km2)")+geom_polygon(data=bt100,aes(lon,lat),fill=NA,color= NA )+geom_polygon(data=map1,aes(long,lat,group=group),fill="lightgoldenrodyellow",color="black")+ggtitle("HOLOFOR zeroinf")

#x11()
#savePlot(filename = paste("Solea abb_e",sep=""),type = "tiff",device = dev.cur(),restoreConsole = T)
summary(adr$predfin);summary(data$TOT); summary(abb$TOT)



