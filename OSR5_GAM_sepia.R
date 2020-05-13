#### Mappe distribuzione Solemon con GAM OSR5

library(mgcv)
library(maps)
library(mapdata)
library(mapplots)
library(shapefiles)
library(dplyr)
library(readr)
library(ROCR)
library(caret)
library(ggplot2)
library(readxl)
library(tiff)
library(raster)
library(rgeos)
library(reshape2)
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")

setwd("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/Progetto SEAFAIR/OSR5")
data <- read_delim("sepia_withEff.csv", ";", escape_double = FALSE, trim_ws = TRUE)
View(data)
summary(data)
colSums(is.na(data))
plot(data$MEAN_LONGITUDE_DEC)
plot(data$MEAN_LATITUDE_DEC)

data <- data %>% dplyr::select(country, area, year.x, month,haul_number.x, shooting_time, hauling_duration,MEAN_DEPTH,MEAN_LATITUDE_DEC, MEAN_LONGITUDE_DEC,SWEPT_AREA, STRATUM_CODE,CHL,TMP_sst,TMP_bot,dox.ss,dox.bot,nit, pho,sal,ph,poc,zc,grain,eff,genus,N_km2,kg_km2) %>% dplyr::rename("X" = "MEAN_LONGITUDE_DEC")%>% dplyr::rename("Y" = "MEAN_LATITUDE_DEC")%>% dplyr::rename("year" = "year.x")%>% dplyr::rename("depth" = "MEAN_DEPTH")%>% na.omit()

#carico il grid adriatico di diego # setto il mese e gli anni della griglia!!
adr<-readRDS("grid.Rdata")%>% dplyr::filter(year > 2004) %>% dplyr::filter(month=="11")  %>% dplyr::filter(depth < 101)%>% dplyr::filter(Y > 40) 
summary(adr)
#adr$X<-round(adr$X,2)
#adr$Y<-round(adr$Y,2)
#adr <- adr%>% distinct(X,Y ,year, .keep_all = TRUE)
#carico il grid adriatico di diego con EFFORT dal 2008 # setto il mese e gli anni della griglia!!
adr_eff<-readRDS("grid_with.eff.rds") %>% dplyr::filter(month=="11")  %>% dplyr::filter(depth < 101)%>% dplyr::filter(Y > 41.5) %>% dplyr::filter(Y < 46) %>% dplyr::filter(X < 16) 
summary(adr_eff)
#adr_eff$X<-round(adr_eff$X,2)
#adr_eff$Y<-round(adr_eff$Y,2)
#adr_eff <- adr_eff%>% distinct(X,Y ,year, .keep_all = TRUE)

# carico la mappa della granulometria
#gr<-'grain_17.tif' 
#gr<-raster(gr)
#grp<-rasterToPoints(gr)
#grp<-as.data.frame(grp)
#grp$X<-round(grp$x,2)
#grp$Y<-round(grp$y,2)
#grp$grain <- grp$grain_17
#grp <- grp%>% distinct(X,Y ,grain ,.keep_all = TRUE)%>% dplyr::select(-x,-y,-grain_17)
# merg adr e grp per avere un grid finale
#adr = merge(adr,grp, by.x=c("X","Y"), by.y=c("X","Y"))
#summary(adr);colSums(is.na(adr))
#plot(adr$X,adr$Y, cex=log(adr$grain), col=3)
#points(data$X, data$Y, cex=0.5, pch = 1,col=1,lty=5,lwd=5)

####
hist(data$kg_km2)
plot(data$kg_km2)
# trasformazione log
data$log.kg_km2 <- log(data$kg_km2+1)
hist(data$log.kg_km2) 
# trasformazione sqrt
data$sqrt.kg_km2 <- sqrt(data$kg_km2)
hist(data$sqrt.kg_km2) 
# creo il database con dati presenza/assenza
data$PA <- rep(0, dim(data)[1])
for (i in 1:(dim(data)[1])) {
  if(data$genus[i] == "SEPI" ) 
  {data$PA[i] <- 1} else
  {data$PA[i] <- 0}
} 
summary(factor(data$PA))

source('C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Mappe distribuzione Solemon/collinearity.R')
###modello gam SOLEA
Z2<-cbind( data$kg_km2, data$sqrt.kg_km2 ,data$log.kg_km2, data$X, data$Y, data$depth, data$grain, as.factor(data$month), data$CHL,data$TMP_sst,data$TMP_bot,data$dox.ss,data$dox.bot, data$nit,data$pho, data$sal,data$ph,data$poc,data$zc, data$eff)
## creo i nomi da visualizzare nel grafico
colnames(Z2) <- c("kg_km2","sqrt.kg_km2","log.kg_km2","X","Y","depth", "grain", "month", "CHL", "SST", "Tbott","Ox", "Oxbott","Nit","Phos","Sal","ph","Poc","Zc", "Effort")
pairs(Z2, lower.panel = panel.smooth2,
      upper.panel = panel.cor, diag.panel = panel.hist)
## decidio di eliminare dalle analisi Poc Zn TMP_sst e dox.ss e CHL
input <- data %>% dplyr::select(-poc, -zc, -TMP_sst, -dox.ss, -ph , -CHL )


####################################
####################################
# passare a script "test50model"   #
####################################
####################################

##########################
# NO EFF data 2005-2018  #
###########################
# best model DELTAGAM     #
###########################

## binomial part
####################
modfin0<- gam(PA~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(pho), family= binomial(),data=input, select=T)
summary(modfin0)
AIC(modfin0)
par(mfrow=c(2,3))
plot(modfin0)
termplot(modfin0)
par(mfrow=c(2,2))
gam.check(modfin0)
## binomial part NO ENV
####################
modfin0_noenv<- gam(PA~ factor(month)+s(X,Y)+s(year)+s(depth), family= binomial(),data=input, select=T)
summary(modfin0_noenv)
# AUC Binomial models were tested for sensitivity by using the area under the receiver operating characteristic curve (AUC). An AUC value of 0.5 indicates that the model performs no better than a random model, whereas a value of 1 indicates that the model is fully capable of distinguishing between occupied and unoccupied sites. AUC values of 0.7-0.9 indicate very good discrimination, while values >0.9 indicate excellent discrimination. (Lauria 2017)
pred.test=predict(modfin0,input,type='response')                  # make prediction for test set
preds.obs=data.frame(pred.test=pred.test,input$PA)           # data frame of preds vs obs 
gam1.eval=prediction(preds.obs$pred.test,preds.obs$input.PA)    # Asses prediction for AUC
attributes(performance(gam1.eval, 'auc'))$y.values[[1]]          # get AUC value
#roc(preds.obs$test.set.PA,preds.obs$pred.test)                   # use gbm package to the AUC
roc <- performance( gam1.eval, "tpr", "fpr")
plot( roc )
#p<-round(p,0);summary(p)
#table(p , data$PA)
#confusionMatrix(table(p, data$PA))

## Gaussuian part 
########################
modfin<- gam(sqrt.kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth)+TMP_bot+s(nit),family= gaussian (link="identity"),data=subset(input,sqrt.kg_km2 >0), select=T)
summary(modfin)
AIC(modfin)
par(mfrow=c(2,3))
plot(modfin)
termplot(modfin)
par(mfrow=c(2,2))
gam.check(modfin)
## Gaussuian part No ENV 
########################
modfin_noenv<- gam(sqrt.kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth),family= gaussian (link="identity"),data=subset(input,sqrt.kg_km2 >0), select=T)
summary(modfin_noenv)

# predict delta gam
preddeltaB<-predict(modfin0, input, type = "response")
summary(preddeltaB)
R2(preddeltaB,input$PA)
MAE(preddeltaB,input$PA)
preddeltaG<-predict(modfin, input , type = "response")
summary(preddeltaG)
R2(preddeltaG,input$sqrt.kg_km2)
MAE(preddeltaG,input$sqrt.kg_km2)

prob.delta= preddeltaB*preddeltaG^2
R2(prob.delta,input$kg_km2)
MAE(prob.delta,input$kg_km2)
summary(prob.delta)
summary(input$kg_km2)
# predict delta gam NO ENV
preddeltaB_noenv<-predict(modfin0_noenv, input, type = "response")
summary(preddeltaB_noenv)
R2(preddeltaB_noenv,input$PA)
MAE(preddeltaB_noenv,input$PA)
preddeltaG_noenv<-predict(modfin_noenv, input , type = "response")
summary(preddeltaG_noenv)
R2(preddeltaG_noenv,input$sqrt.kg_km2)
MAE(preddeltaG_noenv,input$sqrt.kg_km2)

prob.delta_noenv= preddeltaB_noenv*preddeltaG_noenv^2
R2(prob.delta_noenv,input$kg_km2)
MAE(prob.delta_noenv,input$kg_km2)
summary(prob.delta_noenv)
summary(input$kg_km2)
###########################
# best model Tweedie     #
###########################
modtw<- gam(kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth)+(TMP_bot)+s(dox.bot)+(pho)+s(sal),family= tw(),data=input, select=T, gamma=1.4)
summary(modtw)
AIC(modtw)
par(mfrow=c(3,3))
plot(modtw)
termplot(modtw)
par(mfrow=c(2,2))
gam.check(modtw)
#  model Tweedie NO ENV   
###########################
modtw_noenv<- gam(kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth),family= tw(),data=input, select=T, gamma=1.4)
summary(modtw_noenv)
# predict tweedie
predTW<-predict(modtw, input, type = "response")
summary(predTW)
R2(predTW,input$kg_km2)
MAE(predTW,input$kg_km2)
# predict tweedie NO ENV
predTW_noenv<-predict(modtw_noenv, input, type = "response")
summary(predTW_noenv)
R2(predTW_noenv,input$kg_km2)
MAE(predTW_noenv,input$kg_km2)

###########################
# best model GAussian     #
###########################
modGAU<- gam(sqrt.kg_km2~ +s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(nit),family= gaussian (link="identity"),data=input, select=T)
summary(modGAU)
AIC(modGAU)
par(mfrow=c(2,3))
plot(modGAU)
termplot(modGAU)
par(mfrow=c(2,2))
gam.check(modGAU)
# model GAussian  NO ENV   
###########################
modGAU_noenv<- gam(sqrt.kg_km2~ +s(X,Y)+s(year)+s(depth),family= gaussian (link="identity"),data=input, select=T)
summary(modGAU_noenv)
# predict su dati veri gaussian
predGAU2<-(predict(modGAU, input, type = "response"))^2
summary(predGAU2)
R2(predGAU2,input$kg_km2)
MAE(predGAU2,input$kg_km2)
# predict su dati veri gaussian NO ENV
predGAU2_noenv<-(predict(modGAU_noenv, input, type = "response"))^2
summary(predGAU2_noenv)
R2(predGAU2_noenv,input$kg_km2)
MAE(predGAU2_noenv,input$kg_km2)

###########################
# FINAL MODE= DELTA GAM   #
###########################


########################################################
#   PARTIRE DA QUI SE SI HANNO I DATABASE GIA FATTI    #
########################################################
#  MAPPE DISTRIBUZIONE con ggplot
##########################################################

#########################################
# Modello con env. var. da COPERNICUS   #
#########################################
##pred sul grid adriactico 
adr$pred = predict(modfin, adr, type = "response")^2;summary(adr$pred)
adr$pred0 = predict(modfin0, adr, type = "response");summary(adr$pred0)
# DELTA Approach moltiplico le abbondanze predette per la probabilità di PAenza del modello binomiale!!
adr$predfin<-adr$pred*adr$pred0
summary(adr$predfin); summary(input$kg_km2)

#########################################
# Modello Reference NO Env              #
#########################################
##pred sul grid adriactico 
adr$pred_noenv = predict(modfin_noenv, adr, type = "response")^2;summary(adr$pred_noenv )
adr$pred0_noenv  = predict(modfin0_noenv , adr, type = "response");summary(adr$pred0_noenv )
# DELTA Approach moltiplico le abbondanze predette per la probabilità di PAenza del modello binomiale!!
adr$predfin_noenv <-adr$pred_noenv *adr$pred0_noenv 
summary(adr$predfin_noenv ); summary(input$kg_km2)


####################
# MAPPE FINALI     #
####################
# mappa con env. var. da COPERNICUS
ggplot(data=world) +
  geom_sf() +
  coord_sf(xlim = c(11, 20.1875), ylim = c(41,45.9375), expand = FALSE)+
  geom_tile(data=adr,aes(x=X,y=Y, fill=adr$predfin))+
  scale_fill_gradient2(low = "blue" ,high = "red")+#, mid=scales::muted("cyan"),
  #scale_color_gradientn(colours = heat.colors)+
  scale_x_continuous(breaks = c(12, 16, 20))+
  scale_y_continuous(breaks = c(41, 43, 45))+
  theme_bw() +xlab("long")+ylab("lat")+
  theme(legend.title = element_text( size = 10))+
  labs(x=NULL,y=NULL,fill="Biomass (Kg/Km2)")+
  facet_wrap(~year,nrow=4)

# mappa NO ENV
ggplot(data=world) +
  geom_sf() +
  coord_sf(xlim = c(11, 20.1875), ylim = c(41,45.9375), expand = FALSE)+
  geom_tile(data=adr,aes(x=X,y=Y, fill=adr$predfin_noenv))+
  scale_fill_gradient2(low = "blue" ,high = "red")+#, mid=scales::muted("cyan"),
  #scale_color_gradientn(colours = heat.colors)+
  scale_x_continuous(breaks = c(12, 16, 20))+
  scale_y_continuous(breaks = c(41, 43, 45))+
  theme_bw() +xlab("long")+ylab("lat")+
  theme(legend.title = element_text( size = 10))+
  labs(x=NULL,y=NULL,fill="Biomass (Kg/Km2)")+
  facet_wrap(~year,nrow=4)


# diferenza tra env e no_env 
adr$diff <- adr$predfin - adr$predfin_noenv
ggplot(data=world) +
  geom_sf() +
  coord_sf(xlim = c(11, 20.1875), ylim = c(41,45.9375), expand = FALSE)+
  geom_tile(data=adr,aes(x=X,y=Y, fill=adr$diff))+
  scale_fill_gradient2(low = "blue",mid="white" ,high = "red",midpoint = 0)+
  #scale_color_gradientn(colours = heat.colors)+
  scale_x_continuous(breaks = c(12, 16, 20))+
  scale_y_continuous(breaks = c(41, 43, 45))+
  theme(legend.title = element_text( size = 10))+
  labs(x=NULL,y=NULL,fill="Diff in Biomass (Env-NO Env)")+
  theme_bw() +xlab("long")+ylab("lat")+
  facet_wrap(~year,nrow=4)

################################################################################################################

################################################################################################################


##########################
# EFF data 2008-2018  #
###########################
# best model DELTAGAM     #
###########################
input2<-input%>% dplyr::filter(year > 2007)
summary(input2)
## binomial part
####################
modfin0_eff<- gam(PA~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(pho)+(eff), family= binomial(),data=input2, select=T)
summary(modfin0_eff)
AIC(modfin0_eff)
par(mfrow=c(3,3))
plot(modfin0_eff)
termplot(modfin0_eff)
par(mfrow=c(2,2))
gam.check(modfin0_eff)

## Gaussuian part 
########################
modfin_eff<- gam(sqrt.kg_km2~ s(X,Y)+s(year)+s(depth)+TMP_bot+s(nit)+s(eff),family= gaussian (link="identity"),data=subset(input2,sqrt.kg_km2 >0), select=T)
summary(modfin_eff)
AIC(modfin_eff)
par(mfrow=c(2,3))
plot(modfin_eff)
termplot(modfin_eff)
par(mfrow=c(2,2))
gam.check(modfin_eff)

# predict delta gam
preddeltaB_eff<-predict(modfin0_eff, input2, type = "response")
summary(preddeltaB_eff)
R2(preddeltaB_eff,input2$PA)
MAE(preddeltaB_eff,input2$PA)
preddeltaG_eff<-predict(modfin_eff, input2 , type = "response")
summary(preddeltaG_eff)
R2(preddeltaG_eff,input2$sqrt.kg_km2)
MAE(preddeltaG_eff,input2$sqrt.kg_km2)

prob.delta_eff= preddeltaB_eff*preddeltaG_eff^2
R2(prob.delta_eff,input2$kg_km2)
MAE(prob.delta_eff,input2$kg_km2)
summary(prob.delta_eff)
summary(input2$kg_km2)

###########################
# best model Tweedie     #
###########################
modtw_eff<- gam(kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth)+(TMP_bot)+s(dox.bot)+(pho)+s(sal)+(eff),family= tw(),data=input2, select=T, gamma=1.4)
summary(modtw_eff)
AIC(modtw_eff)
par(mfrow=c(3,3))
plot(modtw_eff)
termplot(modtw_eff)
par(mfrow=c(2,2))
gam.check(modtw_eff)
# predict tweedie
predTW_eff<-predict(modtw_eff, input2, type = "response")
summary(predTW_eff)
R2(predTW_eff,input2$kg_km2)
MAE(predTW_eff,input2$kg_km2)


###########################
# best model GAussian     #
###########################
modGAU_eff<- gam(sqrt.kg_km2~ s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(nit)+(eff)+factor(month),family= gaussian (link="identity"),data=input2, select=T)
summary(modGAU_eff)
AIC(modGAU_eff)
par(mfrow=c(3,3))
plot(modGAU_eff)
termplot(modGAU_eff)
par(mfrow=c(2,2))
gam.check(modGAU_eff)
# predict su dati veri gaussian
predGAU2_eff<-(predict(modGAU_eff, input2, type = "response"))^2
summary(predGAU2_eff)
R2(predGAU2_eff,input2$kg_km2)
MAE(predGAU2_eff,input2$kg_km2)

###########################
# FINAL MODE= DELTA GAM   #
###########################


########################################################
#   PARTIRE DA QUI SE SI HANNO I DATABASE GIA FATTI    #
########################################################
#  MAPPE DISTRIBUZIONE con ggplot
##########################################################

#########################################
# Modello con env. var. da COPERNICUS + effort data  #
#########################################
##pred sul grid adriactico adr_eff
adr_eff$pred_eff = predict(modfin_eff, adr_eff, type = "response")^2;summary(adr_eff$pred_eff)
adr_eff$pred0_eff = predict(modfin0_eff, adr_eff, type = "response");summary(adr_eff$pred0_eff)
# DELTA Approach moltiplico le abbondanze predette per la probabilità di PAenza del modello binomiale!!
adr_eff$predfin_eff<-adr_eff$pred_eff*adr_eff$pred0_eff
summary(adr_eff$predfin_eff); summary(input$kg_km2)

# uso il modello no env anche su adr_eff
adr_eff$pred_noenv = predict(modfin_noenv, adr_eff, type = "response")^2;summary(adr_eff$pred_noenv )
adr_eff$pred0_noenv  = predict(modfin0_noenv , adr_eff, type = "response");summary(adr_eff$pred0_noenv )
# DELTA Approach moltiplico le abbondanze predette per la probabilità di PAenza del modello binomiale!!
adr_eff$predfin_noenv <-adr_eff$pred_noenv *adr_eff$pred0_noenv 
summary(adr_eff$predfin_noenv ); summary(input2$kg_km2)


####################
# MAPPE FINALI     #
####################
# mappa con env. var. da COPERNICUS e _eff
ggplot(data=world) +
  geom_sf() +
  coord_sf(xlim = c(11, 20.1875), ylim = c(41,45.9375), expand = FALSE)+
  geom_tile(data=adr_eff,aes(x=X,y=Y, fill=adr_eff$predfin_eff))+
  scale_fill_gradient2(low = "blue" ,high = "red")+#, mid=scales::muted("cyan"),
  #scale_color_gradientn(colours = heat.colors)+
  scale_x_continuous(breaks = c(12, 16, 20))+
  scale_y_continuous(breaks = c(41, 43, 45))+
  theme_bw() +xlab("long")+ylab("lat")+
  theme(legend.title = element_text( size = 10))+
  labs(x=NULL,y=NULL,fill="Biomass (Kg/Km2)")+
  facet_wrap(~year,nrow=4)


# diferenza tra env+_eff e no_env 
adr_eff$diff_eff <- adr_eff$predfin_eff - adr_eff$predfin_noenv
ggplot(data=world) +
  geom_sf() +
  coord_sf(xlim = c(11, 20.1875), ylim = c(41,45.9375), expand = FALSE)+
  geom_tile(data=adr_eff,aes(x=X,y=Y, fill=adr_eff$diff_eff))+
  scale_fill_gradient2(low = "blue",mid="white" ,high = "red",midpoint = 0)+
  #scale_color_gradientn(colours = heat.colors)+
  scale_x_continuous(breaks = c(12, 16, 20))+
  scale_y_continuous(breaks = c(41, 43, 45))+
  theme(legend.title = element_text( size = 10))+
  labs(x=NULL,y=NULL,fill="Diff in Biomass (Env+Eff-NO Env)")+
  theme_bw() +xlab("long")+ylab("lat")+
  facet_wrap(~year,nrow=4)
