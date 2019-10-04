#### Mappe distribuzione Solemon con GAM

library(mgcv)
library(maps)
library(mapdata)
library(mapplots)
library(shapefiles)
library(dplyr)
setwd("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Mappe distribuzione Solemon")
library(readr)
Solea_2019 <- read_delim("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Mappe distribuzione Solemon/Solea_with_grain.csv", 
                                    ",", escape_double = FALSE, trim_ws = TRUE)
#View(Solea_2019)

data<-Solea_2019
summary(data)
colSums(is.na(data))

summary(data$n_km2)
summary(data$n_sqrt)

## carico lo shape della GSA
# MAPS 
GSAtable <<- read.table("gsa_coordinates2.csv", sep=";", header=T)
leg_pos <<- as.character(GSAtable[GSAtable$GSA==19,6])
coord_lon <<- as.numeric( GSAtable[GSAtable$GSA==17,c(2,3)] )
coord_lat <<- as.numeric( GSAtable[GSAtable$GSA==17,c(4,5)] )
GSA17<-read.shapefile("Med_Poly - Copia/Med_Poly")

#carico il grid adriatico da Emodnet
library(readr)
adr <- read_csv("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Mappe distribuzione Solemon/GridADR_with_grain.csv")
#View(adr)
plot(adr$lon,adr$lat, cex=0.1)

source('C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Mappe distribuzione Solemon/collinearity.R')
###modello gam SOLEA
Z2<-cbind( data$n_sqrt,  data$X, data$Y, data$depth, data$grain, data$shooting_time)
## creo i nomi da visualizzare nel grafico
colnames(Z2) <- c("solea","lon","lat","depth", "grain", "shootingtime")
pairs(Z2, lower.panel = panel.smooth2,
      upper.panel = panel.cor, diag.panel = panel.hist)


####
hist(data$n_km2)
plot(data$n_km2)
max(data$n_km2);min(data$n_km2)
summary(data$n_km2)
hist(data$n_sqrt) 
summary(data$n_sqrt)
#modfin<- gam(n_sqrt~  s(X,Y)+s(depth), 
#                 family=gaussian,data=data)

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

#####################
##binomial part
####################
mod<-gam(pres~   s(grain),family= binomial(),data=data, select=T)
summary(mod)
plot(mod)
AIC(mod)
gam.check(mod)
#modfin0<- gam(pres ~  s(grain)+ s(depth) +s(X,Y)+s(year, k=13) +te(Y,year)+te(X,year), 
#             family= binomial(),data=data, select=T)
modfin0<- gam(pres ~ s(grain) +s(depth) +s(X,Y) +s(year)+te(Y,year)+te(X,year),
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


p<-predict(modfin0, data,type = "response");summary(p)
p<-round(p,0);summary(p)
table(p , data$pres)
confusionMatrix(table(p, data$pres))
########################
## Gaussuian part
########################

#modfin<- gam(n_sqrt~  + s(depth) +s(X,Y)+s(year, k=13) +s(hour)+te(Y,year)+te(X,year), 
#            family= gaussian (link="identity"),data=abb, select=T)
modfin<- gam(n_sqrt~  s(grain) + s(depth) +s(X,Y)+s(year) +te(Y,year)+te(X,year), 
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

library(caret)
abb$p<-predict(modfin, abb);summary(abb$p)
abb$p<-round(abb$p,0);summary(abb$p)
abb$valid<-round(abb$n_sqrt,0);summary(abb$valid)
table(abb$p , abb$valid)
abb$valid == 2
abb$valid == 10
zzz= abb[-c(7,89),];summary(zzz)
table(zzz$p , zzz$valid)
confusionMatrix(table(zzz$p, zzz$valid))
##confronto con i dati solemon veri per vedere se il modello è buono
#data$prdsol<-exp(predict(modfin, data, type = "response"))-1;summary(data$prdsol)
#n<-exp(data$n_sqrt)-1;summary(n)
#anno<-data$year
#sol1<-data.frame(data$prdsol,n,anno)
#write.csv(sol1, "sol_pred_distr.csv")
#View(sol1)
#pred<-tapply(sol1$data.prdsol, anno, mean)
#reale<-tapply(sol1$n, anno, mean)
#sol2<-data.frame(pred,reale,unique(anno));sol2
#mediaN_per_cala<-tapply(data$n_km2, data$haul_number, mean);summary(mediaN_per_cala)
#mediaN_per_cala.pred<-tapply(data$prdsol, data$haul_number, mean);summary(mediaN_per_cala.pred)

########################################################
#   PARTIRE DA QUI SE SI HANNO I DATABASE GIA FATTI    #
########################################################
#  MAPPE DISTRIBUZIONE con ggplot
##########################################################

##pred sul grid adriactico (indipendentemente dagli anni)
adr$X<-adr$lon
adr$Y<-adr$lat
adr$grain<-adr$grain_mask
# setto l'anno della mappa!!
adr$year<-rep(2005, 39754)
pred = predict(modfin, adr, type = "response")
predd = predict(modfin0, adr, type = "response")
adr$pred = pred;summary(adr$pred)
adr$predd = predd;summary(adr$predd)
#adr$l<-round(adr$pred,0);summary(adr$l)
#unique(adr$l)
#adr$l2<-adr$pred+3;summary(adr$l2)
#adr$p<-round(exp(adr$pred),0);summary(adr$p)
adr$pred2<-exp(adr$pred)
#adr$pred3<-exp(adr$predd)
summary(adr$pred2)
#summary(adr$pred3)
   #par(mfrow=c(1,2)) 
   #hist(adr$pred2);plot(adr$pred2)
#adr$pp<-(adr$pred2/68);summary(adr$pp)

# DELTA Approach moltiplico le abbondanze predette per la probabilità di presenza del modello binomiale!!
adr$predfin<-adr$pred2*adr$predd
summary(adr$predfin); summary(data$n_km2)





# si trasformano i nem_km2 in classi di "presenza"
adr$abbsol <- rep(0, dim(adr)[1])
for (i in 1:(dim(adr)[1])) {
    if(adr$predfin[i] > 0 &  adr$predfin[i] <= 1) 
    {adr$abbsol[i] <- "(M) <1"} else
      if(adr$predfin[i] > 1 &  adr$predfin[i] <= 100) 
      {adr$abbsol[i] <- "(L) 1-100"} else
        if(adr$predfin[i] > 100 &  adr$predfin[i] <= 200) 
        {adr$abbsol[i] <- "(I) 100-200"} else
          if(adr$predfin[i] > 200 &  adr$predfin[i] <= 300) 
          {adr$abbsol[i] <-"(H) 200-300"} else
            if(adr$predfin[i] > 300 &  adr$predfin[i] <= 400) 
            {adr$abbsol[i] <- "(G) 300-400"} else
              if(adr$predfin[i] > 400 &  adr$predfin[i] <= 500) 
              {adr$abbsol[i] <- "(F) 400-500"} else
                if(adr$predfin[i] > 500 &  adr$predfin[i] <= 700) 
                {adr$abbsol[i] <- "(E) 500-700"} else
                  if(adr$predfin[i] > 700 &  adr$predfin[i] <= 1000) 
                  {adr$abbsol[i] <- "(D) 700-1000"} else
                    if(adr$predfin[i] > 1000 &  adr$predfin[i] <= 1500) 
                    {adr$abbsol[i] <- "(C) 1000-1500"} else
                      if(adr$predfin[i] > 1500 &  adr$predfin[i] <= 2500) 
                  {adr$abbsol[i] <- "(B) 1500-2500"} else
          if(adr$predfin[i] > 2500 ) 
     {adr$abbsol[i] <- "(A) >2500"} 
}
str(adr)
summary(adr)
unique(adr$abbsol)
boxplot(adr$predfin~ adr$abbsol)
#
#adr$abbsol3 <- rep(0, dim(adr)[1])
#for (i in 1:(dim(adr)[1])) {
#  if(adr$pred3[i] > 0 &  adr$pred3[i] <= 1) 
#  {adr$abbsol3[i] <- "(M) <1"} else
#    if(adr$pred3[i] > 1 &  adr$pred3[i] <= 100) 
#    {adr$abbsol3[i] <- "(L) 1-100"} else
#      if(adr$pred3[i] > 100 &  adr$pred3[i] <= 200) 
#      {adr$abbsol3[i] <- "(I) 100-200"} else
#        if(adr$pred3[i] > 200 &  adr$pred3[i] <= 300) 
#        {adr$abbsol3[i] <-"(H) 200-300"} else
#          if(adr$pred3[i] > 300 &  adr$pred3[i] <= 400) 
#          {adr$abbsol3[i] <- "(G) 300-400"} else
#            if(adr$pred3[i] > 400 &  adr$pred3[i] <= 500) 
#            {adr$abbsol3[i] <- "(F) 400-500"} else
#              if(adr$pred3[i] > 500 &  adr$pred3[i] <= 700) 
#              {adr$abbsol3[i] <- "(E) 500-700"} else
#                if(adr$pred3[i] > 700 &  adr$pred3[i] <= 1000) 
#                {adr$abbsol3[i] <- "(D) 700-1000"} else
#                  if(adr$pred3[i] > 1000 &  adr$pred3[i] <= 1500) 
#                  {adr$abbsol3[i] <- "(C) 1000-1500"} else
#                    if(adr$pred3[i] > 1500 &  adr$pred3[i] <= 2500) 
#                    {adr$abbsol3[i] <- "(B) 1500-2500"} else
#                      if(adr$pred3[i] > 2500 ) 
#                      {adr$abbsol3[i] <- "(A) >2500"} 
#}
#
adr$abbsol2 <- rep(0, dim(adr)[1])
for (i in 1:(dim(adr)[1])) {
  if(adr$predfin[i] > 0 &  adr$predfin[i] <= 10) 
  {adr$abbsol2[i] <- "(D) <10"} else
    if(adr$predfin[i] > 10 &  adr$predfin[i] <= 150) 
    {adr$abbsol2[i] <- "(C) 1-150"} else
      if(adr$predfin[i] > 150 &  adr$predfin[i] <= 350) 
      {adr$abbsol2[i] <- "(B) 150-350"} else    
      {adr$abbsol2[i] <- "(A) >350"} 
}
str(adr)
summary(adr)
unique(adr$abbsol2)
boxplot(adr$predfin~ adr$abbsol2)
#
# adr$abbsol4 <- rep(0, dim(adr)[1])
# for (i in 1:(dim(adr)[1])) {
#   if(adr$pred3[i] > 0 &  adr$pred3[i] <= 10) 
#   {adr$abbsol4[i] <- "(D) <10"} else
#     if(adr$pred3[i] > 10 &  adr$pred3[i] <= 150) 
#     {adr$abbsol4[i] <- "(C) 1-150"} else
#       if(adr$pred3[i] > 150 &  adr$pred3[i] <= 350) 
#       {adr$abbsol4[i] <- "(B) 150-350"} else    
#       {adr$abbsol4[i] <- "(A) >350"} 
# }
#######
#plot(adr$X, adr$Y, col=adr$abbsol, cex=0.5, main="2018", xlab = "Longitudine", ylab = "Latitudine",xlim = c(12,15.5))
#draw.shape(GSA17, col = "white", add=T)
#draw.shape(bati100, col = "white")
#legend("bottomleft", title="N km2", c("<1", "1-5","5-12","12-33",">33"), bt="n", box.col="black", pch=c(15), col=c(2,3,4,5,6),pt.cex=c(1.6), cex=0.9, y.intersp=1.1, xjust= 0)
#points(data$lon, data$lat , cex=data$n_sqrt/2, col="black")


library(ggplot2)
#library(rgdal)
library(readr)
map1<-map_data("worldHires") %>% filter(long > 11 & long < 18 & lat > 40 & lat < 48)
Adriatic_Grid <- read_csv("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Mappe distribuzione Solemon/Adriatic_Grid2.csv")
bt100<-Adriatic_Grid%>%  filter(depth<2000 & depth>100 & GSA==17 & lon<18&lon>12,lat>41.5& lat<44)

x11()
ggplot(adr,aes(X,Y,fill=adr$abbsol))+theme_light()+geom_raster(interpolate=T,show.legend = T)+scale_fill_grey()+coord_quickmap(xlim=c(12,15.25),ylim=c(41.9,45.8)) +theme(legend.title = element_text( size = 10))+
  labs(x=NULL,y=NULL,fill="Abbondanza (Num/Km2)")+geom_polygon(data=bt100,aes(lon,lat),fill=NA,color= NA )+geom_polygon(data=map1,aes(long,lat,group=group),fill="lightgoldenrodyellow",color="black")+ggtitle("DeltaGam 2005 with grain")
#+coord_fixed()
savePlot(filename = paste("DeltaGam 2005 with grain",sep=""),type = "tiff",device = dev.cur(),restoreConsole = T)

dev.off()









