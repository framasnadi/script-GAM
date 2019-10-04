library(tiff)
library(raster)
library(readr)

Adriatic_Grid <- read_csv("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Mappe distribuzione Solemon/Adriatic_Grid2.csv")
#View(Adriatic_Grid)
#plot(Adriatic_Grid$lon,Adriatic_Grid$lat, cex=0.1)
library(dplyr)
adr<-Adriatic_Grid%>%  filter(depth<100 & depth>1 & GSA==17 & lon<16&lon>12,lat>41.5)
summary(adr)
length(adr$lon) #54950


# carico la mappa della granulometria
gr<-'grain_mask.tif' 
gr<-raster(gr)
plot(gr)
grp<-rasterToPoints(gr)
grp<-as.data.frame(grp)

grp$lon<-round(grp$x,2)
grp$lat<-round(grp$y,2)
grp <- grp%>% distinct(lon, lat,  .keep_all = TRUE)
# cerco di mergare adr e grp per avere un grid finale
######

adrgrid = merge(adr,grp, by.x=c("lat","lon"), by.y=c("lat","lon"))
View(adrgrid)
summary(adrgrid)
write.csv(adrgrid, "GridADR_with_grain.csv")

plot(adrgrid$lon,adrgrid$lat, cex=log(adrgrid$grain_mask)/2)

#####################################################################
# INSERIRE GRAIN SIZE NELLE CALE SOLEMON per SOLEA   ########################

Solea_stnd_20_6_2019 <- read_delim("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Standardization file/Solemon/Solea stnd (20-6-2019).csv", 
                                   ";", escape_double = FALSE, trim_ws = TRUE)
View(Solea_stnd_20_6_2019)

data<-Solea_stnd_20_6_2019
summary(data)
data$depth <- (data$shooting_depth+data$hauling_depth)/2
data$n_sqrt <- log(data$n_km2+1)
#carico il grid adriatico appena creato
library(readr)
adrgrid <- read_csv("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Mappe distribuzione Solemon/GridADR_with_grain.csv")

# riporto le coordinate di solemon a solo 2 cifre decimali per fare merge con adrigrid
plot(data$X,data$Y, cex=data$shooting_quadrant)
draw.shape(GSA17, col = "white")
data$X<-round(data$X,2)
data$Y<-round(data$Y,2)
points(data$X, data$Y , cex=data$shooting_quadrant, col="red")

adrgrid$Y<-adrgrid$lat
adrgrid$X<-adrgrid$lon

adrgrid<- adrgrid[,c(6, 9, 10, 11)]
plot(adrgrid$X,adrgrid$Y, cex=0.35)
#adrgridcut  <- adrgrid %>% distinct(X, Y,  .keep_all = TRUE)
data<-data[,c("year","haul_number", "n_km2" ,"depth","X","Y","n_sqrt","vessel", "shooting_time")]


solea_grain<- merge(data,adrgrid , by.x=c("X","Y"), by.y=c("X","Y"), all.x = T)
colSums(is.na(solea_grain))
solea_grain$depth<-solea_grain$depth.y
solea_grain$grain<-solea_grain$grain_mask

solea_grain2  <- solea_grain %>% dplyr::select(X, Y,year,haul_number, n_km2 ,depth.x,n_sqrt,vessel, shooting_time,grain,depth) %>% na.omit()

write.csv(solea_grain2, "Solea_with_grain.csv")




# INSERIRE GRAIN SIZE NELLE CALE SOLEMON per SQUILLA   ########################

datasquilla <- read_delim("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Standardization file/Solemon/dataMTS.csv", 
                                   ",", escape_double = FALSE, trim_ws = TRUE)
View(datasquilla)

data<-datasquilla
summary(data)
data$depth <- (data$shooting_depth+data$hauling_depth)/2
data$n_sqrt <- log(data$n_km2+1)
#carico il grid adriatico appena creato
library(readr)
adrgrid <- read_csv("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/SOLEMON/Mappe distribuzione Solemon/GridADR_with_grain.csv")

# riporto le coordinate di solemon a solo 2 cifre decimali per fare merge con adrigrid
plot(data$X,data$Y, cex=data$shooting_quadrant)
draw.shape(GSA17, col = "white")
data$X<-round(data$X,2)
data$Y<-round(data$Y,2)
points(data$X, data$Y , cex=data$shooting_quadrant, col="red")

adrgrid$Y<-adrgrid$lat
adrgrid$X<-adrgrid$lon

adrgrid<- adrgrid[,c(6, 9, 10, 11)]
plot(adrgrid$X,adrgrid$Y, cex=0.35)
#adrgridcut  <- adrgrid %>% distinct(X, Y,  .keep_all = TRUE)
data<-data[,c("year" ,"depth","X","Y","country","MELIKER","SOLEVUL", "SQUIMAN", "SEPIOFF")]
data<-data %>% dplyr::filter(country== "ITA17")

sp_grain<- merge(data,adrgrid , by.x=c("X","Y"), by.y=c("X","Y"), all.x = T)
colSums(is.na(sp_grain))
sp_grain$depth<-sp_grain$depth.y
sp_grain$grain<-sp_grain$grain_mask

sp_grain2  <- sp_grain %>% dplyr::select(X, Y,year ,depth,grain,country,MELIKER,SOLEVUL,SQUIMAN,SEPIOFF) %>% na.omit() 
#sp_grain2  <- sp_grain %>% dplyr::select(X, Y,year ,depth,grain,country,MELIKER,SOLEVUL,SQUIMAN,SEPIOFF) %>% na.omit() %>% filter(country== "ITA17")
write.csv(sp_grain2, "specie_with_grain.csv")



