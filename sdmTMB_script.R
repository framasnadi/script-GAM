#install.packages("sdmTMB", dependencies = TRUE)
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#remotes::install_github("pbs-assess/sdmTMBextra")
library(readr)
library(readr)
library(dplyr)
library(ggplot2)
library(sdmTMB)
library(sdmTMBextra) #package for the simulated residual plot
library(sf)
library(gridExtra)
library(rnaturalearth); library(rnaturalearthdata); library(rnaturalearthhires)

#setwd("D:/BIOLOGIA della PESCA/BIOLOGIA della PESCA/CNR/AIS_Prova_SOM")
dati <- read_csv("Ping_All.csv")
str(dati)
dati2 = as.data.frame(dati[,c("sum_NASC", "Distgroup", "Month", "Lat_M" ,"Lon_M")]) 
str(dati2)

dati2$log.NASC <-  log(dati2$sum_NASC+1)
head(dati2)
summary(dati2)
hist(dati2$log.NASC)

#Proceeding with UTM zone 33N; CRS = 32633.
dati2 <- add_utm_columns(dati2, ll_names = c("Lon_M", "Lat_M"), utm_names = c("X", "Y")) 
KNOT <- 50 # number of knot to draw the mesh, 50 is probably too low but faster to test preliminary models 
mesh <- make_mesh(dati2, xy_cols = c("X", "Y"), type = c("kmeans"), n_knots = KNOT)
#mesh <- make_mesh(dati2, xy_cols = c("X", "Y"), type = c("cutoff"), cutoff = 10)
plot(mesh)

# FIT a MODEL!
fit2 <- sdmTMB(
  log.NASC ~  1, #  s(Distgroup),
  data = dati2,
  mesh = mesh,
  family   = tweedie(link = "log"), #quite good, but you have to better explore other families
  spatial = "on",
  time = "Month",
  spatiotemporal = "IID" #default, but you have to better explore this
)
AIC(fit2)
sanity(fit2)
visreg::visreg(fit2, xvar = "Distgroup")
dati2$resids <- residuals(fit2) # randomized quantile residuals
qqnorm(dati2$resids)
qqline(dati2$resids)

#plot residuals
jpeg("resid_mesh500.jpeg",width = 120, height = 120, units = "mm", res = 1000)
ggplot(dati2, aes(X, Y, col = resids)) +
  scale_colour_gradient2() +
  geom_point() +
  facet_wrap(~Month) +
  coord_fixed()
dev.off()

#simulate residual (better option with mixed effect models)
simm <- simulate(fit2, nsim = 300)
simm |> 
  sdmTMBextra::dharma_residuals(fit2)

#simm <- simulate(fit2, nsim = 500)
#pred_fixed <- fit2$family$linkinv(predict(fit2)$est_non_rf)
#res_simm <- DHARMa::createDHARMa(
#  simulatedResponse = simm,
#  observedResponse = dati2$log.NASC,
#  fittedPredictedResponse = pred_fixed
#)
#DHARMa::testSpatialAutocorrelation(res_simm, x = dati2$X, y = dati2$Y)


#expand grid
summary(dati2)
Xx <- 618.4:692.5
Yy <- 6264:6411
#Distgroup <- unique(dati2$Distgroup)
month <- unique(dati2$Month)
ll <- list(X=Xx,Y=Yy , Month=month)
newdbpred <- expand.grid(ll)

#predict on grid
p <- predict(fit2, newdata = newdbpred)
newdbpred$pred  <- p$est
hist(dati2$log.NASC);hist(exp(newdbpred$pred))

#plot function
plot_map <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    coord_fixed()
}
plot_map(newdbpred, exp(pred)) +
  scale_fill_viridis_c(
    trans = "sqrt",
    # trim extreme high values to make spatial variation more visible
    na.value = "red", limits = c(0, quantile(exp(newdbpred$pred), 0.995))
  ) +
  facet_wrap(~Month) +
  ggtitle("Prediction (fixed effects + all random effects)",
          subtitle = paste("maximum estimated biomass density =", round(max(exp(newdbpred$pred))))
  )

###########
## Making pretty maps with sdmTMB output
map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "sweden")
# Crop the polygon for plotting and efficiency:
# st_bbox(map_data) # find the rough coordinates
bc_coast <- suppressWarnings(suppressMessages(
  st_crop(map_data,
          c(xmin = 16, ymin = 46, xmax =19 , ymax = 59.5))))

#transform our map into UTM 33 coordinates, which is the equal-area projection we fit in:
utm_zone33 <-  32633
bc_coast_proj <- sf::st_transform(bc_coast, crs = utm_zone33)
# ggplot(bc_coast_proj) + geom_sf()
#sf::st_boundary(bc_coast_proj) 
# Finally, we will combine our gridded predictions with the base map. We will multiply the X and Y columns by 1000 because we worked in UTM km for model fitting (to avoid computational issues with the scale of the range parameter):

pred_pl <- ggplot(bc_coast_proj)  +
  geom_raster(data = p, aes(x = X * 1000, y = Y * 1000, fill = exp(est))) +
    geom_sf() +
  xlim(600000, 700000 ) +
  ylim(6250000,6430000 ) +
  scale_fill_viridis_c(
  #  trans = "sqrt",
    # trim extreme high values to make spatial variation more visible
   # na.value = "yellow", limits = c(0, quantile(exp(p$est), 0.995))
  ) +
  theme_light() +
  labs(fill = "Predicted") +
  labs(x = "Longitude", y = "Latitude")+ facet_wrap(~Month)+ ggtitle("Prediction by month")

  #plot raw data
raw_pl <- ggplot(bc_coast_proj)  +
  geom_point(data = dati2, aes(x = X * 1000, y = Y * 1000, col = log.NASC)) +
  scale_colour_viridis_c()+
  geom_sf() +
  xlim(600000, 700000 ) +
  ylim(6250000,6430000 ) +
  theme_light() +
  labs(fill = "Predicted") +
  labs(x = "Longitude", y = "Latitude")+ facet_wrap(~Month)+ ggtitle("Raw data by month")

# save final comparison plot
jpeg("finalComparisonMAP.jpeg",width = 420, height = 320, units = "mm", res = 1000)
 grid.arrange( raw_pl,pred_pl, ncol=2)   
dev.off()




# center of gravity ## not working when there is no Year effect
#pp <- predict(fit2, newdata = newdbpred,return_tmb_object = TRUE)
#cog <- get_cog(pp, format = "wide") 
#ggplot(cog, aes(est_x, est_y, colour = Month)) +
#  geom_pointrange(aes(xmin = lwr_x, xmax = upr_x)) +
#  geom_pointrange(aes(ymin = lwr_y, ymax = upr_y)) +
#  scale_colour_viridis_c()
#plot function
plot_map <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    coord_fixed()
}
#expand grid

X <- 61:69
Y <- 62:65
#Distgroup <- unique(dati2$Distgroup)
griddepth <- expand_grid(X,Y)
griddepth$depth <- -1:-36

plot_map(griddepth, (depth)) +
  scale_fill_viridis_c(  )  +
  ggtitle("Prediction (fixed effects + all random effects)"
  )


Year <- 2000:2008
effect <- c("YES","NO")

newdbpred.step1 <- expand_grid( X = unique(griddepth$X) ,Y = unique(griddepth$Y), Year=Year,  effect=effect)

prediction_dataset <- newdbpred.step1 %>%
  left_join(griddepth, by = c("X", "Y"))

plot_map(prediction_dataset, depth) +
  scale_fill_viridis_c(  )  +
  ggtitle("Prediction"
  ) + facet_wrap(~Year)

