###############################################
# TGAM CONT 
###################################################
# remotes::install_github("saskiaotto/INDperform")
setwd("C:/Users/frma6502/Desktop/SU_Postdoc/Herring Formas/Analisi/variables for PLSSEM/FINAL_DBs")
library(pls)
library(dplyr)
library(readxl)
library(gridExtra)
library(imputeTS)
library(INDperform)
setwd("C:/Users/frma6502/Desktop/SU_Postdoc/Herring Formas/Analisi/variables for PLSSEM/FINAL_DBs")

db_PLSR <- read_excel("Her_FORMAS_db_orig.xlsx", 
                      sheet = "LAND_Ba_Fall")%>% filter(Year > 1994) 
#b <- "BaP"
#ss <- "Fall" 
#plotdir <- paste("output/output",b,ss, "/PLSR/TGAM/", sep="")

# data handling
db_PLSR[14,5:41]  <- NA # for All var in BaP
db_PLSR[2,5] <- NA # for PCB_153 in BaP
#db_PLSR[14,5] <- NA # for PCB_153 in BoS Fall
#db_PLSR[4,7] <- NA # for TCDD_EQV in BoS Fall
#db_PLSR[14,7] <- NA # for PCB_153 in BoS Fall
#db_PLSR[14,9] <- NA # for Cd in BoS Fall
#db_PLSR[21,9] <- NA # for Cd in BoS Fall
#db_PLSR[14,21] <- NA # for Wg in BoS Fall
#db_PLSR[21,21] <- NA # for Wg in BoS Fall
#db_PLSR[21,22] <- NA # for Age in BoS Fall
#db_PLSR[23,22] <- NA # for Age in BoS Fall
#db_PLSR[23,21] <- NA # for Age in BoS Spring
#db_plspmNAint$TP_Nielsen <- db_plspmNAint$TP_Chikaraishi # only for BoS Spring!!!
db_PLSR$Phe_proxyEutr <- -db_PLSR$Phe_proxyEutr # incresing in eutrop proxy
# deal with missing value 
db_PLSR <- na_interpolation(db_PLSR) # ACTUNG tutti gli NA diventano interpolation dei due valori prima e dopo nella timeseries!!!
#attach(db_PLSR) # detach(db_PLSR)

y <- "PCB_153"
ind_ex <- db_PLSR %>% dplyr::select(Year, PCB_153) #TCDD_EQV, Cu, Cd, Hg)
press_ex <- db_PLSR %>% dplyr::select(Year, AirEmi_PCB_153,	Age	, TP_Nielsen 	,d13C_range,	d15N_range,	SEAc, Temp, Sal,Phe_proxyEutr,Cyano,FCA,	Phyto,	Monoporeia)

head(ind_ex)
head(press_ex)

m_trend <- model_trend(ind_tbl = ind_ex[ ,-1],
                       time = db_PLSR$Year)
# Model diagnostics
pd <- plot_diagnostics(model_list = m_trend$model)
pd$all_plots[[1]] # first indicator
# check for outliers in all models
grid.arrange(grobs = pd$cooks_dist, ncol = 3) 
# check normality in all models
grid.arrange(grobs = pd$qq_plot, ncol = 3)
# check homogeneity in all models
grid.arrange(grobs = pd$resid_plot, ncol = 3)
# check for autocorrelation in all models
grid.arrange(grobs = pd$acf_plot, ncol = 3)
# check for partial autocorrelation in all models
grid.arrange(grobs = pd$pacf_plot, ncol = 3)

# Save diagnostic plots per indicator in a pdf
ml <- marrangeGrob(grobs = pd$all_plots, ncol = 1, nrow = 1)
ggsave("Trend_diagnostics_PCB_153.pdf", ml, height = 8, width = 12)


# Inspect trends
pt <- plot_trend(m_trend)
pt$PCB_153 # shows trend of TZA indicator
# show all together
ml <- grid.arrange(grobs = pt, ncol = 1)
ml
# save as pdf
ggsave("Trend_plots_PCB_153.pdf", ml, height = 12, width = 15)
# show only significant trends
grid.arrange(
  grobs = pt[which(m_trend$p_val <= 0.05)],
  ncol = 1
)

#Combine IND with pressures and select training and test period –> default is training data = first 90% of obs)
dat_init <- ind_init(ind_tbl = ind_ex[ ,-1],
                     press_tbl = press_ex[ ,-1], time = ind_ex$Year,
                     train = 0.9, random = FALSE)
````

#### B.2a Model responses using simple GAMs (using default settings here)
m_gam <- model_gam(init_tbl = dat_init, k = 5,
                   family = stats::gaussian(), excl_outlier = NULL)

# Model diagnostics
pd <- plot_diagnostics(model_list = m_gam$model)  # (might take a while)
ml <- marrangeGrob(grobs = pd$all_plots, ncol = 1, nrow = 1)
#ggsave("GAM_diagnostics.pdf", ml, height = 8, width = 12)

# Any outlier?
m_gam$pres_outlier %>% purrr::compact(.)
# - get number of models with outliers detected
purrr::map_lgl(m_gam$pres_outlier, ~!is.null(.)) %>% sum()
# - which models and what observations?
m_gam %>%
  dplyr::select(id, ind, press, pres_outlier) %>%
  dplyr::filter(!purrr::map_lgl(m_gam$pres_outlier, .f = is.null)) %>%
  tidyr::unnest(pres_outlier)

# Exclude outlier in models (using the returned model output tibble as selector)
m_gam <- model_gam(init_tbl = dat_init, excl_outlier = m_gam$pres_outlier)
pd <- plot_diagnostics(model_list = m_gam$model)  # (might take a while)
ml <- marrangeGrob(grobs = pd$all_plots, ncol = 1, nrow = 1)
ggsave("GAM_diagnostics_PCB_153.pdf", ml, height = 8, width = 12)

# Any temporal autocorrelation (TAC) 
sum(m_gam$tac)
# - which models
m_gam %>%
  dplyr::select(id, ind, press, tac) %>%
  dplyr::filter(tac)
#apply GAMM
m_gamm <- model_gamm(init_tbl = dat_init,
                     filter = m_gam$tac) # (apply GAMM only to rows where TAC detected)

# Again, any outlier?
purrr::map_lgl(m_gamm$pres_outlier, ~!is.null(.)) %>% sum()
# Select best GAMM from different correlation structures
# (based on AIC)
best_gamm <- select_model(gam_tbl = m_gam, gamm_tbl = m_gamm)
#View(best_gamm)

# Still any temporal autocorrelation?
sum(best_gamm$tac) 
# - which models
best_gamm %>%
  dplyr::select(id, ind, press, tac, corrstruc) %>%
  dplyr::filter(tac) %>%
  print(n = 100)

# Inspect diagnostic plots of best GAMMs
pd <- plot_diagnostics(model_list = best_gamm$model)
ml <- marrangeGrob(grobs = pd$all_plots, ncol = 1, nrow = 1)
ggsave("bestGAMM_diagnostics_PCB_153.pdf", ml, height = 8, width = 12)

#B.2c Merge GAM and GAMMs
m_merged <- merge_models(m_gam[m_gam$tac == FALSE, ], best_gamm)

# View sign. models
filter(m_merged, p_val <= 0.05) %>%
  View() 

#B.3 Calculate derivatives of (significant) non-linear responses
m_calc <- calc_deriv(init_tbl = dat_init, mod_tbl = m_merged,
                     sign_level = 0.05)

#B.4 Test for pressure interactions (in significant models only=default)
it <- select_interaction(mod_tbl = m_calc)
# (creates combinations to test for)

# takes long time!!!
m_all <- test_interaction(init_tbl = dat_init, mod_tbl = m_calc,
                          interactions = it, sign_level = 0.05)

# Inspect diagnostics of threshold models
pd_thresh <- m_all$thresh_models %>%
  # flatten structure of nested threshold GAMs 
  flatten() %>%
  # and remove empty lists (where no threshold GAMs were applied (=NULL) 
  # or no threshold GAMs were better than the simple GAMs (=NA)
  keep(~is.list(.x))  %>%
  plot_diagnostics()

ml <- marrangeGrob(grobs = pd_thresh$all_plots, ncol = 1, nrow = 1)
ggsave("thresholdGAM_diagnostics_PCB_153.pdf", ml, height = 8, width = 10)

# Look at the development of the generalized cross-validation value 
# at different thresholds level: the plot should show that the GCVV
# of the final threshold should be clearly lower than of others 
# (you should see a single sharp trough, not at the edge)
grid.arrange(grobs = pd_thresh$gcvv_plot, ncol = 3)
# -> some models do NOT show this (e.g. #4,6,7): here we could set the
#  $interaction from TRUE to FALSE (for the scoring later) to ignore
# this potentially spurious interaction found:
#m_all <- m_all %>% mutate(
#  interaction = case_when(
#    (ind == "Micro" & press == "Fsprat") |
#      (ind == "Cod" & press == "Fher") | 
#      (ind == "Cod" & press == "Fcod") ~ FALSE,
#    TRUE ~ interaction))

# Show final significant GAM/GAMM (threshold GAM) plots:
sel <- which(m_all$p_val <= 0.05)
pm <- plot_model(init_tbl = dat_init[sel, ], mod_tbl = m_all[sel, ])
# Save all sign. indicator plots
ml <- gridExtra::marrangeGrob(grobs = pm$all_plots, ncol = 1, nrow = 1)
ggplot2::ggsave("Final_model_results_PCB_153.pdf", ml, height = 10, width = 12)

# save result table
m_all.tb <-m_all[, c(1:16, 30,31)]
write.xlsx((m_all.tb) , "PCB_153_restable.xlsx")

######################################
# C. Scoring based on model output
press_type_ex <- read_excel("output/info.xlsx", sheet = "TGAMind")
scores <- scoring(trend_tbl = m_trend, mod_tbl = m_all, 
                  press_type = press_type_ex)

# Runs a shiny app to modify the score for the subcriterion 10.1:
scores <- expect_resp(mod_tbl = m_all, scores_tbl = scores)

# Generate an easy to read summary from the nested scoring tibble
sum_sc <- summary_sc(scores)
sum_sc

#Visualize scores with a spie chart (using the summary object as input)
spie <- plot_spiechart(sum_sc,
                       lab_size = 4, title_size = 4    )
spie$PCB_153 # shows the spiechart of the indicator TZA
# show all:
#grid.arrange(grobs = spie, ncol = 3)
