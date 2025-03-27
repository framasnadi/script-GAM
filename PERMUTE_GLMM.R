# Load required libraries
library(glmmTMB)
library(emmeans)
library(readr)
library(readxl)
library(tidyr)
library(dplyr)
library(car)
library(DHARMa)
library(ggplot2)

dir <- "C:/Users/frma6502/Desktop/SU/DO_STCKBCK_project/Results"
setwd(dir)
dirMODEL2 <- paste0(dir ,"/NewMethod/output2_GLMM" ) 


db_Nod_Biom <- read_excel("C:/Users/frma6502/Desktop/SU/DO_STCKBCK_project/Results/db_indicator_capt.xlsx")

db_Nod_Biom$EROD<- as.numeric(db_Nod_Biom$EROD)
db_Nod_Biom$GR <- as.numeric(db_Nod_Biom$GR)
db_Nod_Biom$Sampling <- as.factor(db_Nod_Biom$Sampling)
db_Nod_Biom$TankN <- as.factor(db_Nod_Biom$TankN)
db_Nod_Biom$Treatment <- factor(db_Nod_Biom$Treatment, levels = c("Control", "NO-TOX", "LO-MIX", "HI-MIX", "TOXIC"))
#
df_long <- pivot_longer(db_Nod_Biom, cols= c(10:19), values_to = "value", names_to = "Endpoint")  
# trasform the endpoint withbefore the "leveling"
mean.End1 <- df_long %>% filter(Sampling == "1") %>% group_by(Treatment,Endpoint)  %>% dplyr::summarise(value_end1 = mean(value, na.rm = T) , sqrtvalue_end1 = mean(sqrt(value), na.rm = T),cubicvalue_end1 = mean((value)^(1/3), na.rm = T)  )
df_long <- left_join(df_long, mean.End1)
df_long$level <- ((df_long$value / df_long$value_end1)-1)*100  # la stessa cosa di quanto sotto preso da Makaras et al 2020
df_long$eq1makaras <- ((df_long$value - df_long$value_end1)/df_long$value_end1)  # eq 1 in makaras 2020

levels(df_long$Sampling)[levels(df_long$Sampling) == "1"] <- "S1"
levels(df_long$Sampling)[levels(df_long$Sampling) == "2"] <- "S2"
levels(df_long$Sampling)[levels(df_long$Sampling) == "3"] <- "S3"
str(df_long)


###########################################
# Permutational GLMM
############################################
TIPO <- "TRT" # "SMPL"
for (i in unique(df_long$Endpoint)   )   {
# Fit the model
model <- glmmTMB(level ~ Treatment*Sampling +(  Sampling   | TankN) ,
                 data=df_long %>% dplyr::filter(Endpoint == i),  REML=F )
jpeg(paste0(dirMODEL2,"/PERMUTE/Simulated_Resid_GLMM_" ,i,  ".jpeg"),width = 300, height = 170, units = "mm", res = 600)
residuals_sim <- simulateResiduals(fittedModel = model,n = 1000) %>% plot()
dev.off()

# Observed LSMeans and pairwise comparisons
#observed_lsmeans <- emmeans(model, ~ Sampling |Treatment )    # SMPL
observed_lsmeans <- emmeans(model, ~ Treatment  |Sampling )   # TRT
observed_pairs <- pairs(observed_lsmeans)

# Function for grouped permutation
permute_glmmTMB <- function(data, formula, family, group_var, response_col) {
  # Permute the response variable within groups
  data[[response_col]] <- unlist(lapply(split(data[[response_col]], data[[group_var]]), sample))
  
  # Refit the model with permuted data
  perm_model <- glmmTMB(formula, family = family, data = data,  REML=F)
  
  # Extract LSMeans (marginal means) tukey HDS
# emmeans(perm_model, ~ Sampling | Treatment)# post-hoc by Treatment in Sampling; SMPL
  emmeans(perm_model, ~ Treatment | Sampling) # post-hoc by Treatment in Sampling; TRT
}


# Number of permutations
n_possible_permutations <- factorial(length(unique(df_long$TankN)))  # Calculate total unique permutations
# Choose number of permutations (e.g., all possible or a random subset)
n_perm <- ifelse(n_possible_permutations <= 1000, n_possible_permutations, 1000)
#n_perm <- 24

# Perform permutations and extract pairwise differences
set.seed(123)
perm_lsmeans <- replicate(n_perm, permute_glmmTMB(
  data = df_long %>% dplyr::filter(Endpoint == i),
  formula = level ~ Treatment*Sampling + (Sampling | TankN),
  family = gaussian,
  group_var = "TankN",
  response_col = "level"
), simplify = FALSE)

# Extract pairwise differences for each permutation
#perm_diffs <- sapply(perm_lsmeans, function(x) pairs(x)$contrast)
perm_diffs <- do.call(rbind, lapply(perm_lsmeans, function(perm_result) {
  pairwise_result <- pairs(perm_result)  # Get pairwise comparisons
  summary(pairwise_result)$estimate  # Extract the pairwise differences
}))

# Calculate permutation-based p-values for each observed difference
observed_diffs <- summary(observed_pairs)$estimate  # Observed pairwise differences
p_values <- sapply(seq_along(observed_diffs), function(I) {
  mean(abs(perm_diffs[, I]) >= abs(observed_diffs[I]))
})

# Combine observed pairwise differences and p-values
posthoc_results <- data.frame(
  Sampling =      observed_pairs@grid[["Sampling"]],
#  Treatment  =  observed_pairs@grid[["Treatment"]],
  Comparison = summary(observed_pairs)$contrast,  # Pairwise comparisons
  Observed_Difference = observed_diffs,           # Observed differences
  Permutation_p_value = p_values                  # Permutation-based p-values
)
# Print results
# save GLMM result
sink(paste0(dirMODEL2,"/PERMUTE/PERMUTE_output_GLMM_",i,"_",TIPO,".txt", sep=""))
print(summary(model))
print(Anova(model))
print(posthoc_results)
sink()
}


