
#General PSA: This is very much a work in progress - there are likely improper or missed steps here...especially towards the end

# load required packages
library(glmmTMB)
library(performance)
library(MuMIn)
library(tidyverse)

# Part 1 count data 

#simulate zero-inflated count data
set.seed(101)
n <- 500
elevation <- scale(rnorm(n, mean = 1000, sd = 100))
soil_moisture <- scale(runif(n, 0, 1))
treatment <- factor(sample(c("weeded", "unweeded"), n, replace = TRUE))

# Model abundance with environmental factors + treatment
lambda <- exp(1 + 0.6 * elevation - 0.5 * soil_moisture + 
                ifelse(treatment == "weeded", 0.7, 0))

# Simulate counts with excess zeros
zero_prob <- 0.5
abundance <- rpois(n, lambda)
abundance[rbinom(n, 1, zero_prob) == 1] <- 0

rare_plant_df <- data.frame(
  abundance,
  elevation,
  soil_moisture,
  treatment
)

View(rare_plant_df)

#very right skewed! many zeros!
hist(rare_plant_df$abundance)

#question: how does the abundance of our rare plant change with elevation, soil moisture and weeding treatment?

#model 1: poisson
p_model <- glmmTMB(abundance ~ elevation + soil_moisture + treatment,
                     family = poisson(),
                     data = rare_plant_df)
summary(p_model)

#check our model using check_mopdel from performance package

check_model(p_model, check = "pp_check") # predictions don't match observed 
check_model(p_model, check = "overdispersion") # we are overdispersed
check_model(p_model, check = "homogeneity") # heteroscedasticity 
check_model(p_model, check = "qq") # way off


# model 2: ZIP
zip_model <- glmmTMB(abundance ~ elevation + soil_moisture + treatment,
                     ziformula = ~ elevation + soil_moisture + treatment,
                     family = poisson(), 
                     data = rare_plant_df)
summary(zip_model)

check_model(zip_model, check = "pp_check") 
check_model(zip_model, check = "overdispersion") # still overdispersed
check_model(zip_model, check = "homogeneity") # still heteroscedastic
check_model(zip_model, check = "qq") 


# Model 3: ZIP - let's try again using a hurdle
hurdle_model <- glmmTMB(abundance ~ elevation + soil_moisture + treatment,
                        ziformula = ~ elevation + soil_moisture + treatment,
                        family = truncated_poisson(), #truncated poisson
                        data = rare_plant_df)
summary(hurdle_model)

check_model(hurdle_model, check = "pp_check")
check_model(hurdle_model, check = "overdispersion") # still overdispersed
check_model(hurdle_model, check = "homogeneity") # still heteroscedastic
check_model(hurdle_model, check = "qq") 


#hurdle and ZIP do a better job than basic poisson model. they have basically identical AIC. 
AIC(zip_model, hurdle_model, p_model)

#key distinction between hurdle and ZIP are how they handle stuctural vs sampling zeros. 
#in our case we have minimal structural zeros, hence similar performance

#why is the zero inflation model significant in ther hurdle model but not zip?

# we are doing better but still have overdispersion issue... let's try neg binom

#model 4: for consistency sake we can try basic negbinom but we know zeros will likely be an issue
nb_model <- glmmTMB(abundance ~ elevation + soil_moisture + treatment,
                      family = nbinom2, #nbinom2 uses quadratic relationship between mean and variance
                      data = rare_plant_df)
summary(nb_model)

check_model(nb_model, check = "pp_check") # better than basic poisson but still iffy near 0
check_model(nb_model, check = "overdispersion") # major overdispersion
check_model(nb_model, check = "homogeneity") #  heteroscedastic
check_model(nb_model, check = "qq") # bad


#model 5: zi nb
zinb_model <- glmmTMB(abundance ~ elevation + soil_moisture + treatment,
                      ziformula = ~elevation + soil_moisture + treatment,
                      family = nbinom2(),
                      data = rare_plant_df)
summary(zinb_model)

check_model(zinb_model, check = "pp_check")
check_model(zinb_model, check = "overdispersion") 
check_model(zinb_model, check = "homogeneity") 
check_model(zinb_model, check = "qq")


#model 6: hurdle nb
hurdle_nb_model <- glmmTMB(abundance ~ elevation + soil_moisture + treatment,
                      ziformula = ~elevation + soil_moisture + treatment,
                      family = truncated_nbinom2(),
                      data = rare_plant_df)
summary(hurdle_nb_model)

check_model(hurdle_nb_model, check = "pp_check") 
check_model(hurdle_nb_model, check = "overdispersion") 
check_model(hurdle_nb_model, check = "homogeneity") 
check_model(hurdle_nb_model, check = "qq")

#could try tweedie distribution - flexible but likely not better than other options
# tweedie_model <- glmmTMB(abundance ~ elevation + soil_moisture + treatment,
#                            family = tweedie(),
#                            data = rare_plant_df)
# summary(tweedie_model)
# 
# check_model(tweedie_model, check = "homogeneity") 
# check_model(tweedie_model, check = "qq")
# check_model(tweedie_model, check = "pp_check") 

AIC_table<- AIC(zip_model, hurdle_model, hurdle_nb_model, zinb_model)
View(AIC_table)

# which model?

#
#
#
# Part 2 continuous data

set.seed(202)
n <- 500
elevation2 <- scale(rnorm(n, mean = 1000, sd = 100))
canopy_cover <- scale(runif(n, 0, 1))
treatment2 <- factor(sample(c("burned", "unburned"), n, replace = TRUE))

# Linear predictor for mean (log-link for Gamma)
linpred <- 1 + 0.7 * elevation2 - 0.5 * canopy_cover + 
  ifelse(treatment2 == "burned", 0.4, 0)

mu <- exp(linpred)  # ensure positivity
phi <- 1.5          # shape parameter (higher = less variance)

# Simulate from Gamma
flower_cover <- rgamma(n, shape = phi, rate = phi / mu)

# Add zero inflation
flower_cover[rbinom(n, 1, 0.3) == 1] <- 0

# Make dataframe
cover_df <- data.frame(
  flower_cover,
  elevation2,
  canopy_cover,
  treatment2
)

View(cover_df)

#very right skewed! many zeros!
hist(cover_df$flower_cover)


#question: how does flower cover change with elevation, canopy cover and burning treatment?

# model 1: 
norm_model <- glmmTMB(flower_cover ~ elevation2 + canopy_cover + treatment2,
                             family = gaussian(),
                             data = cover_df)
summary(norm_model)

check_model(norm_model, check = "pp_check") 
check_model(norm_model, check = "homogeneity") 
check_model(norm_model, check = "qq")

#as expected things look funky! not appropriate distribution

# model 2 : tweedie - can handle zeros and flexible to skew 
tweedie_model2 <- glmmTMB(flower_cover ~ elevation2 + canopy_cover + treatment2,
                      family = tweedie(),
                      data = cover_df)
summary(tweedie_model2)

check_model(tweedie_model2, check = "pp_check")
check_model(tweedie_model2, check = "homogeneity")  
check_model(tweedie_model2, check = "qq")
 
#still doesn't look right... lets try a hurdle

#check out distribution of non zeros...
cover_df_emit_0s <- cover_df %>%
  filter(flower_cover > 0)
hist(cover_df_emit_0s$flower_cover)

# model 3: hurdle with normal distribution
hurdle_norm_model <- glmmTMB(flower_cover ~ elevation2 + canopy_cover + treatment2,
                             ziformula = ~ elevation2 + canopy_cover + treatment2,
                             family = gaussian(),
                             data = cover_df)
summary(hurdle_norm_model)

check_model(hurdle_norm_model, check = "pp_check") 
check_model(hurdle_norm_model, check = "homogeneity")  
check_model(hurdle_norm_model, check = "qq")

#still missing a lot 

# model 3: gamma might make more sense but how do we deal with zeros? gamma must be positive... hurdle
hurdle_gamma_model <- glmmTMB(flower_cover ~ elevation2 + canopy_cover + treatment2,
                             ziformula = ~ elevation2 + canopy_cover + treatment2,
                             family = ziGamma(link = "log"), # default link is "inverse" - change to "log"
                             data = cover_df)
summary(hurdle_gamma_model)

check_model(hurdle_gamma_model, check = "pp_check") 
check_model(hurdle_gamma_model, check = "homogeneity")  
check_model(hurdle_gamma_model, check = "qq")

AIC_table2 <- AIC(tweedie_model2, hurdle_norm_model, hurdle_gamma_model)
View(AIC_table2)

# we still have issues with our model... ideas? maybe a bigger question about heteroscedacity... nonlinear model?

# other option...data transformations?

cover_df2 <- cover_df %>%
  mutate(scaled_flower_cover = scale(flower_cover + 1)) %>%
  mutate(pos_scaled_flower_cover = scaled_flower_cover + -min(scaled_flower_cover)) %>%
  mutate(log_flower_cover = log(flower_cover+1))

#scaled data and added min value to make positive... maybe not okay??
hist(cover_df2$pos_scaled_flower_cover)
hist(cover_df2$log_flower_cover)

hurdle_gamma_model2 <- glmmTMB(pos_scaled_flower_cover ~ elevation2 + canopy_cover + treatment2,
                               ziformula = ~ elevation2 + canopy_cover + treatment2,
                               family = ziGamma(link = "log"), # default link is "inverse" - change to "log"
                               data = cover_df2)
summary(hurdle_gamma_model2)

check_model(hurdle_gamma_model2, check = "pp_check") 
check_model(hurdle_gamma_model2, check = "homogeneity")  
check_model(hurdle_gamma_model2, check = "qq")

#log transformed should be more normal
hurdle_norm_model2 <- glmmTMB(log_flower_cover ~ elevation2 + canopy_cover + treatment2,
                             ziformula = ~ elevation2 + canopy_cover + treatment2,
                             family = gaussian(),
                             data = cover_df2)
summary(hurdle_norm_model2)

check_model(hurdle_norm_model2, check = "pp_check") 
check_model(hurdle_norm_model2, check = "homogeneity")  
check_model(hurdle_norm_model2, check = "qq")

AIC_table3 <- AIC(tweedie_model2, hurdle_norm_model, hurdle_gamma_model, hurdle_gamma_model2, hurdle_norm_model2)
View(AIC_table3)
