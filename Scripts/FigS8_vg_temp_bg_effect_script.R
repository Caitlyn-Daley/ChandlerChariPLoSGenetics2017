### Script to plot the effects of Genetic Background X Temperature on vg[1] mutation ###

#libraries
library(ggplot2)
library(car)
library(effects)

# location of file
setwd("../data")

# read in file
vg_wing_dat <- read.csv("../Data/JH_2017_vg_temp_background_wing_length_long.csv", h = T)

# quick check of data structure.
str(vg_wing_dat)


# making wild type the default level.
vg_wing_dat$Genotype <- relevel(vg_wing_dat$Genotype, "wt")

# Quick boxplot of data
ggplot(vg_wing_dat, aes(y = length/1000, x = as.factor(Temp), col = Genotype:Background:Sex)) +  
    geom_boxplot(outlier.shape = NA) + 
    labs(y = "Wing Length (mm)", x = "Rearing Temperature", size = 35) + 
    theme(axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20))



# For this study we want to treat temperature as a factor as we are interested in temperature specific effects on mutations, not a general trend from low to high temperature.

vg_wing_dat$TempFactor <- as.factor(vg_wing_dat$Temp)

# fit linear model
mod2 <- lm(length ~ TempFactor*Genotype*Sex*Background, data = vg_wing_dat)
Anova(mod2)
summary(mod2)
plot(allEffects(mod2))

mod2_effects <- effect("TempFactor:Genotype:Sex:Background", mod2)
dat_mod2_effects <- as.data.frame(mod2_effects)


p <- ggplot(dat_mod2_effects, 
    aes(y = fit/1000, x = TempFactor, color = Sex:Background, shape = Genotype)) + 
    geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000), size = 1.3, width = 0.8,
       position = position_dodge(width = 0.5) ) +
    geom_point(size = 5, position = position_dodge(width = 0.5)) + 
    labs(y = "Wing Length (mm)", x = "Rearing Temperature", col="Sex:Background", size = 35) + 
    theme(axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20)) 
p 