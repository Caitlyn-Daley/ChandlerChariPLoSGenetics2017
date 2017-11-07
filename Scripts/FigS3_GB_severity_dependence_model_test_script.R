#### Script to analyze and plot experiment in which sd mutant flies were crossed to DGRP lines ####

# libraries
library(ggplot2)
library(sciplot)
library(lsmeans)
# read in and check data
wing_dat <- read.csv("../Data/GS_Wings_RawArea_Oct2015.csv", h=T)
summary(wing_dat)
str(wing_dat)

# make line (DGRP strain) and allele (sd allele) factors, and relevel
wing_dat$Line <- as.factor(wing_dat$Line)
wing_dat$allele <- relevel(wing_dat$allele, ref="SAM")
wing_dat$LogTotalArea <- log(wing_dat$TotalArea)


# Gayatri noted that line 362 was imaged, despite it being potentially contaminated during the cross. So remove these.
wing_dat <- wing_dat[wing_dat$Line != "362",]
wing_dat$Line <- droplevels(wing_dat$Line)

with(wing_dat, table(Line, allele)) # May need to remove line 341 as sdE3 has only two observations. In the end we retained.

# Some quick plots to check data.
hist(wing_dat$TotalArea)
lattice::bwplot(TotalArea ~ Line|allele, data=wing_dat, 
    ylab = " Wing Size (in pixels)", xlab = "Wild type lines")


# generate line means
wing_means <- with(wing_dat,
                     aggregate(TotalArea, by=list(Line, allele), mean))

log_wing_means <- with(wing_dat,
                     aggregate(LogTotalArea, by=list(Line, allele), mean))

# generate line means for semi-quantitative measure of wing perturbation
wing_means_sq <- with(wing_dat,
                     aggregate(Semiqt_Num, by=list(Line, allele), mean))

wing_means$means_sq <- wing_means_sq[,3]
wing_means$log_mean_wing <- log_wing_means[,3]

colnames(wing_means) <- c("line", "allele", "mean_wing_area", "mean_sq", "log_wing_means")


# wing area and semi-quantitative measure are highly correlated, except at close to wild type phenotypes (where the semi-quantitative is sensitive to perturbation, but wing area has a low signal to noise ratio)


with(wing_means, plot(y = log_wing_means, x = mean_sq))
with(wing_means, plot(y = mean_wing_area, x = mean_sq))



ggplot(wing_means, aes(x = mean_sq, y = mean_wing_area, color = allele)) + geom_point()

ggplot(wing_means, aes(x = mean_sq, y = log_wing_means, color = allele)) + geom_point()
  


# order allele types
levels(wing_means$allele)
ordered(wing_means$allele)

wing_means$Mutation <- factor(wing_means$allele, levels = levels(wing_means$allele)[c(1,2,4,3)]) 
#Re-level the Mutation factor so it is ordered from weakest to strongest mutation

# plot for supplementary figure.
par(mfrow = c(3,1), mar=c(4, 5, 1.1, 1.1))

# untransformed data
lineplot.CI(x.factor = Mutation, response = mean_wing_area, 
    group = line, data = wing_means,
    legend = FALSE, ylab = "Mean Wing Area",
    xaxt = "n", xlab = "", main = NULL)

x.axis.labels <- c(expression(italic(paste(sd^'+'/Y, sep = ''))),
				  expression(italic(paste(sd^'1'/Y, sep=''))),
				  expression(italic(paste(sd^'E3'/Y, sep=''))),
				  expression(italic(paste(sd^'58d'/Y, sep='')))
				)

axis(side=1, at=1:4, labels=x.axis.labels, tick=F, line=-0.3)

# log transformed data
lineplot.CI(x.factor = Mutation, response = log_wing_means, 
    group = line, data = wing_means,
    legend = FALSE, ylab = "Mean Wing Area (log)",
    xaxt = "n", xlab = "", main = NULL)

x.axis.labels <- c(expression(italic(paste(sd^'+'/Y, sep = ''))),
				  expression(italic(paste(sd^'1'/Y, sep=''))),
				  expression(italic(paste(sd^'E3'/Y, sep=''))),
				  expression(italic(paste(sd^'58d'/Y, sep='')))
				)

axis(side=1, at=1:4, labels=x.axis.labels, tick=F, line=-0.3)

# Semi-quantitative measure
lineplot.CI(x.factor = Mutation, response = mean_sq, 
    group = line, data = wing_means,
    legend = FALSE, 
    ylab = "Semiquantitative wing measure", ylim = c(1, 10),
    xlab = "Genotype", xaxt = "n")
axis(side=1, at=1:4, labels=x.axis.labels, tick=F, line=-0.3)


levene_med <- function(x) {
	med.x <- median(x)
	lev.stat <- abs(x - med.x)
	mean(lev.stat)
}


log_levene_med <- function(x) {
	med.x <- median(log(x))
	lev.stat <- abs(log(x) - med.x)
	mean(lev.stat)
}

# Allelic means
(group_means <- with(wing_means,
    tapply(mean_wing_area, allele, mean)))

(group_means_log <- with(wing_means,
    tapply(log_wing_means, allele, mean)))

(group_means_sq <- with(wing_means,
    tapply(mean_sq, allele, mean)))   
     
# For standard deviation
group_sd <- with(wing_means,
    tapply(mean_wing_area, allele, sd))
    
with(wing_means,
    tapply(mean_sq, allele, sd))

#look and relationship between mean and standard deviation
plot(group_sd ~ group_means)

# not a particularly linear relationship, so CV is not going to be a good measure.
    
    
# For Levene Stat
with(wing_means,
    tapply(mean_wing_area, allele, levene_med)) 

    
# For Log Levene Stat
(group_lls <- with(wing_means,
    tapply(mean_wing_area, allele, log_levene_med))
    )

# For statistical modeling
LeveneDeviates <- function(y, group, med=TRUE, log_trans=TRUE) {
    
    #log transform data?
    if (log_trans==TRUE)
        y = log(y)
    
    # allows for use of mean or median as measure of central tendency
    if (med==TRUE)
        meds <- tapply(y, group, median, na.rm=TRUE)
    else 
        meds <- tapply(y, group, mean, na.rm=TRUE) 
    
    # calculates deviates for each observation from a measure of central tendency for a "group"
    abs(y - meds[group])}
    
lev_stat <- with(wing_means, 
    LeveneDeviates(y = mean_wing_area, group = allele, med=TRUE, log_trans=TRUE))

ls.lm <- lm(lev_stat ~ allele, data=wing_means)
summary(ls.lm)
lsmeans(ls.lm, "allele")
cld(lsmeans(ls.lm, "allele"))

