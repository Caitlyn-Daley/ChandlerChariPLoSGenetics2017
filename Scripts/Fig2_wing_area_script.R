#### Script for plotting the alleleic series for sd and vg mutations in Samarkand and Oregon-R genetic backgrounds at 18 and 24 degC using wing area measured in ImageJ ####

#Initial setup

### Linplots with CIs ####
require(sciplot)

### Plotting Wing Area vs Mutation ####

pdf(file="Figure2_all_allelic_series_area.pdf", height=2*4, width=2.1*5)

par(mfrow=c(2,2))


set.margins <- function() {
	par(mar=c(5,5,2,1), mgp=c(2.5,0.75,0))
}

ci.fun <- function(x) {
	x <- na.omit(x)
	std.err <- sd(x)/sqrt(length(x))
	lower.ci <- mean(x) + std.err*(qt(p=0.05/2, df=length(x)-1))
	upper.ci <- mean(x) + std.err*(qt(p=1-(0.05/2), df=length(x)-1))
	return(c(lower.ci, upper.ci))
}


y.axis.lim <- c(0,2.2)

##############################################################################################################
##############################################################################################################
##############################################################################################################

full.data <- read.csv("../Data/Combined_Final_Data_05_14.csv")

#Extract the allelic series data: in which maternal and paternal alleles are the same, and non-hybrid offspring
as.data <- full.data[as.character(full.data$Maternal_Allele) == as.character(full.data$Paternal_Allele),]
as.data$Maternal_Allele <- factor(as.data$Maternal_Allele)
as.data <- as.data[as.data$F1_Background != "HYBRID",]

##############################################################################################################
##############################################################################################################
##############################################################################################################


#sd at 24C

#Subset the data we want -- 24C sd allelic series
sd.wa.data <- as.data[as.data$Temperature==24,] #Get 24C data
sd.rows <- grep("sd", sd.wa.data$Maternal_Allele) #We only want sd mutant and WT flies
wt.rows <- grep("wt", sd.wa.data$Maternal_Allele)
sd.wa.data <- sd.wa.data[c(sd.rows, wt.rows),]

sd.wa.data$Background <- sd.wa.data$Maternal_Background 

sd.wa.data$Mutation <- factor(sd.wa.data$Maternal_Allele) #Re-level the Mutation factor so it is ordered from weakest to strongest mutation
sd.wa.data$Mutation <- factor(sd.wa.data$Mutation, levels(sd.wa.data$Mutation)[c(6,1,4,3,2,5)])

sd.wa.F.data <- sd.wa.data[sd.wa.data$Sex=="F",] ##subsetting data by sex
sd.wa.M.data <- sd.wa.data[sd.wa.data$Sex=="M",]


sd.wa.F.means <- with (sd.wa.F.data, tapply(Wing_Area, INDEX=list(Background,Mutation), FUN=mean)) 
sd.wa.F.means <- t(sd.wa.F.means) ### Means for the female data

sd.wa.M.means <- with (sd.wa.M.data, tapply(Wing_Area, INDEX=list(Background,Mutation), FUN=mean))
sd.wa.M.means <- t(sd.wa.M.means) ### Means for the male data

sd.semiqt.F.means <- with (sd.wa.F.data, tapply(Semiqt_Num, INDEX=list(Background,Mutation), FUN=mean)) 
sd.semiqt.F.means <- t(sd.semiqt.F.means) ### Means for the semiquantitative female data

sd.semiqt.M.means <- with (sd.wa.M.data, tapply(Semiqt_Num, INDEX=list(Background,Mutation), FUN=mean))
sd.semiqt.M.means <- t(sd.semiqt.M.means) ### Means for the semiquantitative female data


set.margins()
lineplot.CI(sd.wa.data$Mutation, sd.wa.data$Wing_Area, group=sd.wa.data$Background, fixed=T,legend=F, ci.fun= ci.fun, xlab= "Genotype", ylab="Wing Area", main="sd at 24°C", col=c("red","blue"), pch=c(19,19), lty=c(1,1), lwd=2, xaxt="n", ylim=y.axis.lim)
legend (x="topright", xjust=0, legend= c("Oregon-R", "Samarkand"), col=c("red","blue"), pch= c(19,19), lty=c(1,1), bty="n", cex=1.0, lwd=2)

x.axis.labels <- c( expression(paste(WT)),
				  expression(italic(paste(sd^'1',sep=''))),
				  expression(italic(paste(sd^'ETX4',sep=''))),
				  expression(italic(paste(sd^'E3',sep=''))),
				  expression(italic(paste(sd^'58d',sep=''))),
				  expression(italic(paste(sd^'G0309',sep='')))
					)

axis(side=1, at=1:6, labels=x.axis.labels, tick=F, line=-0.3)



##############################################################################################################
##############################################################################################################
##############################################################################################################


#sd at 18C

#### Analyses for Scalloped Allelic Series at 18 degC & High Nutrition, Using and Comparing 6-landmark Wing centroid size from tpsDIG, Wing Area from ImageJ And Semiquantitative Wing size measure ####

#Subset the data we want -- 18C sd allelic series
sd.wa.data <- as.data[as.data$Temperature==18,] #Get 18C data
sd.rows <- grep("sd", sd.wa.data$Maternal_Allele) #We only want sd mutant and WT flies
wt.rows <- grep("wt", sd.wa.data$Maternal_Allele)
sd.wa.data <- sd.wa.data[c(sd.rows, wt.rows),]

sd.wa.data$Background <- sd.wa.data$Maternal_Background 

sd.wa.data$Mutation <- factor(sd.wa.data$Maternal_Allele) #Re-level the Mutation factor so it is ordered from weakest to strongest mutation
sd.wa.data$Mutation <- factor(sd.wa.data$Mutation, levels(sd.wa.data$Mutation)[c(5,1,4,3,2)])

sd.wa.F.data <- sd.wa.data[sd.wa.data$Sex=="F",] ##subsetting data by sex
sd.wa.M.data <- sd.wa.data[sd.wa.data$Sex=="M",]

sd.wa.F.means <- with (sd.wa.F.data, tapply(Wing_Area, INDEX=list(Background,Mutation), FUN=mean)) 
sd.wa.F.means <- t(sd.wa.F.means) ### Means for the female data

sd.wa.M.means <- with (sd.wa.M.data, tapply(Wing_Area, INDEX=list(Background,Mutation), FUN=mean))
sd.wa.M.means <- t(sd.wa.M.means) ### Means for the male data

sd.semiqt.F.means <- with (sd.wa.F.data, tapply(Semiqt_Num, INDEX=list(Background,Mutation), FUN=mean)) 
sd.semiqt.F.means <- t(sd.semiqt.F.means) ### Means for the semiquantitative female data

sd.semiqt.M.means <- with (sd.wa.M.data, tapply(Semiqt_Num, INDEX=list(Background,Mutation), FUN=mean))
sd.semiqt.M.means <- t(sd.semiqt.M.means) ### Means for the semiquantitative female data

set.margins()
lineplot.CI(sd.wa.data$Mutation, sd.wa.data$Wing_Area, group=sd.wa.data$Background, fixed=T,legend=F, ci.fun= ci.fun, xlab= "Genotype", ylab="Wing Area", main="sd at 18°C", col=c("red","blue"), pch=c(19,19), lty=c(1,1), lwd=2, xaxt="n", ylim=y.axis.lim)
#legend (x="topright", xjust=0, legend= c("Oregon-R", "Samarkand"), col=c("red","blue"), pch= c(16,16), lty=c(1,1), bty="n", cex=0.8)


x.axis.labels <- c( expression(paste(WT)),
				  expression(italic(paste(sd^'1',sep=''))),
				  expression(italic(paste(sd^'ETX4',sep=''))),
				  expression(italic(paste(sd^'E3',sep=''))),
				  expression(italic(paste(sd^'58d',sep='')))
					)

axis(side=1, at=1:5, labels=x.axis.labels, tick=F, line=-0.3)



##############################################################################################################
##############################################################################################################
##############################################################################################################

# vg at 24C

#### Analyses for Vestigial Allelic Series at 24 degC & High Nutrition Using and Comparing 6-landmark Wing centroid size from tpvgIG, Wing Area from ImageJ and Semiquantitative Wing Size Measure ####

# Entereing the Wing Area data- from ImageJ Macro #
vg.wa.data <- as.data[as.data$Temperature==24,] #Get 24C data
vg.rows <- grep("vg", vg.wa.data$Maternal_Allele) #We only want vg mutant and WT flies
wt.rows <- grep("wt", vg.wa.data$Maternal_Allele)
vg.wa.data <- vg.wa.data[c(vg.rows, wt.rows),]

vg.wa.data$Background <- vg.wa.data$Maternal_Background 

vg.wa.data$Mutation <- factor(vg.wa.data$Maternal_Allele) #Re-level the Mutation factor so it is ordered from weakest to strongest mutation
vg.wa.data$Mutation <- factor(vg.wa.data$Mutation, levels(vg.wa.data$Mutation)[c(7,3,6,2,1,4,5)])

vg.wa.F.data <- vg.wa.data[vg.wa.data$Sex=="F",] ##subsetting data by sex
vg.wa.M.data <- vg.wa.data[vg.wa.data$Sex=="M",]

vg.wa.F.means <- with (vg.wa.F.data, tapply(Wing_Area, INDEX=list(Background,Mutation), FUN=mean)) 
vg.wa.F.means <- t(vg.wa.F.means) ### Means for the female data

vg.wa.M.means <- with (vg.wa.M.data, tapply(Wing_Area, INDEX=list(Background,Mutation), FUN=mean))
vg.wa.M.means <- t(vg.wa.M.means) ### Means for the male data

vg.semiqt.F.means <- with (vg.wa.F.data, tapply(Semiqt_Num, INDEX=list(Background,Mutation), FUN=mean)) 
vg.semiqt.F.means <- t(vg.semiqt.F.means) ### Means for the semiquantitative female data

vg.semiqt.M.means <- with (vg.wa.M.data, tapply(Semiqt_Num, INDEX=list(Background,Mutation), FUN=mean))
vg.semiqt.M.means <- t(vg.semiqt.M.means) ### Means for the semiquantitative male data


### Plotting Wing Area vs Mutation ####

set.margins()
lineplot.CI(vg.wa.data$Mutation, vg.wa.data$Wing_Area, group=vg.wa.data$Background, fixed=T,legend=F, ci.fun= ci.fun, xlab= "Genotype", ylab="Wing Area", main="vg at 24°C",col=c("red","blue"), pch=c(19,19), lty=c(1,1), lwd=2, xaxt="n", ylim=y.axis.lim)
#legend (x="topright", xjust=0, legend= c("Oregon-R", "Samarkand"), col=c("red","blue"), pch= c(16,16), lty=c(1,1), bty="n", cex=0.8)

x.axis.labels <- c( expression(paste(WT)),
				  expression(italic(paste(vg^'2a33',sep=''))),
				  expression(italic(paste(vg^'RNAi',sep=''))),
				  expression(italic(paste(vg^'21-3',sep=''))),
				  expression(italic(paste(vg^'1',sep=''))),
				  expression(italic(paste(vg^'83b27',sep=''))),
				  expression(italic(paste(vg^'f02736',sep='')))
					)

axis(side=1, at=1:7, labels=x.axis.labels, tick=F, line=-0.3)



##############################################################################################################
##############################################################################################################
##############################################################################################################

# vg at 18C

# Entereing the Wing Area data- from ImageJ Macro #
vg.wa.data <- as.data[as.data$Temperature==18,] #Get 18C data
vg.rows <- grep("vg", vg.wa.data$Maternal_Allele) #We only want vg mutant and WT flies
wt.rows <- grep("wt", vg.wa.data$Maternal_Allele)
vg.wa.data <- vg.wa.data[c(vg.rows, wt.rows),]

vg.wa.data$Background <- vg.wa.data$Maternal_Background 

vg.wa.data$Mutation <- factor(vg.wa.data$Maternal_Allele) #Re-level the Mutation factor so it is ordered from weakest to strongest mutation
vg.wa.data$Mutation <- factor(vg.wa.data$Mutation, levels(vg.wa.data$Mutation)[c(6,5,3,2,1,4)])

vg.wa.F.data <- vg.wa.data[vg.wa.data$Sex=="F",] ##subsetting data by sex
vg.wa.M.data <- vg.wa.data[vg.wa.data$Sex=="M",]

vg.wa.F.means <- with (vg.wa.F.data, tapply(Wing_Area, INDEX=list(Background,Mutation), FUN=mean)) 
vg.wa.F.means <- t(vg.wa.F.means) ### Means for the female data

vg.wa.M.means <- with (vg.wa.M.data, tapply(Wing_Area, INDEX=list(Background,Mutation), FUN=mean))
vg.wa.M.means <- t(vg.wa.M.means) ### Means for the male data

vg.semiqt.F.means <- with (vg.wa.F.data, tapply(Semiqt_Num, INDEX=list(Background,Mutation), FUN=mean)) 
vg.semiqt.F.means <- t(vg.semiqt.F.means) ### Means for the semiquantitative female data

vg.semiqt.M.means <- with (vg.wa.M.data, tapply(Semiqt_Num, INDEX=list(Background,Mutation), FUN=mean))
vg.semiqt.M.means <- t(vg.semiqt.M.means) ### Means for the semiquantitative male data

set.margins()
lineplot.CI(vg.wa.data$Mutation, vg.wa.data$Wing_Area, group=vg.wa.data$Background, fixed=T,legend=F, ci.fun= ci.fun, xlab= "Genotype", ylab="Wing Area", main="vg at 18°C",col=c("red","blue"), pch=c(19,19), lty=c(1,1), lwd=2, xaxt="n", ylim=y.axis.lim)
#legend (x="topright", xjust=0, legend= c("Oregon-R", "Samarkand"), col=c("red","blue"), pch= c(16,16), lty=c(1,1), bty="n", cex=0.8)


x.axis.labels <- c( expression(paste(WT)),
				  expression(italic(paste(vg^'RNAi',sep=''))),
				  expression(italic(paste(vg^'2a33',sep=''))),
				  expression(italic(paste(vg^'21-3',sep=''))),
				  expression(italic(paste(vg^'1',sep=''))),
				  expression(italic(paste(vg^'f02736',sep='')))
					)

axis(side=1, at=1:6, labels=x.axis.labels, tick=F, line=-0.3)

##############################################################################################################
##############################################################################################################
##############################################################################################################

# Finishing up

dev.off()