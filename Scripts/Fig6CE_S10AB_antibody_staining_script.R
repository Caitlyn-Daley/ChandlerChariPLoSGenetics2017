#### Analyses for Scalloped allelic series at 24 degC & High Nutrition using Wing Imaginal Disc stained fo cell proliferation (wg_PH3 stained) and ImageJ ####

##### Entering and Checking data #######
#require(effects)
#require(car)
require(sciplot)

# Entering the Wing Area data- from ImageJ Macro #
sd.id.data <- read.csv("../Data/Fig6C_S10A_sd_wg_PH3_reformatted_by_CHC.csv")
head(sd.id.data)
tail(sd.id.data)
str(sd.id.data)
levels(sd.id.data$Mutation)
sd.id.data$Mutation <- factor(sd.id.data$Mutation, levels(sd.id.data$Mutation)[c(5,1,4,3,2)]) ### reordering the mutation levels in increasing order of severity i.e. wildtype 1st and sd58d last

pdf(file="wing_pouch_plots.pdf", width=8, height=8)
par(mfrow=c(2,2))

cellprolif.ylim <- c(3, 25)
discsize.ylim <- c(2e4, 1.15e5)

#scalloped plots
lineplot.CI(	x.factor=sd.id.data$Mutation, 
				response=sd.id.data$logStdized_Count,
				group=sd.id.data$Background,
				fixed=T,
				legend=F,
				ci.fun= function(x) c(mean(x)-1.97*se(x), mean(x)+1.97*se(x)),
				xlab= "Mutation",
				ylab="Log cell count / wing disc area",
				main="Wing imaginal disc cell proliferation",
				col=c("red","blue"),
				pch=c(16,16), lty=c(1,1),
				cex.main=1.0,
				lwd=2,
				xaxt="n",
				ylim=cellprolif.ylim
			)				
legend (x="topright", xjust=0, legend= c("Oregon-R", "Samarkand"), col=c("red","blue"), pch= c(16,16), lty=c(1,1), bty="n",cex=0.8)
x.axis.labels <- c( expression(paste(WT)),
				  expression(italic(paste(sd^'1',sep=''))),
				  expression(italic(paste(sd^'ETX4',sep=''))),
				  expression(italic(paste(sd^'E3',sep=''))),
				  expression(italic(paste(sd^'58d',sep='')))
					)
axis(side=1, at=1:5, labels=x.axis.labels, tick=F, line=-0.3)


#scalloped
lineplot.CI(	x.factor=sd.id.data$Mutation,
				response=sd.id.data$Wingdisc_Area,
				group=sd.id.data$Background,
				fixed=T,
				legend=F,
				ci.fun= function(x) c(mean(x)-1.97*se(x), mean(x)+1.97*se(x)),
				xlab= "Mutation",
				ylab="Wing disc area",
				main="Wing disc size",
				col=c("red","blue"),
				pch=c(16,16),
				lty=c(1,1), cex.main=1.0,
				lwd=2,
				xaxt="n",
				ylim=discsize.ylim
			)
#legend (x="topright", xjust=0, legend= c("Oregon-R", "Samarkand"), col=c("red","blue"), pch= c(16,16), lty=c(1,1), bty="n",cex=0.8)
x.axis.labels <- c( expression(paste(WT)),
				  expression(italic(paste(sd^'1',sep=''))),
				  expression(italic(paste(sd^'ETX4',sep=''))),
				  expression(italic(paste(sd^'E3',sep=''))),
				  expression(italic(paste(sd^'58d',sep='')))
					)
axis(side=1, at=1:5, labels=x.axis.labels, tick=F, line=-0.3)



##############vestigial allelic series for wing pouches

# Entering the Wing Area data- from ImageJ Macro #
vg.id.data <- read.csv("../Data/Fig6E_S10B_vg_wg_PH3_reformatted_by_CHC.csv")
head(vg.id.data)
tail(vg.id.data)
str(vg.id.data)
levels(vg.id.data$Mutation)
vg.id.data$Mutation <- factor(vg.id.data$Mutation,levels(vg.id.data$Mutation)[c(4,3,2,1)]) ### reordering the mutation levels in increasing order of severity i.e. wildtype 1st and vg1 last

#vestigial
lineplot.CI(	x.factor=vg.id.data$Mutation,
				vg.id.data$logStdized_Count,
				group=vg.id.data$Background,
				fixed=T,
				legend=F, 
				ci.fun= function(x) c(mean(x)-1.97*se(x), mean(x)+1.97*se(x)), 
				xlab= "Mutation", 
				ylab="Log cell count / wing disc area", 
				main="Wing imaginal disc cell proliferation", 
				col=c("red","blue"), 
				pch=c(16,16), 
				lty=c(1,1), 
				cex.main=1.0,
				lwd=2,
				xaxt="n",
				ylim=cellprolif.ylim
			)
			
#legend (x="topright", xjust=0, legend= c("Oregon-R", "Samarkand"), col=c("red","blue"), pch= c(16,16), lty=c(1,1), bty="n",cex=0.8)
x.axis.labels <- c( expression(paste(WT)),
				  expression(italic(paste(vg^'2a33',sep=''))),
				  expression(italic(paste(vg^'21-3',sep=''))),
				  expression(italic(paste(vg^'1',sep='')))
					)
axis(side=1, at=1:4, labels=x.axis.labels, tick=F, line=-0.3)


lineplot.CI(	x.factor=vg.id.data$Mutation,
				vg.id.data$Wingdisc_Area, 
				group=vg.id.data$Background, 
				fixed=T,legend=F, 
				ci.fun= function(x) c(mean(x)-1.97*se(x), mean(x)+1.97*se(x)), 
				xlab= "Mutation", 
				ylab="Wing disc area", 
				main="Wing disc size", col=c("red","blue"), 
				pch=c(16,16), 
				lty=c(1,1), 
				cex.main=1.0,
				lwd=2,
				xaxt="n",
				ylim=discsize.ylim
			)

#legend (x="topright", xjust=0, legend= c("Oregon-R", "Samarkand"), col=c("red","blue"), pch= c(16,16), lty=c(1,1), bty="n",cex=0.8)
x.axis.labels <- c( expression(paste(WT)),
				  expression(italic(paste(vg^'2a33',sep=''))),
				  expression(italic(paste(vg^'21-3',sep=''))),
				  expression(italic(paste(vg^'1',sep='')))
					)
axis(side=1, at=1:4, labels=x.axis.labels, tick=F, line=-0.3)


#############################

dev.off()