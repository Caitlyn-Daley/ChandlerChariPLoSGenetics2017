#Read in the raw data
full.data <- read.csv("../Data/Combined_Final_Data_05_14.csv")

#Exclude any "hybrid" offspring
full.data <- full.data[full.data$F1_Background != "HYBRID",]
full.data$F1_Background <- factor(full.data$F1_Background)

#Sort out the complemenation data we want
wt.rows.f <- grep(pattern="wt", x=as.character(full.data$Maternal_Allele))
wt.rows.m <- grep(pattern="wt", x=as.character(full.data$Paternal_Allele))
all(wt.rows.f == wt.rows.m) #Make sure there are no mut x wt crosses

sd.rows.f <- grep(pattern="sd", x=as.character(full.data$Maternal_Allele))
sd.rows.m <- grep(pattern="sd", x=as.character(full.data$Paternal_Allele))
all(sd.rows.f == sd.rows.m) #Make sure that there are no sd x wt or sd x vg crosses in the dataset

vg.rows.f <- grep(pattern="vg", x=as.character(full.data$Maternal_Allele))
vg.rows.m <- grep(pattern="vg", x=as.character(full.data$Paternal_Allele))
all(vg.rows.f == vg.rows.m) #Make sure that there are no vg x wt or vg x sd crosses in the dataset

#Get gene-specific datasets and re-level the allele factors to avoid downstream problems
sd.data <- full.data[c(wt.rows.f, sd.rows.f),]
sd.data$Maternal_Allele <- factor(sd.data$Maternal_Allele)
sd.data$Paternal_Allele <- factor(sd.data$Paternal_Allele)
vg.data <- full.data[c(wt.rows.f, vg.rows.f),]
vg.data$Maternal_Allele <- factor(vg.data$Maternal_Allele)
vg.data$Paternal_Allele <- factor(vg.data$Paternal_Allele)




library(sciplot)

ci.function <- function(x) {
	x <- na.omit(x)
	std.err <- sd(x)/sqrt(length(x))
	lower.ci <- mean(x) + std.err*(qt(p=0.05/2, df=length(x)-1))
	upper.ci <- mean(x) + std.err*(qt(p=1-(0.05/2), df=length(x)-1))
	return(c(lower.ci, upper.ci))
}

#This function generates a single panel based on semi-quantitative wing scores
make.single.plot.semiqt <- function(dataset, temp, allele.1, allele.order) {
	#Subset the data
	cur.data <- dataset[((as.character(dataset$Maternal_Allele) == allele.1) | (as.character(dataset$Paternal_Allele) == allele.1)) & (dataset$Temperature == temp),]
	
	#Re-level the factor so that the alleles get plotted in the right order
	cur.data$Paternal_Allele <- factor(cur.data$Paternal_Allele, levels=allele.order)
	cur.data$Maternal_Allele <- factor(cur.data$Maternal_Allele, levels=allele.order)
	
	cur.data$allele.1 <- ifelse(cur.data$Paternal_Allele == allele.1, as.character(cur.data$Paternal_Allele), as.character(cur.data$Maternal_Allele))
	cur.data$allele.1 <- factor(cur.data$allele.1, levels=allele.order)
	cur.data$allele.2 <- ifelse(cur.data$Paternal_Allele == allele.1, as.character(cur.data$Maternal_Allele), as.character(cur.data$Paternal_Allele))
	cur.data$allele.2 <- factor(cur.data$allele.2, levels=allele.order)
	
	par(mar=c(5, 6.5, 2, 1))
	
	lineplot.CI(x.factor=cur.data$allele.2, response=cur.data$Semiqt_Num, group=cur.data$F1_Background, fixed=T, legend=T, xlab="Allele 2", ylab="Semiquantitative wing measure", col=c("red", "blue"), pch=c(19,19), lty=c(1,1), xlim=c(1,length(allele.order)), lwd=2, cex.leg=1, ci.fun=ci.function, leg.lab=c("Oregon-R", "Samarkand"), ylim=c(0,10))

	
	mtext(paste("Allele 1", allele.1, sep=""), side=2, line=4)
	mtext(sprintf("Semiquantitative wing score, %d°C", temp), side=3, line=0.1)
}

#This function generates a single panel based on wing area
make.single.plot.area <- function(dataset, temp, allele.1, allele.order) {
	#Subset the data
	cur.data <- dataset[((as.character(dataset$Maternal_Allele) == allele.1) | (as.character(dataset$Paternal_Allele) == allele.1)) & (dataset$Temperature == temp),]
	
	#Re-level the factor so that the alleles get plotted in the right order
	cur.data$Paternal_Allele <- factor(cur.data$Paternal_Allele, levels=allele.order)
	cur.data$Maternal_Allele <- factor(cur.data$Maternal_Allele, levels=allele.order)
	
	cur.data$allele.1 <- ifelse(cur.data$Paternal_Allele == allele.1, as.character(cur.data$Paternal_Allele), as.character(cur.data$Maternal_Allele))
	cur.data$allele.1 <- factor(cur.data$allele.1, levels=allele.order)
	cur.data$allele.2 <- ifelse(cur.data$Paternal_Allele == allele.1, as.character(cur.data$Maternal_Allele), as.character(cur.data$Paternal_Allele))
	cur.data$allele.2 <- factor(cur.data$allele.2, levels=allele.order)
	
	par(mar=c(5, 6.5, 2, 1))
	
	lineplot.CI(x.factor=cur.data$allele.2, response=cur.data$Wing_Area, group=cur.data$F1_Background, fixed=T, legend=T, xlab="Allele 2", ylab="Wing area", col=c("red", "blue"), pch=c(19,19), lty=c(1,1), xlim=c(1,length(allele.order)), lwd=2, cex.leg=1, ci.fun=ci.function, leg.lab=c("Oregon-R", "Samarkand"), ylim=c(0, 2.3))
	
	mtext(paste("Allele 1 ", allele.1, sep=""), side=2, line=4)
	mtext(sprintf("Wing area, %d°C", temp), side=3, line=0.1)
}

alleles.sd <- c("sd[1]", "sd[ETX4]", "sd[E3]", "sd[G0309]", "sd[G0315b]", "sd[G0239]", "sd[G0483]")

# alleles.vg <- c("vg[2a33]", "vg[1]", "vg[83b27]", "vg[21-3]", "vg[f02736]")
alleles.vg <- c("vg[2a33]", "vg[1]", "vg[21-3]", "vg[f02736]")

temps <- c(18, 24)
panel.width <- 6
panel.height <- 4

	
	#sd 18 and 24, wing area & semiqt
	cur.file <- ("FigureS4.pdf")
	pdf(file=cur.file, height=panel.height*length(alleles.sd), width=panel.width*4)
	par(mfcol=c(length(alleles.sd), 4))
	for (cur.temp in temps) {
		for (i in 1:length(alleles.sd)) {
			cur.allele.1 <- alleles.sd[i]
			try(make.single.plot.area(sd.data, cur.temp, cur.allele.1, alleles.sd))
		}
		for (i in 1:length(alleles.sd)) {
			cur.allele.1 <- alleles.sd[i]
			try(make.single.plot.semiqt(sd.data, cur.temp, cur.allele.1, alleles.sd))
		}	
	}
	dev.off()
	
	#vg 18 and 24, wing area & semiqt
	cur.file <- ("FigureS5.pdf")
	pdf(file=cur.file, height=panel.height*length(alleles.vg), width=panel.width*4)
	par(mfcol=c(length(alleles.vg), 4))
	for (cur.temp in temps) {
		for (i in 1:length(alleles.vg)) {
			cur.allele.1 <- alleles.vg[i]
			try(make.single.plot.area(vg.data, cur.temp, cur.allele.1, alleles.vg))
		}
		for (i in 1:length(alleles.vg)) {
			cur.allele.1 <- alleles.vg[i]
			try(make.single.plot.semiqt(vg.data, cur.temp, cur.allele.1, alleles.vg))
		}	
	}
	dev.off()
	
	
	
	