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


############################################################################################
############################################################################################
# SET UP THE PDF
pdf(file="Figure3.pdf", width=5, height=12)
par(mfrow=c(3,1))
par(mar=c(4,4,1.5,1))

temperature <- 24

############################################################################################
############################################################################################



############################################################################################
############################################################################################

# PANEL A: sdE3 x sd58d (wing area) [regular non-complementation]

# dataset, allele.1, allele.2, temperature, axis.labels
dataset <- full.data 
allele.1 <- "sd[E3]"
allele.2 <- "sd[58d]"


	#Subset the data & get the non-sd[1] allele
	cur.data <- full.data[full.data$Temperature == temperature,]
	cur.data <- cur.data[(cur.data$Maternal_Allele == allele.1) | (cur.data$Paternal_Allele == allele.1),]
	cur.data$allele.2 <- ifelse(cur.data$Maternal_Allele==allele.1, as.character(cur.data$Paternal_Allele), as.character(cur.data$Maternal_Allele))
	
	order.sd <- c("sd[1]", "sd[ETX4]", "sd[E3]", "sd[58d]", "sd[G0309]", "sd[G0315b]", "sd[G0239]", "sd[G0483]")
	cur.data$allele.2 <- factor(cur.data$allele.2, levels=order.sd)

	lineplot.CI(x.factor=cur.data$allele.2, response=cur.data$Wing_Area, group=cur.data$F1_Background, fixed=T, legend=T, xlab="Second allele", ylab="Wing area", col=c("red", "blue"), pch=c(19,19), lty=c(1,1), cex.leg=0.75, ci.fun=ci.function, xaxt="n", lwd=2, leg.lab=c("Oregon-R", "Samarkand"), ylim=c(0,2.1), y.leg=2.1, x.leg=6, main=expression(paste(italic(sd^'E3'), ' × other ', italic('sd '), 'alleles', sep='')))
	axis.labels <- c( expression(italic(paste(sd^'1',sep=''))),
				  expression(italic(paste(sd^'ETX4',sep=''))),
				  expression(italic(paste(sd^'E3',sep=''))),
				  expression(italic(paste(sd^'58d',sep=''))),
				  expression(italic(paste(sd^'G0309',sep=''))),
				  expression(italic(paste(sd^'G0315b',sep=''))),
				  expression(italic(paste(sd^'G0239',sep=''))),
				  expression(italic(paste(sd^'G0483',sep='')))
				)
	axis(side=1, at=1:8, labels=axis.labels, tick=F, line=-0.3)
	
	#Add horizontal bands for WT and sdG0309 homozygous confidence intervals
	wt.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele=="wt") & (full.data$Paternal_Allele=="wt"),]
	wt.ci <- ci.function(wt.data$Wing_Area)
	rect(xleft=-1, ybottom=wt.ci[1], xright=12, ytop=wt.ci[2], col="#00000031", border="#00000099", lty=2, lwd=0.5)
	
	allele2.ore.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele==allele.2) & (full.data$Paternal_Allele==allele.2) & (full.data$F1_Background=="ORE"),]
	allele2.sam.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele==allele.2) & (full.data$Paternal_Allele==allele.2) & (full.data$F1_Background=="SAM"),]
	allele2.ore.ci <- ci.function(allele2.ore.data$Wing_Area)
	allele2.sam.ci <- ci.function(allele2.sam.data$Wing_Area)
	rect(xleft=-1, ybottom=allele2.sam.ci[1], xright=12, ytop=allele2.sam.ci[2], col="#0000FF31", border="#0000FF99", lty=2, lwd=0.5)
	rect(xleft=-1, ybottom=allele2.ore.ci[1], xright=12, ytop=allele2.ore.ci[2], col="#FF000031", border="#FF000099", lty=2, lwd=0.5)
	
	text("Wild-type, 95% confidence interval", x=3.5, y=mean(wt.ci), pos=4, cex=0.6)
	text(expression(paste(italic(sd^'58d'), " Oregon-R homozygotes", sep="")), x=0.5, y=mean(allele2.ore.ci), pos=4, cex=0.6)
	text(expression(paste(italic(sd^'58d'), " Samarkand homozygotes", sep="")), x=0.5, y=mean(allele2.sam.ci), pos=4, cex=0.6)
	text("A", x=1, y=2.0, cex=1.5)



############################################################################################
############################################################################################

# PANEL B: sdETX4 x sdE3 (wing area) [background-dependent complementation]

# dataset, allele.1, allele.2, temperature, axis.labels
dataset <- full.data 
allele.1 <- "sd[ETX4]"
allele.2 <- "sd[E3]"


	#Subset the data & get the non-sd[1] allele
	cur.data <- full.data[full.data$Temperature == temperature,]
	cur.data <- cur.data[(cur.data$Maternal_Allele == allele.1) | (cur.data$Paternal_Allele == allele.1),]
	cur.data$allele.2 <- ifelse(cur.data$Maternal_Allele==allele.1, as.character(cur.data$Paternal_Allele), as.character(cur.data$Maternal_Allele))
	
	order.sd <- c("sd[1]", "sd[ETX4]", "sd[E3]", "sd[58d]", "sd[G0309]", "sd[G0315b]", "sd[G0239]", "sd[G0483]")
	cur.data$allele.2 <- factor(cur.data$allele.2, levels=order.sd)

	lineplot.CI(x.factor=cur.data$allele.2, response=cur.data$Wing_Area, group=cur.data$F1_Background, fixed=T, legend=T, xlab="", ylab="Wing area", col=c("red", "blue"), pch=c(19,19), lty=c(1,1), cex.leg=0.75, ci.fun=ci.function, xaxt="n", lwd=2, leg.lab=c("Oregon-R", "Samarkand"), ylim=c(0,2.1), y.leg=2.1, x.leg=6, main=expression(paste(italic(sd^'ETX4'), ' × other ', italic('sd '), 'alleles', sep='')))
	axis.labels <- c( expression(italic(paste(sd^'1',sep=''))),
				  expression(italic(paste(sd^'ETX4',sep=''))),
				  expression(italic(paste(sd^'E3',sep=''))),
				  expression(italic(paste(sd^'58d',sep=''))),
				  expression(italic(paste(sd^'G0309',sep=''))),
				  expression(italic(paste(sd^'G0315b',sep=''))),
				  expression(italic(paste(sd^'G0239',sep=''))),
				  expression(italic(paste(sd^'G0483',sep='')))
				)
	axis(side=1, at=1:8, labels=axis.labels, tick=F, line=-0.3)
	
	#Add horizontal bands for WT and sdG0309 homozygous confidence intervals
	wt.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele=="wt") & (full.data$Paternal_Allele=="wt"),]
	wt.ci <- ci.function(wt.data$Wing_Area)
	rect(xleft=-1, ybottom=wt.ci[1], xright=12, ytop=wt.ci[2], col="#00000031", border="#00000099", lty=2, lwd=0.5)
	
	allele2.ore.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele==allele.2) & (full.data$Paternal_Allele==allele.2) & (full.data$F1_Background=="ORE"),]
	allele2.sam.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele==allele.2) & (full.data$Paternal_Allele==allele.2) & (full.data$F1_Background=="SAM"),]
	allele2.ore.ci <- ci.function(allele2.ore.data$Wing_Area)
	allele2.sam.ci <- ci.function(allele2.sam.data$Wing_Area)
	rect(xleft=-1, ybottom=allele2.sam.ci[1], xright=12, ytop=allele2.sam.ci[2], col="#0000FF31", border="#0000FF99", lty=2, lwd=0.5)
	rect(xleft=-1, ybottom=allele2.ore.ci[1], xright=12, ytop=allele2.ore.ci[2], col="#FF000031", border="#FF000099", lty=2, lwd=0.5)
	
	text("Wild-type, 95% confidence interval", x=3.5, y=mean(wt.ci), pos=4, cex=0.6)
	text(expression(paste(italic(sd^'E3'), " Oregon-R homozygotes", sep="")), x=0.5, y=mean(allele2.ore.ci), pos=4, cex=0.6)
	text(expression(paste(italic(sd^'E3'), " Samarkand homozygotes", sep="")), x=0.5, y=mean(allele2.sam.ci), pos=4, cex=0.6)
	text("B", x=1, y=2.0, cex=1.5)





############################################################################################
############################################################################################

# PANEL C: sd1 x sdG0309 (wing area) [background dependence in heterozygote, not in corresponding homozygotes]

# dataset, allele.1, allele.2, temperature, axis.labels
dataset <- full.data 
allele.1 <- "sd[1]"
allele.2 <- "sd[G0309]"


	#Subset the data & get the non-sd[1] allele
	cur.data <- full.data[full.data$Temperature == temperature,]
	cur.data <- cur.data[(cur.data$Maternal_Allele == allele.1) | (cur.data$Paternal_Allele == allele.1),]
	cur.data$allele.2 <- ifelse(cur.data$Maternal_Allele==allele.1, as.character(cur.data$Paternal_Allele), as.character(cur.data$Maternal_Allele))
	
	order.sd <- c("sd[1]", "sd[ETX4]", "sd[E3]", "sd[58d]", "sd[G0309]", "sd[G0315b]", "sd[G0239]", "sd[G0483]")
	cur.data$allele.2 <- factor(cur.data$allele.2, levels=order.sd)

	lineplot.CI(x.factor=cur.data$allele.2, response=cur.data$Wing_Area, group=cur.data$F1_Background, fixed=T, legend=T, xlab="", ylab="Wing area", col=c("red", "blue"), pch=c(19,19), lty=c(1,1), cex.leg=0.75, ci.fun=ci.function, xaxt="n", lwd=2, leg.lab=c("Oregon-R", "Samarkand"), ylim=c(0,2.1), y.leg=2.1, x.leg=6, main=expression(paste(italic(sd^'1'), ' × other ', italic('sd '), 'alleles', sep='')))
	axis.labels <- c( expression(italic(paste(sd^'1',sep=''))),
				  expression(italic(paste(sd^'ETX4',sep=''))),
				  expression(italic(paste(sd^'E3',sep=''))),
				  expression(italic(paste(sd^'58d',sep=''))),
				  expression(italic(paste(sd^'G0309',sep=''))),
				  expression(italic(paste(sd^'G0315b',sep=''))),
				  expression(italic(paste(sd^'G0239',sep=''))),
				  expression(italic(paste(sd^'G0483',sep='')))
				)
	axis(side=1, at=1:8, labels=axis.labels, tick=F, line=-0.3)
	
	#Add horizontal bands for WT and sdG0309 homozygous confidence intervals
	wt.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele=="wt") & (full.data$Paternal_Allele=="wt"),]
	wt.ci <- ci.function(wt.data$Wing_Area)
	rect(xleft=-1, ybottom=wt.ci[1], xright=12, ytop=wt.ci[2], col="#00000031", border="#00000099", lty=2, lwd=0.5)
	
	allele2.ore.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele==allele.2) & (full.data$Paternal_Allele==allele.2) & (full.data$F1_Background=="ORE"),]
	allele2.sam.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele==allele.2) & (full.data$Paternal_Allele==allele.2) & (full.data$F1_Background=="SAM"),]
	allele2.ore.ci <- ci.function(allele2.ore.data$Wing_Area)
	allele2.sam.ci <- ci.function(allele2.sam.data$Wing_Area)
	rect(xleft=-1, ybottom=allele2.sam.ci[1], xright=12, ytop=allele2.sam.ci[2], col="#0000FF31", border="#0000FF99", lty=2, lwd=0.5)
	rect(xleft=-1, ybottom=allele2.ore.ci[1], xright=12, ytop=allele2.ore.ci[2], col="#FF000031", border="#FF000099", lty=2, lwd=0.5)
	
	text("Wild-type, 95% confidence interval", x=4, y=mean(wt.ci), pos=4, cex=0.6)
	text(expression(paste(italic(sd^'G0309'), " Oregon-R homozygotes", sep="")), x=0.5, y=mean(allele2.ore.ci)+0.05, pos=4, cex=0.6)
	text(expression(paste(italic(sd^'G0309'), " Samarkand homozygotes", sep="")), x=0.5, y=mean(allele2.sam.ci-0.05), pos=4, cex=0.6)
	text("C", x=1, y=2.0, cex=1.5)






############################################################################################
############################################################################################
# CLOSE THE PDF 
dev.off()
############################################################################################
############################################################################################









############################################################################################
############################################################################################
# SET UP THE PDF
pdf(file="Figure4.pdf", width=5, height=12)
par(mfrow=c(3,1))
par(mar=c(4,4,1.5,1))

############################################################################################
############################################################################################


############################################################################################
############################################################################################
# PANEL A: vg83b27 x vgf02736 (wing area)

# dataset, allele.1, allele.2, temperature, axis.labels
dataset <- full.data 
allele.1 <- "vg[83b27]"
allele.2 <- "vg[f02736]"


	#Subset the data & get the non-sd[1] allele
	cur.data <- full.data[full.data$Temperature == temperature,]
	cur.data <- cur.data[(cur.data$Maternal_Allele == allele.1) | (cur.data$Paternal_Allele == allele.1),]
	cur.data$allele.2 <- ifelse(cur.data$Maternal_Allele==allele.1, as.character(cur.data$Paternal_Allele), as.character(cur.data$Maternal_Allele))
	
	order.vg <- c("vg[2a33]", "vg[1]", "vg[83b27]", "vg[21-3]", "vg[f02736]")
	cur.data$allele.2 <- factor(cur.data$allele.2, levels=order.vg)

	lineplot.CI(x.factor=cur.data$allele.2, response=cur.data$Wing_Area, group=cur.data$F1_Background, fixed=T, legend=T, xlab="", ylab="Wing area", col=c("red", "blue"), pch=c(19,19), lty=c(1,1), cex.leg=0.75, ci.fun=ci.function, xaxt="n", lwd=2, leg.lab=c("Oregon-R", "Samarkand"), ylim=c(0,2.1), y.leg=2.1, x.leg=4.5, main=expression(paste(italic(vg^'83b27'), ' × other ', italic('vg '), 'alleles', sep='')))
	axis.labels <- c( expression(italic(paste(vg^'2a33',sep=''))),
				  expression(italic(paste(vg^'1',sep=''))),
				  expression(italic(paste(vg^'83b27',sep=''))),
				  expression(italic(paste(vg^'21-3',sep=''))),
				  expression(italic(paste(vg^'f02736',sep='')))
				)
	axis(side=1, at=1:5, labels=axis.labels, tick=F, line=-0.3)
	
	#Add horizontal bands for WT and sdG0309 homozygous confidence intervals
	wt.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele=="wt") & (full.data$Paternal_Allele=="wt"),]
	wt.ci <- ci.function(wt.data$Wing_Area)
	rect(xleft=-1, ybottom=wt.ci[1], xright=12, ytop=wt.ci[2], col="#00000031", border="#00000099", lty=2, lwd=0.5)
	
	allele2.ore.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele==allele.2) & (full.data$Paternal_Allele==allele.2) & (full.data$F1_Background=="ORE"),]
	allele2.sam.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele==allele.2) & (full.data$Paternal_Allele==allele.2) & (full.data$F1_Background=="SAM"),]
	allele2.ore.ci <- ci.function(allele2.ore.data$Wing_Area)
	allele2.sam.ci <- ci.function(allele2.sam.data$Wing_Area)
	rect(xleft=-1, ybottom=allele2.sam.ci[1], xright=12, ytop=allele2.sam.ci[2], col="#0000FF31", border="#0000FF99", lty=2, lwd=0.5)
	rect(xleft=-1, ybottom=allele2.ore.ci[1], xright=12, ytop=allele2.ore.ci[2], col="#FF000031", border="#FF000099", lty=2, lwd=0.5)
	
	text("Wild-type, 95% confidence interval", x=2.5, y=mean(wt.ci), pos=4, cex=0.6)
	text(expression(paste(italic(vg^'f02736'), " Oregon-R homozygotes", sep="")), x=0.8, y=mean(allele2.ore.ci)+0.05, pos=4, cex=0.6)
	text(expression(paste(italic(vg^'f02736'), " Samarkand homozygotes", sep="")), x=0.8, y=mean(allele2.sam.ci)-0.05, pos=4, cex=0.6)
	text("A", x=1, y=2.0, cex=1.5)




############################################################################################
############################################################################################


# PANEL B: vg2a33 x vg21-3 (wing area)

# dataset, allele.1, allele.2, temperature, axis.labels
dataset <- full.data 
allele.1 <- "vg[2a33]"
allele.2 <- "vg[21-3]"


	#Subset the data & get the non-sd[1] allele
	cur.data <- full.data[full.data$Temperature == temperature,]
	cur.data <- cur.data[(cur.data$Maternal_Allele == allele.1) | (cur.data$Paternal_Allele == allele.1),]
	cur.data$allele.2 <- ifelse(cur.data$Maternal_Allele==allele.1, as.character(cur.data$Paternal_Allele), as.character(cur.data$Maternal_Allele))
	
	order.vg <- c("vg[2a33]", "vg[1]", "vg[83b27]", "vg[21-3]", "vg[f02736]")
	cur.data$allele.2 <- factor(cur.data$allele.2, levels=order.vg)

	lineplot.CI(x.factor=cur.data$allele.2, response=cur.data$Wing_Area, group=cur.data$F1_Background, fixed=T, legend=T, xlab="", ylab="Wing area", col=c("red", "blue"), pch=c(19,19), lty=c(1,1), cex.leg=0.75, ci.fun=ci.function, xaxt="n", lwd=2, leg.lab=c("Oregon-R", "Samarkand"), ylim=c(0,2.1), y.leg=2.1, x.leg=4.5, main=expression(paste(italic(vg^'2a33'), ' × other ', italic('vg '), 'alleles', sep='')))
	axis.labels <- c( expression(italic(paste(vg^'2a33',sep=''))),
				  expression(italic(paste(vg^'1',sep=''))),
				  expression(italic(paste(vg^'83b27',sep=''))),
				  expression(italic(paste(vg^'21-3',sep=''))),
				  expression(italic(paste(vg^'f02736',sep='')))
				)
	axis(side=1, at=1:5, labels=axis.labels, tick=F, line=-0.3)
	
	#Add horizontal bands for WT and sdG0309 homozygous confidence intervals
	wt.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele=="wt") & (full.data$Paternal_Allele=="wt"),]
	wt.ci <- ci.function(wt.data$Wing_Area)
	rect(xleft=-1, ybottom=wt.ci[1], xright=12, ytop=wt.ci[2], col="#00000031", border="#00000099", lty=2, lwd=0.5)
	
	allele2.ore.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele==allele.2) & (full.data$Paternal_Allele==allele.2) & (full.data$F1_Background=="ORE"),]
	allele2.sam.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele==allele.2) & (full.data$Paternal_Allele==allele.2) & (full.data$F1_Background=="SAM"),]
	allele2.ore.ci <- ci.function(allele2.ore.data$Wing_Area)
	allele2.sam.ci <- ci.function(allele2.sam.data$Wing_Area)
	rect(xleft=-1, ybottom=allele2.sam.ci[1], xright=12, ytop=allele2.sam.ci[2], col="#0000FF31", border="#0000FF99", lty=2, lwd=0.5)
	rect(xleft=-1, ybottom=allele2.ore.ci[1], xright=12, ytop=allele2.ore.ci[2], col="#FF000031", border="#FF000099", lty=2, lwd=0.5)
	
	text("Wild-type, 95% confidence interval", x=3.5, y=mean(wt.ci), pos=4, cex=0.6)
	text(expression(paste(italic(vg^'21-3'), " Oregon-R homozygotes", sep="")), x=0.8, y=mean(allele2.ore.ci)-0.05, pos=4, cex=0.6)
	text(expression(paste(italic(vg^'21-3'), " Samarkand homozygotes", sep="")), x=0.8, y=mean(allele2.sam.ci)+0.05, pos=4, cex=0.6)
	text("B", x=1, y=2.0, cex=1.5)



############################################################################################
############################################################################################


# PANEL C: vg1 x vg2a33 (wing area)

# dataset, allele.1, allele.2, temperature, axis.labels
dataset <- full.data 
allele.1 <- "vg[1]"
allele.2 <- "vg[21-3]"


	#Subset the data & get the non-sd[1] allele
	cur.data <- full.data[full.data$Temperature == temperature,]
	cur.data <- cur.data[(cur.data$Maternal_Allele == allele.1) | (cur.data$Paternal_Allele == allele.1),]
	cur.data$allele.2 <- ifelse(cur.data$Maternal_Allele==allele.1, as.character(cur.data$Paternal_Allele), as.character(cur.data$Maternal_Allele))
	
	order.vg <- c("vg[2a33]", "vg[1]", "vg[83b27]", "vg[21-3]", "vg[f02736]")
	cur.data$allele.2 <- factor(cur.data$allele.2, levels=order.vg)

	lineplot.CI(x.factor=cur.data$allele.2, response=cur.data$Wing_Area, group=cur.data$F1_Background, fixed=T, legend=T, xlab="Second allele", ylab="Wing area", col=c("red", "blue"), pch=c(19,19), lty=c(1,1), cex.leg=0.75, ci.fun=ci.function, xaxt="n", lwd=2, leg.lab=c("Oregon-R", "Samarkand"), ylim=c(0,2.1), y.leg=2.1, x.leg=4.5, main=expression(paste(italic(vg^'1'), ' × other ', italic('vg '), 'alleles', sep='')))
	axis.labels <- c( expression(italic(paste(vg^'2a33',sep=''))),
				  expression(italic(paste(vg^'1',sep=''))),
				  expression(italic(paste(vg^'83b27',sep=''))),
				  expression(italic(paste(vg^'21-3',sep=''))),
				  expression(italic(paste(vg^'f02736',sep='')))
				)
	axis(side=1, at=1:5, labels=axis.labels, tick=F, line=-0.3)
	
	#Add horizontal bands for WT and sdG0309 homozygous confidence intervals
	wt.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele=="wt") & (full.data$Paternal_Allele=="wt"),]
	wt.ci <- ci.function(wt.data$Wing_Area)
	rect(xleft=-1, ybottom=wt.ci[1], xright=12, ytop=wt.ci[2], col="#00000031", border="#00000099", lty=2, lwd=0.5)
	
	allele2.ore.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele==allele.2) & (full.data$Paternal_Allele==allele.2) & (full.data$F1_Background=="ORE"),]
	allele2.sam.data <- full.data[(full.data$Temperature == temperature) & (full.data$Maternal_Allele==allele.2) & (full.data$Paternal_Allele==allele.2) & (full.data$F1_Background=="SAM"),]
	allele2.ore.ci <- ci.function(allele2.ore.data$Wing_Area)
	allele2.sam.ci <- ci.function(allele2.sam.data$Wing_Area)
	rect(xleft=-1, ybottom=allele2.sam.ci[1], xright=12, ytop=allele2.sam.ci[2], col="#0000FF31", border="#0000FF99", lty=2, lwd=0.5)
	rect(xleft=-1, ybottom=allele2.ore.ci[1], xright=12, ytop=allele2.ore.ci[2], col="#FF000031", border="#FF000099", lty=2, lwd=0.5)
	
	text("Wild-type, 95% confidence interval", x=3.5, y=mean(wt.ci), pos=4, cex=0.6)
	text(expression(paste(italic(vg^'21-3'), " Oregon-R homozygotes", sep="")), x=2.3, y=mean(allele2.ore.ci)-0.05, pos=4, cex=0.6)
	text(expression(paste(italic(vg^'21-3'), " Samarkand homozygotes", sep="")), x=2.3, y=mean(allele2.sam.ci)+0.05, pos=4, cex=0.6)
	text("C", x=1, y=2.0, cex=1.5)









############################################################################################
############################################################################################
# CLOSE THE PDF 
dev.off()
############################################################################################
############################################################################################