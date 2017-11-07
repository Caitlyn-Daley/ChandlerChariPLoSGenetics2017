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
	result <- c(mean(x) - 2*se(x), mean(x) + 2*se(x))
}

order.sd.paternal <- c("wt", "sd[1]", "sd[ETX4]", "sd[E3]", "sd[58d]", "sd[G0309]")
order.sd.maternal <- c("wt", "sd[1]", "sd[ETX4]", "sd[E3]", "sd[58d]", "sd[G0309]", "sd[G0315b]", "sd[G0239]", "sd[G0483]")
mutations.sd <- c("WT", "1", "ETX4", "E3", "G0309", "G0315b", "G0239", "G0483")

order.vg <- c("wt", "vg[2a33]", "vg[RNAi]", "vg[1]", "vg[83b27]", "vg[21-3]", "vg[f02736]")
mutations.vg <- c("WT", "2a33", "RNAi", "1", "83b27", "21-3", "f02736")

all.alleles <- c(order.sd.maternal, order.vg[2:length(order.vg)])
all.mutations <- c(mutations.sd, mutations.vg[2:length(mutations.vg)])

make.single.plot <- function(cur.data, cur.temp, allele.1, allele.2, response.trait="Semiqt_Num") {
	#Subset the data to just the two alleles of interest, and throw away the homozygous flies
	keepers.1 <- cur.data$Maternal_Allele == allele.1 & cur.data$Paternal_Allele == allele.2
	keepers.2 <- cur.data$Maternal_Allele == allele.2 & cur.data$Paternal_Allele == allele.1
	cur.data <- cur.data[keepers.1 | keepers.2,]
				
	cur.data$cross.direction <- paste(as.character(cur.data$Maternal_Allele), " x\n", as.character(cur.data$Paternal_Allele), sep="")	
	
	if (length(grep(pattern="sd", allele.1)) >= 1) {
		cur.gene <- "sd"
	} else {
		cur.gene <- "vg"
	}
	
	par(mar=c(4,4,1,1))
	
	x.lab <- paste("Maternal allele x Paternal allele, ", cur.temp, "C", sep="")
	
	if (response.trait == "Semiqt_Num") {
		y.axis.lim <- c(1.5, 10)
		y.lab <- "Semiquantitive wing score"
	} else {
		y.axis.lim <- c(0.1, 2.3)
		y.lab <- "Wing area"
	}
	
	lineplot.CI(x.factor=cur.data$cross.direction, response=cur.data[,response.trait], group=cur.data$F1_Background, xlab=x.lab, ylab=y.lab, col=c("red", "blue"), pch=c(16,16), lty=c(1,1), main="", fixed=T, leg.lab=c("Oregon-R", "Samarkand"), lwd=2, x.leg=1.6, y.leg=y.axis.lim[2], ylim=y.axis.lim)

	#axis(side=1, at=1:2, labels=c(axis.label.1.e, axis.label.2.e), tick=F)
}

temps <- c(18, 24)

#Do SD
order.sd <- order.sd.paternal #There are fewer paternal sd alleles, so if we want reciprocal combinations, we have to stick to those present in the father
pdf("FigureS6.pdf", height=6*4, width=2*4)
par(mfrow=c(6,2))
for (cur.temp in temps) {
	#Loop through allele combinations and pick out which ones we do have good data for
	n.alleles <- length(order.sd)
	for (i in 1:(n.alleles-1)) {
		allele.1 <- order.sd[i]
		for (j in (i+1):n.alleles) {
			allele.2 <- order.sd[j]

			#We only want to keep the current combination if we have reciprocal crosses for both SAM and ORE
			cur.data <- sd.data[sd.data$Temperature==cur.temp,]
			sam.data <- cur.data[cur.data$F1_Background=="SAM",]
			ore.data <- cur.data[cur.data$F1_Background=="ORE",]
			reciprocal.sam <- any(sam.data$Maternal_Allele == allele.1 & sam.data$Paternal_Allele == allele.2) & any(sam.data$Maternal_Allele == allele.2 & sam.data$Paternal_Allele == allele.1)
			reciprocal.ore <- any(ore.data$Maternal_Allele == allele.1 & ore.data$Paternal_Allele == allele.2) & any(ore.data$Maternal_Allele == allele.2 & ore.data$Paternal_Allele == allele.1)
			if (reciprocal.sam & reciprocal.ore) {
				#We have reciprocal cross data for both SAM & ORE
				#Now make the plot
				#	y-axis: area/score
				#	x-axis: two reciprocal crosses
				#	multiple lines for SAM & ORE
				
				make.single.plot(cur.data, cur.temp, allele.1, allele.2, "Semiqt_Num")
				make.single.plot(cur.data, cur.temp, allele.1, allele.2, "Wing_Area")
				
				print(c(allele.1, allele.2))
			}
		}
	}
	
}
dev.off()



#Do VG
pdf("FigureS7.pdf", height=5*4, width=2*4)
par(mfrow=c(7,2))
for (cur.temp in temps) {
	#Loop through allele combinations and pick out which ones we do have good data for
	n.alleles <- length(order.vg)
	for (i in 1:(n.alleles-1)) {
		allele.1 <- order.vg[i]
		for (j in (i+1):n.alleles) {
			allele.2 <- order.vg[j]

			#We only want to keep the current combination if we have reciprocal crosses for both SAM and ORE
			cur.data <- vg.data[vg.data$Temperature==cur.temp,]
			sam.data <- cur.data[cur.data$F1_Background=="SAM",]
			ore.data <- cur.data[cur.data$F1_Background=="ORE",]
			reciprocal.sam <- any(sam.data$Maternal_Allele == allele.1 & sam.data$Paternal_Allele == allele.2) & any(sam.data$Maternal_Allele == allele.2 & sam.data$Paternal_Allele == allele.1)
			reciprocal.ore <- any(ore.data$Maternal_Allele == allele.1 & ore.data$Paternal_Allele == allele.2) & any(ore.data$Maternal_Allele == allele.2 & ore.data$Paternal_Allele == allele.1)
			if (reciprocal.sam & reciprocal.ore) {
				#We have reciprocal cross data for both SAM & ORE
				#Now make the plot
				#	y-axis: area/score
				#	x-axis: two reciprocal crosses
				#	multiple lines for SAM & ORE
				
				make.single.plot(cur.data, cur.temp, allele.1, allele.2, "Semiqt_Num")
				make.single.plot(cur.data, cur.temp, allele.1, allele.2, "Wing_Area")
				
				print(c(allele.1, allele.2))
			}
		}
	}
	
}
dev.off()
