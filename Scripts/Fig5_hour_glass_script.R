### Script for plotting the theoretical as well as empirical hour glass model ###

####################################

#Read in the raw complementation data

#Read in the raw data
full.data <- read.csv("../Data/Combined_Final_Data_05_14.csv")

#Exclude any "hybrid" offspring
full.data <- full.data[full.data$F1_Background != "HYBRID",]
full.data$F1_Background <- factor(full.data$F1_Background)



#####################################

# Plotting functions




#Calculates the marginal effect of the allele and the variance in background effects in the simplest way possible:
#   Marginal effect = mean(mutant) - mean(wt)
#   Variance = effect of allele in ORE - effect of allele in SAM
#...while pooling across sex, environment, etc.

simple.effects <- function(my.data) {

	alleles <- unique(c(as.character(my.data$Maternal_Allele), as.character(my.data$Paternal_Allele)))
	alleles <- alleles[alleles != "wt"] #We don't want to include WT in the analysis... WT is the reference allele

	wt.data <- my.data[my.data$Maternal_Allele=="wt" & my.data$Paternal_Allele=="wt",]
	wt.data.ore <- wt.data[wt.data$F1_Background=="ORE",]
	wt.data.sam <- wt.data[wt.data$F1_Background=="SAM",]

	blank.col <- rep(NA, length(alleles)^2)
	results <- data.frame(allele=blank.col, effect=blank.col, variance=blank.col, ore.effect=blank.col, sam.effect=blank.col)

	result.row <- 1
	#Loop across all allelic combinations
	for (i in 1:length(alleles)) {
		cur.allele.1 <- alleles[i]
		for (j in i:length(alleles)) {
			cur.allele.2 <- alleles[j]
			cur.data <- my.data[((my.data$Maternal_Allele == cur.allele.1) & (my.data$Paternal_Allele == cur.allele.2)) | ((my.data$Maternal_Allele == cur.allele.2) & (my.data$Paternal_Allele == cur.allele.1)),]
			
			if (nrow(cur.data) >= 10) {
				if ((sum(cur.data$F1_Background=="ORE") >= 5) & (sum(cur.data$F1_Background=="SAM") >= 5)) {
					marginal.effect <- -1*( mean(cur.data$Wing_Area) - mean(wt.data$Wing_Area) )
					
					cur.data.ore <- cur.data[cur.data$F1_Background=="ORE",]
					cur.data.sam <- cur.data[cur.data$F1_Background=="SAM",]
					
					ore.effect <- -1* ( mean(cur.data.ore$Wing_Area) - mean(wt.data.ore$Wing_Area) )
					sam.effect <- -1* ( mean(cur.data.sam$Wing_Area) - mean(wt.data.sam$Wing_Area) )
					
					if (cur.allele.1 == cur.allele.2) {
						results$allele[result.row] <- cur.allele.1
					} else {
						results$allele[result.row] <- paste(as.character(cur.allele.1), as.character(cur.allele.2), sep=".")
					}
					results$effect[result.row] <- marginal.effect
					results$ore.effect[result.row] <- ore.effect
					results$sam.effect[result.row] <- sam.effect
					results$variance[result.row] <- (ore.effect - sam.effect)^2
					result.row <- result.row + 1
				}
			}
		}
	}
	
	results <- results[1:result.row - 1,]
	results <- results[order(results$effect, decreasing=TRUE),]
	
	return(results)
}

make.plot.range.3 <- function(my.data, jitter=T) {

	plot.data <- simple.effects(my.data)
    plot.data 
	ylim <- range(c((plot.data$ore.effect - plot.data$effect ), (plot.data$sam.effect -plot.data$effect)))

	par(mar=c(4,5,1,2))
	plot(x=NULL, y=NULL, ylim=ylim, xlim=range(plot.data$effect),
	    xlab="Average magnitude of genotypic effect (reduction in wing area)", 
	    ylab="Range in genotypic\neffects across backgrounds", 
	    #xaxt="n"
	    )
	
	#axis(side=1, at=1:nrow(plot.data), labels=plot.data$allele)
	
	het.color <- "#33AA00AA"
	hom.color <- "#882288AA"
	
	bar.colors <- rep(hom.color, nrow(plot.data)) #Use different colors for homozygous...
	bar.colors[grep(pattern=".", fixed=TRUE, x=as.character(plot.data$allele))] <- het.color # ... and heterozygous genotypes [complementation]
	
	jittered.x <- jitter(plot.data$effect, factor=25)
	
	segments(x0=jittered.x, x1=jittered.x, 
	    y0=(plot.data$ore.effect - plot.data$effect ),
	    y1= (plot.data$sam.effect- plot.data$effect), lwd=4, col=bar.colors)
	    
	legend(x=1.15, y=0.6, legend=c("Homozygous genotypes", "Heterozygous genotypes"), col=c(hom.color, het.color), lwd=4, bty="n", cex=0.8)
}


######################################

pdf(file="Figure5_full_hourglass_wing_area.pdf", width=14, height=5)
par(mfrow=c(1,2))

######################################
#Do panel A
make.plot.range.3(full.data)
text(x=0.035, y=0.65, "A", cex=2)

######################################
#Do panel B
#Set up the plot window
allelic.effects <- seq(from=-pi, to=pi, by=0.05)
par(mar=c(4,5,1,2))
plot(x=NULL, y=NULL, xlim=range(allelic.effects), ylim=c(-2, 2), xaxt="n", yaxt="n", xlab="", ylab="")

#Draw the hourglass figure
variation.1 <- cos(allelic.effects)+1 #Curve representing the top half of the inverted hourglass
variation.2 <- -cos(allelic.effects)-1 #Curve representing the bottom half of the inverted hourglass
polygon(c(allelic.effects, rev(allelic.effects)), c(variation.1, rev(variation.2)), col="#000000AA", border=NA) #Can change this to include a border if we want

#Add the x-axis
axis(side=1, at=c(-pi+0.5, pi-0.5), labels=c("Weak", "Severe"), tick="F", line=-0.5)
mtext("Allelic effects", side=1, line=3)
arrow.y <- -2.4
arrows(x0=-pi+1, y0=arrow.y, x1=pi-1, y1=arrow.y, xpd=TRUE, code=2, length=0.1)

#Add the y-axis label
mtext(text="Phenotypic variation due to\ngenetic background", side=2, line=1, adj=NA)

text(x=-pi+0.05, y=1.8, "B", cex=2)


######################################
dev.off()