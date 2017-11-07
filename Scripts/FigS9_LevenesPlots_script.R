### Script the perform the analysis and plotting of Levene's statistic for the experiment ###

####################################

#Read in the raw complementation data

#Read in the raw data
full.data <- read.csv("../Data/Combined_Final_Data_05_14.csv")

#Exclude any "hybrid" offspring
full.data <- full.data[full.data$F1_Background != "HYBRID",]
full.data$F1_Background <- factor(full.data$F1_Background)



#####################################

cv <- function(cv.data) {
	data.mean <- mean(cv.data, na.rm=T)
	data.sd <- sd(cv.data, na.rm=T)
	cv.result <- data.sd / data.mean
	return(cv.result)
}


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


subset.data <- function(allele.1, allele.2, background=NULL, temp=NULL, sex=NULL) {
	subset.filter <- (full.data$Maternal_Allele == allele.1) & (full.data$Paternal_Allele == allele.2)
	subset.filter <- subset.filter | ((full.data$Maternal_Allele == allele.2) & (full.data$Paternal_Allele == allele.1))
	if (!is.null(background)) {
		subset.filter <- subset.filter & (full.data$F1_Background == background)
	}
	if (!is.null(temp)) {
		subset.filter <- subset.filter & (full.data$Temperature == temp)
	}
	if (!is.null(sex)) {
		subset.filter <- subset.filter & (full.data$Sex == sex)
	}
	
	result <- full.data[subset.filter,]
	
	return(result)
}

get.effects <- function() {
	test.alleles <- c("wt", "sd[1]", "sd[58d]", "sd[E3]", "sd[ETX4]", "sd[G0309]", "vg[1]", "vg[21-3]", "vg[2a33]", "vg[83b27]", "vg[f02736]", "vg[RNAi]")
	blank.col <- rep(NA, length(test.alleles)^2)
	effects.table <- data.frame(allele.1=blank.col, allele.2=blank.col, background=blank.col, area=blank.col, effect=blank.col, cv=blank.col, levene.med=blank.col, log.levene.med=blank.col, n=blank.col, area.wt=blank.col, cv.wt=blank.col, levene.med.wt=blank.col, log.levene.med.wt=blank.col, n.wt=blank.col)
	
	cur.row <- 1
	
	for (cur.allele.1.index in 1:(length(test.alleles)-1)) {
		cur.allele.1 <- test.alleles[cur.allele.1.index]
		for (cur.allele.2.index in (cur.allele.1.index):(length(test.alleles))) {
			cur.allele.2 <- test.alleles[cur.allele.2.index]
			for (cur.background in c("ORE", "SAM")) {
				
				cur.data <- subset.data(allele.1=cur.allele.1, allele.2=cur.allele.2, background=cur.background)
				cur.area <- NA
				cur.effect <- NA
				cur.cv <- NA
				cur.levene.med <- NA
				cur.log.levene.med <- NA
				cur.n <- nrow(cur.data)
				
				cur.data.wt <- subset.data(allele.1="wt", allele.2="wt", background=cur.background)
				cur.area.wt <- NA
				cur.cv.wt <- NA
				cur.levene.med.wt <- NA
				cur.log.levene.med.wt <- NA
				cur.n.wt <- nrow(cur.data.wt)
				
				if ((cur.n >= 5) & (cur.n.wt >= 5)) {
					cur.area <- mean(cur.data$Wing_Area)
					cur.effect <- mean(cur.data$Wing_Area) - mean(cur.data.wt$Wing_Area)
					cur.cv <- cv(cur.data$Wing_Area)
					cur.levene.med <- levene_med(cur.data$Wing_Area)
					cur.log.levene.med <- log_levene_med(cur.data$Wing_Area)
					cur.area.wt <- mean(cur.data.wt$Wing_Area)
					cur.cv.wt <- cv(cur.data.wt$Wing_Area)
					cur.levene.med.wt <- levene_med(cur.data.wt$Wing_Area)
					cur.log.levene.med.wt <- log_levene_med(cur.data.wt$Wing_Area)
				}
				
				effects.table[cur.row,] <- c(cur.allele.1, cur.allele.2, cur.background, cur.area, cur.effect, cur.cv, cur.levene.med, cur.log.levene.med, cur.n, cur.area.wt, cur.cv.wt, cur.levene.med.wt, cur.log.levene.med.wt, cur.n.wt)
				
				cur.row <- cur.row + 1
			}
		}
	}
	
	effects.table <- effects.table[!is.na(effects.table$effect),]
	
	return(effects.table)
}

effects.table <- get.effects()
effects.table$plot.color <- ifelse(effects.table$background=="ORE", "#FF0000", "#0000FF")
effects.table$plot.symbol <- 15 #square -- will be used for sd mutants
effects.table$plot.symbol[(effects.table$allele.1 == "wt") & (effects.table$allele.2 == "wt")] <- 19 #circle -- for WT
effects.table$plot.symbol[grep(pattern="vg", x=effects.table$allele.1)] <- 17 #triangle -- for vg mutants

pdf(file="FigureS9_Levenes_plots.pdf", width=10, height=10)

par(mfrow=c(2,2))

par(mar=c(4,4,1,1))

plot(x=effects.table$area, y=effects.table$levene.med, pch=effects.table$plot.symbol, xlab="Wing area", ylab="Levene's statistic (median)", col=effects.table$plot.color)

panel.id.cex <- 1.5
text("A", x=0.05, y=0.36, cex=panel.id.cex)

legend(x=0.2, y=0.375, legend=c("Oregon-R", "Samarkand"), col=c("red", "blue"), lty=1, pch=16, box.lty=0)

plot(x=effects.table$area, y=effects.table$log.levene.med, pch=effects.table$plot.symbol, xlab="Wing area", ylab="Levene's statistic (log)", col=effects.table$plot.color)


legend.labels <- c("Wild-type",
				  expression(paste(italic(vg), " mutants", sep="")),
				  expression(paste(italic(sd), " mutants", sep=""))
				)



legend(x=1.35, y=0.56, legend=legend.labels, col="black", lty=0, pch=c(19, 17, 15), box.lty=0)
text("B", x=0.05, y=0.535, cex=panel.id.cex)


effects.table$background.effect <- NA
for (i in 1:nrow(effects.table)) {
	cur.allele.1 <- effects.table$allele.1[i]
	cur.allele.2 <- effects.table$allele.2[i]
	cur.background <- effects.table$background[i]
	alternate.row <- which((effects.table$allele.1 == cur.allele.1) & (effects.table$allele.2 == cur.allele.2) & (effects.table$background != cur.background))[1]
	effects.table$background.effect[i] <- abs(as.numeric(effects.table$effect[i]) - as.numeric(effects.table$effect[alternate.row]))
}

plot(x=effects.table$background.effect, y=effects.table$levene.med, pch=effects.table$plot.symbol, xlab="Difference in effects across backgrounds", ylab="Levene's statistic (median)", col=effects.table$plot.color)
text("C", x=0.0015, y=0.36, cex=panel.id.cex)

plot(x=effects.table$background.effect, y=effects.table$log.levene.med, pch=effects.table$plot.symbol, xlab="Difference in effects across backgrounds", ylab="Levene's statistic (log)", col=effects.table$plot.color)
text("D", x=0.0015, y=0.535, cex=panel.id.cex)

dev.off()


cor.test(as.numeric(effects.table[,17]), as.numeric(effects.table[,7]))
cor.test(as.numeric(effects.table[,17]), as.numeric(effects.table[,8]))

mod1 <- lm( as.numeric(effects.table[,8]) ~ as.factor(effects.table[,3]))