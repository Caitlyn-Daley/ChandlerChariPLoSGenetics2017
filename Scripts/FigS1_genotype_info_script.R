### Script for plotting the information from genotyping data for all the mutations in Samarkand and Oregon-R genetic backgrounds ###

wide.data <- read.csv(file="../Data/BackgroundGenotypes.csv", skip=7, na.strings=c("?"), colClasses=rep("character", 102))

#Reorganize and clean up the data table for easier plottings
tall.data <- as.data.frame(t(wide.data[,2:ncol(wide.data)]))
colnames(tall.data) <- as.character(wide.data[,1])
for (c in 1:ncol(tall.data)) {
	tall.data[,c] <- as.character(tall.data[,c])
}

rownames(tall.data) <- sub(pattern="X", replacement="", x=rownames(tall.data))
raw.chroms <- substr(rownames(tall.data), start=1, stop=2)
chroms <- ifelse(raw.chroms == "10", "X",
		  ifelse(raw.chroms == "21", "2L",
		  ifelse(raw.chroms == "22", "2R",
		  ifelse(raw.chroms == "31", "3L",
		  ifelse(raw.chroms == "32", "3R", NA)))))
snp.positions <- as.numeric(substr(rownames(tall.data), start=3, stop=nchar(rownames(tall.data))))

tall.data$chrom <- chroms
tall.data$pos <- snp.positions

tall.data <- tall.data[order(tall.data$pos),]
tall.data <- tall.data[order(tall.data$chrom),]

#Make the plot -- color code by base
allele.colors <- c("#FF0066", "#6600FF", "#00FF66", "#CC9911") #A C G T plotting colors
alleles <- c("A", "C", "G", "T")
chrom.sizes <- c(20e6, 20e6, 21e6, 22e6, 27e6) #X, 2L, 2R, 3L, 3R total sizes
white.space <- 3e6

chrom.starts <- c(2e7)
chrom.ends <- c(chrom.starts[1] + chrom.sizes[1])

for (i in 2:5) {
	chrom.starts <- c(chrom.starts, chrom.ends[i-1] + white.space)
	chrom.ends <- c(chrom.ends, chrom.starts[i] + chrom.sizes[i])
}


tall.data$plot.position <- tall.data$pos
tall.data$plot.position <- ifelse(tall.data$chrom == "X", tall.data$plot.position + chrom.starts[1], tall.data$plot.position)
tall.data$plot.position <- ifelse(tall.data$chrom == "2L", tall.data$plot.position + chrom.starts[2], tall.data$plot.position)
tall.data$plot.position <- ifelse(tall.data$chrom == "2R", tall.data$plot.position + chrom.starts[3], tall.data$plot.position)
tall.data$plot.position <- ifelse(tall.data$chrom == "3L", tall.data$plot.position + chrom.starts[4], tall.data$plot.position)
tall.data$plot.position <- ifelse(tall.data$chrom == "3R", tall.data$plot.position + chrom.starts[5], tall.data$plot.position)



plot.genotypes <- function(column, y.value, label) {
	rect.x.side <- 200000
	rect.y.side <- 0.2

	cur.data <- tall.data[,column]
	allele.1 <- substr(cur.data, start=1, stop=1)
	allele.2 <- substr(cur.data, start=3, stop=3)
	plot.colors.1 <- allele.colors[match(allele.1, table=alleles)]
	plot.colors.2 <- allele.colors[match(allele.2, table=alleles)]
	#Draw the chromosomes
	for (c in 1:5) {
		lines(x=c(chrom.starts[c], chrom.ends[c]), y=c(y.value, y.value))
	}
	y.vector <- rep(y.value, nrow(tall.data))
	rect(xleft=tall.data$plot.position - rect.x.side, xright=tall.data$plot.position + rect.x.side, ybottom=y.vector, ytop=y.vector + rect.y.side, col=plot.colors.1, border=NA)
	rect(xleft=tall.data$plot.position - rect.x.side, xright=tall.data$plot.position + rect.x.side, ybottom=y.vector - rect.y.side, ytop=y.vector, col=plot.colors.2, border=NA)
	#points(x=tall.data$plot.position, y=rep(y.value, nrow(tall.data)), col=plot.colors.1, cex=point.cex, pch=15)
	text(x=chrom.starts[1], y=y.value, labels=paste(label, "   ", sep=""), adj=c(1, NA), cex=0.75)

}


pdf(file="FigureS1_BackgroundGenotypePlot.pdf", width=12, height=6)

par(mar=c(0,0,0,0))
plot(x=NULL, y=NULL, xlim=range(c(0, chrom.ends)), ylim=c(0, 18), bty="n", xaxt="n", yaxt="n")
chrom.names <- c("X", "2L", "2R", "3L", "3R")
for (c in 1:5) {
	text(x=(chrom.starts[c] + chrom.ends[c])/2, y=0, label=chrom.names[c])
}
for (c in 1:18) {
	plot.genotypes(column=c, y.value=c, label=colnames(tall.data)[c])
}

legend(x=0, y=18, legend=alleles, pch=15, col=allele.colors, bty="n", cex=0.75)

dev.off()