### Script for analysing and plotting the wing margin proportion for sd mutants from antibody staining data###
# Library and Data
require(sciplot)

sd.margin.data <- read.csv("../Data/Fig6D_sd_margin_length_formatted.csv",h=T)
#sd.margin.data <- read.csv("~/Dropbox/ChandlerChariEtAl2013/ChandlerChariPLoSGenetics2017_DataScripts/Data/Fig6D_sd_margin_length_formatted.csv",h=T)

str(sd.margin.data)

sd.margin.data$Mutation <- factor(sd.margin.data$Mutation, 
    levels(sd.margin.data$Mutation)[c(5,1,4,3,2)]) # reorder data so that mutations are ordered in the order of increasing severity

# Custom CI function. But used the native one from sciplot for the paper (see below)
ci.fun <- function(x) {
	x <- na.omit(x)
	std.err <- sd(x)/sqrt(length(x))
	lower.ci <- mean(x) + std.err*(qt(p=0.05/2, df = length(x) - 1))
	upper.ci <- mean(x) + std.err*(qt(p=1-(0.05/2), df = length(x) - 1))
	return(c(lower.ci, upper.ci))
}


lineplot.CI(sd.margin.data$Mutation, sd.margin.data$Margin_Proportion,
    group=sd.margin.data$Background, fixed=T, legend=F, 
    ci.fun = ci.fun,
    #ci.fun= function(x) c(mean(x) - 1.97*se(x), mean(x) + 1.97*se(x)),
    xaxt="n", xlab= "Mutation", ylab="Margin Proportion", 
    col=c("red","blue"), pch=c(16,16), lty=c(1,1), 
    cex.main=1.0, lwd=c(2,2), cex.lab=1.1)
    
legend (x="topright", xjust=0, legend= c("Oregon-R", "Samarkand"),
    col=c("red","blue"), pch= c(16,16), lty=c(1,1), lwd=c(2,2), bty="n", cex=1)

xaxis.labels.sd <- c( expression(paste(WT)),
				  expression(italic(paste(sd^'1',sep=''))),
				  expression(italic(paste(sd^'ETX4',sep=''))),
				  expression(italic(paste(sd^'E3',sep=''))),
				  expression(italic(paste(sd^'58d',sep='')))
				  )

axis(side=1, at=1:5, labels=xaxis.labels.sd, tick=F, line=-0.3)


#### Script for analysing and plotting the wing margin proportion for vg mutants from antibody staining data ###

vg.margin.data <- read.csv("../Data/Fig6F_vg_margin_length_formatted.csv", h = T)
str(vg.margin.data) 
vg.margin.data$Mutation <- factor(vg.margin.data$Mutation,
    levels(vg.margin.data$Mutation)[c(4,3,2,1)]) #reorder data so that mutations are ordered in the order of increasing severity

lineplot.CI(vg.margin.data$Mutation, vg.margin.data$Margin_Proportion,
    group=vg.margin.data$Background, 
    fixed=T, legend=F, 
    #ci.fun= function(x) c(mean(x) - 1.97*se(x), mean(x) + 1.97*se(x)), 
    ci.fun = ci.fun, 
    xaxt="n", xlab= "Mutation", ylab="Margin Proportion", 
    col=c("red","blue"), pch=c(16,16), lty=c(1,1), 
    cex.main=1.0, lwd=c(2,2), cex.lab=1.1)
    
legend (x="topright", xjust=0, legend= c("Oregon-R", "Samarkand"), col=c("red","blue"), pch= c(16,16), lty=c(1,1), bty="n",cex=1)

x.axis.labels.vg <- c( expression(paste(WT)),
				  expression(italic(paste(vg^'2a33',sep=''))),
				  expression(italic(paste(vg^'21-3',sep=''))),
				  expression(italic(paste(vg^'1',sep='')))
				  )

axis(side=1, at=1:4, labels=x.axis.labels.vg, tick=F, line=-0.3)


