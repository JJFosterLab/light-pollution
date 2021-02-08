rm(list = ls())
graphics.off()
#################################################################
#	Useful Functions											#
#################################################################
source(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LP_PlotFunctions', '.R'))

#####################################################################
#	Install and Load Packages										#
#####################################################################

Instalload(c('circular', 'beeswarm'))
# lp.all <- read.table(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPallOrganised20200221', '.txt') , header = T)#, sep  = '\t')
lp.all <- read.table(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPdataOrganised20200505', '.txt') , header = T)#, sep  = '\t')
head(lp.all)
lp.all$Beetle <- as.factor(lp.all$Beetle)
lp.all$Condition <- relevel(lp.all$Condition, 'EyesCoveredWallsUrban')
levels(lp.all$Condition)	

#####################################################################
#	Plot Individually												#
#####################################################################
#	Plot exits														#

for(cnd in levels(lp.all$Condition)){
	ttt <- subset(lp.all, Condition == cnd)
	bbb <- length(unique(ttt$Beetle))
	dev.new(height = 3/2 *ceiling(bbb/5), width = 6.5)
	par(mfrow = c(ceiling(bbb/5),5), mai = c(0,0,0,0))
	for(btl in unique(ttt$Beetle)){
		#while I still have missing data for one beetle
		#TEMPORARY FIX
		if(sum(!is.na(subset(ttt, Beetle == btl)$Heading))){
			Cplot2(subset(ttt, Beetle == btl)$Heading, cex = 1.5, sp  = 0.08, col = c('gray20', 'orange4')[subset(ttt, Beetle == btl)$Direct+1], bg = c('gray', 'orange')[subset(ttt, Beetle == btl)$Indirect+1], pch = 21, lwd = 1.5)
		}#if(sum(!is.na(subset(ttt, Beetle == btl)$Heading)))
	}#	for(btl in unique(ttt$Beetle))
	mtext(unique(ttt$ShortName), line = -6.25*ceiling(bbb/5), outer = T, cex = 0.6)
	PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/RawPlots/'), Experiment = 'LP', PlotName = paste0("Cplot2_", cnd))
	graphics.off()
}#	for(cnd in levels(lp.all$Condition))

#####################################################################
#	Plot Summary													#
#####################################################################
#	Plot mean vectors for each experiment						#

# rhoPlot(lp.all, cnd, 'Beetle', col = 'mediumaquamarine', ci.col = 'darkgreen')
rho.col = #rgb((1+118)/256, (1+194)/256, (1+167)/256, 0.90)#cyan
			rgb((1+234)/256,(1+079)/256, (1+069)/256, 0.90)#red-pink
rho.ci.col = #'darkgreen'
			rgb((1+ 249)/256,(1+ 178)/256, (1+051)/256, 1.00)#red-pink
dev.new(width = 3, height = 3)
par(mai = c(0,0,0,0))
for(cnd in levels(lp.all$Condition)){
	ttt <- subset(lp.all, Condition == cnd)
	rhoPlot(lp.all, cnd, 'Beetle', col = rho.col, ci.col = rho.ci.col)
	mtext(cnd, line = -0.9, outer = T, cex = 0.4)
	lines(sqrt(-log(0.05)/10)*cos(seq(-pi,pi, length.out = 10^3)), 
			sqrt(-log(0.05)/10)*sin(seq(-pi,pi, length.out = 10^3)),
			col = 'gray', lty = 3, lwd = 2.5, lend = 'butt') 
	lines(1*cos(seq(-pi,pi, length.out = 10^3)), 
			1*sin(seq(-pi,pi, length.out = 10^3)),
			col = 'white', lwd = 3, lend = 'butt') 
	PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/RawPlots/'), Experiment = 'LP',  PlotName = paste0("rhoPlot_red_", cnd))
}#	for(cnd in levels(lp.all$Condition))

#####################################################################
#	Plot Mean Vectors Beeswarm Style								#
#####################################################################
#	Plot mean vectors for each condition as boxplot and beeswarm	#

#load the saved dataset
# lp.all.rho <- read.table(file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPallRho20200221', '.txt'))
lp.all.rho <- read.table(file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPallRho20200505', '.txt'))
head(lp.all.rho)
lp.all.rho$Condition <- relevel(lp.all.rho$Condition, 'EyesCoveredWallsUrban')
levels(lp.all.rho$Condition)	

#	All data
dev.new(width =10); par(mai = c(0,0.8, 0, 0), lend = 'butt')
#plot on the appropriate scale transformed	
boxplot(rho~Condition, data = lp.all.rho,
		cex = 0.5, outline = F, border = rgb(0,0,0,0.5),
		pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5),
		axes = F, ylim =(c(0,1)),# xlim = c(0.5, 4.5),
		ylab = 'Mean Vector Length', xlab = '')
polygon(c(0,2+rep(length(levels(lp.all.rho$Condition)),2),0), c(0,0, sqrt(-log(0.05)/10),sqrt(-log(0.05)/10)), col = rgb(1,0,0,0.05), border = NA)
legend('bottom', inset = sqrt(-log(0.05)/10), legend = '   Rayleigh test     \n       p<0.05     \n            \n       p>0.05     ', cex = 0.5, bty = 'n')
beeswarm(rho~Condition, data = lp.all.rho,
		pch = 21,  cex = 1.8/2, #method = 'center',
		pwcol = c('gray20', 'orange4')[lp.all.rho$Direct+1], 
		pwbg = c('gray', 'orange')[lp.all.rho$Indirect+1],
		add = T)
#	
	axis(2)#
	mtext(levels(lp.all.rho$Condition), side = 1, at = 1:length(levels(lp.all.rho$Condition))-0.5, line = -0.25-24*with(lp.all.rho, aggregate(rho, by = list(Condition = Condition), function(x) quantile(x, 0.75, na.rm=T)))$x, las = 2, cex = 0.75 )
	abline(h =c(0,1), lwd = 0.25)
	
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/RawPlots/'), Experiment = 'LP',  PlotName = paste0("Beeswarm_", 'all'))

#	Data that will be used in the paper
# lp.data.rho <- read.table(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPdataRho20200225', '.txt') , header = T)#, sep  = '\t')
lp.data.rho <- read.table(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPdataRho20200505', '.txt') , header = T)#, sep  = '\t')
head(lp.data.rho)
lp.data.rho$Beetle <- as.factor(lp.data.rho$Beetle)
lp.data.rho$Condition <- relevel(lp.data.rho$Condition, 'EyesCoveredWallsUrban')
cbind(levels(lp.data.rho$Condition))
#attempt reorganisation
lp.data.rho$Condition <- factor(lp.data.rho$Condition, levels(lp.data.rho$Condition)[c(2:3,10:11,5:6,14,4,8:9,1,7,12,13)])	

dev.new(width =10); par(mai = c(0,0.8, 0, 0), lend = 'butt')
#plot on the appropriate scale transformed	
boxplot(rho~Condition, data = lp.data.rho,
		cex = 0.5, outline = F, border = rgb(0,0,0,0.5),
		pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5),
		axes = F, ylim =(c(0,1)),# xlim = c(0.5, 4.5),
		ylab = 'Mean Vector Length', xlab = '')
polygon(c(0,2+rep(length(levels(lp.data.rho$Condition)),2),0), c(0,0, sqrt(-log(0.05)/10),sqrt(-log(0.05)/10)), col = rgb(1,0,0,0.05), border = NA)
legend('bottom', inset = sqrt(-log(0.05)/10), legend = '   Rayleigh test     \n       p<0.05     \n            \n       p>0.05     ', cex = 0.5, bty = 'n')
beeswarm(rho~Condition, data = lp.data.rho,
		pch = 21,  cex = 1.8/2, #method = 'center',
		pwcol = c('gray20', 'orange4')[lp.data.rho$Direct+1], 
		pwbg = c('gray', 'orange')[lp.data.rho$Indirect+1],
		add = T)
#	
	axis(2)#
	segments(rep(9,10), 
	subset(lp.data.rho, Condition == 'SecurityLightOFFRural')$rho,
	rep(10,10),
	subset(lp.data.rho, Condition == 'SecurityLightONRural')$rho, 
	lty = 3, rgb(0,0,1,0.3))
	mtext(levels(lp.data.rho$Condition), side = 1, at = 1:length(levels(lp.data.rho$Condition))-0.5, line = -0.25-24*with(lp.data.rho, aggregate(rho, by = list(Condition = Condition), function(x) quantile(x, 0.75, na.rm=T)))$x, las = 2, cex = 0.75 )
	abline(h =c(0,1), lwd = 0.25)
	
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/RawPlots/'), Experiment = 'LP',  PlotName = paste0("Beeswarm_", 'allpaper'))




#	Just walls
lp.walls.rho <- subset(lp.all.rho, grepl('Walls', Condition))
lp.walls.rho$Condition <- as.character(lp.walls.rho$Condition)
lp.walls.rho$Condition <- as.factor(lp.walls.rho$Condition)
lp.walls.rho$Condition <- relevel(lp.walls.rho$Condition, 'EyesCoveredWallsUrban')
dev.new(width =length(unique(lp.walls.rho$Condition))/2	+0.5); par(mai = c(0,0.8, 0, 0), lend = 'butt')
#plot on the appropriate scale transformed	
boxplot(rho~Condition, data = lp.walls.rho,
		cex = 0.5, outline = F, border = rgb(0,0,0,0.5),
		pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5),
		axes = F, ylim =(c(0,1)),# xlim = c(0.5, 4.5),
		ylab = 'Mean Vector Length', xlab = '')
polygon(c(0,2+rep(length(levels(lp.walls.rho$Condition)),2),0), c(0,0, sqrt(-log(0.05)/10),sqrt(-log(0.05)/10)), col = rgb(1,0,0,0.05), border = NA)
# legend('bottom', inset = sqrt(-log(0.05)/10), legend = '   Rayleigh test     \n       p<0.05     \n            \n       p>0.05     ', cex = 0.5, bty = 'n')
beeswarm(rho~Condition, data = lp.walls.rho,
		pch = 21,  cex = 1.8/2, #method = 'center',
		pwcol = c('gray20', 'orange4')[lp.walls.rho$Direct+1], 
		pwbg = c('gray', 'orange')[lp.walls.rho$Indirect+1],
		add = T)
#	
	axis(2)#
	mtext(levels(lp.walls.rho$Condition), side = 1, at = 1:length(levels(lp.walls.rho$Condition))-0.5, line = -0.25-24*with(lp.walls.rho, aggregate(rho, by = list(Condition = Condition), function(x) quantile(x, 0.75, na.rm=T)))$x, las = 2, cex = 0.75 )
	abline(h =c(0,1), lwd = 0.25)
	
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/RawPlots/'), Experiment = 'LP',  PlotName = paste0("Beeswarm_", 'Walls'))

#	Just walls for paper
lp.walls.paper.rho <- subset(lp.data.rho, grepl('Walls', Condition))
lp.walls.paper.rho$Condition <- as.character(lp.walls.paper.rho$Condition)
lp.walls.paper.rho$Condition <- as.factor(lp.walls.paper.rho$Condition)
lp.walls.paper.rho$Condition <- relevel(lp.walls.paper.rho$Condition, 'EyesCoveredWallsUrban')
dev.new(width =length(unique(lp.walls.paper.rho$Condition))/2	+1); par(mai = c(0,0.8, 0, 0), lend = 'butt')
#plot on the appropriate scale transformed	
boxplot(rho~Condition, data = lp.walls.paper.rho,
		cex = 0.5, outline = F, border = rgb(0,0,0,0.5),
		pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5),
		axes = F, ylim =(c(0,1)),# xlim = c(0.5, 4.5),
		ylab = 'Mean Vector Length', xlab = '')
polygon(c(0,2+rep(length(levels(lp.walls.paper.rho$Condition)),2),0), c(0,0, sqrt(-log(0.05)/10),sqrt(-log(0.05)/10)), col = rgb(1,0,0,0.05), border = NA)
# legend('bottom', inset = sqrt(-log(0.05)/10), legend = '   Rayleigh test     \n       p<0.05     \n            \n       p>0.05     ', cex = 0.5, bty = 'n')
beeswarm(rho~Condition, data = lp.walls.paper.rho,
		pch = 21,  cex = 1.8/2, #method = 'center',
		pwcol = c('gray20', 'orange4')[lp.walls.paper.rho$Direct+1], 
		pwbg = c('gray', 'orange')[lp.walls.paper.rho$Indirect+1],
		add = T)
#	
	axis(2)#
	mtext(levels(lp.walls.paper.rho$Condition), side = 1, at = 1:length(levels(lp.walls.paper.rho$Condition))-0.5, line = -0.25-24*with(lp.walls.paper.rho, aggregate(rho, by = list(Condition = Condition), function(x) quantile(x, 0.75, na.rm=T)))$x, las = 2, cex = 0.75 )
	abline(h =c(0,1), lwd = 0.25)
	
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/RawPlots/'), Experiment = 'LP',  PlotName = paste0("Beeswarm_", 'paperWalls'))

#	Just starry
lp.stars.rho <- subset(lp.all.rho, grepl('Star', Condition) | grepl('MW', Condition))
lp.stars.rho$Condition <- as.character(lp.stars.rho$Condition)
lp.stars.rho$Condition <- as.factor(lp.stars.rho$Condition)
lp.stars.rho$Condition <- relevel(lp.stars.rho$Condition, 'StarsRural')
dev.new(width =length(unique(lp.stars.rho$Condition))/2	+0.5); par(mai = c(0,0.8, 0, 0), lend = 'butt')
#plot on the appropriate scale transformed	
boxplot(rho~Condition, data = lp.stars.rho,
		cex = 0.5, outline = F, border = rgb(0,0,0,0.5),
		pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5),
		axes = F, ylim =(c(0,1)),# xlim = c(0.5, 4.5),
		ylab = 'Mean Vector Length', xlab = '')
polygon(c(0,2+rep(length(levels(lp.stars.rho$Condition)),2),0), c(0,0, sqrt(-log(0.05)/10),sqrt(-log(0.05)/10)), col = rgb(1,0,0,0.05), border = NA)
# legend('bottom', inset = sqrt(-log(0.05)/10), legend = '   Rayleigh test     \n       p<0.05     \n            \n       p>0.05     ', cex = 0.5, bty = 'n')
beeswarm(rho~Condition, data = lp.stars.rho,
		pch = 21,  cex = 1.8/2, #method = 'center',
		pwcol = c('gray20', 'orange4')[lp.stars.rho$Direct+1], 
		pwbg = c('gray', 'orange')[lp.stars.rho$Indirect+1],
		add = T)
#	
	axis(2)#
	mtext(levels(lp.stars.rho$Condition), side = 1, at = 1:length(levels(lp.stars.rho$Condition))-0.5, line = -0.25-24*with(lp.stars.rho, aggregate(rho, by = list(Condition = Condition), function(x) quantile(x, 0.75, na.rm=T)))$x, las = 2, cex = 0.75 )
	abline(h =c(0,1), lwd = 0.25)
	
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/RawPlots/'), Experiment = 'LP',  PlotName = paste0("Beeswarm_", 'Stars'))

#	Just starry
lp.stars.paper.rho <- subset(lp.data.rho, grepl('Star', Condition) )
lp.stars.paper.rho$Condition <- as.character(lp.stars.paper.rho$Condition)
lp.stars.paper.rho$Condition <- as.factor(lp.stars.paper.rho$Condition)
lp.stars.paper.rho$Condition <- relevel(lp.stars.paper.rho$Condition, 'StarsRural')
dev.new(width =length(unique(lp.stars.paper.rho$Condition))/2	+1); par(mai = c(0,0.8, 0, 0), lend = 'butt')
#plot on the appropriate scale transformed	
boxplot(rho~Condition, data = lp.stars.paper.rho,
		cex = 0.5, outline = F, border = rgb(0,0,0,0.5),
		pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5),
		axes = F, ylim =(c(0,1)),# xlim = c(0.5, 4.5),
		ylab = 'Mean Vector Length', xlab = '')
polygon(c(0,2+rep(length(levels(lp.stars.paper.rho$Condition)),2),0), c(0,0, sqrt(-log(0.05)/10),sqrt(-log(0.05)/10)), col = rgb(1,0,0,0.05), border = NA)
# legend('bottom', inset = sqrt(-log(0.05)/10), legend = '   Rayleigh test     \n       p<0.05     \n            \n       p>0.05     ', cex = 0.5, bty = 'n')
beeswarm(rho~Condition, data = lp.stars.paper.rho,
		pch = 21,  cex = 1.8/2, #method = 'center',
		pwcol = c('gray20', 'orange4')[lp.stars.paper.rho$Direct+1], 
		pwbg = c('gray', 'orange')[lp.stars.paper.rho$Indirect+1],
		add = T)
#	
	axis(2)#
	mtext(levels(lp.stars.paper.rho$Condition), side = 1, at = 1:length(levels(lp.stars.paper.rho$Condition))-0.5, line = -0.25-24*with(lp.stars.paper.rho, aggregate(rho, by = list(Condition = Condition), function(x) quantile(x, 0.75, na.rm=T)))$x, las = 2, cex = 0.75 )
	abline(h =c(0,1), lwd = 0.25)
	
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/RawPlots/'), Experiment = 'LP',  PlotName = paste0("Beeswarm_", 'paperStars'))

#	Just cloudy
lp.clouds.rho <- subset(lp.all.rho, grepl('Cloud', Condition) | grepl('Overcast', Condition))
lp.clouds.rho$Condition <- as.character(lp.clouds.rho$Condition)
lp.clouds.rho$Condition <- as.factor(lp.clouds.rho$Condition)
lp.clouds.rho$Condition <- relevel(lp.clouds.rho$Condition, 'OvercastRural')
dev.new(width =length(unique(lp.clouds.rho$Condition))/2	+0.5); par(mai = c(0,0.8, 0, 0), lend = 'butt')
#plot on the appropriate scale transformed	
boxplot(rho~Condition, data = lp.clouds.rho,
		cex = 0.5, outline = F, border = rgb(0,0,0,0.5),
		pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5),
		axes = F, ylim =(c(0,1)),# xlim = c(0.5, 4.5),
		ylab = 'Mean Vector Length', xlab = '')
polygon(c(0,2+rep(length(levels(lp.clouds.rho$Condition)),2),0), c(0,0, sqrt(-log(0.05)/10),sqrt(-log(0.05)/10)), col = rgb(1,0,0,0.05), border = NA)
# legend('bottom', inset = sqrt(-log(0.05)/10), legend = '   Rayleigh test     \n       p<0.05     \n            \n       p>0.05     ', cex = 0.5, bty = 'n')
beeswarm(rho~Condition, data = lp.clouds.rho,
		pch = 21,  cex = 1.8/2, #method = 'center',
		pwcol = c('gray20', 'orange4')[lp.clouds.rho$Direct+1], 
		pwbg = c('gray', 'orange')[lp.clouds.rho$Indirect+1],
		add = T)
#	
	axis(2)#
	mtext(levels(lp.clouds.rho$Condition), side = 1, at = 1:length(levels(lp.clouds.rho$Condition))-0.5, line = -0.25-24*with(lp.clouds.rho, aggregate(rho, by = list(Condition = Condition), function(x) quantile(x, 0.75, na.rm=T)))$x, las = 2, cex = 0.75 )
	abline(h =c(0,1), lwd = 0.25)
	
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/RawPlots/'), Experiment = 'LP',  PlotName = paste0("Beeswarm_", 'Clouds'))

#	Just cloudy for paper
lp.clouds.paper.rho <- subset(lp.data.rho, grepl('Cloud', Condition) | grepl('Overcast', Condition))
lp.clouds.paper.rho$Condition <- as.character(lp.clouds.paper.rho$Condition)
lp.clouds.paper.rho$Condition <- as.factor(lp.clouds.paper.rho$Condition)
lp.clouds.paper.rho$Condition <- relevel(lp.clouds.paper.rho$Condition, 'OvercastRural')
dev.new(width =length(unique(lp.clouds.paper.rho$Condition))/2	+1); par(mai = c(0,0.8, 0, 0), lend = 'butt')
#plot on the appropriate scale transformed	
boxplot(rho~Condition, data = lp.clouds.paper.rho,
		cex = 0.5, outline = F, border = rgb(0,0,0,0.5),
		pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5),
		axes = F, ylim =(c(0,1)),# xlim = c(0.5, 4.5),
		ylab = 'Mean Vector Length', xlab = '')
polygon(c(0,2+rep(length(levels(lp.clouds.paper.rho$Condition)),2),0), c(0,0, sqrt(-log(0.05)/10),sqrt(-log(0.05)/10)), col = rgb(1,0,0,0.05), border = NA)
# legend('bottom', inset = sqrt(-log(0.05)/10), legend = '   Rayleigh test     \n       p<0.05     \n            \n       p>0.05     ', cex = 0.5, bty = 'n')
beeswarm(rho~Condition, data = lp.clouds.paper.rho,
		pch = 21,  cex = 1.8/2, #method = 'center',
		pwcol = c('gray20', 'orange4')[lp.clouds.paper.rho$Direct+1], 
		pwbg = c('gray', 'orange')[lp.clouds.paper.rho$Indirect+1],
		add = T)
#	
	axis(2)#
	mtext(levels(lp.clouds.paper.rho$Condition), side = 1, at = 1:length(levels(lp.clouds.paper.rho$Condition))-0.5, line = -0.25-24*with(lp.clouds.paper.rho, aggregate(rho, by = list(Condition = Condition), function(x) quantile(x, 0.75, na.rm=T)))$x, las = 2, cex = 0.75 )
	abline(h =c(0,1), lwd = 0.25)
	
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/RawPlots/'), Experiment = 'LP',  PlotName = paste0("Beeswarm_", 'paperClouds'))


#	Just moon
lp.moon.rho <- subset(lp.all.rho, grepl('Moon', Condition) | grepl('SecurityLightOFFRural', Condition))
lp.moon.rho$Condition <- as.character(lp.moon.rho$Condition)
lp.moon.rho$Condition <- as.factor(lp.moon.rho$Condition)
# lp.moon.rho$Condition <- relevel(lp.moon.rho$Condition, 'OvercastRural')
dev.new(width =length(unique(lp.moon.rho$Condition))/2	+0.5); par(mai = c(0,0.8, 0, 0), lend = 'butt')
#plot on the appropriate scale transformed	
boxplot(rho~Condition, data = lp.moon.rho,
		cex = 0.5, outline = F, border = rgb(0,0,0,0.5),
		pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5),
		axes = F, ylim =(c(0,1)),# xlim = c(0.5, 4.5),
		ylab = 'Mean Vector Length', xlab = '')
polygon(c(0,2+rep(length(levels(lp.moon.rho$Condition)),2),0), c(0,0, sqrt(-log(0.05)/10),sqrt(-log(0.05)/10)), col = rgb(1,0,0,0.05), border = NA)
# legend('bottom', inset = sqrt(-log(0.05)/10), legend = '   Rayleigh test     \n       p<0.05     \n            \n       p>0.05     ', cex = 0.5, bty = 'n')
beeswarm(rho~Condition, data = lp.moon.rho,
		pch = 21,  cex = 1.8/2, #method = 'center',
		pwcol = c('gray20', 'orange4')[lp.moon.rho$Direct+1], 
		pwbg = c('gray', 'orange')[lp.moon.rho$Indirect+1],
		add = T)
#	
	axis(2)#
	mtext(levels(lp.moon.rho$Condition), side = 1, at = 1:length(levels(lp.moon.rho$Condition))-0.5, line = -0.25-24*with(lp.moon.rho, aggregate(rho, by = list(Condition = Condition), function(x) quantile(x, 0.75, na.rm=T)))$x, las = 2, cex = 0.75 )
	abline(h =c(0,1), lwd = 0.25)
	
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/RawPlots/'), Experiment = 'LP',  PlotName = paste0("Beeswarm_", 'Moon'))