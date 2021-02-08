rm(list = ls())
graphics.off()
#################################################################
#	Useful Functions											#
#################################################################
#run these functions to add them to your workspace
#(highlight, cmd + enter)
#Check for installed packages, choose Swedish CRAN mirror, install packages, load packages from library #New and improved
Instalload <- function(Required){
	if(sum(rownames(installed.packages()) %in% Required, na.rm = T)<length(Required)){
	mirrors <- getCRANmirrors()
	chooseCRANmirror(graphics=FALSE, ind = which(mirrors$Country == 'Sweden'))
	install.packages(Required)
 }
 for(i in 1:length(Required)){
 	library( paste0(Required[i]) , character.only = T)}
###################	END OF FUNCTION	###########################
 }

#	Save current plot with dimensions displayed on screen and label accordingly
PDFsave <- function(Directory, Experiment, Species, PlotName, Dim){
	if(missing(Directory)){Directory <- paste0(Sys.getenv('HOME'),'/Dropbox/Plotmania/')}
	if(missing(Experiment)){Experiment = '_'}
	if(missing(Species)){Species = '_'}
	if(missing(PlotName)){PlotName = date()}
	if(missing(Dim)){Dim <- par("din")}
	#Copy to invisible device to save
	suppressWarnings(
	dev.copy(pdf, paste0(Directory, Experiment, Species, PlotName, '.pdf'),
	width = Dim[1], height = Dim[2], useDingbats = F)
	)
	#Save the image actual size
	dev.off()#Doesn't save until invisible copy is closed
	dev.set(dev.prev())
###################	END OF FUNCTION	###########################
}

AN <- as.numeric	#shorten 'as.numeric'
degrad <- pi/180#ratio of radians to degrees

#	Functions for Circular Data									#
# convert to bearings in "circular" package's preferred format
#i.e. clockwise starting at the graph's top
 mycirc <- function(angles, clock){
	if(missing(clock)){clock <- T}
	if(clock){
	return(
	as.circular(angles,
		units = 'degrees',
		type = 'angles', #don't set this to directions, apparently very different
		modulo = '2pi',
		zero = pi/2,
		rotation = 'clock',
		template = 'none')
		)
		}else{
		as.circular(angles,
		units = 'degrees',
		type = 'angles', #don't set this to directions, apparently very different
		modulo = '2pi',
		zero = pi/2,
		rotation = 'counter',
		template = 'none')
		}#if(clock)
###################	END OF FUNCTION	###########################
}

# plot a circular histogram with stacked points, a mean vector and ±1 s.d. error bars
Cplot <- function(headings, sp, bt, ax, ...){
	#fit mean axis, fits mean direction unless otherwise specified
	if(missing(ax)){ax <- F}
	#spacing of stacked points, now automatically stacks towards centre unless otherwise specified
	if(missing(sp) & missing(bt)){sp <- 0.04}
	#bt specifies the stacking by a multipicative factor, 1 = stacked, 2 = 1 point's space between, 0.5 = half overlapping
	if( missing(sp) & !(missing(bt)) ){sp <- bt*.04}
	#	Get functions and packages
	if(sum(rownames(installed.packages()) %in% c('CircStats', 'circular'), na.rm = T)<2){
		install.packages(c('CircStats', 'circular'))
	}#if(sum(rownames(installed.packages()) %in% c('CircStats', 'circular'), na.rm = T)<2)
	if(!(	sum('CircCI'%in% ls())	)){
		CircCI <- function(mn, lci, uci, out, zro, drc, lng, ...){
			if(missing(lng)){lng<-10*360/5};	if(missing(drc)){drc<-'clock'}
			if(missing(zro)){zro <- pi/2};if(missing(out)){out <- 0.05}
			if(missing(uci)){uci <- lci}
			lwr <- mn - lci;	upr <- mn + uci
			circ.pos <- ( ((drc == 'clock')-1)*2 +1) *
				-seq( pi*lwr/180, pi*upr/180, length.out = lng) + zro
			circ.x <- cos(circ.pos)*(1+out);	circ.y <- sin(circ.pos)*(1+out)
			lines(circ.x, circ.y, ...)
			lines.circular( as.circular(rep(lwr,2),units = 'degrees',
				type = 'angles', modulo = '2pi', zero = zro,
				rotation = drc, template = 'none'),
				out*c(0.5, 1.5), modulo = '2pi',
				zero = zro, rotation = drc, ...)
			lines.circular(as.circular(rep(upr,2),units = 'degrees',
				type = 'angles', modulo = '2pi', zero = zro,
				rotation = drc, template = 'none'),
			 	out*c(0.5, 1.5), modulo = '2pi', zero = zro,
			 	rotation = drc, ...)
		}#CircCI <- function(mn, lci, uci, out, zro, drc, lng, ...)
	}#if(!(	sum('CircCI'%in% ls())	))
	if(!(	sum('mycirc'%in% ls())	)){
		mycirc <- function(angles, clock){
			if(missing(clock)){clock <- T}
			if(clock){
			return(		as.circular(angles,units='degrees',type='angles',
			modulo='2pi',zero=pi/2,rotation='clock',	template='none')	)
				}else{
				as.circular(angles,units='degrees',type='angles',
				 modulo='2pi',zero=pi/2,rotation='counter',template='none')
				}#if(clock)
			}#mycirc <- function(angles, clock)
	}#if(!(	sum('mycirc'%in% ls())	))
	#circular plot settings
	increments <- 5 #degrees
	zr <- pi/2 #start at top of screen (pi*	90	/180)
	bn <- 10*10*360/5 #bins
	degrad <- 180/pi #conversion from radians to degrees
	tcl <- rgb(1,1,1,0)#transparent colour
	pcl <- rgb(.3,.1,.1,.5)#point colour
	#plot characters
	lw <- 0.5 #line width
	pnt <- 2.5 #point size
	arw <- 10 #arrowhead angle
	arl <- 0.1 #arrowhead length
	#	set up input variables
	hd <- mycirc(headings)
	sm <- summary(hd)
	sv <- degrad*sd.circular(hd, na.rm=T)
	lbl <- 90*(1:4-1)
	plot(hd, col=tcl, main="", zero=zr, axes=F, shrink=1,tol=0.075)
	axis.circular(1, at = mycirc(lbl), labels = paste0(lbl, 'º'))
	par(new=T)
	plot.circular(hd, col=tcl,main="",zero=zr,axes=F,shrink=1.05,tol=0.075)
	points(hd,stack=T,bin=bn,sep=-sp,zero=zr,...)
	if(!(ax)){
		arrows.circular( mycirc(sm['Mean']),zero=zr,col='red4',lwd=3,
		 length=arl,angle=arw,shrink = sm['Rho'])
		 CircCI(sm['Mean'], sv, out = 0.15, zro=zr, drc='clock',col='red4',lwd=1)	}else{
		 sm2 <- summary(mycirc(hd*2))
		 sv2 <- degrad*sd.circular(hd*2, na.rm=T)/2
		 arrows.circular( mycirc(sm2['Mean']/2),zero=zr,col='red4',lwd=3,
		 length=arl,angle=arw,shrink = sm2['Rho'])
		 arrows.circular( mycirc(180+sm2['Mean']/2),zero=zr,col='red4',lwd=3,
		 length=arl,angle=arw,shrink = sm2['Rho'])
		 CircCI(sm2['Mean']/2, sv2, out = 0.15, zro=zr, drc='clock',col='red4',lwd=1)
		 CircCI(180+sm2['Mean']/2, sv2, out = 0.15, zro=zr, drc='clock',col='red4',lwd=1)
	 }#if(!(ax))
###################	END OF FUNCTION	###########################
}
#####################################################################
#	Plot Settings											#
#####################################################################

#axis size multiplication factor
mf <- 2.559055# plots are 65mm tall 2.375 #plots are 60.325 mm tall
pl <- T#F

#####################################################################
#	Install and Load Packages										#
#####################################################################

Instalload(c('circular', 'beeswarm'))
lp.all <- read.table(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/SecurityLightBergsig', '.txt') , header = T, sep  = '\t')

lp.all.rho <- data.frame(Experiment = sort(rep(levels(lp.all$Experiment), 10)),
		Beetle = rep(1:10, length(levels(lp.all$Experiment))) )

for(epr in levels(lp.all$Experiment)){
	for(btl in 1:10){
	lp.all.rho$rho[lp.all.rho$Experiment == epr & lp.all.rho$Beetle == btl] <-
	 rho.circular( mycirc(subset(lp.all, Experiment == epr & Beetle == btl)$Heading), na.rm = T )
	}#for(btl in 1:10)
}#for(exp in levels(lp.all$Experiment))

lp.all.rho$ShortName <- lp.all.rho$Experiment
levels(lp.all.rho$ShortName) <- c('Stars\nRural','Lamp\nRural')

write.table(lp.all.rho, file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPallRho-SecurityLightBergsig', '.txt'))

lp.all.rho <- read.table(file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPallRho-SecurityLightBergsig', '.txt'))

dev.new(height = 3, width = 6.5)
par(mfrow = c(2,5), mai = c(0,0,0,0))
for(btl in 1:10){
Cplot(subset(lp.all, Experiment == 'LPexper.20191117.SecurityLightBergsigOFF' & Beetle == btl)$Heading, cex = 2, bt  = 2, col = 'gray')
}#for(btl in 1:10)
PDFsave(Directory = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/'), Experiment = 'LP', PlotName = "Stars.SecurityLightBergsig")

# dev.new(height = 7, width = 14)
par(mfrow = c(2,5), mai = c(0,0,0,0))
for(btl in 1:10){
Cplot(subset(lp.all, Experiment == 'LPexper.20191117.SecurityLightBergsigON' & Beetle == btl)$Heading, cex = 2, bt = 2, col = 'orange4')
}#for(btl in 1:10)
# mtext('Gibbous Moon Wits', outer = T, line = -3.2, cex = 3)
PDFsave(Directory = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/'), Experiment = 'LP', PlotName = "Lamp.SecurityLightBergsig")

dev.new(width =3); par(mai = c(0.8,0.8, 0, 0))
#plot on the appropriate scale transformed
boxplot(rho~ShortName, data = lp.all.rho, ylim =(c(0,1)), xlim = c(0.5, 2.5),
cex = 0.5, outline = F, border = rgb(0,0,0,1),
		pars = list(boxwex = 0.3/2, staplewex = 0.5, outwex = 0.5),
#colour code by sex, turn axes off
	axes = F, col = c('salmon', 'slateblue'),
	ylab = 'Mean Vector Length', xlab = '')
#log scale axes with original data values
	axis(2)#, at =log10( pretty( range(soil.dt2$CO2_Flux) *2 )/2 ), labels = pretty( range(soil.dt2$CO2_Flux) *2 )/2)
	mtext(levels(lp.all.rho$ShortName), side = 3, at = 1:length(levels(lp.all.rho$ShortName)), line = -12 )
	abline(h =c(0,1), lwd = 0.25)

beeswarm(rho~ShortName, data = lp.all.rho,
		pch = 20, method = 'center', cex = 1.8,
		pwcol = c('orange4','gray')[lp.all.rho$ShortName],
		add = T, axes = F)

PDFsave(Directory = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/'), Experiment = 'LP', PlotName = "SecurityLightBergsig-Beeswarm")

boxplot(rho~ShortName, data = lp.all.rho, ylim =(c(0,1)), xlim = c(0.5, 2.5),
cex = 0.5, outline = F, border = rgb(0,0,0,1),
		pars = list(boxwex = 0.3/2, staplewex = 0.5, outwex = 0.5),
#colour code by sex, turn axes off
	axes = F, col = c('salmon', 'slateblue'),
	ylab = 'Mean Vector Length', xlab = '')
#log scale axes with original data values
	axis(2)#, at =log10( pretty( range(soil.dt2$CO2_Flux) *2 )/2 ), labels = pretty( range(soil.dt2$CO2_Flux) *2 )/2)
	mtext(levels(lp.all.rho$ShortName), side = 3, at = 1:length(levels(lp.all.rho$ShortName)), line = -12 )
	abline(h =c(0,1), lwd = 0.25)

for(btl in unique(lp.all.rho$Beetle)){
	lines(2:1, lp.all.rho$rho[lp.all.rho$Beetle == btl], type = 'o', pch = 21, bg = 'white', col = as.numeric(btl))
}#for(btl in levels(wd$Beetle))

PDFsave(Directory = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/'), Experiment = 'LP', PlotName = "SecurityLightBergsig-Paired")

lp.all.mu <- data.frame(Experiment = sort(rep(levels(lp.all$Experiment), 10)),
		Beetle = rep(1:10, length(levels(lp.all$Experiment))) )

for(epr in levels(lp.all$Experiment)){
	for(btl in 1:10){
	lp.all.mu$mu[lp.all.mu$Experiment == epr & lp.all.mu$Beetle == btl] <- mean.circular( mycirc(subset(lp.all, Experiment == epr & Beetle == btl)$Heading), na.rm = T )
	}#for(btl in 1:10)
}#for(exp in levels(lp.all$Experiment))

lp.all.mu$ShortName <- lp.all.mu$Experiment
levels(lp.all.mu$ShortName) <- c('Stars\nRural','Lamp\nRural')

dev.new(height = 3, width = 6.5)
par(mfrow = c(1,3), mai = c(0,0,0,0))
Cplot(subset(lp.all.mu, Experiment == 'LPexper.20191117.SecurityLightBergsigON')$mu, cex = 2, bt  = 2, col = 'orange4')
mtext('Mean Headings\nSecurity Light On', line = -3)
Cplot(subset(lp.all.mu, Experiment == 'LPexper.20191117.SecurityLightBergsigOFF')$mu, cex = 2, bt  = 2, col = 'gray')
mtext('Mean Headings\nSecurity Light Off', line = -3)
Cplot(subset(lp.all.mu, Experiment == 'LPexper.20191117.SecurityLightBergsigOFF')$mu - subset(lp.all.mu, Experiment == 'LPexper.20191117.SecurityLightBergsigON')$mu, cex = 2, bt  = 2, col = 'salmon')
mtext('Mean On - Mean Off\nIndividual Changes', line = -3)

PDFsave(Directory = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/'), Experiment = 'LP', PlotName = "Mean Headings.SecurityLightBergsig")
