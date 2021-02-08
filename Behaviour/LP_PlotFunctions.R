#run these functions to add them to your workspace
#(highlight, cmd + enter)
#Check for installed packages, choose local CRAN mirror, install packages, load packages from library #New and improved
Instalload <- function(Required, Country = 'Sweden'){
	if(sum(rownames(installed.packages()) %in% Required, na.rm = T)<length(Required)){
	mirrors <- getCRANmirrors()
	chooseCRANmirror(graphics=FALSE, ind = which(mirrors$Country == Country))
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

#	Update of Cplot												#
# plot points stacked towards the centre of a circle
# add a mean vector	
# add mean vector's CI in the correct direction
# add a circular confidence interval
# accomodate axial data 
Cplot2 <- function(x, alpha, ax, sp , rho.col, out.by, kappa.ci, ...){
	if(missing(ax)){ax <- F}#fit mean axis (not direction)
	if(missing(rho.col)){rho.col <- 'deeppink4'}#fit mean axis (not direction)
	if(missing(out.by)){out.by = 0.05} # draw outside (or inside) by
	if(missing(alpha)){alpha = 0.05} #proportion NOT to plot across
		#spacing of stacked points, now automatically stacks towards centre unless otherwise specified
	if(missing(kappa.ci)){kappa.ci = F} #proportion NOT to plot across
		#spacing of stacked points, now automatically stacks towards centre unless otherwise specified
	if(missing(sp)){sp <- 0.04}
	if(sum(rownames(installed.packages()) %in% c('CircStats', 'circular'), na.rm = T)<2){
		install.packages(c('CircStats', 'circular'))
	}#if(sum(rownames(installed.packages()) %in% c('CircStats', 'circular'), na.rm = T)<2)
	if(sum(c('CircStats', 'circular')  %in% (.packages()), na.rm = T)<2){
		library(circular)
		library(CircStats)
	}#if(sum(rownames(installed.packages()) %in% c('CircStats', 'circular'), na.rm = T)<2)
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
	hd <- as.circular(x,units='degrees',type='angles',
			modulo='2pi',zero=pi/2,rotation='clock',	template='none')
	sm <- summary(hd)
	sv <- degrad*sd.circular(hd, na.rm=T)
	lbl <- 90*(1:4-1)
	plot(hd, col=tcl, main="", zero=zr, axes=F, shrink=1,tol=0.075)
	axis.circular(1, at = mycirc(lbl), labels = paste0(lbl, 'º'))
	par(new=T)
	plot.circular(hd, col=tcl,main="",zero=zr,axes=F,shrink=1.05,tol=0.075)
	points(hd,stack=T,bin=bn,sep=-sp,zero=zr,...)
	par(new=T)
	plot(hd, col=tcl, main="", zero=zr, axes=F, shrink=1,tol=0.075)
	if(ax == F){
		if(kappa.ci){
		lines.circular(mycirc(rep(mean(hd, na.rm =T),2)), y = -(1-A1(mle.vonmises.bootstrap.ci(hd, alpha = 0.05, bias = T, reps = 10^4)$kappa.ci)), col = rgb(1,0,0,0.2), lwd = 7, lend = 'butt')
		}#if(kappa.ci)
		arrows.circular( mean(hd, na.rm =T),zero=zr, col = rho.col,lwd=3,
			 length=arl,angle=arw,shrink = rho.circular(hd,na.rm =T))
		#find the lower and upper confidence interval bounds
		ci <- mle.vonmises.bootstrap.ci(hd, alpha = alpha, bias = T, 
	          reps = 10^4)$mu.ci
		# if(abs(diff(ci)) > 180){#check if the default angle is "little"
			# #otherwise, plot in the reverse direction
			# sq = seq( from = as.numeric(min(ci)),
						# to = -as.numeric(360 - max(ci)), length.out = 10^3)}else{ #if it is "little" plot from min to max
			# sq = seq( from = as.numeric(min(ci)),
						# to = as.numeric(max(ci)), length.out = 10^3)
		# }#if(abs(diff(qt)) > pi)
		
		if(sign(max(ci) - mean(hd, na.rm =T)) == -1 | sign(mean(hd, , na.rm =T) - min(ci)) == -1 ){#check if mean is right or left of max CI
			#otherwise, plot in the reverse direction
			sq = seq( from = as.numeric(min(ci)),
						to = -as.numeric(360 - max(ci)), length.out = 10^3)}else{ #if it is "little" plot from min to max
			sq = seq( from = as.numeric(min(ci)),
						to = as.numeric(max(ci)), length.out = 10^3)
		}#if(abs(diff(qt)) > pi)
		
		#draw a line in coordinates of rcos(x) rsin(y)
		#any other line arguments (i.e. col or lwd) input as "..."
		lines.circular(mycirc(sq), rep(out.by,10^3), col = rho.col, lwd = 3)
	}else{
			if(!(sum(rownames(installed.packages()) %in% c('CircMLE'), na.rm = T))){
		install.packages(c('CircMLE'))
	}#if(sum(rownames(installed.packages()) %in% c('CircStats', 'circular'), na.rm = T)<2)
		if(	!('CircMLE' %in% (.packages()))	){
			library(CircMLE)
		}#if(!('CircMLE' %in% .packages))
		warning('Axial confidence intervals require more computing time...')
		qM4A <- function(x){mm <- M4A(x, 10^9, 4, 'BFGS', 1000, 0.25); 							return(mm$par)}
		mean.bs <- boot(data = hd, statistic = qM4A, R = 10^3, stype="i", sim = 'parametric', parallel = c("multicore"), ncpus = parallel::detectCores()-1)
		mu.estimates <- -mean.bs$t[,1]*180/pi -90
		kappa.estimates <- mean.bs$t[,2]
		m1 <- median(mycirc(mu.estimates), na.rm = T)[1]
		k1 <- median(kappa.estimates, na.rm = T)[1]
		m1.dist <- mu.estimates[
					mycirc(mu.estimates) < mycirc(m1 + 90) &
		 			mycirc(mu.estimates) > mycirc(m1 - 90) ]
		m2.dist <- mu.estimates[
					!(mu.estimates %in%  m1.dist)]#the rest
		m1.q <- quantile(mycirc(m1.dist), c(0.025, 0.975))
		m2.q <- quantile(mycirc(m2.dist), c(0.025, 0.975))
		k.q <- quantile(kappa.estimates, c(0.025, 0.975))
		if(kappa.ci){
			lines.circular(mycirc(rep(m1,2)), y = -(1-A1(k.q)), col = rgb(1,0,0,0.2), lwd = 7, lend = 'butt')
			lines.circular(mycirc(rep(m1 +180,2)), y = -(1-A1(k.q)), col = rgb(1,0,0,0.2), lwd = 7, lend = 'butt')
		}#if(kappa.ci)
		arrows.circular( mycirc(m1),zero=zr, col = rho.col,lwd=3,
			 length=arl,angle=arw,shrink = A1(k1))
		arrows.circular( mycirc(m1+180),zero=zr, col = rho.col,lwd=3,
			 length=arl,angle=arw,shrink = A1(k1))
		if(abs(diff(m1.q)) > 180){#check if the default angle is "little"
		#otherwise, plot in the reverse direction
			sq1 = seq( from = as.numeric(min(m1.q)),
					to = -as.numeric(360 - max(m1.q)), length.out = 10^3)}else{
			sq1 = seq( from = as.numeric(min(m1.q)),
						to = as.numeric(max(m1.q)), length.out = 10^3)	}
		if(abs(diff(m2.q)) > 180){#check if the default angle is "little"
		#otherwise, plot in the reverse direction
			sq2 = seq( from = as.numeric(min(m2.q)),
					to = -as.numeric(360 - max(m2.q)), length.out = 10^3)}else{
			sq2 = seq( from = as.numeric(min(m2.q)),
						to = as.numeric(max(m2.q)), length.out = 10^3)	}
		
		lines.circular(mycirc(sq1), rep(out.by,10^3), col = rho.col, lwd = 3)
		lines.circular(mycirc(sq2), rep(out.by,10^3), col = rho.col, lwd = 3)
			 
	}#if(ax == F)
}###################	END OF FUNCTION	###########################

#	Function to plot mean vectors and bootstrapped CI				#
rhoPlot <- function(xx, Exp, IDfactor, ci.col = 'black', ...){
	plot.circular(NA, axes = F, , stack = T, bins = 360)
	for(exp in Exp){
		dt <- subset(xx, Condition == exp)
		# dev.new(width = 3, height = 3)
		# par(mai = c(0,0,0,0))
		for(id in unique(dt[,IDfactor])){
			mv <- dt$Heading[dt[,IDfactor] == id]
			if(sum(!is.na(mv))>0){

			mns <- mean(mycirc(mv), na.rm = T)
			arrows.circular(mean(circular(mv,  units = 'degrees', zero = pi/2, rotation = 'clock'), na.rm = T),
			 length = 0, shrink = rho.circular(circular(mv,  units = 'degrees', zero = pi/2, rotation = 'clock'), na.rm = T), 
			  lwd = 5, lend = 'butt',...)
			 }#if(!is.na(mv))
		}#for(id in unique(subset(xx, Experiment == Experiment)[IDfactor]))
		for(id in unique(dt[,IDfactor])){
			mv <- dt$Heading[dt[,IDfactor] == id]
			if(sum(!is.na(mv))>0){
				
				mns <- mean(mycirc(mv), na.rm = T)
				bts <- mle.vonmises.bootstrap.ci(mycirc(mv), alpha = 0.05, 
					bias = T, reps = 10^4)$kappa.ci
				arrows(A1(bts[1])*sin(as.numeric(mns)*pi/180),
				A1(bts[1])*cos(as.numeric(mns)*pi/180),
				A1(bts[2])*sin(as.numeric(mns)*pi/180), 
				A1(bts[2])*cos(as.numeric(mns)*pi/180), 
				code = 3, angle = 90, length = 0.01, 
				col = ci.col, lwd = 2.5)
			 }#if(!is.na(mv))
		}#for(id in unique(subset(xx, Experiment == Experiment)[IDfactor]))
	}#	for(exp in Exp){
###################	END OF FUNCTION	###########################	
}
MeanVec <- function(x){rho.circular(mycirc(x), na.rm = T)}

#####################################################################
#	Plot Settings											#
#####################################################################

#axis size multiplication factor
mf <- 2.559055# plots are 65mm tall# 2.375 #plots are 60.325 mm tall
pl <- T#F