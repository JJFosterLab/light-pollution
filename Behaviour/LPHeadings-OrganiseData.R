#################################################################
#Some experiments from different nights can be combined
#These amalgamations will be called "Conditions"
#Some data may need to be excluded
#Beetles within each Condition may need renaming
#Add variables decribing direct and indirect light pollution
#################################################################

rm(list = ls())
graphics.off()
#################################################################
#	Useful Functions											#
#################################################################
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
#plotting functions excluded, plotting will be performed in another script

#####################################################################
#	Plot Settings											#
#####################################################################

#axis size multiplication factor
mf <- 2.559055# plots are 65mm tall# 2.375 #plots are 60.325 mm tall
pl <- T#F

#####################################################################
#	Install and Load Packages										#
#####################################################################

Instalload(c('circular', 'beeswarm'))
lp.all <- read.table(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LP Headings all 20200505', '.txt') , header = T, sep  = '\t')

lp.all$Beetle <- as.factor(lp.all$Beetle)
#####################################################################
#	Add details to describe conditions								#
#####################################################################
#give experiments a shorter name
lp.all$ShortName <- factor(lp.all$Experiment)
levels(lp.all$Experiment)#these are the levels we have to work with
		# [1] "LPexper.20171116.LPexper.lowMWayStonehenge"        
		 # [2] "LPexper.20171121.lowMWayThornwood"                 
		 # [3] "LPexper.20171121.lowMWayWits"                      
		 # [4] "LPexper.20171121.LPexper.OvercastThornwood"        
		 # [5] "LPexper.20171121.LPexper.OvercastWits"             
		 # [6] "LPexper.20171121.LPexper.VerandaThornwood"         
		 # [7] "LPexper.20181122.Fullmoon.Thornwood"               
		 # [8] "LPexper.20181122.WaxGibbous.Wits"                  
		 # [9] "LPexper.20181123.Fullmoon.overcast.Thornwood"      
		# [10] "LPexper.20181123.FullMoonHazy.Wits"                
		# [11] "LPexper.20181127.LowMWay.Thornwood"                
		# [12] "LPexper.20181128.LowMWay.Wits"                     
		# [13] "LPexper.20181128.LowMWayandMushrooms.Wits"         
		# [14] "LPexper.20181128.LowMWayPartlyCloudy.Thornwood"    
		# [15] "LPexper.20191117.SecurityLightBergsigOFF"          
		# [16] "LPexper.20191117.SecurityLightBergsigON"           
		# [17] "LPexper.20191120.OvercastMarcusWallsBergsig"       
		# [18] "LPexper.20191122.ClearMarcusWallsWits"             
		# [19] "LPexper.20191122.CloudsMarcusWallsWits"            
		# [20] "LPexper.20191123.ClearGOGGLESMarcusWallsWits"      
		# [21] "LPexper.20191123.ClearMarcusWallsWits"             
		# [22] "LPexper.20191123.CloudsMarcusWallsWits"            
		# [23] "LPexper.20191125.ClearMarcusWallsBergsig"          
		# [24] "LPexper.20200214.MostlyClearMarcusWalls.WitsCLEAR" 
		# [25] "LPexper.20200214.MostlyClearMarcusWalls.WitsCLOUDY"
		# [26] "LPexper.20200215.ClearMarcusWalls.Wits"            
		# [27] "LPexper.20200216.WhiteTackMarcusWalls.Wits"           
levels(lp.all$ShortName)[1:length(levels(lp.all$ShortName))] <- 
				c('noMW\nRural', #1
				'Stars\nRural', #2
				'Stars\nUrban', #3
				'Overcast\nRural', #4
				'Overcast\nUrban', #5
				'"Veranda"\nRural', #6
				'Moon\nRural', #7
				'Moon\nUrban', #8
				'Hazy Moon\nRural', #9
				'Hazy Moon\n Urban',#10
				'Stars 2\nRural', #11
				'Stars 2\nUrban', #12
				'"Mushrooms"\nUrban', #13
				'Cloudy\nRural', #14
				'Security Light OFF\nRural',  #15
				'Security Light ON\nRural',  #16
				'Overcast (Walls)\nRural', #17
				'Some Stars (Walls)\nUrban', #18
				'Cloudy (Walls)\nUrban', #19
				'Stars GOGGLES (Walls)\nUrban', #20
				'Stars 2 (Walls)\nUrban', #21
				'Cloudy 2 (Walls)\nUrban', #22
				'Stars (Walls)\nRural', #23
				'Stars 3 (Walls)\nUrban', #24
				'Cloudy 3 (Walls)\nUrban', #25
				'Stars (Walls)\nUrban', #26
				'White Tack (Walls)\nUrban' #27
				)[1:length(levels(lp.all$ShortName))]#keeps levels in the right order
				
#by default, there was no light pollution
lp.all$Indirect <- FALSE
#in urban environments, there is always indirect light pollution
lp.all$Indirect[grepl('Urban', lp.all$ShortName)] <- TRUE
#by default, there was no light pollution
lp.all$Direct <- FALSE
#in urban environements, there was usually direct light pollution
lp.all$Direct[grepl('Urban', lp.all$ShortName)] <- TRUE
#except when there were walls blocking the lower 60° of elevation
lp.all$Direct[grepl('Walls', lp.all$ShortName)] <- FALSE
#in two rural experiments, there were security lights
lp.all$Direct[lp.all$ShortName == '"Veranda"\nRural'] <- TRUE
lp.all$Direct[lp.all$ShortName == 'Security Light ON\nRural'] <- TRUE
levels(lp.all$ShortName)
#ShortName is good for plotting, but Condition can be simpler		
lp.all$Condition <- gsub('\n', '', lp.all$ShortName)		
lp.all$Condition <- gsub(' ', '', lp.all$Condition)		
lp.all$Condition <- gsub('\\(', '', lp.all$Condition)
lp.all$Condition <- gsub('\\)', '', lp.all$Condition)
lp.all$Condition <- gsub('\\"', '', lp.all$Condition)
unique(lp.all$Condition)			
#combine the "eyes covered" conditions
lp.all$Condition[ 
	lp.all$ShortName=='Stars GOGGLES (Walls)\nUrban'| lp.all$ShortName=='White Tack (Walls)\nUrban' ] <- 'EyesCoveredWallsUrban'
#combine the urban wall experiments (that I trust)
lp.all$Condition[ 
	lp.all$ShortName=='Stars 3 (Walls)\nUrban'] <- 'StarsWallsUrban'
lp.all$Condition[ 
	lp.all$ShortName=='Cloudy 3 (Walls)\nUrban'] <- 'CloudsWallsUrban'
#and the ones I'm less sure about
lp.all$Condition[ 
	lp.all$ShortName=='Some Stars (Walls)\nUrban' |lp.all$ShortName=='Stars 2 (Walls)\nUrban'] <- 'SomeStarsWallsUrban'
lp.all$Condition[ 
	lp.all$ShortName=='Cloudy (Walls)\nUrban' |lp.all$ShortName=='Cloudy 2 (Walls)\nUrban'] <- 'SomeCloudsWallsUrban'
#combine the clear sky experiments (including ones Marie doesn't trust)
lp.all$Condition[#Marie is suspicious of heading bias in StarsRural
	lp.all$ShortName=='Stars 2\nRural'] <- 'StarsRural'
lp.all$Condition[
	lp.all$ShortName=='Stars 2\nUrban'] <- 'StarsUrban'	
#now organise and order Condition	
lp.all$Condition <- as.factor(lp.all$Condition)
lp.all$Condition <- relevel(lp.all$Condition, 'EyesCoveredWallsUrban')
levels(lp.all$Condition)	
#some conditions now have beetle numbers twice, and most beetle numbers are shared between conditions without truly being the same beetle
#only in the "SecurityLight" experiment were beetles rolled twice
lp.all$Beetle <- as.numeric(lp.all$Beetle)
lp.all$Beetle[#find all conditions except SecurityLight
	!(grepl('SecurityLight', lp.all$Condition))] <- 
	sort(#beetle numbers will be ordered	
	rep(1:(length(lp.all$Heading[#give each beetle a unique number
		!(grepl('SecurityLight', lp.all$Condition))]) / 
		10), 10) )#knowing each beetle rolled 10 times, we can divide the number or headings by 10 to get the number of beetles, and repeat each beetle's number 10 times
lp.all$Beetle[#find all conditions with SecurityLight
	grepl('SecurityLight', lp.all$Condition)] <- 
	as.numeric(lp.all$Beetle[#find all conditions with SecurityLight
	grepl('SecurityLight', lp.all$Condition)]) + 
	(length(lp.all$Heading[#add the highest number used
		!(grepl('SecurityLight', lp.all$Condition))]) / 
		10) #knowing each beetle rolled 10 times, we can divide the number or headings by 10 to get the total number of beetles, and add it to the security light beetle numbers
lp.all$Beetle <- as.factor(lp.all$Beetle)
lp.all[c(-1,-6,-7)]
#save this dataset
write.table(lp.all, file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPallOrganised20200505', '.txt'))
#load the saved dataset
lp.all <- read.table(file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPallOrganised20200505', '.txt'))

#N.B. still missing data from the first beetle of StarsRural
subset(lp.all, Beetle == 41)#UPDATE not anymore

#####################################################################
#	Calculate Mean Vectors for Each Beetle							#
#####################################################################
MeanVec <- function(x){rho.circular(mycirc(x), na.rm = T)}
lp.all.rho <- with(lp.all,aggregate(Heading, by = list(ShortName = ShortName, Condition = Condition, Beetle = Beetle, Direct = Direct, Indirect = Indirect), MeanVec))		
names(lp.all.rho)[length(names(lp.all.rho))] <- 'rho'
head(lp.all.rho)
write.table(lp.all.rho, file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPallRho20200505', '.txt'))
#load the saved dataset
lp.all.rho <- read.table(file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPallRho20200505', '.txt'))

MeanAng <- function(x){mean.circular(mycirc(x), na.rm = T)}
lp.all.mu <- with(lp.all,aggregate(Heading, by = list(ShortName = ShortName, Condition = Condition, Beetle = Beetle, Direct = Direct, Indirect = Indirect), MeanAng))
names(lp.all.mu)[length(names(lp.all.mu))] <- 'mu'
head(lp.all.mu)

write.table(lp.all.mu, file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPallMu20200528', '.txt'))
#load the saved dataset
lp.all.mu <- read.table(file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPallMu20200528', '.txt'))

#####################################################################
#	Remove Conditions that Won't be Considered						#
#####################################################################
# We don't have any interesting findings from:
# CloudsWallsUrban	#Could be interesting, not enough data.
# HazyMoonRural	#Good data, no interesting hypotheses.
# HazyMoonUrban	#Good data, no interesting hypotheses.	
# noMWrural		#Good data, no interesting hypotheses.
# Cloudyrural		#Good data, no interesting hypotheses.
# SomeCloudsWallsUrban	#Strange data, would rather not use.
# SomeStarsWallsUrban	#Strange data, would rather not use.
lp.data <- lp.all[#find and exclude these conditions
	!with(lp.all, 
	grepl('CloudsWalls', Condition)	|
	grepl('Cloudy', Condition)	|
	grepl('Hazy', Condition)	|
	grepl('Some', Condition)	|
	grepl('noMW', Condition)	),	] 
	
#make a "sky" variable for sky types
# lp.data$sky <- NA
# lp.data$sky[with(lp.data, 
	# grepl('Star', Condition))] <- 'starry'
# lp.data$sky[with(lp.data, 
	# grepl('Overcast', Condition))] <- 'overcast'
# lp.data$sky[with(lp.data, 
	# grepl('Moon', Condition))] <- 'moonlit'
# lp.data$sky[with(lp.data, 
	# grepl('Moon', Condition))] <- 'moonlit'
#save this dataset
write.table(lp.data, file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPdataOrganised20200505', '.txt'))

lp.data.rho <- lp.all.rho[#find and exclude these conditions
	!with(lp.all.rho, 
	grepl('CloudsWalls', Condition)	|
	grepl('Cloudy', Condition)	|
	grepl('Hazy', Condition)	|
	grepl('Some', Condition)	|
	grepl('noMW', Condition)	),	] 
#save this dataset
write.table(lp.data.rho, file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPdataRho20200505', '.txt'))


lp.data.mu <- lp.all.mu[#find and exclude these conditions
	!with(lp.all.mu, 
	grepl('CloudsWalls', Condition)	|
	grepl('Cloudy', Condition)	|
	grepl('Hazy', Condition)	|
	grepl('Some', Condition)	|
	grepl('noMW', Condition)	),	] 
#save this dataset
write.table(lp.data.mu, file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPdataRho20200528', '.txt'))

#####################################################################
#	Plot Summary													#
#####################################################################




#####################################################################
#	Experiments with walls											#
#####################################################################

#	select all the "walls" experiments
lp.walls.rho <- subset(lp.all.rho, Experiment == 'LPexper.20200214.MostlyClearMarcusWalls.WitsCLEAR' | Experiment == 'LPexper.20200215.ClearMarcusWalls.Wits' | Experiment == 'LPexper.20191125.ClearMarcusWallsBergsig'| Experiment == 'LPexper.20200216.WhiteTackMarcusWalls.Wits' | Experiment == 'LPexper.20191123.ClearGOGGLESMarcusWallsWits' |  Experiment == 'LPexper.20191120.OvercastMarcusWallsBergsig')
unique(lp.walls.rho$ShortName)
lp.walls.rho$ShortName <- as.character(lp.walls.rho$ShortName)
#Starry skies are the same on the 14th and 15th
lp.walls.rho$ShortName[ 
	lp.walls.rho$ShortName=='Stars 3 (Walls)\nUrban'] <- 'Stars (Walls)\nUrban'
#Goggles and white tack have the same effect
lp.walls.rho$ShortName[ 
	lp.walls.rho$ShortName=='Stars GOGGLES (Walls)\nUrban'| lp.walls.rho$ShortName=='White Tack (Walls)\nUrban' ] <- 'Eyes Covered (Walls)\nUrban'
unique(lp.walls.rho$ShortName)
#reorder and remove the "(Walls)" part of the name since they all have that in common
lp.walls.rho$ShortName <- factor(lp.walls.rho$ShortName, unique(lp.walls.rho$ShortName)[c(1,3,4,2)])
levels(lp.walls.rho$ShortName) <- gsub(" \\(Walls\\)", '', levels(lp.walls.rho$ShortName))
lp.walls.rho[,-1]
# remove the confusing escape characters for modelling
lp.walls.rho$condition <- gsub('\n', '', lp.walls.rho$ShortName)





#	plot the "walls" experiments
dev.new(width =5); par(mai = c(0.8,0.8, 0, 0))
#plot on the appropriate scale transformed
boxplot(rho~ShortName, data = lp.walls.rho, ylim =(c(0,1)), xlim = c(0.5, 4.5),
cex = 0.5, outline = F, border = rgb(0,0,0,1),
		pars = list(boxwex = 0.3/2, staplewex = 0.5, outwex = 0.5),
#colour code by sex, turn axes off
	axes = F, col = rgb(0,0,0,0.1),
	ylab = 'Mean Vector Length', xlab = '')
#log scale axes with original data values	
	axis(2)#, at =log10( pretty( range(soil.dt2$CO2_Flux) *2 )/2 ), labels = pretty( range(soil.dt2$CO2_Flux) *2 )/2)
	mtext(levels(lp.walls.rho$ShortName), side = 1, at = 1:length(levels(lp.walls.rho$ShortName)), line = 0.5 )
	abline(h =c(0,1), lwd = 0.25)
	# abline(h = 0.36503, col = 'green', lty = 2)
polygon(c(0,5,5,0), c(0,0, 0.36503,0.36503), col = rgb(1,0,0,0.05), border = NA)
legend('bottom', inset = 0.36503, legend = '   Rayleigh test     \n       p<0.05     \n            \n       p>0.05     ', cex = 0.65, bty = 'n')
	
beeswarm(rho~ShortName, data = lp.walls.rho,
		pch = 20,  cex = 1.8,#method = 'center',
		pwcol = c('gray50','darkblue','orange2', 'salmon')[lp.walls.rho$ShortName] , 
		add = T, axes = F)
		
PDFsave(Directory = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/'), Experiment = 'LP', PlotName = "WallsAll-Beeswarm")




#	Non-parametric comparisons										#

#	Location change (median), Kruskal-Wallis & Wilcoxon
kruskal.test(rho~condition, data = lp.walls.rho)
# Kruskal-Wallis chi-squared = 6.8847, df = 3, p-value = 0.07567
#slightly annoying, I shouldn't really do Wilcoxon post-hocs like this
#whole dataset for reference
kruskal.test(rho~ShortName, data = lp.all.rho)
#Kruskal-Wallis chi-squared = 124.28, df = 24, p-value = 1.665e-15
#can just quote this instead
#need to fix conditions in the dataset as a whole first though...
#..can do that later
levels(lp.walls.rho$ShortName)#testing will be in this order
# "Overcast\nRural" "Stars\nRural" "Stars\nUrban" "Eyes Covered\nUrban"
#1st level is tested against the 2nd (i.e. 1st - 2nd)
 # OvercastRural - Eyes CoveredUrban
wilcox.test(rho ~ ShortName, data = subset(lp.walls.rho, condition == 'Eyes CoveredUrban' | condition == 'OvercastRural'), alternative = 'less')
		# W = 36, p-value = 0.5566

# StarsRural Eyes - CoveredUrban
wilcox.test(rho ~ ShortName, data = subset(lp.walls.rho, condition == 'Eyes CoveredUrban' | condition == 'StarsRural'), alternative = 'greater')
		# W = 55, p-value = 0.02766

# StarsUrban - Eyes CoveredUrban
wilcox.test(rho ~ ShortName, data = subset(lp.walls.rho, condition == 'Eyes CoveredUrban' | condition == 'StarsUrban'), alternative = 'greater')
		# W = 58, p-value = 0.1753

# StarsRural - OvercastRural
wilcox.test(rho ~ ShortName, data = subset(lp.walls.rho, condition == 'OvercastRural' | condition == 'StarsRural'), alternative = 'less')
		# W = 23, p-value = 0.02163

# StarsUrban - OvercastRural
wilcox.test(rho ~ ShortName, data = subset(lp.walls.rho, condition == 'OvercastRural' | condition == 'StarsUrban'), alternative = 'less')
		# W = 46, p-value = 0.1284

# StarsRural - StarsUrban 
wilcox.test(rho ~ ShortName, data = subset(lp.walls.rho, condition == 'StarsRural' | condition == 'StarsUrban'), alternative = 'greater')
		# W = 94, p-value = 0.03843
 
     
    
        
     
#	Scale change (variance), Brown-Forsythe tests
Instalload(c('onewaytests'))
bf.test(rho ~ condition, data = lp.walls.rho)
		  # statistic  : 3.612992 
		  # num df     : 3 
		  # denom df   : 23.54427 
		  # p.value    : 0.02804121 
  
# Eyes CoveredUrban - OvercastRural
bf.test(rho ~ ShortName, data = subset(lp.walls.rho, condition == 'Eyes CoveredUrban' | condition == 'OvercastRural'))
		  # statistic  : 0.001463094 
		  # num df     : 1 
		  # denom df   : 9.417792 
		  # p.value    : 0.9702866 

# Eyes CoveredUrban - StarsRural  
bf.test(rho ~ ShortName, data = subset(lp.walls.rho, condition == 'Eyes CoveredUrban' | condition == 'StarsRural'))
		  # statistic  : 5.217001 
		  # num df     : 1 
		  # denom df   : 14.77199 
		  # p.value    : 0.03760092 

# Eyes CoveredUrban - StarsUrban
bf.test(rho ~ ShortName, data = subset(lp.walls.rho, condition == 'Eyes CoveredUrban' | condition == 'StarsUrban'))
		  # statistic  : 0.8949162 
		  # num df     : 1 
		  # denom df   : 12.47293 
		  # p.value    : 0.3621155 

# OvercastRural - StarsRural
bf.test(rho ~ ShortName, data = subset(lp.walls.rho, condition == 'OvercastRural' | condition == 'StarsRural'))
		  # statistic  : 7.034754 
		  # num df     : 1 
		  # denom df   : 11.59236 
		  # p.value    : 0.02162861
   
# OvercastRural - StarsUrban
bf.test(rho ~ ShortName, data = subset(lp.walls.rho, condition == 'OvercastRural' | condition == 'StarsUrban'))
		  # statistic  : 1.787053 
		  # num df     : 1 
		  # denom df   : 20.4126 
		  # p.value    : 0.1959948 

# StarsRural - StarsUrban
bf.test(rho ~ ShortName, data = subset(lp.walls.rho, condition == 'StarsRural' | condition == 'StarsUrban'))
		  # statistic  : 3.023938 
		  # num df     : 1 
		  # denom df   : 13.92372 
		  # p.value    : 0.1040941 





#	Model with Beta distribution									#
Instalload(c('betareg','lmtest','emmeans','multcomp'))#for the sake of argument, let's try a Beta model

breg.walls.1<- betareg(rho~condition, data = lp.walls.rho)
#B-F tests suggest variance also changes as a function of condition
breg.walls.2<- betareg(rho~condition|condition, data = lp.walls.rho)
breg.walls.0 <- betareg(rho~1, data = lp.walls.rho)
lrtest(breg.walls.1, breg.walls.0)
		# Model 1: rho ~ condition
		# Model 2: rho ~ 1
		  # #Df LogLik Df  Chisq Pr(>Chisq)  
		# 1   5 14.524                       
		# 2   2 10.313 -3 8.4222    0.03805 *
lrtest(breg.walls.2, breg.walls.1)		
AIC(breg.walls.2)#totally worth it
AIC(breg.walls.1)
AIC(breg.walls.0)

summary(breg.walls.2)

#marginal means
emm.breg.walls <- emmeans(breg.walls.2, list(pairwise ~ condition))
test(emm.breg.walls)
# $`emmeans of condition`
 # condition         emmean     SE  df z.ratio p.value
 # Eyes CoveredUrban  0.261 0.0565 Inf 4.616   <.0001 
 # OvercastRural      0.259 0.0310 Inf 8.348   <.0001 
 # StarsRural         0.500 0.0778 Inf 6.428   <.0001 
 # StarsUrban         0.341 0.0442 Inf 7.703   <.0001 


# $`pairwise differences of condition`
 # contrast                          estimate     SE  df z.ratio p.value
 # Eyes CoveredUrban - OvercastRural   0.0022 0.0645 Inf  0.034  1.0000 
 # Eyes CoveredUrban - StarsRural     -0.2395 0.0962 Inf -2.490  0.0615 
 # Eyes CoveredUrban - StarsUrban     -0.0799 0.0718 Inf -1.113  0.6817 
 # OvercastRural - StarsRural         -0.2417 0.0838 Inf -2.885  0.0205 
 # OvercastRural - StarsUrban         -0.0821 0.0540 Inf -1.519  0.4258 
 # StarsRural - StarsUrban             0.1596 0.0895 Inf  1.783  0.2815
 
 #marginal precision
 phi.emm.breg.walls <- emmeans(breg.walls.2, list(pairwise ~ condition), mode = 'phi.link')
 test(phi.emm.breg.walls)
 # $`emmeans of condition`
 # condition         emmean    SE  df z.ratio p.value
 # Eyes CoveredUrban   2.02 0.514 Inf 3.934   0.0001 
 # OvercastRural       2.94 0.439 Inf 6.702   <.0001 
 # StarsRural          1.07 0.390 Inf 2.739   0.0062 
 # StarsUrban          2.05 0.372 Inf 5.499   <.0001 

# Results are given on the log (not the response) scale. 

# $`pairwise differences of condition`
 # contrast                          estimate    SE  df z.ratio p.value
 # Eyes CoveredUrban - OvercastRural  -0.9208 0.675 Inf -1.363  0.5225 
 # Eyes CoveredUrban - StarsRural      0.9508 0.645 Inf  1.474  0.4533 
 # Eyes CoveredUrban - StarsUrban     -0.0272 0.634 Inf -0.043  1.0000 
 # OvercastRural - StarsRural          1.8716 0.587 Inf  3.186  0.0079 
 # OvercastRural - StarsUrban          0.8937 0.575 Inf  1.553  0.4058 
 # StarsRural - StarsUrban            -0.9780 0.539 Inf -1.813  0.2672 

 
 betapred.lcl <- summary(emm.breg.walls)$`emmeans of condition`['asymp.LCL']
 betapred.ucl <- summary(emm.breg.walls)$`emmeans of condition`['asymp.UCL']

betapred <- predict(breg.walls, newdata = data.frame(condition = unique(lp.walls.rho$condition)), type = c("response"))
betapred.q <- predict(breg.walls, newdata = data.frame(condition = unique(lp.walls.rho$condition)), type = "quantile", at = c(0.025, 0.975))

dev.new(width =5); par(mai = c(0.8,0.8, 0, 0), lend = 'butt')
#plot on the appropriate scale transformed	
beeswarm(rho~ShortName, data = lp.walls.rho,
		pch = 20,  cex = 1.8, #method = 'center',
		pwcol = c('gray50','darkblue','orange2', 'salmon')[lp.walls.rho$ShortName] , 
			axes = F, ylim =(c(0,1)), xlim = c(0.5, 4.5),
	ylab = 'Mean Vector Length', xlab = '')
#log scale axes with original data values	
	axis(2)#, at =log10( pretty( range(soil.dt2$CO2_Flux) *2 )/2 ), labels = pretty( range(soil.dt2$CO2_Flux) *2 )/2)
	mtext(levels(lp.walls.rho$ShortName), side = 1, at = 1:length(levels(lp.walls.rho$ShortName)), line = 0.5 )
	abline(h =c(0,1), lwd = 0.25)
	# abline(h = 0.36503, col = 'green', lty = 2)
	polygon(c(0,5,5,0), c(0,0, 0.36503,0.36503), col = rgb(1,0,0,0.05), border = NA)
	# legend('bottom', inset = 0.36503, legend = 'Rayleigh test     \n    p<0.05     \n         \n    p>0.05     ', cex = 0.65, bty = 'n')
	
points(1:length(unique(lp.walls.rho$condition)), betapred[c(2,3,4,1)], pch = 10, col = 'black', cex = 2, lwd = 3 )
arrows(1:length(unique(lp.walls.rho$condition)),
		betapred.q[c(2,3,4,1),'q_0.025'],
		1:length(unique(lp.walls.rho$condition)),
		betapred.q[c(2,3,4,1),'q_0.975'],
		angle = 90, length = 0, code = 3, lwd = 0.5,
		col = 'black'
		)
arrows(1:length(unique(lp.walls.rho$condition)),
		betapred.lcl[c(2,3,4,1),1],
		1:length(unique(lp.walls.rho$condition)),
		betapred.ucl[c(2,3,4,1),1],
		angle = 90, length = 0.05, code = 3, lwd = 4, lend = 'butt',
		col = 'black'
		)
		
legend('topright',
		legend = c('Fitted mean and 95% CI', 'Prediction interval'),
		pch = c(10,NA),
		lwd = c(4, 0.5),
		lty = c(1,1),
		 col = c('black'), inset = .05, cex = 0.75, pt.cex = 2)

		
PDFsave(Directory = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/'), Experiment = 'LP', PlotName = "WallsAll-BetaregModel")







#####################################################################
#	Other Experiments												#
#####################################################################

minmax <- cbind(subset(moon, ShortName == 'Moon\nUrban')$rho, subset(overcast, ShortName == 'Overcast\nRural')$rho)
bf.test(rho ~ ShortName, data = minmax)

moon <- subset(lp.all.rho, ShortName %in% c('Moon\nRural', 'Moon\nUrban'))
wilcox.test(subset(moon, ShortName == 'Moon\nUrban')$rho, subset(moon, ShortName == 'Moon\nRural')$rho, alternative = 'greater')
#W = 66, p-value = 0.1237
bf.test(rho ~ ShortName, data = moon)
  # statistic  : 1.732272 
  # num df     : 1 
  # denom df   : 10.2582 
  # p.value    : 0.2167715 

overcast <- subset(lp.all.rho, ShortName %in% c('Overcast\nRural', 'Overcast\nUrban'))
boxplot(rho~ShortName, data = overcast)
wilcox.test(subset(overcast, ShortName == 'Overcast\nUrban')$rho, subset(overcast, ShortName == 'Overcast\nRural')$rho, alternative = 'greater')
# W = 89, p-value = 0.001045
bf.test(rho ~ ShortName, data = overcast)
  # statistic  : 16.87455 
  # num df     : 1 
  # denom df   : 12.52775 
  # p.value    : 0.001330328 
  
stars <- subset(lp.all.rho, ShortName %in% c('Stars\nRural', 'Stars\nUrban'))
boxplot(rho~ShortName, data = stars)
wilcox.test(subset(stars, ShortName == 'Stars\nUrban')$rho, subset(stars, ShortName == 'Stars\nRural')$rho, alternative = 'greater')
#W = 59, p-value = 0.1388
bf.test(rho ~ ShortName, data = stars)
  # statistic  : 0.5009182 
  # num df     : 1 
  # denom df   : 16.88239 
  # p.value    : 0.4887536 

lamps <- subset(lp.all.rho, ShortName %in% c('"Veranda"\nUrban', '"Mushrooms"\nUrban'))
boxplot(rho~ShortName, data = lamps)
wilcox.test(subset(lamps, ShortName == '"Mushrooms"\nUrban')$rho, subset(lamps, ShortName == '"Veranda"\nUrban')$rho, alternative = 'less')
#W = 28, p-value = 0.05256
bf.test(rho ~ ShortName, data = lamps)
  # statistic  : 1.61044 
  # num df     : 1 
  # denom df   : 15.54667 
  # p.value    : 0.2231008 

lamps1 <- subset(lp.all.rho, ShortName %in% c('Stars\nRural','"Veranda"\nUrban'))
boxplot(rho~ShortName, data = lamps1)
wilcox.test(subset(lamps1, ShortName == 'Stars\nRural')$rho, subset(lamps1, ShortName == '"Veranda"\nUrban')$rho, alternative = 'less')
#W = 31, p-value = 0.1388
bf.test(rho ~ ShortName, data = lamps1)
  # statistic  : 0.284861 
  # num df     : 1 
  # denom df   : 16.98931 
  # p.value    : 0.6004468 

lamps2 <- subset(lp.all.rho, ShortName %in% c('Stars\nRural','"Mushrooms"\nUrban'))
boxplot(rho~ShortName, data = lamps2)
wilcox.test(subset(lamps2, ShortName == 'Stars\nRural')$rho, subset(lamps2, ShortName == '"Mushrooms"\nUrban')$rho, alternative = 'greater')
#W = 56, p-value = 0.2001
bf.test(rho ~ ShortName, data = lamps2)
  # statistic  : 0.7957844 
  # num df     : 1 
  # denom df   : 14.57924 
  # p.value    : 0.3868436 