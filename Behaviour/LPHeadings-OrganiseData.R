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