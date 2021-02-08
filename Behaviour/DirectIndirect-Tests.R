rm(list = ls())
graphics.off()
#################################################################
#	Useful Functions											#
#################################################################
source(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LP_PlotFunctions', '.R'))

#####################################################################
#	Install and Load Packages										#
#####################################################################

Instalload(c('circular', 'beeswarm', 'onewaytests', 'muStat'))
#	Heading data that will be used in the paper
lp.data <- read.table(file = paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPdataOrganised20200505', '.txt'))
head(lp.data)
lp.data$Beetle <- as.factor(lp.data$Beetle)
#reorganise conditions
lp.data$Condition <- relevel(lp.data$Condition, 'EyesCoveredWallsUrban')
lp.data$Condition <- factor(lp.data$Condition, levels(lp.data$Condition)[c(1,2:3,10:11,5:6,14,4,8:9,7,12,13)])	
cbind(levels(lp.data$Condition))#is this what you want?

#	Mean vector length data that will be used in the paper
lp.data.rho <- read.table(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPdataRho20200505', '.txt') , header = T)#, sep  = '\t')
head(lp.data.rho)
lp.data.rho$Beetle <- as.factor(lp.data.rho$Beetle)
#reorganise conditions
lp.data.rho$Condition <- relevel(lp.data.rho$Condition, 'EyesCoveredWallsUrban')
lp.data.rho$Condition <- factor(lp.data.rho$Condition, levels(lp.data.rho$Condition)[c(1,2:3,10:11,5:6,14,4,8:9,7,12,13)])	
cbind(levels(lp.data.rho$Condition))#is this what you want?	

#	Mean vector length data that will be used in the paper
lp.data.mu <- read.table(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPdataRho20200528', '.txt') , header = T)#, sep  = '\t')
head(lp.data.mu)
lp.data.mu$Beetle <- as.factor(lp.data.mu$Beetle)
lp.data.mu$Condition <- as.factor(lp.data.mu$Condition)
#reorganise conditions
lp.data.mu$Condition <- relevel(lp.data.mu$Condition, 'EyesCoveredWallsUrban')
lp.data.mu$Condition <- factor(lp.data.mu$Condition, levels(lp.data.mu$Condition)[c(1,2:3,10:11,5:6,14,4,8:9,7,12,13)])	
cbind(levels(lp.data.mu$Condition))#is this what you want?	

#####################################################################
#	Focus on Urban verus Rural										#
#####################################################################

#Choose data that appears in Figure
DirectIndirect <- subset(lp.data.rho,
							( grepl('SecurityLight', Condition) |
							grepl('Veranda', Condition) | 
							grepl('Mushrooms', Condition)|
							grepl('StarsWallsRural', Condition)|
							grepl('StarsWallsUrban', Condition)|
							grepl('OvercastWallsRural', Condition)|
						grepl('EyesCoveredWallsUrban', Condition) ) )
#Make sky an accurate reflection of sky consitions
DirectIndirect$sky <- with(DirectIndirect, gsub('Rural','',Condition))
DirectIndirect$sky <- with(DirectIndirect, gsub('Urban','',sky))
DirectIndirect$sky <- with(DirectIndirect, gsub('Walls','',sky))
DirectIndirect$sky <- with(DirectIndirect, gsub('Mushrooms','Stars',sky))
DirectIndirect$sky <- with(DirectIndirect, gsub('SecurityLight','Stars',sky))
DirectIndirect$sky <- with(DirectIndirect, gsub('Veranda','Stars',sky))
DirectIndirect$sky <- with(DirectIndirect, gsub('ON','',sky))
DirectIndirect$sky <- with(DirectIndirect, gsub('OFF','',sky))
DirectIndirect$sky <- factor(DirectIndirect$sky)
DirectIndirect$Condition <- factor(DirectIndirect$Condition)#, 
									# c('SecurityLightOFF',
									# 'SecurityLightON',
									# 'VerandaRural',
									# 'MushroomsUrban'))
									
head(DirectIndirect)

#####################################################################
#	Location Tests													#
#####################################################################
kt.DirectIndirect <- kruskal.test(rho~Condition, data = DirectIndirect)
# Kruskal-Wallis chi-squared = 53.969, df = 7, p-value = 2.386e-09

#	Post-hoc tests													#
#remember condition order is important here
 cbind(levels(DirectIndirect$Condition))
 
#	compare to reference											#

#	Direct Light Pollution
 # Floodlight ON - Floodlight OFF 
dt.FloodlightON.OFF.Direct <- subset(DirectIndirect,
								Condition == 'SecurityLightONRural' | 
								Condition == 'SecurityLightOFFRural')
dt.FloodlightON.OFF.Direct$Condition <- factor(
								dt.FloodlightON.OFF.Direct$Condition,
			levels = c('SecurityLightONRural','SecurityLightOFFRural'))
								#ON is alternative
wc.FloodlightON.OFF.Direct <- wilcox.test(rho ~ Condition,
									data = dt.FloodlightON.OFF.Direct,
									alternative = 'greater') 
# print(FloodlightON.OFF.Direct)
		# W = 75, p-value = 0.03151
#Use paired!
dt.ON <- subset(dt.FloodlightON.OFF.Direct, grepl('ON', Condition))
dt.OFF <- subset(dt.FloodlightON.OFF.Direct, grepl('OFF', Condition))
dt.ON <- dt.ON[order(dt.ON$Beetle),]
dt.OFF <- dt.OFF[order(dt.OFF$Beetle),]
#Wilcoxon signed-ranks
wc.sr.FloodlightON.OFF.Direct <- wilcox.test(dt.ON$rho,dt.OFF$rho,
											paired =T,
											alternative = 'greater')
#V = 45, p-value = 0.04199

 # VerandaRural - Floodlight OFF 
dt.Veranda.OFF.Direct <- subset(DirectIndirect,
								Condition == 'VerandaRural' |
								Condition == 'SecurityLightOFFRural')
dt.Veranda.OFF.Direct$Condition <- factor(
										dt.Veranda.OFF.Direct$Condition,
					levels = c('VerandaRural','SecurityLightOFFRural'))
							#VerandaRural is alternative					
wc.Veranda.OFF.Direct <- wilcox.test(rho ~ Condition,
										data = dt.Veranda.OFF.Direct,
										 alternative = 'greater')
# print(wc.Veranda.OFF.Direct)
		# W = 63, p-value = 0.1763
		
 # MushroomsUrban - Floodlight OFF 
dt.Mushrooms.OFF.Direct <- subset(DirectIndirect,
								Condition == 'MushroomsUrban' |
								Condition == 'SecurityLightOFFRural')
dt.Mushrooms.OFF.Direct$Condition <- factor(dt.Mushrooms.OFF.Direct$Condition,
				levels = c('MushroomsUrban','SecurityLightOFFRural'))
						#MushroomsUrban is the alternative
wc.Mushrooms.OFF.Direct <- wilcox.test(rho ~ Condition,
										data = dt.Mushrooms.OFF.Direct,
											alternative = 'greater')
 # print(wc.Mushrooms.OFF.Direct)
		# W = 36, p-value = 0.8601 

#	Indirect Light Pollution
 # StarsWallsRural - OvercastWallsRural 
dt.Stars.Cloud.Indirect <- subset(DirectIndirect,
									 Condition == 'StarsWallsRural' |
									 Condition == 'OvercastWallsRural')
dt.Stars.Cloud.Indirect$Condition <- 
							factor(dt.Stars.Cloud.Indirect$Condition,
					levels = c('StarsWallsRural', 'OvercastWallsRural'))
					#OvercastWallsRural is the alternative
wc.Stars.Cloud.Indirect <- wilcox.test(rho ~ Condition,
										data = dt.Stars.Cloud.Indirect,
											alternative = 'greater')
 # print(wc.Stars.Cloud.Indirect)
		# W = 77, p-value = 0.02163
		
 # StarsWallsUrban - OvercastWallsRural 
dt.Urban.Cloud.Indirect <- subset(DirectIndirect,
					Condition == 'StarsWallsUrban' |
					Condition == 'OvercastWallsRural')
dt.Urban.Cloud.Indirect$Condition <- 
							factor(dt.Urban.Cloud.Indirect$Condition,
					levels = c('StarsWallsUrban','OvercastWallsRural'))
wc.Urban.Cloud.Indirect <- wilcox.test(rho ~ Condition,
										data = dt.Urban.Cloud.Indirect,
										alternative = 'greater')
 # print(wc.Urban.Cloud.Indirect)
		# W = 84, p-value = 0.1284 	
		
#Urban Versus Rural
 # StarsRural - StarsUrban 
dt.Stars.Rural.Urban <- subset(DirectIndirect,
										Condition == 'StarsWallsRural' |
										Condition == 'StarsWallsUrban')
dt.Stars.Rural.Urban$Condition <- factor(dt.Stars.Rural.Urban$Condition,
						levels = c('StarsWallsRural','StarsWallsUrban'))
wc.Stars.Rural.Urban <- wilcox.test(rho ~ Condition,
											data = dt.Stars.Rural.Urban,
											alternative = 'greater')
 # print(wc.Stars.Rural.Urban)
		# W = 94, p-value = 0.03843
		
#Eyes covered versus everything
 # SecurityON - EyesCovered 
dt.SecurityON.EyesCovered <- subset(DirectIndirect,
								Condition == 'SecurityLightONRural' |
								Condition == 'EyesCoveredWallsUrban')
dt.SecurityON.EyesCovered$Condition <- factor(dt.SecurityON.EyesCovered$Condition,
			levels = c('SecurityLightONRural','EyesCoveredWallsUrban'))
wc.SecurityON.EyesCovered <- wilcox.test(rho ~ Condition,
									data = dt.SecurityON.EyesCovered,
									alternative = 'greater')
 # print(wc.SecurityON.EyesCovered)
		# W = 70, p-value = 5.142e-05

 # SecurityOFF - EyesCovered 
dt.SecurityOFF.EyesCovered <- subset(DirectIndirect,
								Condition == 'SecurityLightOFFRural' |
								Condition == 'EyesCoveredWallsUrban')
dt.SecurityOFF.EyesCovered$Condition <- factor(dt.SecurityOFF.EyesCovered$Condition,
			levels = c('SecurityLightOFFRural','EyesCoveredWallsUrban'))
wc.SecurityOFF.EyesCovered <- wilcox.test(rho ~ Condition,
									data = dt.SecurityOFF.EyesCovered,
									alternative = 'greater')
 # print(wc.SecurityON.EyesCovered)
		# W = 70, p-value = 5.142e-05
		
 # Veranda - EyesCovered 
dt.Veranda.EyesCovered <- subset(DirectIndirect,
								Condition == 'VerandaRural' |
								Condition == 'EyesCoveredWallsUrban')
dt.Veranda.EyesCovered$Condition <- factor(dt.Veranda.EyesCovered$Condition,
			levels = c('VerandaRural','EyesCoveredWallsUrban'))
wc.Veranda.EyesCovered <- wilcox.test(rho ~ Condition,
									data = dt.Veranda.EyesCovered,
									alternative = 'greater')
 # print(wc.SecurityON.EyesCovered)
		# W = 69, p-value = 0.0001028		
		
 # MushroomsUrban - EyesCovered 
dt.MushroomsUrban.EyesCovered <- subset(DirectIndirect,
								Condition == 'MushroomsUrban' |
								Condition == 'EyesCoveredWallsUrban')
dt.MushroomsUrban.EyesCovered$Condition <- factor(dt.MushroomsUrban.EyesCovered$Condition,
			levels = c('MushroomsUrban','EyesCoveredWallsUrban'))
wc.MushroomsUrban.EyesCovered <- wilcox.test(rho ~ Condition,
								data = dt.MushroomsUrban.EyesCovered,
									alternative = 'greater')
 # print(wc.SecurityON.EyesCovered)
		# W = 65, p-value = 0.000977		
		
 # StarsRural - EyesCovered 
dt.Rural.EyesCovered <- subset(DirectIndirect,
								Condition == 'StarsWallsRural' |
								Condition == 'EyesCoveredWallsUrban')
dt.Rural.EyesCovered$Condition <- factor(dt.Rural.EyesCovered$Condition,
			levels = c('StarsWallsRural','EyesCoveredWallsUrban'))
wc.Rural.EyesCovered <- wilcox.test(rho ~ Condition,
								data = dt.Rural.EyesCovered,
									alternative = 'greater')
 # print(wc.SecurityON.EyesCovered)
		# W = 55, p-value = 0.02766

 # StarsUrban - EyesCovered 
dt.Urban.EyesCovered <- subset(DirectIndirect,
								Condition == 'StarsWallsUrban' |
								Condition == 'EyesCoveredWallsUrban')
dt.Urban.EyesCovered$Condition <- factor(dt.Urban.EyesCovered$Condition,
			levels = c('StarsWallsUrban','EyesCoveredWallsUrban'))
wc.Urban.EyesCovered <- wilcox.test(rho ~ Condition,
									data = dt.Urban.EyesCovered,
										alternative = 'greater')
 # print(wc.SecurityON.EyesCovered)
		# W = 58, p-value = 0.1753
		
 # OvercastRural - EyesCovered 
dt.Clouds.EyesCovered <- subset(DirectIndirect,
								Condition == 'OvercastWallsRural' |
								Condition == 'EyesCoveredWallsUrban')
dt.Clouds.EyesCovered$Condition <- 
								factor(dt.Clouds.EyesCovered$Condition,
			levels = c('OvercastWallsRural','EyesCoveredWallsUrban'))
wc.Clouds.EyesCovered <- wilcox.test(rho ~ Condition,
									data = dt.Clouds.EyesCovered,
										alternative = 'greater')
 # print(wc.SecurityON.EyesCovered)
		# W = 36, p-value = 0.4811
#EDIT 20200712, in the paper this is now:
wc.Clouds.EyesCovered <- wilcox.test(rho ~ Condition,
									data = dt.Clouds.EyesCovered,
										alternative = 'less')
# data:  rho by Condition
# W = 36, p-value = 0.5566

# OvercastRural - EyesCovered 
dt.Veranda.Mushrooms <- subset(DirectIndirect,
								Condition == 'VerandaRural' |
								Condition == 'MushroomsUrban')
dt.Veranda.Mushrooms$Condition <- 
								factor(dt.Veranda.Mushrooms$Condition,
			levels = c('VerandaRural','MushroomsUrban'))
wc.Veranda.Mushrooms <- wilcox.test(rho ~ Condition,
									data = dt.Veranda.Mushrooms,
										alternative = 'greater')
 # print(wc.SecurityON.EyesCovered)
		# W = 72, p-value = 0.05256

median(dt.Veranda.Mushrooms$rho)#[1] 0.8630919

DirectIndirect.wc.p <- c(wc.sr.FloodlightON.OFF.Direct$p.value, 
			wc.Veranda.OFF.Direct$p.value, 
			wc.Mushrooms.OFF.Direct$p.value,
			wc.Stars.Cloud.Indirect$p.value,#not really relevant
			wc.Urban.Cloud.Indirect$p.value, #not really relevant
			wc.Stars.Rural.Urban$p.value,  
			wc.SecurityON.EyesCovered$p.value,
			wc.SecurityOFF.EyesCovered$p.value,
			wc.Veranda.EyesCovered$p.value,
			wc.MushroomsUrban.EyesCovered$p.value,
			wc.Rural.EyesCovered$p.value,
			wc.Urban.EyesCovered$p.value,
			wc.Clouds.EyesCovered$p.value)
names(DirectIndirect.wc.p) <- c('wc.sr.FloodlightON.OFF.Direct', 
			'wc.Veranda.OFF.Direct',
			'wc.Mushrooms.OFF.Direct',
			'wc.Stars.Cloud.Indirect',#not really relevant
			'wc.Urban.Cloud.Indirect',#not really relevant
			'wc.Stars.Rural.Urban', 
			'wc.SecurityON.EyesCovered',  
			'wc.SecurityOFF.EyesCovered',
			'wc.Veranda.EyesCovered',
			'wc.MushroomsUrban.EyesCovered',
			'wc.Rural.EyesCovered',
			'wc.Urban.EyesCovered',
			'wc.Clouds.EyesCovered')
cbind(round(p.adjust(DirectIndirect.wc.p, 'BH'),5))
# wc.sr.FloodlightON.OFF.Direct 0.06824 #does not stand up to adjustment
# wc.Veranda.OFF.Direct         0.20840
# wc.Mushrooms.OFF.Direct       0.86007
# wc.Stars.Cloud.Indirect       0.05623
# wc.Urban.Cloud.Indirect       0.18552
# wc.Stars.Rural.Urban          0.06824 #does not stand up to adjustment
# wc.SecurityON.EyesCovered     0.00033 #Floodlight contributes
# wc.SecurityOFF.EyesCovered    0.00033 #Rural surroundings contribute
# wc.Veranda.EyesCovered        0.00045 #Veranda lights contribute
# wc.MushroomsUrban.EyesCovered 0.00318 #Mushrooms contribute
# wc.Rural.EyesCovered          0.05994 #does not stand up to adjustment
# wc.Urban.EyesCovered          0.20840 
# wc.Clouds.EyesCovered         0.52122

DirectIndirect.prent <- with(DirectIndirect, prentice.test(
									y = rho,
								groups = Direct,
								blocks = sky,
								alternative = 'greater'))#asking specifically if condition light pollution conditions were more oriented))
# statistic: chi-square = 18.436, df = 1, p-value = 8.784e-06
#Excluding eyes covered!
DirectIndirect.Eyes <- subset(DirectIndirect, sky != "EyesCovered")
DirectIndirect.Eyes.prent <- with(DirectIndirect.Eyes, prentice.test(
									y = rho,
								groups = Direct,
								blocks = sky,
								alternative = 'greater'))#asking 
# statistic: chi-square = 18.436, df = 1, p-value = 8.784e-06
#Actually, the grouping is doing nothing here, better to compare the starry skies
#Perhaps ignore conditions for this test
wilcox.test(rho~!Direct, data = subset(DirectIndirect, sky == 'Stars'), alternative = 'greater')#somehow TRUE is the reference condition!
#W = 807, p-value = 3.745e-06

# grouping asks a stupid question
# first only the direct ones
Direct.wc.p <- c(wc.sr.FloodlightON.OFF.Direct$p.value, 
			wc.Veranda.OFF.Direct$p.value, 
			wc.Mushrooms.OFF.Direct$p.value)
names(Direct.wc.p) <- c('wc.sr.FloodlightON.OFF.Direct', 
			'wc.Veranda.OFF.Direct',
			'wc.Mushrooms.OFF.Direct')
cbind(round(p.adjust(Direct.wc.p, 'BH'),5))
# wc.sr.FloodlightON.OFF.Direct 0.12598 #Pretty much all the same
# wc.Veranda.OFF.Direct         0.26451
# wc.Mushrooms.OFF.Direct       0.86007

# then only the indirect ones
Indirect.wc.p <- c(#wc.Stars.Cloud.Indirect$p.value,#not relevant
			#wc.Urban.Cloud.Indirect$p.value, #not relevant
			wc.Stars.Rural.Urban$p.value,
			wc.Stars.Cloud.Indirect$p.value,  
			wc.Rural.EyesCovered$p.value,
			wc.Urban.EyesCovered$p.value,
			wc.Clouds.EyesCovered$p.value)
names(Indirect.wc.p) <- c(#'wc.Stars.Cloud.Indirect',#not relevant
			#'wc.Urban.Cloud.Indirect',#not relevant
			'wc.Stars.Rural.Urban', 
			'wc.Stars.Cloud.Indirect',
			'wc.Rural.EyesCovered',
			'wc.Urban.EyesCovered',
			'wc.Clouds.EyesCovered')
cbind(round(p.adjust(Indirect.wc.p, 'BH'),5))
# wc.Stars.Rural.Urban    0.06405 ???
# wc.Stars.Cloud.Indirect 0.06405 ???
# wc.Rural.EyesCovered    0.06405 ???
# wc.Urban.EyesCovered    0.21919
# wc.Clouds.EyesCovered   0.48113


# DirectIndirect.prent <- with(DirectIndirect, prentice.test(
									# y = rho,
								# groups = Direct,
								# blocks = Indirect,
								# alternative = 'greater'))#asking specifically if direct light pollution conditions were more oriented))
# # statistic: chi-square = 30.136, df = 1, p-value = 2.014e-08
# #On average, orientation precision was higher in the presence of light pollution
bf.test(rho ~ Condition, data = DirectIndirect)
  # statistic  : 28.01555 
  # num df     : 7 
  # denom df   : 43.82327 
  # p.value    : 3.353887e-14 
#certainly some differences in variance, let's ignore these for now

bf.test(rho ~ Condition, data = dt.Stars.Rural.Urban)
  # statistic  : 3.023938 #dissappointing
  # num df     : 1 
  # denom df   : 13.92372 
  # p.value    : 0.1040941 
#####################################################################
#	Heading Choice													#
#####################################################################
source(paste0(Sys.getenv('HOME'), '/Dropbox/My Papers/Light Pollution/',
'MMRayleigh.test.R'))

for(cnd in levels(DirectIndirect$Condition)){
	aaa <- subset(lp.data.mu, Condition == cnd)$mu
	rrr <- subset(lp.data.rho, Condition == cnd)$rho
	print(cnd)
	rslt <- MMRayleigh.test(aaa,rrr)
	print(data.frame(rslt))
	print(rao.spacing.test(mycirc(aaa)))								
}
# [1] "EyesCoveredWallsUrban"
        # α.        R       R. n  p.value
# 1 56.84123 19.38251 1.046557 7 0.057414

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 163.9528 
# P-value > 0.10 
 
# [1] "VerandaRural"
        # α.       R        R.  n  p.value
# 1 173.4142 13.2287 0.4183283 10 0.655285

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 106.2432 
# P-value > 0.10 
 
# [1] "MushroomsUrban"
         # α.        R        R.  n  p.value
# 1 -75.51058 23.21957 0.7342673 10 0.257123

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 103.4185 
# P-value > 0.10 
 
# [1] "SecurityLightOFFRural"
        # α.        R         R.  n  p.value
# 1 106.8135 3.030379 0.09582899 10 0.978612

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 93.3948 
# P-value > 0.10 
 
# [1] "SecurityLightONRural"
         # α.        R       R.  n p.value
# 1 -149.2838 52.49815 1.660137 10       0

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 242.3485 
# P-value < 0.001 
 
# [1] "OvercastWallsRural"
        # α.       R        R.  n  p.value
# 1 10.14601 19.2112 0.6075116 10 0.402267

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 137.3724 
# P-value > 0.10 
 
# [1] "StarsWallsRural"
         # α.       R        R.  n  p.value
# 1 -82.43595 19.0986 0.6039508 10 0.406665

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 85.3382 
# P-value > 0.10 
 
# [1] "StarsWallsUrban"
        # α.        R        R.  n  p.value
# 1 127.7187 27.00809 0.5762075 13 0.425267

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 145.8593 
# P-value > 0.10 

range(mycirc(subset(lp.data.mu, Condition == 'SecurityLightONRural')$mu))
# 81.65148°
range(mycirc(subset(lp.data.mu, Condition == 'SecurityLightOFFRural')$mu))
# 258.3639°

IQR(mycirc(subset(lp.data.mu, Condition == 'SecurityLightONRural')$mu))
#5/10 within 44°
diff(quantile(mycirc(subset(lp.data.mu, Condition == 'SecurityLightONRural')$mu), c(1, 8)/10 ))
#7/10 within just 34°
#comp with light off
diff(quantile(mycirc(subset(lp.data.mu, Condition == 'SecurityLightOFFRural')$mu), c(1, 8)/10 ))
#203°
#####################################################################
#	Plot Data														#
#####################################################################

cbind(levels(DirectIndirect$Condition))
#Re-order for plotting
DirectIndirect$Condition <- factor(DirectIndirect$Condition,
			levels = c('SecurityLightOFFRural','SecurityLightONRural',
					'VerandaRural','MushroomsUrban',
					'StarsWallsRural','StarsWallsUrban',
					'OvercastWallsRural','EyesCoveredWallsUrban') )

#	Open plot
dev.new(width =5); par(mai = c(0,0.8, 0, 0), lend = 'butt')
#plot on the appropriate scale transformed	
boxplot(rho~Condition, data = DirectIndirect,
		cex = 0.5, outline = F, border = rgb(0,0,0,0.5),
		pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5),
		axes = F, ylim =(c(0,1)),# xlim = c(0.5, 4.5),
		ylab = 'Mean Vector Length', xlab = '')
polygon(c(0,2+rep(length(levels(DirectIndirect$Condition)),2),0), c(0,0, sqrt(-log(0.05)/10),sqrt(-log(0.05)/10)), col = rgb(1,0,0,0.05), border = NA)
legend('bottom', inset = sqrt(-log(0.05)/10), legend = '   Rayleigh test     \n       p<0.05     \n            \n       p>0.05     ', cex = 0.5, bty = 'n')
beeswarm(rho~Condition, data = DirectIndirect,
		pch = 21,  cex = 1.8/2, #method = 'center',
		pwcol = c('gray20', 'orange4')[DirectIndirect$Direct+1], 
		pwbg = c('gray', 'orange')[DirectIndirect$Indirect+1],
		add = T)
#	
	axis(2)#
	mtext(levels(DirectIndirect$Condition), side = 1, at = 1:length(levels(DirectIndirect$Condition))-0.5, line = -0.25-24*with(DirectIndirect, aggregate(rho, by = list(Condition = Condition), function(x) quantile(x, 0.75, na.rm=T)))$x, las = 2, cex = 0.75 )
	abline(h =c(0,1), lwd = 0.25)
	segments(4.35,0,4.35,1, lwd = 0.25, lty = 1)
	#Light ON to OFF lines
	segments(rep(1,10), 
	subset(lp.data.rho, Condition == 'SecurityLightOFFRural')$rho,
	rep(2,10),
	subset(lp.data.rho, Condition == 'SecurityLightONRural')$rho, 
	lty = 1, rgb(0,0,0,0.3))
	
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/'), Experiment = 'LP',  PlotName = paste0("Beeswarm_", 'DirectIndirect'))

#####################################################################
#	Add significant differences										#
#####################################################################
cbind(names(Direct.wc.p[p.adjust(Direct.wc.p, 'BH')<=0.05] ))
cbind(names(Indirect.wc.p[p.adjust(Indirect.wc.p, 'BH')<=0.05] ))
cbind(names(Indirect.wc.p[p.adjust(Indirect.wc.p, 'BH')<=0.1] ))
# [1,] "wc.Stars.Rural.Urban"
# [2,] "wc.Rural.EyesCovered"

highest <- 0.9
spc <- 0.03#spacing
#also from Light OFF to Light OFF
lines(c(1,2), rep(0.7-spc*0,2), col = 'darkred', lwd = 2)
text(1.5,0.7-spc*0.75, '*', cex = 1.5)

# So from StarsRural to EyesCovered
lines(c(5,8), rep(highest-spc*0,2), col = 'darkred', lwd = 0.5)
text(6.5,highest-spc*-0.5, '.', cex = 1.5)
# from StarsRural to StarsUrban
lines(c(5,6), rep(highest-spc*1,2), col = 'darkred', lwd = 0.5)
text(5.5,highest-spc*0.5, '.', cex = 1.5)
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/'), Experiment = 'LP',  PlotName = paste0("Hypotheses_", 'DirectIndirect'))