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
lp.data$Condition <- as.factor(lp.data$Condition)
#reorganise conditions
lp.data$Condition <- relevel(lp.data$Condition, 'EyesCoveredWallsUrban')
lp.data$Condition <- factor(lp.data$Condition, levels(lp.data$Condition)[c(1,2:3,10:11,5:6,14,4,8:9,7,12,13)])	
cbind(levels(lp.data$Condition))#is this what you want?

#	Mean vector length data that will be used in the paper
lp.data.rho <- read.table(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/LPdataRho20200505', '.txt') , header = T)#, sep  = '\t')
head(lp.data.rho)
lp.data.rho$Beetle <- as.factor(lp.data.rho$Beetle)
lp.data.rho$Condition <- as.factor(lp.data.rho$Condition)
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

#Choose data that appears in Figure 1
UrbanRural <- subset(lp.data.rho, !grepl('Walls', Condition)& 
								(grepl('Moon', Condition) |
								grepl('Stars', Condition) | 
								grepl('Overcast', Condition)) )
UrbanRural$sky <- with(UrbanRural, gsub('Rural','',Condition))
UrbanRural$sky <- with(UrbanRural, gsub('Urban','',sky))
UrbanRural$Condition <- factor(UrbanRural$Condition)


#####################################################################
#	Location Tests													#
#####################################################################
kt.UrbanRural <- kruskal.test(rho~Condition, data = UrbanRural)
# Kruskal-Wallis chi-squared = 20.701, df = 5, p-value = 0.0009224

#	Post-hoc tests													#
#remember condition order is important here
 cbind(levels(UrbanRural$Condition))
 
#	compare to reference											#

#	Rural
 # MoonRural - OvercastRural 
dt.Moon.Cloud.Rural <- subset(UrbanRural,
										Condition == 'MoonRural' | 
										Condition == 'OvercastRural')
dt.Moon.Cloud.Rural$Condition <- factor(
										dt.Moon.Cloud.Rural$Condition,
								levels = c('MoonRural','OvercastRural'))
								#MoonRural is alternative
wc.Moon.Cloud.Rural <- wilcox.test(rho ~ Condition,
									data = dt.Moon.Cloud.Rural,
									alternative = 'greater') 
# print(wc.Moon.Cloud.Rural)
		# W = 79, p-value = 0.0144 
 # StarsRural - OvercastRural 
dt.Stars.Cloud.Rural <- subset(UrbanRural,
										Condition == 'StarsRural' |
										Condition == 'OvercastRural')
dt.Stars.Cloud.Rural$Condition <- factor(
										dt.Stars.Cloud.Rural$Condition,
							levels = c('StarsRural','OvercastRural'))
							#StarsRural is alternative					
wc.Stars.Cloud.Rural <- wilcox.test(rho ~ Condition,
										data = dt.Stars.Cloud.Rural,
										 alternative = 'greater')
# print(wc.Stars.Cloud.Rural)
		# W = 156, p-value = 0.006365
		
#	Urban
 # MoonUrban - OvercastUrban
dt.Moon.Cloud.Urban <- subset(UrbanRural,
										Condition == 'MoonUrban' |
										Condition == 'OvercastUrban')
dt.Moon.Cloud.Urban$Condition <- factor(dt.Moon.Cloud.Urban$Condition,
								levels = c('MoonUrban','OvercastUrban'))
wc.Moon.Cloud.Urban <- wilcox.test(rho ~ Condition,
											data = dt.Moon.Cloud.Urban,
											alternative = 'greater')
 # print(wc.Moon.Cloud.Urban)
		# W = 63, p-value = 0.1763 

 # StarsUrban - OvercastUrban
dt.Stars.Cloud.Urban <- subset(UrbanRural,
										 Condition == 'StarsUrban' |
										 Condition == 'OvercastUrban')
dt.Stars.Cloud.Urban$Condition <- factor(dt.Stars.Cloud.Urban$Condition,
							levels = c('StarsUrban', 'OvercastUrban'))
wc.Stars.Cloud.Urban <- wilcox.test(rho ~ Condition,
											data = dt.Stars.Cloud.Urban,
											alternative = 'greater')
 # print(wc.Stars.Cloud.Urban)
		# W = 78, p-value = 0.8359
		
 # MoonRural - MoonUrban 
dt.Moon.Rural.Urban <- subset(UrbanRural, Condition == 'MoonRural' | Condition == 'MoonUrban')
dt.Moon.Rural.Urban$Condition <- factor(dt.Moon.Rural.Urban$Condition,
									levels = c('MoonUrban','MoonRural'))
wc.Moon.Rural.Urban <- wilcox.test(rho ~ Condition,
										data = dt.Moon.Rural.Urban,
										alternative = 'greater')
 # print(wc.Moon.Rural.Urban)
		# W = 34, p-value = 0.1237 	
		
 # StarsRural - StarsUrban 
dt.Stars.Rural.Urban <- subset(UrbanRural,
											Condition == 'StarsRural' |
											Condition == 'StarsUrban')
dt.Stars.Rural.Urban$Condition <- factor(dt.Stars.Rural.Urban$Condition,
								levels = c('StarsUrban','StarsRural'))
wc.Stars.Rural.Urban <- wilcox.test(rho ~ Condition,
											data = dt.Stars.Rural.Urban,
											alternative = 'greater')
 # print(wc.Stars.Rural.Urban)
		# W = 136, p-value = 0.04296
		
 # OvercastRural - OvercastUrban 
dt.Cloud.Rural.Urban <- subset(UrbanRural,
										Condition == 'OvercastRural' |
										Condition == 'OvercastUrban')
dt.Cloud.Rural.Urban$Condition <- factor(dt.Cloud.Rural.Urban$Condition,
						levels = c('OvercastUrban','OvercastRural'))
wc.Cloud.Rural.Urban <- wilcox.test(rho ~ Condition,
									data = dt.Cloud.Rural.Urban,
									alternative = 'greater')
 # print(wc.Cloud.Rural.Urban)
		# W = 11, p-value = 0.001045
		
 # MoonUrban - OvercastRural
dt.Moon.Cloud.Urban.Rural <- subset(UrbanRural,
										 Condition == 'MoonUrban' |
										 Condition == 'OvercastRural')
dt.Moon.Cloud.Urban.Rural$Condition <- 
							factor(dt.Moon.Cloud.Urban.Rural$Condition,
							levels = c('MoonUrban', 'OvercastRural'))
wc.Moon.Cloud.Urban.Rural <- wilcox.test(rho ~ Condition,
									data = dt.Moon.Cloud.Urban.Rural,
											alternative = 'greater')
 # print(wc.Moon.Cloud.Urban.Rural)
		# W = 91, p-value = 0.000525

 # StarsUrban - OvercastRural
dt.Stars.Cloud.Urban.Rural <- subset(UrbanRural,
										 Condition == 'StarsUrban' |
										 Condition == 'OvercastRural')
dt.Stars.Cloud.Urban.Rural$Condition <- 
							factor(dt.Stars.Cloud.Urban.Rural$Condition,
							levels = c('StarsUrban', 'OvercastRural'))
wc.Stars.Cloud.Urban.Rural <- wilcox.test(rho ~ Condition,
									data = dt.Stars.Cloud.Urban.Rural,
											alternative = 'greater')
 # print(wc.Stars.Cloud.Urban.Rural)
		# W = 162, p-value = 0.002643
#?
 # MoonRural - StarsRural 
dt.Moon.Stars.Rural <- subset(UrbanRural,
										Condition == 'MoonRural' | 
										Condition == 'StarsRural')
dt.Moon.Stars.Rural$Condition <- factor(
										dt.Moon.Stars.Rural$Condition,
								levels = c('MoonRural','StarsRural'))
								#MoonRural is alternative
wc.Moon.Stars.Rural <- wilcox.test(rho ~ Condition,
									data = dt.Moon.Stars.Rural,
									alternative = 'greater') 
#W = 139, p-value = 0.04529



UrbanRural.wc.p <- c(wc.Moon.Cloud.Rural$p.value, 
			wc.Stars.Cloud.Rural$p.value, 
			wc.Moon.Cloud.Urban$p.value,
			wc.Stars.Cloud.Urban$p.value,
			wc.Moon.Rural.Urban$p.value, 
			wc.Stars.Rural.Urban$p.value,  
			wc.Cloud.Rural.Urban$p.value, 
			wc.Moon.Cloud.Urban.Rural$p.value,
			wc.Stars.Cloud.Urban.Rural$p.value)
names(UrbanRural.wc.p) <- c('wc.Moon.Cloud.Rural', 
			'wc.Stars.Cloud.Rural',
			'wc.Moon.Cloud.Urban',
			'wc.Stars.Cloud.Urban',
			'wc.Moon.Rural.Urban',
			'wc.Stars.Rural.Urban', 
			'wc.Cloud.Rural.Urban', 
			'wc.Moon.Cloud.Urban.Rural',
			'wc.Stars.Cloud.Urban.Rural')
cbind(round(p.adjust(UrbanRural.wc.p, 'BH'),5))
# wc.Moon.Cloud.Rural  0.03361 #Moon contributes
# wc.Stars.Cloud.Rural 0.02228 #Stars contribute
# wc.Moon.Cloud.Urban  0.20573 #
# wc.Stars.Cloud.Urban 0.83591 #
# wc.Moon.Rural.Urban  0.17322 #
# wc.Stars.Rural.Urban 0.07517 #Urban contributes
# wc.Cloud.Rural.Urban 0.00731 #Urban contributes
# wc.Moon.Cloud.Urban.Rural  0.00470 #Urban contributes
# wc.Stars.Cloud.Urban.Rural 0.00793 #Urban contributes
UrbanRural.prent <- with(UrbanRural, prentice.test(
									y = rho,
								groups = Indirect,
								blocks = sky,
								alternative = 'greater'))#asking specifically if condition light pollution conditions were more oriented))
# statistic: chi-square = 8.2491, df = 1, p-value = 0.002039
median(subset(UrbanRural, Indirect == F)$rho)#no light pollution
# 0.8673096
median(subset(UrbanRural, Indirect == T)$rho)#light pollution
# 0.9534928
#On average, orientation precision was higher in the presence of light pollution
bf.test(rho ~ Condition, data = UrbanRural)
  # statistic  : 7.489311 
  # num df     : 5 
  # denom df   : 39.36693 
  # p.value    : 5.10413e-05 
#certainly some differences in variance, let's ignore these for now
#	Rural

 # MoonRural - StarsRural 
# dt.Moon.Stars.Rural <- subset(UrbanRural,
										# Condition == 'MoonRural' | 
										# Condition == 'StarsRural')
# dt.Moon.Stars.Rural$Condition <- factor(
										# dt.Moon.Stars.Rural$Condition,
								# levels = c('MoonRural','StarsRural'))
								# #MoonRural is alternative
bf.Moon.Stars.Rural <- bf.test(rho ~ Condition,
									data = dt.Moon.Stars.Rural) 
#####################################################################
#	Heading Choice													#
#####################################################################									
MMRayleigh.test <- function(mang, mvec, method = 'simulation'){
	mrstar <- function(ang,vec){
		nn <- length(ang)#number of vectors
		mv.rank <- rank(vec)#rank of each vector
		ma.rad <- ang*pi/180#angles in radians
		x.rank <- sum(mv.rank*cos(ma.rad))#x projection of angles by rank
		y.rank <- sum(mv.rank*sin(ma.rad))#y projection of angles by rank
		MM.R <- sqrt(x.rank^2+y.rank^2)#Resultant vector
		MM.mv <- MM.R/nn#Normalised resultant vector (not used)
		MM.rstar <- MM.R/(nn^(3/2))#Tranformation used as test statistic
		return(MM.rstar)
	}
	nn <- length(mang)#number of vectors
	mv.rank <- rank(mvec)#rank of each vector
	ma.rad <- mang*pi/180#angles in radians
	x.rank <- sum(mv.rank*cos(ma.rad))#x projection of angles by rank
	y.rank <- sum(mv.rank*sin(ma.rad))#y projection of angles by rank
	MM.R <- sqrt(x.rank^2+y.rank^2)#Resultant vector
	MM.mv <- MM.R/nn#Normalised resultant vector (not used)
	MM.rstar <- MM.R/(nn^(3/2))#Tranformation used as test statistic
	MM.ma <- atan2(y.rank,x.rank)*180/pi#Resultant angle
	sig.sq <- nn*(nn+1)*(2*nn+1)/12#estimated variance of rstar distribution
	#P value calculation is a long step
	#inspiration taken from
	#https://github.com/ruthcfong/Family-of-Rayleigh-Statistics
	#In Matlab this takes a permutation approach which is quite conservative
	#focussing on the question of whether vectors are longer near pop-mean
	#To calculate answers closer to Moore, I use a semi-empirical,
	#simulation method, reassigning vector lengths to
	#a uniform distribution on the circle
    num.samples = 1e6;
    perm.rstar <- rep(NA, num.samples)
    for(i in 1:num.samples){
        if(method == 'permutation'){
        # perm.rstar[i] <- mrstar(mang,sample(mvec))
        }else{
        perm.rstar[i] <- mrstar(runif(nn)*360-.Machine$double.eps,mvec)
        }
    }
    # #faster as a table?
    # permi <- sort(rep(1:num.samples,nn))
    # dtf <- data.frame(i  = permi,
				    # aa = runif(nn* num.samples)*360-.Machine$double.eps,
				    # mv = rep(mvec,num.samples))    
	p.val <- sum(perm.rstar>MM.rstar)/sum(!is.na(perm.rstar))
	# p.val <- pchisq(MM.rstar^2/sig.sq, nn)#probability of > rstar
	result <- list(`α*` = MM.ma, R = MM.R, `R*` = MM.rstar, n = nn, p.value = p.val)
	return(result)
}#MMRayleigh.test

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
}###################	END OF FUNCTION	###########################	

for(cnd in levels(UrbanRural$Condition)){
	aaa <- subset(lp.data.mu, Condition == cnd)$mu
	rrr <- subset(lp.data.rho, Condition == cnd)$rho
	print(cnd)
	rslt <- MMRayleigh.test(aaa,rrr)
	print(data.frame(rslt))
	print(rao.spacing.test(mycirc(aaa)))								
}
# [1] "MoonRural"
        # α.        R        R.  n  p.value
# 1 32.71862 18.98214 0.6002679 10 0.411248

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 157.9346 
# P-value > 0.10 
 
# [1] "MoonUrban"
       # α.        R       R.  n  p.value
# 1 73.5945 30.09391 0.951653 10 0.091417

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 191.7851 
# 0.01 < P-value < 0.05 
 
# [1] "StarsRural"
        # α.        R      R.  n  p.value
# 1 79.21921 141.9599 1.58716 20 0.000285

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 180.558 
# 0.001 < P-value < 0.01 
 
# [1] "StarsUrban"
         # α.        R       R.  n  p.value
# 1 -172.1731 70.43802 0.787521 20 0.180238

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 131.3955 
# P-value > 0.10 
 
# [1] "OvercastRural"
         # α.        R        R.  n  p.value
# 1 -152.6841 15.39424 0.4868087 10 0.562372

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 185.6164 
# 0.01 < P-value < 0.05 
 
# [1] "OvercastUrban"
      # α.        R       R.  n  p.value
# 1 -99.87 22.50903 0.711798 10 0.280059

       # Rao's Spacing Test of Uniformity 
 
# Test Statistic = 153.9682 
# P-value > 0.10 


#Pairwise comparisons
#Moon
watson.wheeler.test(mycirc(mu)~ Indirect,
					 subset(lp.data.mu,
						 grepl('Moon', Condition) &
						 !grepl('Walls', Condition)
						 )
					 )
# Watson-Wheeler test for homogeneity of angles

# data:  mycirc(mu) by Indirect
# W = 1.5754, df = 2, p-value = 0.4549

#Stars
watson.wheeler.test(mycirc(mu)~Indirect,
					 subset(lp.data.mu,
						 grepl('Star', Condition) &
						 !grepl('Walls', Condition)
						 )
					 )
# Watson-Wheeler test for homogeneity of angles

# data:  mycirc(mu) by Indirect
# W = 10.724, df = 2, p-value = 0.004691
			 
#Clouds
watson.wheeler.test(mycirc(mu)~Indirect,
					 subset(lp.data.mu,
						 grepl('Overcast', Condition) &
						 !grepl('Walls', Condition)
						 )
					 )

# Watson-Wheeler test for homogeneity of angles

# data:  mycirc(mu) by Indirect
# W = 2.0706, df = 2, p-value = 0.3551
#####################################################################
#	Plot Data														#
#####################################################################

cbind(levels(UrbanRural$Condition))
#Marie's suggestion for ordering
UrbanRural$Condition <- factor(UrbanRural$Condition,
					levels = c('MoonRural','StarsRural','OvercastRural',
					'MoonUrban','StarsUrban','OvercastUrban') )

#	Open plot
dev.new(width =5); par(mai = c(0,0.8, 0, 0), lend = 'butt')
#plot on the appropriate scale transformed	
boxplot(rho~Condition, data = UrbanRural,
		cex = 0.5, outline = F, border = rgb(0,0,0,0.5),
		pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5),
		axes = F, ylim =(c(0,1)),# xlim = c(0.5, 4.5),
		ylab = 'Mean Vector Length', xlab = '')
polygon(c(0,2+rep(length(levels(UrbanRural$Condition)),2),0), c(0,0, sqrt(-log(0.05)/10),sqrt(-log(0.05)/10)), col = rgb(1,0,0,0.05), border = NA)
legend('bottom', inset = sqrt(-log(0.05)/10), legend = '   Rayleigh test     \n       p<0.05     \n            \n       p>0.05     ', cex = 0.5, bty = 'n')
beeswarm(rho~Condition, data = UrbanRural,
		pch = 21,  cex = 1.8/2, #method = 'center',
		pwcol = c('gray20', 'orange4')[UrbanRural$Direct+1], 
		pwbg = c('gray', 'orange')[UrbanRural$Indirect+1],
		add = T)
#	
	axis(2)#
	mtext(levels(UrbanRural$Condition), side = 1, at = 1:length(levels(UrbanRural$Condition))-0.5, line = -0.25-24*with(UrbanRural, aggregate(rho, by = list(Condition = Condition), function(x) quantile(x, 0.75, na.rm=T)))$x, las = 2, cex = 0.75 )
	abline(h =c(0,1), lwd = 0.25)
	segments(3.35,0,3.35,1, lwd = 0.25, lty = 2)
	
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/'), Experiment = 'LP',  PlotName = paste0("Beeswarm_", 'UrbanRural'))

#####################################################################
#	Add significant differences										#
#####################################################################
cbind(names(UrbanRural.wc.p[p.adjust(UrbanRural.wc.p, 'BH')<=0.05] ))
# [1,] "wc.Moon.Cloud.Rural" 
# [2,] "wc.Stars.Cloud.Rural"
# [3,] "wc.Cloud.Rural.Urban"
# So from MoonRural to CloudRural
# from StarsRural to CloudRural
# and from CloudRural to CloudUrban

lines(c(1,3), rep(0.4-0.03*0,2), col = 'darkred', lwd = 2)
text(1.5,0.4-0.03*0.75, '*', cex = 1.5)
lines(c(2,3), rep(0.4-0.03*1,2), col = 'darkred', lwd = 2)
text(2.5,0.4-0.03*1.75, '*', cex = 1.5)
lines(c(3,6), rep(0.4-0.03*2,2), col = 'darkred', lwd = 2)
text(4.5,0.4-0.03*2.75, '**', cex = 1.5)
PDFsave(paste0(Sys.getenv('HOME'),'/Dropbox/My Papers/Light Pollution/'), Experiment = 'LP',  PlotName = paste0("Hypotheses_", 'UrbanRural'))