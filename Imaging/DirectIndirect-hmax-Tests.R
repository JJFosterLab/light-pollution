cnd <- c('floodlightON','floodlightOFF',
			'veranda', 'mushrooms',
			'starsWallsUrbanFullview','starsWallsUrban', 
			'starsWallsRural', 'cloudsWallsRural')
di.01.90 <- c(12,20,
				35,16,
				15,82,
				22,44)
names(di.01.90) <- cnd
di.01.55 <- c(12,20,
				33,16,
				24,82,
				22,49)
names(di.01.55) <- cnd
di.55.90 <- c(24,71,
				48,119,
				32,91,
				55,49)
names(di.55.90) <- cnd

# wilcox.test(urban.01.90, rural.01.90,alternative = 'less')# paired = T, 
# #W = 3, p-value = 0.35
# wilcox.test(urban.55.90, rural.55.90,alternative = 'less')# paired = T, 
# #W = 0, p-value = 0.05

snap.index <- function(x){round(180/x,1)}

snap.index(di.01.90)
snap.index(di.01.55)
snap.index(di.55.90)

t(rbind(
		snap.index(di.01.90),
		snap.index(di.01.55),
		snap.index(di.55.90)
		))
                        # [,1] [,2] [,3]
# # floodlightON          15.0 15.0  7.5
# floodlightOFF            9.0  9.0  2.5
# veranda                  5.1  5.5  3.8
# mushrooms               11.2 11.2  1.5
# starsWallsUrbanFullview 12.0  7.5  5.6
# starsWallsUrban          2.2  2.2  2.0
# starsWallsRural          8.2  8.2  3.3
# cloudsWallsRural         4.1  3.7  3.7