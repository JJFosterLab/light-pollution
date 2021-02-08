cnd <- c('moon','stars','cloud')
urban.01.90 <- c(14,13,28); names(urban.01.90) <- cnd
urban.01.55 <- c(16,26,35); names(urban.01.55) <- cnd
urban.55.90 <- c(13,15,24); names(urban.55.90) <- cnd
rural.01.90 <- c(17,23,26); names(rural.01.90) <- cnd
rural.01.55 <- c(17,23,24); names(rural.01.55) <- cnd
rural.55.90 <- c(115,48,46); names(rural.55.90) <- cnd


wilcox.test(urban.01.90, rural.01.90,alternative = 'less')# paired = T, 
#W = 3, p-value = 0.35
wilcox.test(urban.55.90, rural.55.90,alternative = 'less')# paired = T, 
#W = 0, p-value = 0.05

snap.index <- function(x){round(180/x,1)}
						# moon stars cloud 
snap.index(urban.01.90)# 12.9  13.8   6.4 
snap.index(urban.01.55)# 11.2   6.9   5.1 
snap.index(urban.55.90)# 13.8  12.0   7.5
snap.index(rural.01.90)# 10.6   7.8   6.9 
snap.index(rural.01.55)# 10.6   7.8   7.5
snap.index(rural.55.90)#  1.6   3.8   3.9 
