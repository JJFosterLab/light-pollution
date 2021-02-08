#Let's write a Moore's Modified Rayleigh test

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
}
