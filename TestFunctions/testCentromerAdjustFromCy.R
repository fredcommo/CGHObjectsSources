
#centromerAdjust <- function(object){

	arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
	hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')

#	cnSet <- getCNset(object1)
	cy3 <- log2(cnSet$gMedianSignal)
	cy5 <- log2(cnSet$rMedianSignal) 
#	pos <- cnSet$ChrStart
#	chr <- cnSet$ChrNum

#	centro <- newLR <- test <- lDist <- rDist <- c()

	nprobes = 100
	change = 1e-6
	K = 5

#	if(getInfo(object, 'platform') == "Affymetrix") {K = K*6; nprobes = nprobes*6}
#	Rmed <- runmed(lr, k = K)
#	medlr <- median(Rmed, na.rm = T)
#	for(i in 1:24){
		
	Values <- lapply(seq(1:24), function(i){
		cat('\nchr:', i)
		output <- c()
		cL <- hg19$centromerStart[i]
		cR <- hg19$centromerEnd[i]
		index <- which(chr == i)
		tmpCy3 <- cy3[index]
		tmpCy5 <- cy5[index]
#		tmp.rmed <- Rmed[index]
#		tmp.pos <- pos[index]
		Start <- cnSet$ChrStart[index]
		DL <- cL - Start[Start<=cL]
		DR <- Start[Start>=cR] - cR
			if(length(DL)>0){
				last <- length(which(Start<=cL))
				first <- last - nprobes + 1
				first <- ifelse(first>0, first, 1)
				cat('\nLeft' ,'\tlast', last, '\tfirst', first, '\tlen:', last-first)
				cy3L <- tmpCy3[first:last]
				cy5L <- tmpCy5[first:last]
				posL <- Start[first:last]
				Left <- cbind(chr = rep(i, length(cy3L)), pos = posL, cy3 = cy3L, cy5 = cy5L, loc = rep('left', length(cy3L)))
				output <- rbind(output, Left)
				}
			if(length(DR)>=nprobes){
				first <- which(Start >= cR)[1]
				last <- first + nprobes
				cat('\nRight' ,'\tlast', last, '\tfirst', first, '\tlen:', last-first)
				cy3R <- tmpCy3[first:last]
				cy5R <- tmpCy5[first:last]
				posR <- Start[first:last]
				Right <- cbind(chr = rep(i, length(cy3R)), pos = posR, cy3 = cy3R, cy5 = cy5R, loc = rep('right', length(cy3R)))
				output <- rbind(output, Right)
				}
		output
		}
	)
		
cent <- do.call(rbind.data.frame, Values)
for(i in 1:4) cent[,i] <- as.numeric(as.character(cent[,i]))
head(cent)
plot(density(cent$cy3, na.rm = TRUE))
lines(density(cent$cy5, na.rm = TRUE))

# Adjust on peak distribution
par(mfrow = c(1, 2))
for(loc in c('left', 'right')){
	loc = 'right'
	X <- cent$cy5 #[cent$loc == loc]
	if(any(is.na(X))) X <- X[!is.na(X)]
	leftModel <- buildEMmodel(X, by = 1, G = 1:6)
	n <- length(X)
	m <- leftModel$m
	p <- leftModel$p
	s <- leftModel$s
	n; m; p; s
	#mergedC <- mergePeaks(centroModel$nG, m, s, p, MergeVal = 0.03)
	#m <- mergedC$m
	#p <- mergedC$p
	#s <- mergedC$s
	simD <- computeDensities(n, m, p, s)
	bestC <- chooseBestPeak(simD$peaks, m, 0.9)
	plotEMmodel(X, simD$dList, m, bestC, cut = c(5, 15), Title = paste0('Centromer: Cy3 distribution'))
	}


		probeslo <- probeshi <- rnorm(nprobes, tukey.biweight(tmp.rmed), sd(tmp.rmed))
	
		if(length(DL) >= nprobes){
			last <- length(which(Start<=cL))
			first <- last - nprobes + 1
			cy3L <- tmpCy3[first:last]
			cy5L <- tmpCy5[first:last]
			posL <- Start[first:last]

			plot(posL, cy3L, cex = 0.5, col = 'blue')
			points(posL, cy5L, cex = 0.5, col = 'red')
			plot(cy3L, cy5L, asp = 1)
			abline(0, 1)
			test <- lm(cy5L~cy3L)
		summary(test)

			centro <- rbind(centro, cbind(val = probesL, pos = posL, Med = rep(median(probesL, na.rm = T), length(probesL)), loc = rep("left", length(probesL))))
			d <- DL[first]
			lDist <- rbind(lDist, cbind(chr = i, d = d))
			}

		if(length(DR) >= nprobes){
			first <- which(Start >= cR)[1]
			last <- first + nprobes
			probesR <- tmp.rmed[first:last]
			posR <- tmp.pos[first:last]
			centro <- rbind(centro, cbind(val = probesR, pos = posR, Med = rep(median(probesR, na.rm = T), length(probesR)), loc = rep("right", length(probesR))))
			d <- DR[nprobes]
			rDist <- rbind(rDist, cbind(chr = i, d = d))
			}
	
	#p <- wilcox.test(probesL, probesR)$p.value
	p <- t.test(probesL, probesR)$p.value
	test <- rbind(test, c(chr = i, left = mean(probesL, na.rm = TRUE), right = mean(probesR, na.rm = TRUE), pvalue = p, adjusted = ifelse(p>change, '*', '')))
	if(p > change)
		tmpLR <- tmpLR - tukey.biweight(c(probesL, probesR))
		#tmpLR <- tmpLR - mean(c(probesL, probesR))
	newLR <- c(newLR, tmpLR)
	}
	centro <- as.data.frame(centro)
	for(Col in 1:3) centro[,Col] <- as.numeric(as.character(centro[,Col]))
	return(list(test = as.data.frame(test), centroValues = centro, adjLr = newLR))
}

adjC <- centromerAdjust(object5)
centroValues <- adjC$centroValues

# Adjust on peak distribution
X <- rep(centroValues$Med, 5)
centroModel <- buildEMmodel(X, by = 1, G = 1:6)
n <- length(centroValues$Med)
m <- centroModel$m
p <- centroModel$p
s <- centroModel$s
n; m; p; s
mergedC <- mergePeaks(centroModel$nG, m, s, p, MergeVal = 0.03)
m <- mergedC$m
p <- mergedC$p
s <- mergedC$s
simD <- computeDensities(n, m, p, s)
bestC <- chooseBestPeak(simD$peaks, m, 0.9)
plotEMmodel(X, simD$dList, m, bestC, cut = c(-0.5, 0.5), Title = 'Centromer distribution')

x <- log(seq(1, length(m))/2)
model <- lm(m ~ x)
b = coef(model)[2]
a = coef(model)[1]
r <- summary(model)$r.squared
m2 = m + (x - (a+b*m))*r
plot(x, m, ylim = range(-1, 1), col = 'red', asp = 1)
abline(model)
abline(0, 1)
points(x, m2, col = 'blue')

adjLR <- adjC$adjLr + (scale(adjC$adjLr)-(a+b*adjC$adjLr))*r
