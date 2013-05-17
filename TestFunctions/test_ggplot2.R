
locateChr <- function(y){
	colText = 'grey40'
	colLines = 'grey80'
	arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
	hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
	cumLen = cumsum(as.numeric(hg19$length))
	cumCentr <- 1/2*cumLen[1]
	for(chr in 2:length(cumLen)) cumCentr = c(cumCentr, cumLen[chr-1] + 1/2*hg19$length[chr])
	+ geom_hline(yintercept = 0)
	+ geom_vline(xintercept = cumLen[1:23], color = 'grey30', linetype = 2, size = 0.1)
	# annotate(0, y, labels = "chr", cex = 0.75, col = colText)
	# ltext(cumLen[1]/2, y, labels = 1, cex = 0.75, col = colText)
	# for(i in 2:length(cumCentr)){
		# x <- (hg19$length[i]/2 + cumLen[i-1])
		# ltext(x, y, labels = i, cex = 0.75, col = colText)
		# }
	}

createProfile <- function(object, gain = log2(2.25/2), loss = log2(1.80/2)){
	require(ggplot2, quietly = TRUE)
	wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")

	hg19 <- read.csv(paste0(arrayInfoPath, 'human.chrom.info.hg19.FC.txt'), header = TRUE, sep = '\t')
	cumLen = cumsum(as.numeric(hg19$length))
	cumCentr <- 1/2*cumLen[1]
	for(chr in 2:length(cumLen)) cumCentr = c(cumCentr, cumLen[chr-1] + 1/2*hg19$length[chr])

	cnSet = getCNset(object)
	segTable = getSegTable(object)
	Title = paste(getInfo(object, 'sampleId'), '-', getInfo(object, 'analyseDate'),
						'\nGain threshold: ', round(gain, 3), ' Loss threshold:', round(loss, 3))

	if(any(is.na(cnSet))){
		NAs <- c()
		for(i in 1:ncol(cnSet))
			NAs <- c(NAs, which(is.na(cnSet[,i])))
		NAs <- unique(NAs)
		cnSet <- cnSet[-NAs,]
	}

	Samp = seq(1, nrow(cnSet), len = 15e3)
	cnSet = cnSet[Samp,]
	cnSet <- cbind.data.frame(cnSet, rMed = runmed(cnSet$Log2Ratio, k=71))

	gPlot <- ggplot(data = cnSet, aes(x = genomicPos, y = Log2Ratio)) +
					geom_point(pch = 19, cex = 0.25, col = 'grey80') +
					geom_hline(yintercept = 0) +
					geom_vline(xintercept = cumLen[1:23], color = 'grey30', linetype = 2, size = 0.1) +
					annotate('text', x = c(0, cumCentr), y = rep(1.4, 25), label = c("Chr", seq(1, 24)), size = 2.5, colour = 'grey30') +
					geom_point(aes(x = genomicPos, y = rMed), cex = 0.25, type = 'l') +
					geom_point(aes(x = genomicPos, y = Segm), cex = 1,
									col = ifelse(cnSet$Segm<= loss, 'red3', ifelse(cnSet$Segm>= gain, 'dodgerblue2', 'black'))) +
					coord_cartesian(ylim = range(-1.5, 1.5))+
					ggtitle(Title) +
					theme_bw() +
					theme(	panel.grid.major = element_blank(),
								panel.grid.minor = element_blank(),
								plot.title = element_text(lineheight=.8, face="bold"))

	object@gProfile <- gPlot
	return(object)
	}


tagGene <- function(object, gain = log2(2.25/2), loss = log2(1.80/2), tag = NULL){
	if(!'geneDB' %in% ls()){
		arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
		geneDB <- readRDS(paste0(arrayInfoPath, 'myGeneDB.rds'))
	}
	
	if(!is.null(tag)){
		segTable = getSegTable(object)
		tmp <- geneDB[which(geneDB$NomenclatureSymbol == tag),]
		LR <- segTable$seg.med[which(segTable$loc.start<= tmp$genomicStart & segTable$loc.end >= tmp$genomicStop)]
		Col <- ifelse(LR<= loss, 'red3', ifelse(LR>= gain, 'dodgerblue2', 'black'))
		dat <- data.frame(xstart = c(tmp$genomicStart-5e7, tmp$genomicStart),
									xend = c(tmp$genomicStart, tmp$genomicStart),
									ystart = c(1, 1),
									yend = c(1, LR),
									Col = rep(Col, 2))
		getProfile(object) +
		annotate("text", x = tmp$genomicStart-1.5e8, y = 1, label = paste0(tmp$NomenclatureSymbol, '\n(Log2R = ', round(LR, 3), ')'), cex = 5) +
		geom_segment(data = dat, aes(x = xstart, y = ystart, xend = xend, yend = yend), colour = Col)
	}
	else getProfile(object)
}

tagGene(object6, tag = 'EGFR')


myPlot 
annotate("text", 2e9, 1, label = 'another gene', cex = 5)