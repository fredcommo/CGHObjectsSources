# Accessors :Functions to have access to attributes

setGeneric("getInfo", function(object, arg = NULL,...) standardGeneric("getInfo"))
setMethod('getInfo', signature('cghObj'), function(object, arg = NULL,...) {
											out = data.frame(info = object@info)
											if (is.null(arg))
												return (out)
											else
												if (arg %in% rownames(out)) return (as.character(out$info[rownames(out) == arg]))
												else cat(paste('\'', arg, '\'', sep = ''), 'is not an available item.\n')
												}
											)

setGeneric('getCNset', function(object) standardGeneric('getCNset'))
setMethod('getCNset', 'cghObj', function(object) object@cnSet)

setGeneric('getSNPset', function(object) standardGeneric('getSNPset'))
setMethod('getSNPset', 'cghObj', function(object) object@snpSet)

setGeneric('getParam', function(object) standardGeneric('getParam'))
setMethod('getParam', 'cghObj', function(object) object@param)

setGeneric('getSegTable', function(object) standardGeneric('getSegTable'))
setMethod('getSegTable', 'cghObj', function(object) object@segTable)

setGeneric('getDensity', function(object) standardGeneric('getDensity'))
setMethod('getDensity', 'cghObj', function(object) print(object@probesDensity))

setGeneric('getProfile', function(object,...) standardGeneric('getProfile'))
#setMethod('getProfile', 'cghObj', function(object,...) return(object@gProfile))
setMethod('getProfile', 'cghObj', function(object, ylim = NULL,...) {gPlot <- object@gProfile
																							if(!is.null(ylim)){
																								hg19 <- read.csv(paste0(arrayInfoPath, 'human.chrom.info.hg19.FC.txt'),
																															header = TRUE, sep = '\t')
																								cumLen = cumsum(as.numeric(hg19$length))
																								cumCentr <- 1/2*cumLen[1]
																								for(chr in 2:length(cumLen)) cumCentr = c(cumCentr, cumLen[chr-1] + 1/2*hg19$length[chr])
																								gPlot <- gPlot+
																											coord_cartesian(ylim = ylim)+
																											annotate('text', x = c(0, cumCentr[1:23]), y = rep(max(ylim)-0.1, 24),
																															label = c("Chr", seq(1, 23)), size = 2.5, colour = 'grey30')
																								}
																							return(gPlot)
																							})

setGeneric('tagMyGene', function(object,...) standardGeneric('tagMyGene'))
setMethod('tagMyGene', 'cghObj', function(object, tag = NULL, ylim = range(-1.5, 1.5), gain = log2(2.25/2), loss = log2(1.80/2)){
	myBlue <- rgb(0, 0.45, 1, 1)
	if(!'geneDB' %in% ls()){
		arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
		geneDB <- readRDS(paste0(arrayInfoPath, 'myGeneDB.rds'))
		}
	if(is.null(tag)) gPlot
	else if(sum(is.element(tag, geneDB$NomenclatureSymbol)) == 0){
		cat('No valid gene symbol(s)\n')
		gPlot
	}
	else {
		valid <- is.element(tag, geneDB$NomenclatureSymbol)
		if(any(!valid)) cat('Not a valid gene symbol:\t', tag[!valid], '\n')
		tmp <- geneDB[which(geneDB$NomenclatureSymbol %in% tag),]
		segTable = getSegTable(object)
		gPlot <- getProfile(object)
#		LR <- segTable$seg.med[which(segTable$loc.start<= tmp$genomicStart & segTable$loc.end >= tmp$genomicStop)]
		LR <- lapply(1:nrow(tmp), function(i){
					startAt <- tmp$genomicStart[i]
					stopAt <- tmp$genomicStop[i]
					segTable$seg.med[which(segTable$loc.start<= startAt & segTable$loc.end >= stopAt)]
					})
		LR <- do.call(c, LR)
		Col <- ifelse(LR<= loss, 'red3', ifelse(LR>= gain, myBlue, 'grey40'))
		dat <- data.frame(xstart = c(tmp$genomicStart-2.5e7, tmp$genomicStart),
								xend = c(tmp$genomicStart, tmp$genomicStart),
								#ystart = c(rep(1, length(LR)), rep(1, length(LR))),
								#yend = c(rep(1, length(LR)), LR),
								ystart = c(LR/abs(LR)*1.2, LR/abs(LR)*1.2),
								yend = c(LR/abs(LR)*1.2, LR),
								Col = rep(Col, 2))
		gPlot +
		coord_cartesian(ylim = ylim)+
		annotate("text", x = tmp$genomicStart - 2.2e8, y = LR/abs(LR)*1.2,
						label = paste0(tmp$NomenclatureSymbol, '\n(Log2R = ', round(LR, 3), ')'), cex = 4) +
		geom_segment(data = dat, aes(x = xstart, y = ystart, xend = xend, yend = yend), colour = as.character(dat$Col), size = 0.75)
		}
	}
)

setGeneric('zoomIn', function(object,...) standardGeneric('zoomIn'))
setMethod('zoomIn', 'cghObj', function(object, chr = NA){
	if(!is.na(chr)){
		X11(type = 'dbcairo')
		cnSet = getCNset(object)
		index <- which(cnSet$ChrNum == chr)
		minY <- min(min(cnSet$Segm[index], na.rm = TRUE), -1)
		maxY <- max(max(cnSet$Segm[index], na.rm = TRUE), -1)
		xRange = range(cnSet$genomicPos[index], na.rm = TRUE)
		gPlot <- getProfile(object)
		gPlot + 	coord_cartesian(xlim = xRange, ylim = range(minY, 1.5)) +
		annotate('rect', xmin = min(xRange), xmax = max(xRange), ymin = 1.32, ymax = 1.5, fill = 'grey75')+
		annotate('text', x = median(xRange), y = 1.4, label = paste0("Chr", chr), size = 7.5, colour = 'grey15')+
		xlab('Chromosome position (Kb)')+
		scale_x_continuous(breaks = seq(min(xRange), max(xRange), len = 10),	
									labels = signif(seq(min(cnSet$ChrStart, na.rm = TRUE), max(cnSet$ChrStart, na.rm = TRUE), len = 10)/1e3, 3))
		}
	}
)

