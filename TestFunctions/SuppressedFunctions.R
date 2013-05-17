# Unused or suppressed functions


geneRequest <- function(geneList){
	# Suppressed because biomaRt doesn't return up-to-date information
	require(biomaRt)
	human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	output = getBM(attributes=c('hgnc_symbol', 'description', 'chromosome_name', 'band', 'start_position','end_position', 'entrezgene'),
			filters = 'hgnc_symbol', values = geneList, mart = human)
	output$description =  gsub(' \\[.*\\]', '', output$description)
	output = as.data.frame(output)
	# Add genomic positions
	arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
	hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
	cumLen = cumsum(as.numeric(hg19$length))
	cumLen = c(0, cumLen[-length(cumLen)])
	output = cbind.data.frame(output[,-ncol(output)], genomicStart = rep(NA, nrow(output)), genomicEnd = rep(NA, nrow(output)), entrezgene = output[,ncol(output)])
	foreach(i = 1:nrow(output)) %do%{
		chr = as.numeric(output$chromosome_name[i])
		output$genomicStart[i] = output$start_position[i] + cumLen[chr]
		output$genomicEnd[i] = output$end_position[i] + cumLen[chr]
		}
	return(output)
}


#################
# Plot functions
#################

locateChr <- function(y){
	colText = 'grey40'
	colLines = 'grey80'
	arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
	hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
	cumLen = cumsum(as.numeric(hg19$length))
	cumCentr <- 1/2*cumLen[1]
	for(chr in 2:length(cumLen)) cumCentr = c(cumCentr, cumLen[chr-1] + 1/2*hg19$length[chr])
	panel.abline(h = 0)#, lty = 3)
	panel.abline(v = cumLen[1:23], col = colLines, lty = 2)
	ltext(0, y, labels = "chr", cex = 0.75, col = colText)
	ltext(cumLen[1]/2, y, labels = 1, cex = 0.75, col = colText)
	for(i in 2:length(cumCentr)){
		x <- (hg19$length[i]/2 + cumLen[i-1])
		ltext(x, y, labels = i, cex = 0.75, col = colText)
		}
	}

addRunMed <- function(genomicPos, LogRatio){
	if(any(is.na(LogRatio))) LogRatio = LogRatio[!is.na(LogRatio)]
	Samp = seq(1, length(LogRatio), len = length(LogRatio)/10)
	g = genomicPos[Samp]; rmed = runmed(LogRatio[Samp], k = length(LogRatio)/1500)
	llines(g, rmed, col = 'black', lwd = 0.5)
	}

locateProbes <- function(genomicPos, LogRatio){
	if(any(is.na(LogRatio))) LogRatio = LogRatio[!is.na(LogRatio)]
	Samp = seq(1, length(LogRatio), len = length(LogRatio)/20)
	g = genomicPos[Samp]; LR = LogRatio[Samp]
	lpoints(g, LR, pch = 19, cex = 0.1, col = 'grey80')
	#rgb(0.8, 0.8, 0.8, 0.5)
	}

locateSegments <- function(genomicPos, Segments, thresh){
	Samp = seq(1, length(Segments), len = length(Segments)/10)
	g = genomicPos[Samp]; s = Segments[Samp]
	lpoints(g, s, pch = 19, cex = 0.25, col = ifelse(s>thresh, 'dodgerblue3', ifelse(s < -thresh , 'red3', 'grey50')))
	}

setMethod('createProfile', 'cghObj',
function(object, gene = NULL, Thresh = 0.1){
	cnSet = getCNset(object)
	lr = cnSet$Log2Ratio
	g = cnSet$genomicPos
	s = cnSet$Segm
	Title = getInfo(object, 'sampleId')
	subTitle = paste(getInfo(object, 'platform'), getInfo(object, 'gridGenomicBuild'), getInfo(object, 'analyseDate'), sep = ' - ')
	xRange = range(g, na.rm = TRUE); yRange = c(-1.5, 1.5) 
	object@gProfile <-xyplot(yRange~xRange, type = 'n',
											xlab = 'Genomic Position', ylab = 'Log2Ratio',
											main = Title, sub = subTitle,
											panel = function(){panel.xyplot(xRange, yRange, type = 'n')
																		locateChr(max(yRange))								# helper function
																		locateProbes(g, lr)										# helper function
																		addRunMed(g, lr)										# helper function
																		locateSegments(g, s, Thresh)						# helper function
																		if(!is.null(gene)){
																			geneInfo <- geneOfInt(object, gene)
																			x = geneInfo$genomicStart
																			if(!is.na(x)){
																				y0 = geneInfo$Log2Ratio
																				y1 = min(1.5, max(y0, 0.80))
																				label = as.character(geneInfo$NomenclatureSymbol)
																				ltext(x, y1, labels = label, cex = 1, font = 2)
																				lsegments(x0 = x, y0 = y0, x1 = x, y1 = y1*0.85, lwd = 2)
																				lsegments(x0 = x-1e8, y0 = y1*0.9, x1 = x+1e8, y1 = y1*0.9,
													 											lty = 3, lwd = 2, col = ifelse(y0>Thresh, 'dodgerblue3', ifelse(y0<(-Thresh), 'red3', 'grey')))
																				}
																			else cat('No location for', gene, '\n')
																			}
																		}
			)
	return(object)
	}
)

# zoom in a from_to location
setGeneric('zoomIn', function(object,...) standardGeneric('zoomIn'))
setMethod('zoomIn', 'cghObj', function(object, from = NA, to = NA){
	X11(type = 'dbcairo')
	segTable = getSegTable(object)
	index <- which(segTable$loc.start>=from & segTable$loc.end <= to)
	minY <- min(min(segTable$seg.med[index], na.rm = TRUE), -1)
	maxY <- max(max(segTable$seg.med[index], na.rm = TRUE), -1)
	gPlot <- getProfile(object)
	gPlot + 	coord_cartesian(xlim = range(from, to), ylim = range(minY, 1.5))
	}
)

