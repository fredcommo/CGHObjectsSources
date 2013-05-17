
myPlot <- function(object, gene = NULL, Thresh = 0.1){
	cnSet = getCNset(object5)
	lr = cnSet$Log2Ratio
	g = cnSet$genomicPos
	s = cnSet$Segm
	Title = getInfo(object5, 'sampleId')
	subTitle = paste(getInfo(object5, 'platform'), getInfo(object5, 'gridGenomicBuild'), getInfo(object5, 'analyseDate'), sep = ' - ')
	xRange = range(g, na.rm = TRUE); yRange = c(-1.5, 1.5) 
	xyplot(yRange~xRange, type = 'n',
			xlab = 'Genomic Position', ylab = 'Log2Ratio',
			main = Title, sub = subTitle,
			panel = function(){panel.xyplot(xRange, yRange, type = 'n')
										locateChr(max(yRange))
										locateProbes(g, lr)
										addRunMed(g, lr)
										locateSegments(g, s, Thresh)
										if(!is.null(gene)){
											geneInfo <- geneOfInt(object5, gene)
											x = geneInfo$genomicStart
											if(!is.na(x)){
												y0 = geneInfo$Log2Ratio
												y1 = min(1.5, max(y0, 0.75))
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
}


myPlot(object5, 'PIK3CB')
myPlot(object5, 'rototo')
myPlot(object5, 'EGFR')
myPlot(object5, 'CCND1')
myPlot(object5, 'MYC')




