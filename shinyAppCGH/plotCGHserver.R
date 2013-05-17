# library(shiny)

cnSet = getCNset(object5)
segTable = getSegTable(object5)
# call synId here

locateChr <- function(y){
	colText = 'grey40'
	colLines = 'grey80'
	arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
	hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
	cumLen = cumsum(as.numeric(hg19$length))
	cumCentr <- 1/2*cumLen[1]
	for(chr in 2:length(cumLen)) cumCentr = c(cumCentr, cumLen[chr-1] + 1/2*cumLen[chr])
	abline(h = 0)#, lty = 3)
	abline(v = cumLen[1:23], col = colLines, lwd = 3, lty = 2)
	text(0, y, labels = "chr", cex = 1.25, col = colText)
	text(cumLen[1]/2, y, labels = 1, cex = 1.25, col = colText)
	for(i in 2:23){
		x <- (hg19$length[i]/2 + cumLen[i-1])
		text(x, y, labels = i, cex = 1.25, col = colText)
		}
	}

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {

  # Function that generates a plot of the genomic profile and provide the information about specified genes.
  # It's updated each time a new gene symbol is called.
 	# createPlot is the reactive function called by renderPlot() to generate the main plot.
  	# createTable is the reactive function called by renderTable()
 	
 	createPlot <- reactive({
		Thresh = log2(2.2/2)
		x <- cnSet$genomicPos[cnSet$ChrNum != 24]
		y <- cnSet$Log2Ratio[cnSet$ChrNum != 24]
  		n = length(y)
  		Samp = seq(1, n, len = n/20)
  		plot(x[Samp], runmed(y[Samp], k = 11),
  			ylim = range(-1.5, 1.5), cex = 0.1, col = 'grey30',
  			cex.axis = 1, cex.lab = 1.5, las = 1, mar = c(10, 10, 10, 10), mgp = c(3, 1, 0),
  			xlab = 'Genomic position', ylab = 'Log2Ratio')#, main = 'Genomic profile')
  		segCols <- ifelse(segTable$seg.med <= -Thresh, 'red3', ifelse(segTable$seg.med >= Thresh, 	myBlue <- rgb(0, 0.45, 1, 1), 'grey75'))
  		lines(x[Samp], runmed(y[Samp], k = 501))
  		segments(x0 = segTable$loc.start, y0 = segTable$seg.med, x1 = segTable$loc.end, lwd = 5, col = segCols)
  		locateChr(1.5)
  		gene <- input$PredefSymb
  		if(input$OtherSymb != '') gene <- input$OtherSymb
  		if(gene != 'NONE' & gene != ''){
  			tmp <- try(geneOfInt(object5, gene), silent = TRUE)
  				if(class(tmp) != 'try-error'){
	  			segments(x0 = tmp$genomicStart, y0 = tmp$Log2Ratio, y1 = 1, lwd = 3)
	  			text(max(2.5e8, tmp$genomicStart), 1.1, labels = paste0(gene, ' (Log2Ratio = ', round(tmp$Log2Ratio, 3), ')'), cex = 1.5, font = 2)
	  			}
     		}
	})
	output$Profile <- renderPlot({createPlot()})
	
	createTable <- reactive({
		gene <- input$PredefSymb
  		if(input$OtherSymb != '') gene = input$OtherSymb
  		if(gene != 'NONE' & gene != ''){
			tmp <- try(geneOfInt(object5, gene), silent = TRUE)
  			if(class(tmp) != 'try-error')
				data.frame(Symbol = tmp$Symbol, entrezgene = tmp$entrezgene, Description = tmp$Description, Cytoband = tmp$Cytoband, Log2Ratio = round(tmp$Log2Ratio, 3))
			else data.frame(Symbol = 'Not a valid symbol')
			}
		})
		
	output$geneTable <- renderTable({createTable()})
	
	# geneList = c("CCND1", "ALK", "MDM2", "FRS2", "MET", "RPTOR", "ESR1", "PGR", "FGFR1", "FGFR2",
					# "MYC", "FGF4", "FGF9", "EGFR", "ERBB2", "TOP2A", "IGF1", "IGF1R", "BRCA1", "BRCA2",
					# "NOTCH4", "VEGFA", "PTEN", "PIK3CB", "PAK1")
	# output$summaryTable <- renderTable({geneOfInt(object5, geneList)[,-c(1, 5:6, 8:12)]})
	}
)


#setwd('/Users/fredcommo/Documents/Stats/Docs R/shinyApps/')
#runApp('./shinyAppCGH/')