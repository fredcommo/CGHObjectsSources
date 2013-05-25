# library(shiny)

require(synapseClient)
# Load class def & accessors
cat('Loading functions...')
e1  <- loadEntity('syn1864348')
source(file.path(e1$cacheDir, e1$files), echo = FALSE, verbose = FALSE)
e2  <- loadEntity('syn1864353')
source(file.path(e2$cacheDir, e2$files), echo = FALSE, verbose = FALSE)
e3  <- loadEntity('syn1864359')
source(file.path(e3$cacheDir, e3$files), echo = FALSE, verbose = FALSE)

# Load annot tables
e <- synGet('syn1877556')
HG19 <- read.csv(e@filePath, header = TRUE, sep = '\t')
e <- synGet('syn1877601')
geneDB <- readRDS(e@filePath)

# Load data
cat('\n\nLoading data...')
cgh  <- loadEntity('syn1864342')
load(file.path(cgh$cacheDir, cgh$files))
cnSet = getCNset(cghProfile)
segTable = getSegTable(cghProfile)
cat('\n\n')

###########################
# Helper functions
locateChr <- function(y, hg19 = HG19){
	colText = 'grey40'
	colLines = 'grey80'
#	arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
#	hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
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

geneOfInt <- function(object, geneList, DB = geneDB){
            segTable = getSegTable(object)
            output <- lapply(geneList, function(gene){
              tmp <- DB[which(DB$Symbol == gene),]
              startLoc <- tmp$genomicStart
              stopLoc <- tmp$genomicStop
              containStart <- which(segTable$loc.start<=startLoc & segTable$loc.end>=startLoc)
              containStop <- which(segTable$loc.start<=stopLoc & segTable$loc.end>=stopLoc)
              if( is.na(startLoc) | is.na(stopLoc) )
                tmp <- cbind(tmp, Log2Ratio = NA)
              else if(containStart == containStop)
                tmp <- cbind(tmp, Log2Ratio = segTable$seg.med[containStart])
              else{
                Log2Ratio <- segTable$seg.med[union(containStart, containStop)]
                tmp <- cbind(tmp, Log2Ratio)
              }
              return(tmp)
            })
            output <- do.call(rbind, output)
            rownames(output) = seq(1, nrow(output))
            return(output)
          }
# End helper functions
###########################

#
shinyServer(function(input, output) {

  # Function that generates a plot of the genomic profile and provide the information about specified genes.
  # It's updated each time a new gene symbol is called.
 	# createPlot is the reactive function called by renderPlot() to generate the main plot.
  	# createTable is the reactive function called by renderTable()

#  	access <- reactive({
#      if(input$userId != 'test') stop('Access denied')
#  	  })
#   output$Access <- renderText({access()})
     
 	createPlot <- reactive({
 	  myBlue <- rgb(0, 0.475, 1, 1)
    Up = log2(1.15) # gain = 25%
 	  Lo = log2(0.90)  # lost = 20%
 	  x <- cnSet$genomicPos[cnSet$ChrNum != 24]
		y <- cnSet$Log2Ratio[cnSet$ChrNum != 24]
  	n <- length(y)
  	Samp <- seq(1, n, len = n/40)
  	plot(x[Samp], runmed(y[Samp], k = 11),
         ylim = range(-1.5, 1.5), cex = 0.1, col = 'grey10',
         cex.axis = 1, cex.lab = 1.5, las = 1, mar = c(10, 10, 10, 10), mgp = c(3, 1, 0), cex.main = 1.5,
         xlab = 'Genomic position', ylab = 'Log2Ratio',
         main = paste(getInfo(cghProfile, 'sampleId'), '\n',
                     'Gain threshold:', round(Up, 3), '- Lost threshold:', round(Lo, 3)))
  	segCols <- ifelse(segTable$seg.med <= Lo, 'red3',
                        ifelse(segTable$seg.med >= Up, myBlue, 'grey60'))
  	lines(x[Samp], runmed(y[Samp], k = 201))
  	segments(x0 = segTable$loc.start, y0 = segTable$seg.med, x1 = segTable$loc.end, lwd = 5, col = segCols)
  	locateChr(1.5)
  	gene <- toupper(input$PredefSymb)
  	if(input$OtherSymb != '') gene <- toupper(input$OtherSymb)
  	if(gene != 'NONE' & gene != ''){
  		tmp <- try(geneOfInt(cghProfile, gene), silent = TRUE)
  		if(class(tmp) != 'try-error'){
          lr <- tmp$Log2Ratio
          Col <- ifelse(lr <= Lo, 'red3', ifelse(lr >= Up, myBlue, 'grey30'))
  	  		segments(x0 = tmp$genomicStart, y0 = lr, y1 = 0.95*lr/abs(lr), lwd = 3, col = Col)
	    		text(max(2.5e8, tmp$genomicStart), 1.15*lr/abs(lr), labels = paste0(gene, '\n(Log2Ratio = ', round(tmp$Log2Ratio, 3), ')'), cex = 1.5, font = 2)
	  	  	}
     		}
	})
  output$Profile <- renderPlot({createPlot()})
	
	createTable <- reactive({
		gene <- toupper(input$PredefSymb)
  	if(input$OtherSymb != '') gene = toupper(input$OtherSymb)
		if(gene == 'NONE' | gene == '')
		  data.frame(Symbol = 'Select a gene in the "Predefined list" or enter a valid symbol in "Other symbol",
                 \n then click on "Update view"')
    else{
			tmp <- try(geneOfInt(cghProfile, gene), silent = TRUE)
  		if(class(tmp) != 'try-error')
				data.frame(Symbol = tmp$Symbol, entrezgene = tmp$entrezgene, Description = tmp$Description, Cytoband = tmp$Cytoband, Log2Ratio = round(tmp$Log2Ratio, 3))
			else data.frame(Symbol = paste0('Too bad!! "', gene,'" is not a valid symbol'))
			}
		})
 	output$geneTable <- renderTable({createTable()})
	}
)


#setwd('/Users/fredcommo/Documents/Stats/Docs R/shinyApps/')
#runApp('./shinyAppCGH/')

# geneList = c("CCND1", "ALK", "MDM2", "FRS2", "MET", "RPTOR", "ESR1", "PGR", "FGFR1", "FGFR2",
# "MYC", "FGF4", "FGF9", "EGFR", "ERBB2", "TOP2A", "IGF1", "IGF1R", "BRCA1", "BRCA2",
# "NOTCH4", "VEGFA", "PTEN", "PIK3CB", "PAK1")
# output$summaryTable <- renderTable({geneOfInt(object5, geneList)[,-c(1, 5:6, 8:12)]})
