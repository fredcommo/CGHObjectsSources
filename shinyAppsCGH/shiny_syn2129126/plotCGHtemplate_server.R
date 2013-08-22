require(synapseClient)
cgh  <- synGet('syn2129126')

require(synapseClient)
cgh  <- synGet('syn2129126')

require(synapseClient)
cgh  <- synGet('syn2129126')

# library(shiny)

###########################
# Helper functions
.locateChr <- function(y, hg19 = HG19){
  colText = 'grey40'
  colLines = 'grey80'
  cumLen = cumsum(as.numeric(hg19$length))
  cumCentr <- 1/2*cumLen[1]
  for(chr in 2:length(cumLen)) cumCentr = c(cumCentr, cumLen[chr-1] + 1/2*cumLen[chr])
  abline(h = 0)#, lty = 3)
  abline(v = cumLen[1:23], col = colLines, lwd = 3, lty = 2)
  text(0, y, labels = "chr", cex = 1.1, col = colText)
  text(cumLen[1]/2, y, labels = 1, cex = 1.1, col = colText)
  for(i in 2:23){
    x <- (hg19$length[i]/2 + cumLen[i-1])
    text(x, y, labels = i, cex = 1.1, col = colText)
  }
}

.geneOfInt <- function(segTable, geneList, DB = geneDB){
  output <- lapply(geneList, function(gene){
    gene <- toupper(gene)
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

.mainPlot <- function(gPos, lr, Up, Lo){
  plot(gPos, runmed(lr, k = 11),
       ylim = range(-1.5, 1.5), cex = 0.1, col = 'grey20',
       cex.axis = 1, cex.lab = 1.5, las = 1, mar = c(10, 10, 10, 10), mgp = c(3, 1, 0), cex.main = 1.5,
       xlab = 'Genomic position', ylab = 'Log2Ratio',
       main = paste(synId, '\nGain threshold:', round(Up, 3), '- Lost threshold:', round(Lo, 3)))
  lines(gPos, runmed(lr, k = 81))
}

.addSegments <- function(segTable, Up, Lo, lossCol, normCol, gainCol){
  segCols <- ifelse(segTable$seg.med <= Lo, lossCol,
                    ifelse(segTable$seg.med >= Up, gainCol, normCol))
  segments(x0 = segTable$loc.start, y0 = segTable$seg.med,
           x1 = segTable$loc.end, lwd = 5, col = segCols)
}
  
.addLabel <- function(tmp, Up, Lo, lossCol, normCol, gainCol){
  geneValue <- tmp$Log2Ratio
  symbol <- as.character(tmp$Symbol)
  Col <- ifelse(geneValue <= Lo, lossCol, ifelse(geneValue >= Up, gainCol, normCol))
  x0 <- c(max(1e8, tmp$genomicStart-2.5e7), tmp$genomicStart)
  x1 = c(tmp$genomicStart, tmp$genomicStart)
  y0 = c(geneValue/abs(geneValue)*1.2, geneValue/abs(geneValue)*1.2)
  y1 = c(geneValue/abs(geneValue)*1.2, geneValue)
  Col = rep(Col, 2)
  segments(x0, y0, x1, y1, col = Col, lwd = 3.5)
  text( x = max(2.5e8, tmp$genomicStart - 2.0e8), y = geneValue/abs(geneValue)*1.25,
      #  labels = paste0(tmp$Symbol, '\n(Log2Ratio = ', round(geneValue, 3), ')'), cex = 1.5, font = 2)
        labels = c(symbol, paste0('\n\n(Log2R = ', round(geneValue, 3), ')')), cex = c(1.5, 1.25), font = 2)
}
# End helper functions
###########################


# Load annot tables
e <- synGet('syn1877556')
HG19 <- read.csv(e@filePath, header = TRUE, sep = '\t')
e <- synGet('syn1877601')
geneDB <- readRDS(e@filePath)

# Load data
cat('\nLoading data...')
shinyData <- get(load(cgh@filePath))
segTable = shinyData$segTable
gPos <- shinyData$gPos
lr <- shinyData$L2R
fileName <- shinyData$fileName
synId <- shinyData$synId
cat('\n\n')


#
shinyServer(function(input, output) {

  # Function that display the genomic profile and provide the information about specified genes.
  # Plot and table are updated each time a new gene symbol is called.
 	# createPlot is the reactive function called by renderPlot() to generate the main plot.
  # createTable is the reactive function called by renderTable()

  output$Id <- renderText(paste0(gsub('.txt.bz2', '', fileName), ' - ', synId))
  
 	createPlot <- reactive({
    Up = log2(1.15) # gain = 15%
 	  Lo = log2(0.90)  # lost = 10%
    gainCol <- rgb(0, 0.475, 1, 1)
    lossCol <- 'red3'
    normCol <- 'grey60'
    .mainPlot(gPos, lr, Up, Lo)
    .addSegments(segTable, Up, Lo, lossCol, normCol, gainCol)
    .locateChr(1.5) # yvalue: what heigh to write
  	gene <- toupper(input$PredefSymb)
  	if(input$OtherSymb != '') gene <- toupper(input$OtherSymb)
  	if(gene != 'NONE' & gene != ''){
  		tmp <- try(.geneOfInt(segTable, gene), silent = TRUE)
  		if(class(tmp) != 'try-error')
  		  .addLabel(tmp, Up, Lo, lossCol, 'grey30', gainCol)
  	  }
	})
  output$Profile <- renderPlot({createPlot()})
	
	createTable <- reactive({
		gene <- toupper(input$PredefSymb)
  	if(input$OtherSymb != '') gene = toupper(input$OtherSymb)
		if(gene == 'NONE' | gene == '')
		  data.frame(Symbol = 'Select a gene in the "Predefined list" or enter a valid symbol in "Other symbol",
                 \nthen click on "Update view"')
    else{
			tmp <- try(.geneOfInt(segTable, gene), silent = TRUE)
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
