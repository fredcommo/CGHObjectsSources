# library(shiny)
require(SOAR)
Remove(shinyData, gPos, l2r, segTable, currentId, currentSample)
require(synapseClient)

###########################
# Helper functions

.locateChr <- function(y, hg19 = HG19){
  cat('LocateChr\n')
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
  cat('geneOfInt\n')
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

.mainPlot <- function(gPos, l2r, Up, Lo){
  cat('mainPlot\n')
  plot(gPos, runmed(l2r, k = 11),
       ylim = range(-1.5, 1.5), cex = 0.1, col = 'grey20',
       cex.axis = 1, cex.lab = 1.5, las = 1, mar = c(10, 10, 10, 10), mgp = c(3, 1, 0), cex.main = 1.5,
       xlab = 'Genomic position', ylab = 'Log2Ratio',
       main = paste(currentSample, ' - ', currentId, '\nGain threshold:', round(Up, 3), '- Lost threshold:', round(Lo, 3)))
  lines(gPos, runmed(l2r, k = 81))
}

.addSegments <- function(segTable, Up, Lo, lossCol, normCol, gainCol){
  cat('addSegments\n')
  segCols <- ifelse(segTable$seg.med <= Lo, lossCol,
                    ifelse(segTable$seg.med >= Up, gainCol, normCol))
  segments(x0 = segTable$loc.start, y0 = segTable$seg.med,
           x1 = segTable$loc.end, lwd = 5, col = segCols)
}
  
.addLabel <- function(tmp, Up, Lo, lossCol, normCol, gainCol){
  cat('addLabel\n')
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

.getData <- function(synId){
  cat('getData\n')
  cgh <- synGet(synId)
  shinyData <- get(load(cgh@filePath))
  segTable <- getSegTable(shinyData) #@segTable
  cnSet <- getCNset(shinyData) #@cnSet
  gPos <- cnSet$genomicPos[cnSet$ChrNum %in% 1:23]
  l2r <- cnSet$Log2Ratio[cnSet$ChrNum %in% 1:23]
  currentSample <- gsub('_syn(.)*', '', c(propertyValue(cgh, 'name')))
  Samp <- seq(1, length(l2r), len = 20e3)
  gPos <- gPos[Samp]
  l2r <- l2r[Samp]
  if(any(is.na(l2r))){
    NAs <- which(is.na(l2r))
    l2r <- l2r[-NAs]
    gPos <- gPos[-NAs]
  }
  segTable <- segTable[which(segTable$chrom %in% 1:23),]
  currentId <- synId
  StoreData(cgh, gPos, l2r, segTable, currentSample, currentId)
}

# End helper functions
###########################

# Load annot tables
e <- synGet('syn1877556')
HG19 <- read.csv(e@filePath, header = TRUE, sep = '\t')
e <- synGet('syn1877601')
geneDB <- readRDS(e@filePath)

# Load workflow
if(!exists('cghObj')){
  workflow <- synGet('syn2128342')
  source(workflow@filePath)
  }

#
shinyServer(function(input, output) {

  # Function that display the genomic profile and provide the information about specified genes.
  # Plot and table are updated each time a new gene symbol is called.
 	# createPlot is the reactive function called by renderPlot() to generate the main plot.
  # createTable is the reactive function called by renderTable()
  
 	createPlot <- reactive({
 	  if(!exists('currentId') | (exists('currentId') && currentId != input$Samp)){
       cat('new data\n')
       #Remove(shinyData, gPos, l2r, segTable, currentId, currentSample)
       .getData(input$Samp)
 	      }
 	  cat(currentSample, '\t', currentId, '\n')
 	  output$Id <- renderText(paste(currentSample, input$Samp, sep = ' - '))
 	  
    Up = log2(1.15); Lo = log2(0.90) # gain = 15%, lost = 10%
    gainCol <- rgb(0, 0.475, 1, 1); lossCol <- 'red3'; normCol <- 'grey60'
    .mainPlot(gPos, l2r, Up, Lo)
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

  output$Profile <- renderPlot({createPlot()})
  output$geneTable <- renderTable({createTable()})
}
)


#setwd('/Users/fredcommo/Documents/Stats/Docs R/shinyApps/')
#runApp('./shinyAppCGH/')

# geneList = c("CCND1", "ALK", "MDM2", "FRS2", "MET", "RPTOR", "ESR1", "PGR", "FGFR1", "FGFR2",
# "MYC", "FGF4", "FGF9", "EGFR", "ERBB2", "TOP2A", "IGF1", "IGF1R", "BRCA1", "BRCA2",
# "NOTCH4", "VEGFA", "PTEN", "PIK3CB", "PAK1")
# output$summaryTable <- renderTable({geneOfInt(object5, geneList)[,-c(1, 5:6, 8:12)]})
