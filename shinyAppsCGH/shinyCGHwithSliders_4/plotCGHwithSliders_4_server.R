# library(shiny)
# require(SOAR)
# Remove(shinyData, gPos, l2r, segTable, currentId, currentSample)
require(synapseClient)

###########################
# Helper functions
.getData <- function(synId){
  cat('getData\n')
  cgh <- synGet(synId)
  shinyData <- read.csv(cgh@filePath, header = TRUE, sep = '\t')
  chrSelect <- which(shinyData$chrNum %in% 1:23)
  shinyData <- shinyData[chrSelect,]
  return(shinyData[chrSelect,])
}
.mainPlot <- function(gPos, l2r, inputSamp){
  cat('mainPlot...\n')
  set.seed(12345)
  Samp <- sort(sample(1:length(gPos), 15e3))
  gPos <- gPos[Samp]
  l2r <- l2r[Samp]
  plot(gPos, runmed(l2r, k = 5),
       ylim = range(-1.5, 1.5), cex = 0.1, col = 'grey20',
       cex.axis = 1, cex.lab = 1.5, las = 1, mar = c(10, 10, 10, 10), mgp = c(3, 1, 0), cex.main = 1.5,
       xlab = 'Genomic position', ylab = 'Log2Ratio',
       main = gsub('-', ' - ',inputSamp), cex.main = 2)
  lines(gPos, runmed(l2r, k = 81))
  cat('Any NA:', any(is.na(gPos)), '\t', any(is.na(l2r)), '\n')
}
.addSegments <- function(segTable, Up, Lo, lossCol, normCol, gainCol){
  cat('addSegments...\n')
  cat('Up:', Up, 'Lo:', Lo, '\n')
  segCols <- ifelse(segTable$seg <= Lo, lossCol, ifelse(segTable$seg >= Up, gainCol, normCol))
  segments(x0 = segTable$Start, y0 = segTable$seg, x1 = segTable$Stop, lwd = 5, col = segCols)
}
.makeSegTable <- function(loc, seg){
  cat('Making seg table...\n')
  .diff <- diff(seg, 1)
  idx <- c(1, which(.diff != 0)+1, length(seg))
  idx <- embed(idx, 2)
  Start <- loc[idx[,2]]
  seg <- seg[idx[,2]]
  Stop <- loc[idx[,1]-1]
  out <- cbind.data.frame(Start = Start, Stop = Stop, seg = seg)
  return(cbind.data.frame(Start = Start, Stop = Stop, seg = seg))
}
# # Fix errors in .makeSegTable
# tmp <- try(.geneOfInt(Values$seg, 'esr1'), silent = TRUE)
# s1 <- 1214550000
# s2 <- 1214970000
# Table <- Input$Table
# loc <- Table$genomicPos
# seg <- Table$Segm1
# idxloc <- which(loc>s1 & loc<s2)
# head(Table[idxloc, 7:15], n = 20)
# tmp

.locateChr <- function(y, hg19 = HG19){
  cat('LocateChr...\n')
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
.geneOfInt <- function(segTable, gene, DB = geneDB){
  cat('geneOfInt...\n')
    gene <- toupper(gene)
    tmp <- DB[which(DB$Symbol == gene),]
    geneStart <- tmp$genomicStart
    geneStop <- tmp$genomicStop
    bound1 <- max(which(segTable$Start <= geneStart))
    bound2 <- min(which(segTable$Stop >= geneStop))
    bounds <- unique(c(bound1, bound2))
    if(is.na(geneStart) | is.na(geneStop))
      tmp <- cbind(tmp, Log2Ratio = NA, segment.len = NA)
    else if(length(bounds) > 0){
      geneLR <- unique(segTable$seg[bounds])
      tmp <- cbind(tmp, Log2Ratio = geneLR, segment.len = abs(segTable$Stop[bounds] - segTable$Start[bounds]))
      }
    else{
      tmp <- cbind(tmp, Log2Ratio = NA, , segment.len = NA)
      }
  return(tmp)
}
.addLabel <- function(tmp, Up, Lo, lossCol, normCol, gainCol){
  cat('addLabel...:', nrow(tmp),'rows\n')
  lapply(1:nrow(tmp), function(i){
    geneValue <- tmp$Log2Ratio[i]
    symbol <- as.character(tmp$Symbol[i])
    Col <- ifelse(geneValue <= Lo, lossCol, ifelse(geneValue >= Up, gainCol, normCol))
    x0 <- c(max(1e8, tmp$genomicStart[i]-2.5e7), tmp$genomicStart[i])
    x1 = c(tmp$genomicStart[i], tmp$genomicStart[i])
    y0 = c(geneValue/abs(geneValue)*1.2, geneValue/abs(geneValue)*1.2)
    y1 = c(geneValue/abs(geneValue)*1.2, geneValue)
    Col = rep(Col, 2)
    segments(x0, y0, x1, y1, col = Col, lwd = 3.5)
    text( x = max(2.5e8, tmp$genomicStart - 2.0e8), y = geneValue/abs(geneValue)*1.25,
        labels = c(symbol, paste0('\n\n(Log2R = ', round(geneValue, 3), ')')), cex = c(1.5, 1.25), font = 2)
    }
  )
}
.dlrs <- function(x){
  nx <- length(x)
  if (nx<3) {
    stop("Vector length>2 needed for computation")
  }
  tmp <- embed(x,2)
  diffs <- tmp[,2]-tmp[,1]
  dlrs <- IQR(diffs, na.rm = TRUE)/(sqrt(2)*1.34)
  return(dlrs)
}
.defStatus <- function(dLRs){
  status <- ifelse(dLRs<=.1, 'Excellent (dLRs < 0.1)',
                   ifelse(dLRs>.1 & dLRs<=.2, 'Good (0.1 < dLRs < 0.2)',
                          ifelse(dLRs>.2 & dLRs<=.3, 'Poor (0.2 < dLRs < 0.3)', 'Bad (dLRs > .3)')))
  return(status)
}

# End helper functions
###########################

# Load annot tables
e <- synGet('syn1877556')
HG19 <- read.csv(e@filePath, header = TRUE, sep = '\t')
e <- synGet('syn1877601')
geneDB <- readRDS(e@filePath)

# Colors
gainCol <- rgb(0, 0.475, 1, 1); lossCol <- 'red3'; normCol <- 'grey60'

#
shinyServer(function(input, output) {
  # Function that display the genomic profile and provide the information about specified genes.
  # Plot and table are updated each time a new gene symbol is called.
 	# createPlot is the reactive function called by renderPlot() to generate the main plot.
  # createTable is the reactive function called by renderTable()
  
 	createPlot <- reactive({
    .mainPlot(Values$gPos, Values$l2r, input$Samp)
    .addSegments(Values$seg, input$gain, input$loss, lossCol, normCol, gainCol)
    .locateChr(1.5) # yvalue: what heigh to write
 	  if(Gene$current != 'NONE'){
  		tmp <- try(.geneOfInt(Values$seg, Gene$current), silent = TRUE)
  		if(class(tmp) != 'try-error')
  		  .addLabel(tmp, input$gain, input$loss, lossCol, 'grey30', gainCol)
  	  }
	})
   
	createTable <- reactive({		
	  if(Gene$current == 'NONE' | Gene$current == '')
		  data.frame(Symbol = 'Select a gene in the "Predefined gene list" or enter a valid symbol in "Other symbol",
                 \nthen click on "Update view"')
    else{
#      segTable <- cbind.data.frame(loc = Values$gPos, seg = Values$seg)
			tmp <- try(.geneOfInt(Values$seg, Gene$current), silent = TRUE)
  		if(class(tmp) != 'try-error')
				data.frame(Symbol = tmp$Symbol, entrezgene = tmp$entrezgene,
                   Description = tmp$Description, Cytoband = tmp$Cytoband,
				          Segment.length = tmp$segment.len, Log2Ratio = round(tmp$Log2Ratio, 3))
			else data.frame(Symbol = paste0('Too bad!! "', Gene$current,'" is not a valid symbol'))
			}
		})

  Input <- reactiveValues(synId = character(),
                          Table = data.frame()
                          )
  observe({
    Input$synId <- gsub('(.)*-', '', input$Samp)
    Input$Table <- .getData(Input$synId)
    })

  Values <- reactiveValues(
    Samp = numeric(),
    l2r = numeric(),
    gPos = numeric(),
    seg = data.frame(),
    dlrs = numeric(),
    MAD = numeric(),
    status = character()
    )
  observe({
    Values$gPos <- Input$Table$genomicPos
    Values$l2r <- switch(input$modtype,
                         ncp = (Input$Table$Log2Ratio1),
                         lcp = (Input$Table$Log2Ratio2),
                         zcp = (Input$Table$Log2Ratio3),
                         ccp = (Input$Table$Log2Ratio4))
    Values$seg <- switch(input$modtype,
                              ncp = .makeSegTable(Input$Table$genomicPos,Input$Table$Segm1),
                              lcp = .makeSegTable(Input$Table$genomicPos,Input$Table$Segm2),
                              zcp = .makeSegTable(Input$Table$genomicPos,Input$Table$Segm3),
                              ccp = .makeSegTable(Input$Table$genomicPos,Input$Table$Segm4)
                              )
    Values$dlrs <- round(.dlrs(Values$l2r), 3)
    Values$MAD <- NA
    Values$status <- .defStatus(Values$dlrs)
    })
   
 	Gene <- reactiveValues(current = character())
 	observe({Gene$current <- toupper(input$PredefSymb)})
 	observe({Gene$current <- toupper(input$OtherSymb)})
 
#  output$Id <- renderText({paste(Input$SampleId, Input$synId, sep = ' - ')})
 	output$Id <- renderText({input$Samp})
 	output$dlrs <- renderText({sprintf('dLRs: %s', Values$dlrs)})
 	output$mad <- renderText({sprintf('MAD: %s', Values$MAD)})
 	output$status <- renderText({sprintf('QC status: %s', Values$status)})
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
