###############################
# helper fuctions
require(ggplot2)
.plotProfile <- function(object, ylim = range(-1.5, 1.5),...) {
  object$gProfile + coord_cartesian(ylim = ylim)
}

.getValue <- function(object, geneDB, tag){
  segTable = object$segTable
  tmp <- geneDB[which(geneDB$Symbol %in% tag),]
  # LR <- segTable$seg.med[which(segTable$loc.start<= tmp$genomicStart & segTable$loc.end >= tmp$genomicStop)]
  LR <- lapply(1:nrow(tmp), function(i){
    startAt <- tmp$genomicStart[i]
    stopAt <- tmp$genomicStop[i]
    segTable$seg.med[which(segTable$loc.start<= startAt & segTable$loc.end >= stopAt)]
  })
  tmp$LR <- do.call(c, LR)
  return(tmp)
}

.tagMyGene <- function(object, geneDB, tag = NULL, ylim = range(-1.5, 1.5),
                      gain = log2(2.25/2), loss = log2(1.80/2)){
  myBlue <- rgb(0, 0.45, 1, 1)
  if(is.null(tag)) .plotProfile(object)
  else if(sum(is.element(tag, geneDB$Symbol)) == 0){
    cat('Not a valid gene symbol(s)\n')
    .plotProfile(object)
  }
  else {
    tmp <- .getValue(object, geneDB, tag)
    Col <- ifelse(tmp$LR<= loss, 'red3', ifelse(tmp$LR>= gain, myBlue, 'grey40'))
    dat <- data.frame(xstart = c(max(1e8, tmp$genomicStart-2.5e7), tmp$genomicStart),
                      xend = c(tmp$genomicStart, tmp$genomicStart),
                      ystart = c(tmp$LR/abs(tmp$LR)*1.2, tmp$LR/abs(tmp$LR)*1.2),
                      yend = c(tmp$LR/abs(tmp$LR)*1.2, tmp$LR),
                      Col = rep(Col, 2))
    object$gProfile +
      coord_cartesian(ylim = ylim)+
      annotate("text", x = max(2.5e8, tmp$genomicStart - 3.2e8), y = tmp$LR/abs(tmp$LR)*1.2,
               label = paste0(tmp$Symbol, '\n(Log2R = ', round(tmp$LR, 3), ')'), cex = 6) +
      geom_segment(data = dat, aes(x = xstart, y = ystart, xend = xend, yend = yend), colour = as.character(dat$Col), size = 1)
  }
}

# End helper functions
###########################
# file <- File('/Users/fredcommo/Documents/Projet Safir/CGHObjectsSources/shinyHelperFunctions.R',
#              parentId = 'syn1864121')
# file <- synStore(file)
