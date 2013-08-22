# Build Nimblegene CGH object from Synapse.


ent <- synGet('syn1680110')
rawDat <- strsplit(readLines(file.path(ent$cacheDir, ent$files), n = 1), '\\t')

buildNumble <- function(synId){
  # Load cghData from a synapse entity and build an AgilentObject.
  entity <- synGet('syn1680110') # synGet(synId)
  
  # Need to create a NimbleObj class!
  object <- NimbleObj(info = c(fileName = propertyValue(entity, 'name'),
                                synapseId = propertyValue(entity, 'id'),
                                platform = 'Nimblegene'))
  object@info <- c(object@info, .readNimbleInfo(entity))
  object@cnSet <- .readNimbleMatrix(entity)
  object <- suppressFlags(object)
  object <- suppressDuplic(object)
  object <- preset(object)
  return (object)
}
.readNimbleInfo(entity) <- function(entity){
  row1 <- strsplit(readLines(file.path(entity$cacheDir, entity$files), n = 1), '\\t')
  row1 <- unlist(row1)
  labId = "Broad"
  barCode = gsub('(.)*/', '', row1[grep('imagefile', row1)])
  gridName = gsub('(.)*=', '', row1[grep('designname', row1)])
  scanDate = gsub('(.)*=', '', row1[grep('date', row1)])
  arrayProtocol = gsub('(.)*=', '', row1[grep('designid', row1)])
  gridGenomicBuild = gsub('(.)*=', '', row1[grep('designname', row1)])
  ref = 'Dual color hybridization'
  #	ScanName = as.character(a.info[2, which(tmpNames == "Scan_ScannerName")])
  #	version = as.character(a.info[2, which(tmpNames == "FeatureExtractor_Version")])
  cat('\tDone.\n')
  return(c(labId = labId, barCode = barCode, gridName = gridName,
           scanDate = scanDate, arrayProtocol = arrayProtocol,
           gridGenomicBuild = gridGenomicBuild, ref = ref,
           analyseDate = format(Sys.Date(), "%m-%d-%Y"))
  )
  
}
.readNimbleMatrix <- function(entity){
  cat('Reading values...')
  cnSet <- read.csv(file.path(entity$cacheDir, entity$files),
                    header = T, skip = 1, sep = "\t", stringsAsFactors = FALSE)
  cat('\tDone.\n')
  cnSet <- .curateNimbleCnSet(cnSet)
  return(cnSet)
}
.curateNimbleCnSet <- function(cnSet){
  colNames <- c(	"ProbeName", "SystematicName",
     "gMedianSignal", "rMedianSignal",
     "gIsSaturated", "rIsSaturated",
     "gIsFeatNonUnifOL", "rIsFeatNonUnifOL",
     "gIsWellAboveBG", "rIsWellAboveBG")
  )															
}
