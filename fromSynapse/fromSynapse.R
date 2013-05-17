# Load entity function reads in the tab separated files
buildAgilent <- function(synId){
  # Load cghData from a synapse entity and build an AgilentObject.
  	entity <- loadEntity(synId)
	object <- AgilentObj(info = c(fileName = propertyValue(entity, 'name'), sampleId = propertyValue(entity, 'id'), platform = 'Agilent'))
	object@info <- c(object@info, .readAgilentInfo(entity))
	object@cnSet <- .readAgilentMatrix(entity)
	object <- suppressFlags(object)
	object <- suppressDuplic(object)
	object <- preset(object)
	return (object)
}

.readAgilentInfo <- function(entity){
	cat('Reading information...')
	a.info <- read.csv(file.path(entity$cacheDir, entity$files), header = F, fill = T, skip = 1, nrows = 8, sep = "\t")
	tmpNames <- as.vector(a.info[1,])
	labId = as.character(a.info[2, which(tmpNames == "FeatureExtractor_UserName")])
	barCode = as.character(a.info[2, which(tmpNames == "FeatureExtractor_Barcode")])
	gridName = as.character(a.info[2, which(tmpNames == "Grid_Name")])
	scanDate = as.character(a.info[2, which(tmpNames == "Scan_Date")])
	scanDate <- unlist(strsplit(scanDate, " "))[1]
	arrayProtocol = as.character(a.info[2, which(tmpNames == "Protocol_Name")])
	gridGenomicBuild = as.character(a.info[2, which(tmpNames == "Grid_GenomicBuild")])
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

.readAgilentMatrix <- function(entity){
	cat('Reading values...')
	cnSet <- read.csv(file.path(entity$cacheDir, entity$files), header = T, skip = 9, sep = "\t", stringsAsFactors = FALSE)
	cat('\tDone.\n')
	cnSet <- .curateAgilentCnSet(cnSet)
	return(cnSet)
}

.curateAgilentCnSet <- function(cnSet){
	keepCol <- which(as.character(colnames(cnSet)) %in%
							c(	"ProbeName", "SystematicName",
								"gMedianSignal", "rMedianSignal",
								"gIsSaturated", "rIsSaturated",
								"gIsFeatNonUnifOL", "rIsFeatNonUnifOL",
								"gIsWellAboveBG", "rIsWellAboveBG")
								)															
		cat('\tFiltering control probes...')
	isChr = grep('^chr[^Mrandom]*$', cnSet$SystematicName)
	cnSet <- cnSet[isChr, keepCol]
		cat('\tDone.\n')
		
		cat('\tAnnotate chrX, chrY...')
	systNames = sapply(cnSet$SystematicName, function(x){unlist(strsplit(x, ':'))})
	systNames = gsub('X', '23', systNames); systNames = gsub('Y', '24', systNames)
		cat('\tDone.\n')

		cat('\tGetting chr nums...')
	chrNum = systNames[seq(1, length(systNames), by = 2)]
	chrNum = as.numeric(gsub('chr', '', chrNum))
		cat('\tDone.\n')

		cat('\tGetting probe positions...')
	positions = systNames[seq(2, length(systNames), by = 2)]
	start = sapply(positions, function(x){unlist(strsplit(x, '-'))[1]})
	cnSet <- cbind.data.frame(ProbeName = cnSet[,1], ChrNum = chrNum, ChrStart = as.numeric(start), cnSet[,-c(1:2)])
	cnSet = cnSet[order(cnSet$ProbeName), ]
		cat('\tDone.\n')

	return(cnSet)
}


# Load workflow
scriptPath = "/Users/fredcommo/Documents/Projet Safir/Safir R Sources/CGHObjects/"
setwd(scriptPath)
source('SourceCodeObj.R')
