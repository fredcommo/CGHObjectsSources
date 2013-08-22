# Load entity function reads in the tab separated files
buildAgilent <- function(synId){
  # Load cghData from a synapse entity and build an AgilentObject.
  entity <- synGet(synId)
	object <- AgilentObj(info = c(fileName = propertyValue(entity, 'name'),
                                synapseId = propertyValue(entity, 'id'),
                                platform = 'Agilent'))
	object@info <- c(object@info, .readAgilentInfo(entity))
	object@cnSet <- .readAgilentMatrix(entity)
	object <- suppressFlags(object)
	object <- suppressDuplic(object)
	object <- preset(object)
	return (object)
}

.readAgilentInfo <- function(entity){
	cat('Reading information...')
	a.info <- read.csv(entity@filePath, header = F, fill = T, skip = 1, nrows = 8, sep = "\t")
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
	cnSet <- read.csv(entity@filePath, header = T, skip = 9, sep = "\t", stringsAsFactors = FALSE)
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

pushToSynapse <- function(cghObject, parentId = 'syn2116884', codeId = 'syn2117484'){
  Name <- paste0(gsub('.txt.bz2', '', getInfo(cghObject, 'fileName')))
  rawId <- getInfo(cghObject, 'synapseId')
  fileName <- sprintf('%s_%s_Profile.rds', Name, rawId)
  save(cghObject, file = file.path(tempdir(), fileName))
  file <- File(file.path(tempdir(), fileName), parentId = parentId, name = fileName)
  file <- synStore(file,
                      activityName = 'Genomic Profile',
                      used = list(list(entity = synGet(rawId, downloadFile = FALSE), wasExecuted = FALSE),
                                  list(entity = synGet(codeId, downloadFile = FALSE), wasExecuted = TRUE)))
  return(file)
}

buildShinyData <- function(cghObject, builtFrom = Profile, parentId = 'syn2116884'){
  # Save a smaller object to be called in shinyApps
  Name <- paste0(gsub('.txt.bz2', '', getInfo(cghObject, 'fileName')))
  rawId <- getInfo(cghObject, 'synapseId')
  fileName <- sprintf('%s_%s_shiny.rds', Name, rawId)
  cnSet = getCNset(cghObject)
  x <- cnSet$genomicPos[cnSet$ChrNum != 24]
  y <- cnSet$Log2Ratio[cnSet$ChrNum != 24]
  segTable <- getSegTable(cghObject)
  Samp <- seq(1, length(y), len = 20e3)
  shinyData <- list(gPos = x[Samp],
                    L2R = y[Samp],
                    segTable = segTable[segTable$chrom %in% 1:23,],
                    fileName = getInfo(cghObject, 'fileName'),
                    synId = getInfo(cghObject, 'synapseId'))
  save(shinyData, file = file.path(tempdir(), fileName))
  shinyData <- File(file.path(tempdir(), fileName), parentId = parentId)
  shinyData <- synStore(shinyData, used = list(list(entity = synGet(propertyValue(builtFrom, 'id'), downloadFile = FALSE),
                                                    wasExecuted = FALSE))
                        )
  .buildShinyApp(propertyValue(shinyData, 'id'))
  return(shinyData)
}

.buildShinyApp <- function(synId, mainPath = '/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/shinyAppsCGH'){
  newDir <- file.path(mainPath, paste0('shiny_', synId))
  dir.create(newDir)
  # Copy files from synapse
  server <- synGet('syn2128542'); file.copy(server@filePath, newDir)
  ui <- synGet('syn2128545'); file.copy(ui@filePath, newDir)
  # Add the target on the server.R
  fConn <- file(file.path(newDir, 'plotCGHtemplate_server.R'), 'r+')
  Lines <- readLines(fConn)
  writeLines(c(sprintf("require(synapseClient)\ncgh  <- synGet('%s')\n", synId), Lines), con = fConn)
  close(fConn)
  # Copy the folder on the Shiny server & add the link in the wiki.
  system(sprintf("scp -r %s toShiny:~/ShinyApps", newDir))
  .buildWiki(synId)
}

.buildWiki <- function(synId){
  # Create the folder wiki uri.
  iframe <- .createIFrame(synId)
  folderWikiUri <- sprintf("/entity/%s/wiki", 'syn2129126')#synId)
  folderWiki <- try(synRestGET(folderWikiUri), silent = TRUE)
  # If already exists, update.
  if(class(folderWiki) != 'try-error')
    .updateIframe(file.path(folderWikiUri, folderWiki$id), folderWiki, iframe)
  # Else: Start a wiki, and initialize the markdown as empty.
  else{
    folderWiki <- list()
    folderWiki$attachmentFileHandleIds <- list()
    folderWiki$markdown <- iframe
    folderWiki <- synRestPOST(folderWikiUri, folderWiki)
  }
}
.createIFrame <- function(synId){
  startFrame <- "${iframe?site=http%3A%2F%2Fshiny%2Esynapse%2Eorg%3A3838%2Fusers%2Ffcommo%2F"
  endFrame <- "%2F&height=1100}"
  return(paste0(startFrame, paste0('shiny_', synId), endFrame))
}
.updateIframe <- function(folderWikiUri, folderWiki, newFrame){
  folderWiki$markdown <- newFrame
  folderWiki <- synRestPUT(folderWikiUri, folderWiki)
}

## Not used!
# updateShiny <- function(synId){
#   # Load the ui.R from synapse
#   # Add the synId in the idsList if not on it yet.
#   # Push to synapse, and push to the synapse.server
#   cat('Reading the current ui.R\n')
#   ui <- synGet('syn2145898'); file.copy(ui@filePath, tempdir())
#   fConn <- file(ui@filePath, 'r+')
#   Lines <- readLines(fConn)
#   if(!grepl(synId, Lines[1])){
#     cat('Adding new id in the list\n')
#     newTag <- sprintf(', \"%s\")\n', synId)
#     L1 <- gsub(')', newTag, Lines[1])
#     writeLines(c(L1, Lines[-1]), con = fConn)
#     }
#   close(fConn)
#   
#   # Update the synapse version: shiny folder is 'syn2145860'
#   cat('Push to synapse\n')
#   file <- File(ui@filePath, parentId = 'syn2145860')
#   file <- synStore(file)
#   
#   # Copy the folder on the Shiny server & add the link in the wiki.
#   cat('Push to the shiny server\n')
#   system(sprintf("scp %s toShiny:~/ShinyApps/shinyCGHselectSample/", ui@filePath))
# }
