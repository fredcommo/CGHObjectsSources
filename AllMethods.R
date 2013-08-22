
# Method do read and store array information from Agilent files.
setMethod('readInfo', 'AgilentObj',
	function(object){
		'
		platform = Agilent
		Open the selected file and read the 8 first lines
		Assign the items and values as a list to object@info
		Return an object of the same class
		'	
		fileName = getInfo(object, 'fileName')
		a.info <- read.csv(fileName, header = F, fill = T, skip = 1, nrows = 8, sep = "\t")
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

		object@info = c(	object@info,
									labId = labId,
									barCode = barCode,
									gridName = gridName,
									scanDate = scanDate, 
									arrayProtocol = arrayProtocol,
									gridGenomicBuild = gridGenomicBuild,
									ref = ref,
									analyseDate = format(Sys.Date(), "%m-%d-%Y")
									)
		return (object)
		}
)


# Method do read and store array information from Affymetrix files.
setMethod('readInfo', 'AffyObj',
	function(object){
		'
		platform = Affymetrix
		Open the selected file and read the 400 first lines
		Assign the items and values as a list to object@info
		Return an object of the same class
		'
		getTagValue <- function(x){
			return(unlist(strsplit(x, '='))[2])
			}

		fileName = getInfo(object, 'fileName')
		# Replace by readLine, split by '=', and find 'ProbeSet'
		arrayInfos <- readLines(fileName, n = 400)
		path = unlist(strsplit(getwd(), '/'))
		labId = path[length(path)]
		fileInfo = unlist(strsplit(fileName, '_'))
		barCode = fileInfo[3]
		arraySet = getTagValue(arrayInfos[grep("#state-array-name", arrayInfos)])
		Date = getTagValue(arrayInfos[grep("#state-time-start", arrayInfos)])
		Date <- unlist(strsplit(Date, ' '))
		scanDate = paste(Date[1], Date[2], Date[3], sep = '-')
		arrayProtocol = getTagValue(arrayInfos[grep("#state-program-name", arrayInfos)])
		grid1 = getTagValue(arrayInfos[grep("#%genome-version-ucsc", arrayInfos)])
		grid2 = getTagValue(arrayInfos[grep("#%genome-version-ncbi", arrayInfos)])
		gridGenomicBuild = paste(grid1, grid2, sep = '/')
		ref = getTagValue(arrayInfos[grep("#state-reference-input", arrayInfos)])
		object@info = c(	object@info,
									labId = labId,
									barCode = barCode,
									arraySet = arraySet,
									scanDate = scanDate, 
									arrayProtocol = arrayProtocol,
									gridGenomicBuild = gridGenomicBuild,
									ref = ref,
									analyseDate = format(Sys.Date(), "%m-%d-%Y")
									)
		return (object)
		}
)

# Method do read and store the Cy3 & Cy5 values from Agilent files.
setMethod('readCN', 'AgilentObj',
	function(object){
		'
		platform = Agilent
		Open the selected file and read after the 8 first line to to get the values
		Select only the columns of interest: 
			- Probes Ids: ProbeName, SystematicName,
			- Probe Signal: gMedianSignal, rMedianSignal,
			- Probes QC: gIsSaturated, rIsSaturated, gIsFeatNonUnifOL, rIsFeatNonUnifOL, gIsWellAboveBG, rIsWellAboveBG
		Assign the matrix to object@cnSet
		Return an object of the same class
		'
		# Read and get some suppl. info
		fileName = getInfo(object, 'fileName')
		fileSize = file.info(fileName)$size
		cat('\n\tsize: ', round(fileSize/1e3), 'Kb\n')
		s = system.time(cnSet <- read.csv(fileName, header = T, skip = 9, sep = "\t", stringsAsFactors = FALSE))
		cat('\ttime:', s[1], 'sec\n')

		# Select columns of interest : probeIds, Cy3 value, Cy5 value, QC columns
		keepCol <- which(as.character(colnames(cnSet)) %in% c(	"ProbeName", "SystematicName",
																								"gMedianSignal", "rMedianSignal",
																								"gIsSaturated", "rIsSaturated",
																								"gIsFeatNonUnifOL", "rIsFeatNonUnifOL",
																								"gIsWellAboveBG", "rIsWellAboveBG")
																								)
															
		# Select QC columns and rows containing probes with Ids of type 'A_xxx'
		# cnSet <- cnSet[grep('^A', cnSet$ProbeName), keepCol]
		# cnSet <- cnSet[order(cnSet$ProbeName), ]
	
		# use grep methods
		cat('\tFiltering control probes...')
		isChr = grep('^chr[^Mrandom]*$', cnSet$SystematicName)
		cnSet <- cnSet[isChr, keepCol]
		# cnSet <- cnSet[order(cnSet$ProbeName), ]
		cat('\tDone.\n')

		cat('\tAnnotate chrX, chrY...')
		systNames = sapply(cnSet$SystematicName, function(x){unlist(strsplit(x, ':'))})
		systNames = gsub('X', '23', systNames); systNames = gsub('Y', '24', systNames)
		cat('\tDone.\n')

		cat('\tGet chr nums...')
		chrNum = systNames[seq(1, length(systNames), by = 2)]
		chrNum = as.numeric(gsub('chr', '', chrNum))
		cat('\tDone.\n')

		cat('\tGet probe positions...')
		positions = systNames[seq(2, length(systNames), by = 2)]
		start = sapply(positions, function(x){unlist(strsplit(x, '-'))[1]})
		cat('\tDone.\n')

		cnSet <- cbind.data.frame(ProbeName = cnSet[,1], ChrNum = chrNum, ChrStart = as.numeric(start), cnSet[,-c(1:2)])
		object@cnSet = cnSet = cnSet[order(cnSet$ChrNum, cnSet$ChrStart, cnSet$ProbeName), ]
		return(object)
		}
)

# Method to get Log values from Affymetrix files.
# To do: add an attribute snpProbes.
setMethod('readCN', 'AffyObj',
	function(object){
		'
		platform = Affymetrix
		Open the selected file and read 1000 rows
		Define where ProbeSet is : starting point to get values.
		Compute the genomic positions according to hg19.
		Filter and keep the CN probes and SNP probes separately.
		Assign the matrix to object@cnSet
		Return an object of the same class
		'
		fileName = getInfo(object, 'fileName')
		a.info <- read.csv(fileName, header = F, fill = T, skip = 0, nrows = 1000, sep = "\t")
		datastart <- which(a.info == "ProbeSet")

		# Read from 'ProbeSet' position.
 		s = system.time(cnSet <- read.csv(fileName, header = T, skip = datastart - 1, sep = "\t"))
		fileSize = file.info(fileName)$size
		cat('\tsize: ', round(fileSize/1e3), 'Kb\n')
		cat('\ttime:', s[1], 'sec\n')
		
		colnames(cnSet)[1:3] <- c("ProbeName", "ChrNum", "ChrStart")
		levels(cnSet$ChrNum)[23] <- 23
		levels(cnSet$ChrNum)[24] <- 24
		cnSet$ChrNum <- as.numeric(as.character(cnSet$ChrNum))
		cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart), 1:6]
		
		# # Genomic postions
		# arrayInfoPath = '/Users/fredcommo/Documents/Projet_Safir/Arrays Infos/'
		# hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
		# cumLen = cumsum(as.numeric(hg19$length))
		# cumLen = c(0, cumLen[-length(cumLen)])
		# chrTable = table(cnSet$ChrNum, useNA = 'ifany')
		# genomic = c()
		# for(i in 1:length(chrTable)){
			# n = chrTable[i]
			# genomic = c(genomic, rep(cumLen[i], n))
			# }
		# genomic = ifelse(!is.na(cnSet$ChrStart), genomic + cnSet$ChrStart, NA)
		# cnSet <- cbind.data.frame(cnSet[,c("ProbeName", "ChrNum", "ChrStart")], genomicPos = genomic, Log2Ratio = cnSet[,'Log2Ratio'])

		# Select the CN probes (probeId starting with CN)
		cnProbes = which(substr(cnSet$ProbeName, 1, 2) == "CN")
		object@cnSet = cnSet[cnProbes,]
		object@snpSet = cnSet[-cnProbes,] 
		return (object)
		}
)

# This method define what spots have to be flagged, according to Agilent criteria. In case of flag, Cy3 and Cy5 values are replaced by NA.
setMethod('suppressFlags', 'cghObj',
	function(object){
		'
		If Agilent:
			Filter the probes according to Agilent QC criteria.
			In case of flags, Cy3 & Cy5 values are replaced by NA
			Assign the filtered CNset to object@cnSet
		Return an object of the same class
		'
		if(getInfo(object, 'platform') == 'Agilent'){
			cat('\tSuppressing flagged probes...')
			cnSet <- getCNset(object)
			flags <- which(	cnSet$gIsSaturated == 1 | 			# 1 = non valid gIsSaturated
									cnSet$rIsSaturated == 1 | 					# 1 = non valid rIsSaturated
									cnSet$gIsFeatNonUnifOL == 1 | 			# 1 = non valid gIsFeatureNonUnifOL
									cnSet$rIsFeatNonUnifOL == 1 | 			# 1 = non valid rIsFeatureNonUnifOL
									cnSet$gIsWellAboveBG == 0 | 				# 0 = non valid gIsWellAboveBG
									cnSet$rIsWellAboveBG == 0)					# 0 = non valid rIsWellAboveBG
			cnSet$gMedianSignal[flags] = cnSet$rMedianSignal[flags] = NA
			flagNames <- c('gIsSaturated', 'rIsSaturated', 'gIsFeatNonUnifOL', 'rIsFeatNonUnifOL', 'gIsWellAboveBG', 'rIsWellAboveBG')
			object@cnSet = cnSet[,-which(colnames(cnSet) %in% flagNames)]
			}
		cat('\tDone.\n')
		return(object)
		}
)


setMethod('suppressDuplic', 'cghObj',
	function(object){
	'
	The column containing the probeIds have to be named $ProbeName
	In case of duplicated probes, compute the median value.
	Suppress the duplicated probes.
	Return an object of the same class.
	'
	cat('\tSuppressing duplicated probes...')
	cnSet <- getCNset(object)[,1:5]
	if (!any(colnames(cnSet) == 'ProbeName')) stop('None of the columns can be identifed as ProbeNames')
	duplicIndex <- which(duplicated(cnSet$ProbeName))
	if(length(duplicIndex>0)) {
		duplicProbes <- as.character(unique(cnSet$ProbeName[duplicIndex]))
		duplicSet = subset(cnSet, cnSet$ProbeName %in% duplicProbes)
		medianSet <- ddply(.data = duplicSet, .variables=.(ProbeName, ChrNum, ChrStart), summarize,
										# ChrNum = ChrNum, ChrStart = ChrStart,
										rMedianSignal = median(rMedianSignal), gMedianSignal = median(gMedianSignal))
		cnSet <- rbind.data.frame(cnSet[-which(cnSet$ProbeName %in% duplicProbes),], medianSet)
		}
	cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart),]
	rownames(cnSet) <- seq(1, nrow(cnSet))
	object@cnSet <- cnSet
	cat('\tDone.\n')
	return(object)
	}
)

# This method defines the default values (paramaters) used in the segmentation step.
setMethod('preset', 'cghObj',
	function(object, Ksmooth = 1000, Kmax = 20, Nmin = Kmax*8, Mwidth = 2, UndoSD = 0.5, Alpha = 1e-10){
		'
		Define the default settings called by the segmentation step.
		'
		cat('Adding CBS presettings...')
		if(getInfo(object, 'platform') == 'Affymetrix') Ksmooth = 5000
		param = list(Ksmooth = Ksmooth, Kmax = Kmax, Nmin = Nmin, Mwidth = Mwidth, UndoSD = UndoSD, Alpha = Alpha)
		object@param = param
		cat('\tDone.\n')
		return (object)
		}
)


setMethod('adjustSignal', 'cghObj', function(object, Cy = TRUE, Fract = 1, GC = TRUE, Centr = TRUE){

# object:				an object of class(cghObj)
# Cy, CG, Centr:	Indicates if theses corrections have to be applied. If is.Agilent = FALSE, Cy = CG = FALSE
# Fract:				Indicates the approximative proportion of tumor cells (do not use before the method has been validated) 

	cnSet <- getCNset(object)
	if (getInfo(object, 'platform') ==  'Agilent'){
		cnSet <- CyAdjust(cnSet, Fract, Centr)			# helper function
		object@param$CyAdjusted = TRUE
		if(GC){
			cnSet <- GCadjust(cnSet, gridName = getInfo(object, 'gridName'))						# helper function
			object@param$GCAdjusted = TRUE
			}
		}
	else{
		object@param$CyAdjusted = FALSE
		object@param$GCAdjusted = FALSE
		}
	# Need to add genomic positions, then order by genomicPos
	cnSet <- .addGenomicPos(cnSet)
	cnSet$Log2Ratio <- .supprOutliers(cnSet$Log2Ratio)
	object@param$dLRs <- .dlrs(cnSet$Log2Ratio)
	object@param$MAD <- .MAD(cnSet$Log2Ratio)
	object@cnSet <- cnSet
	cat('dLRs:', object@param$dLRs, '\tMAD:', object@param$MAD, '\n')
	return(object)
	}
)

##########################
setMethod('EMnormalize', 'cghObj',
	function(object, cut = c(-0.5, 0.5), G = 3:7, peakThresh = 0.9, MergePeaks = FALSE, MergeVal = 0.025, method = "Left",
										Plot = TRUE, Expand = 1.75, useX11 = FALSE, Save = FALSE){

# object: 				an object of class cghObj
# cut:					Quantiles thresholds in EM centralization
# G:						Number of groups to consider in EM centralization
# peakThresh: 		Proportion of the maximum density value to consider a peak as the central reference
# MergePeaks:		Allow to merge two peaks if there distance is lower than Mergeval
# MergeVal:			Minimum distance to consider two peaks as different.
# method:			Define the method to choose the reference peak. Left = left major peak, Zero = major peak closed to Zero, Right = right major peak.
# Plot:	 				Edit a density plot with EM peaks
# Expand:				A graphic parameter to expand the x axis
# Save:				To save automatically the density plot. The current default folder is Root/~/Safir CentralizeProfiles
# Root: 				The root to save the density plot.

	'
	Use a EM algorithm to model the LogRatio density as a mixture of gaussian distributions
	'
	cnSet <- getCNset(object)
	testLR <- cnSet$Log2Ratio
	runLR <- smoothLR(testLR, getInfo(object, 'platform'), cut, K = 301)		# helper function

	cat('Building model...')
	EMmodel <- buildEMmodel(runLR, G, by = 20)										# helper function
	nG <- EMmodel$nG; m <- EMmodel$m; s <- EMmodel$s ; p <- EMmodel$p
	cat('\tDone.\n')
	
	# merge Classes: depending on MergePeaks and MergeVal
	if(MergePeaks){
		cat('Merging the closer peaks...')
		mergedValues <- mergePeaks(nG, m, s, p, MergeVal)								# helper function
		nG <- mergedValues$nG; m <- mergedValues$m; s <- mergedValues$s ;p <- mergedValues$p
		cat('\tDone.\n')
		}

	# compute densities
	cat('Computing densities...')
	denLR <- density(runLR)
	computeD <- computeDensities(length(runLR), m, p, s)								# helper function
	dList <- computeD$dList
	peaks <- computeD$peaks
	kurt <- computeD$kurt
	cat('\tDone.\n')
	
	cat('Gaussian mixture:')
	cat("\nn.peaks = ", nG) 
	bestPeak <- chooseBestPeak(peaks, m, peakThresh)
	correct = m[bestPeak]
	# cat('\nmeans:', m, '\nProp:', p, '\nSigmaSq:', s)
	cat ('\nCVs per group:\n')
	for (grp in 1:length(m)){
		cat('\n')
		cat('\tGrp', grp, '\tprop:', p[grp], ':\tmean:', m[grp],
			'\tSd:', sqrt(s[grp]), '\tCV:', signif(sqrt(s[grp])/abs(m[grp]), 4),
			'\tKurtosis:', kurt[grp])
			}
	cat("\n\nEM correction factor = ", correct, "\n\n")

  sampleId <- gsub('\\.(.*)', '', getInfo(object, 'fileName'))
  synapseId <- getInfo(object, 'synapseId')
	Title = paste(paste(sampleId, synapseId, sep = ' - '), "\nEM centralization: correction.factor =", round(correct, 5))
	EMplot = plotEMmodel(runLR, dList, m, bestPeak, cut, Title)

	# return the full cgh objet containing the Log2Ratios adjusted for Cy and CG, and centralized by EM
	cnSet$Log2Ratio = cnSet$Log2Ratio - correct
	object@cnSet = cnSet
	object@param$EMcentralized = TRUE
	object@param$nPeak = nG
	object@param$PeakValues = as.numeric(m)
	object@param$SigmaSq = s
	object@param$correctionValue = as.numeric(correct)
	object@probesDensity = EMplot
	return(object)
	}
)


setMethod('SegmentCGH', 'cghObj',
	function(object, Smooth = TRUE, UndoSD = NULL){
		cnSet = getCNset(object)
		params = getParam(object)

		# Segmentation parmeters
		Ksmooth <- params$Ksmooth; Kmax <- params$Kmax; Nmin <- params$Nmin; Mwidth <- params$Mwidth; Alpha <- params$Alpha
		if(is.null(UndoSD)) UndoSD <- 0.5 + max(sqrt(params$SigmaSq))*3
		#

		LR <- cnSet$Log2Ratio
		Chr <- cnSet$ChrNum
		genomicPos <- cnSet$genomicPos

		cat('\nBuilding CNA object...')
		cna.obj <- CNA(LR, Chr, genomicPos, data.type = "logratio", sampleid = paste(getInfo(object, 'synapseId'), getInfo(object, 'barCode'), sep = "_"))
		cat('\tDone.')
		
		cat('\nSegmenting...')
		if(Smooth){
			smooth.cna.obj <- smooth.CNA(cna.obj, smooth.region = Ksmooth)
			seg.cna.obj <- segment(smooth.cna.obj, undo.splits = "sdundo", undo.SD = UndoSD, alpha = Alpha, kmax = Kmax, nmin = Nmin, min.width = Mwidth)
			}
		else seg.cna.obj <- segment(cna.obj, undo.splits = "sdundo", undo.SD = UndoSD, alpha = Alpha, kmax = Kmax, nmin = Nmin, min.width = Mwidth)
		cat('\n\tUndoSD:', UndoSD)
		cat('\n\tDone.')
		
		# To do!
		# Ncgh <- QCSegm(seg.cna.obj, Ncgh)		# Segmentation Quality Control

		# Implement these methods in HelperFunctions 
		cat('\nAdding segment values in the main table...')
		seg.cna.obj <- MedSegm(seg.cna.obj, cnSet)								# helper function
		cnSet <- AddSegments(seg.cna.obj, cnSet)								# helper function
		cat('\tDone.\n')
		
		# Update
		# object@info$dLRsd <- as.character(dLRsd(cnSet$Log2Ratio))
		# object@info$ScriptVersion <- SVersion 	! How to add the script-version number ?
    object@param$UndoSD <- UndoSD
		object@segTable = seg.cna.obj$output
		object@cnSet = cnSet
		return(object)
	}
)

setMethod('createProfile', 'cghObj',
function(object, gain = log2(2.25/2), loss = log2(1.80/2), Ylim = NULL){
	require(ggplot2, quietly = TRUE)
	myBlue <- rgb(0, 0.45, 1, 1)
# 	arrayInfoPath = '/Users/fredcommo/Documents/Projet_Safir/Arrays_Infos/'
# 	hg19 <- read.csv(paste0(arrayInfoPath, 'human.chrom.info.hg19.FC.txt'), header = TRUE, sep = '\t')

  # Load hg19 chromosome length table
  if(!exists('hg19')){
    ent <- synGet('syn2141399')
    hg19 <- read.csv(ent@filePath, header = TRUE, sep = '\t')
  }
	cumLen = cumsum(as.numeric(hg19$length))
	cumCentr <- 1/2*cumLen[1]
	for(chr in 2:length(cumLen)) cumCentr = c(cumCentr, cumLen[chr-1] + 1/2*hg19$length[chr])

	cnSet = getCNset(object)
	cnSet <- cnSet[cnSet$ChrNum != 24,]
	segTable = getSegTable(object)
	Title = paste(gsub('.txt.bz2', '', getInfo(object, 'fileName')), '-',
                getInfo(object, 'synapseId'), '-',
                getInfo(object, 'analyseDate'),
						'\nGain threshold: ', round(gain, 3), ' Loss threshold:', round(loss, 3))

	if(any(is.na(cnSet))){
		NAs <- c()
		for(i in 1:ncol(cnSet))
			NAs <- c(NAs, which(is.na(cnSet[,i])))
		NAs <- unique(NAs)
		cnSet <- cnSet[-NAs,]
	}

	Samp = seq(1, nrow(cnSet), len = 20e3)
	cnSet = cnSet[Samp,]
	cnSet <- cbind.data.frame(cnSet, rMed1 = runmed(cnSet$Log2Ratio, k=9), rMed2 = runmed(cnSet$Log2Ratio, k=81))
	minY <- min(min(segTable$seg.med, na.rm = TRUE), -1)
	gPlot <- ggplot(data = cnSet, aes(x = genomicPos, y = rMed1)) +
					geom_point(pch = 19, cex = 0.1, col = 'grey30') +
					geom_hline(yintercept = 0) +
					geom_vline(xintercept = cumLen[1:23], color = 'grey30', linetype = 2, size = 0.2) +
					annotate('text', x = c(0, cumCentr[1:23]), y = rep(1.4, 24), label = c("Chr", seq(1, 23)), size = 2.5, colour = 'grey30') +
					geom_line(aes(x = genomicPos, y = rMed2), linetype = 1, size = 0.15) +
					geom_point(aes(x = genomicPos, y = Segm), cex = 1,
									col = ifelse(cnSet$Segm<= loss, 'red3', ifelse(cnSet$Segm>= gain, myBlue, 'grey40'))) +
					coord_cartesian(ylim = range(minY, 1.5))+
					# coord_cartesian(ifelse(is.null(Ylim), ylim = range(minY, 1.5), ylim = Ylim))+
					ggtitle(Title) +
					xlab('Genomic position (bp)')+
					ylab('Log2(Ratio)') +
					theme_bw() +
					theme(	panel.grid.major = element_blank(),
								panel.grid.minor = element_blank(),
								plot.title = element_text(lineheight=.8, face="bold"))

	object@gProfile <- gPlot
	return(object)
	}
)

# arrayInfoPath = '/Users/fredcommo/Documents/Projet_Safir/Arrays_Infos/'		
# geneDB <- readRDS(paste0(arrayInfoPath, 'myGeneDB_2013_Mar_26.rds'))
ent <- synGet('syn2141368')
geneDB <- readRDS(ent@filePath)

setMethod('geneOfInt', 'cghObj',
		function(object, geneList, DB = geneDB){
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
)

