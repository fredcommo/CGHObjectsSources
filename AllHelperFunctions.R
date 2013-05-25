
# The function activates an interactive selection of file to load
# Collect microarray informations and Cy5/Cy3 intensity values,
# and suppress flags and duplicated probes.
# Root : A letter indicating the hard drive partition to use. May run with other indications (~/home)

##########################################################################################
###########################

# Class are defined in AllClasses.R
# Accessors 	are defined in AllAccessors.R
# Generic methids are defined in AllGenerics.R
# Methods are in AllMethods .R

# validity method : To do !

# Common functions : To do !
	# setGeneric('Norm', function(object) standardGeneric('Norm'))
	# setGeneric('Segment', function(object) standardGeneric('Segment'))
	# setGeneric('Plot', function(object) standardGeneric('Plot'))
	# setGeneric('Save', function(object) standardGeneric('Save'))
	# setGeneric('EditSupTab', function(object) standardGeneric('EditSubTab'))



###########################
###########################

buildCGHObj.v01 <- function(){
	f = try(getFile(), silent = TRUE)								# helper function
	if(class(f) == 'try-error') stop ('No selected file!\n')
	
	cat("\nGetting Array Information...")
	object = getAnnot(f)												# helper function
	object <- readInfo(object)										# Class function
	
	cat('\n\nReading', f, '...')
	object <- readCN(object)											# Class function
	
	cat('\nCleaning data...\n')	
	object <- suppressFlags(object)									# Class function
	object <- suppressDuplic(object)								# Class function
	object <- preset(object)											# Class function
	
	cat("\n\n\t*** Everything is fine. Lucky you :-) ***\n\n")
	return(object)
}


###########################
getFile <- function(){
	
	'
	Called by buildCGHObj()
	Open a dialogue box and return the name of the file selected by the user.
	'

	#fileName <- tclvalue(tkgetOpenFile()) 									# Open current the directory to select the file to load.
	fileName = file.choose()
	if (!nchar(fileName)) 
		tkmessageBox(message = "No file selected!")
	#else 
    #	tkmessageBox(message = paste("Selected File:", fileName))

	fileName = unlist(strsplit(fileName, '/'))
	path = fileName[1]
	for (i in 2:(length(fileName)-1))
		path = paste(path, fileName[i], sep = '/')
	setwd(path)
	return (fileName[length(fileName)])															# add Folder
}

###########################
getAnnot <- function(fileName){
	'
	Called by buildCGHObj()
	Open the selected file and read the first line to identify the type of array (platform = Agilent or Affy)
	Return an object of class Agilent or Affymetrix.
	Assign the fileName (without its path), the sampleId, and the platform to object@info
	'	
	obj = list()
	fileInfo = unlist(strsplit(fileName, '_'))
	sampleId = 	fileInfo[1]
	a <- try(read.csv(fileName, header = F, fill = T, skip = 0, nrows = 1, sep = "\t"), silent = T)
	{
		if(class(a)!="try-error") cat("\nFile summary:", "\n\tpath: ", getwd(), "\n\tfileName:", fileName)
		else stop("\t", fileName, ": No such file on directory (glup!). Please fix it!\n\n")
	}

	# According to 'platform', create a corresponding cghObject
	if (a[1,1] == '#GenomeWideSNP_6.na32.annot.db'){
		platform = 'Affymetrix'
		Obj <- AffyObj(info = c(fileName = fileName, sampleId = sampleId, platform = platform))
		}
	else{
		platform = 'Agilent'
		Obj <- AgilentObj(info = c(fileName = fileName, sampleId = sampleId, platform = platform))
	}
	return (Obj)
}

###########################
# In construction
.computeFract <- function(object){
	cnSet <- getCNset(object)
}

###########################

CyAdjust <- function(cnSet, Fract, Centr){
	'
	Called by adjustSignal(object)
	If Fract (prop of tumor in the sample), adjust the tumor signal (Cy5)
	Compute and adjust the Log2(Cy3/Cy5)
	Return the data.frame with a supplementary column: Log2Ratio
	'
	cat('\nCy effect adjustment...')
	g <- log2(cnSet$gMedianSignal)					# Ref in Cy3
	r <- log2(cnSet$rMedianSignal)					# Test in Cy5

	# Calculates weights to correct the dilution effect (due to tumor cell rate). No effect if Fract = 1. DO NOT USE until validation !
	if(!is.null(Fract)){
		cat('\nCompute Fract adjustment...')
		Q <- quantile(r, Fract, na.rm = TRUE)
		w <- 1/(1+exp(1/sqrt(Fract)*(Q-r)))
		r <- r*(1 + (1-Fract)*w)
		}
	M <- r - g
	A <- (r + g)/2
	Loess <- loessFit(M, A)$fitted
	LR <- M - Loess
	cnSet$Log2Ratio <- LR 	
	if (Centr)
		cnSet$Log2Ratio - median(cnSet$Log2Ratio, na.rm = T)
	cat('\tDone.\n')
	return (cnSet)
	}

#################

GCadjust <- function(cnSet, gridName){
	'
	Called by adjustSignal(object)
	Adjust the GC%
	'
	arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
	cat('GC% adjustment...')
	cnSet <- cnSet[order(cnSet$ProbeName),]
	AgilentDB = readRDS(paste0(arrayInfoPath, gridName, '.rds'))
	AgilentDB <- AgilentDB[order(AgilentDB$ProbeID),]
	AgilentDB <- AgilentDB[which(AgilentDB$ProbeID %in% cnSet$ProbeName),]

	# Check the probeNames
	if(!all(as.character(AgilentDB$ProbeID) == as.character(cnSet$ProbeName)))
		stop('Agilent DB: Probe names do not match.')
		
	lr = cnSet$Log2Ratio
	GC <- AgilentDB$GCpercent
	adjLr <- lr - loessFit(lr, GC)$fitted
	cnSet$Log2Ratio = adjLr
	cat('\tDone.\n')
	cnSet <- cbind.data.frame(cnSet[,c('ProbeName', 'ChrNum', 'ChrStart')], Log2Ratio = cnSet[,'Log2Ratio'])
	return(cnSet)
	}

#################

.addGenomicPos <- function(cnSet){
		cat('\n\nAdding hg19 genomic positions...')
		arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
		hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
		cumLen = cumsum(as.numeric(hg19$length))
		cumLen = c(0, cumLen[-length(cumLen)])
		cnSet = cnSet[order(cnSet$ChrNum, cnSet$ChrStart, cnSet$ProbeName), ]
		chrTable = table(cnSet$ChrNum, useNA = 'ifany')
		genomic = c()
		for(i in 1:length(chrTable)){
			n = chrTable[i]
			genomic = c(genomic, rep(cumLen[i], n))
			}
		genomic = ifelse(!is.na(cnSet$ChrStart), genomic + cnSet$ChrStart, NA)
		cnSet <- cbind.data.frame(cnSet[,c("ProbeName", "ChrNum", "ChrStart")], genomicPos = genomic, Log2Ratio = cnSet[,'Log2Ratio'])
		.checkOverlaps(cnSet)
		return(cnSet)
}

#################

.checkOverlaps <- function(cnSet){
	overlaps <- sapply(seq(2, 24), function(chr){
		last <- max(cnSet$genomicPos[cnSet$ChrNum == chr-1], na.rm = TRUE)
		first <- min(cnSet$genomicPos[cnSet$ChrNum == chr], na.rm = TRUE)
		last > first
	})
	if(any(overlaps)) cat('overlap between', which(overlaps)[-1], "and", which(overlaps), '\n')
	else cat('No overlap\n')
}

#################

.Iqr <- function(x, n){
	'
	Called by .supprOutliers
	If xn is an outlier, then replace xn by the median of the +/-n neighbours (computed without xn)
	'
	xprim = x[-(n+1)]
	q = quantile(xprim, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
	iqr = q[3] - q[1]
	if(is.na(x[n+1]) | (x[n+1] < q[2]-1.5*iqr | x[n+1] > q[2]+1.5*iqr))
		x[n+1] <- q[2]
	return(x)
}

.supprOutliers <- function(x, n = 5){
	'
	Called by adjustSignal(object)
	Replace the outliers by the medians of the +/- n neighbours
	'
	cat('\nSuppressing outliers...')
	p = length(x)
	xnew <- xprim <- c(x[1:n], x, x[(p-n+1):p])
	for(i in (n+1):p){
		tmp <- xprim[(i-n):(i+n)]
		if(sum(is.na(tmp)) <= n) xnew[i] <- .Iqr(tmp, n)[n+1]
		else xnew[i] <- NA
		}
	xnew <- xnew[(n+1):(p+n)]
	cat('\tDone.')
	return(xnew)
}

#################

smoothLR <- function(LR, Platform, cut, K){
	'
	Called by EMnormalize(object)
	Smoothing the LR vector to improve the EM classification.
	'
	if(Platform == 'Affymetrix')
		LR <- LR[seq(1, length(LR), by = 6)]
	if(any(is.na(LR))) LR <- LR[!is.na(LR)]
	runLR <- runmed(LR, k = K)	
	q1 <- cut[1]; q2 <- cut[2]
	runLR = runLR[which(runLR>=q1 & runLR<=q2)]
	return (runLR)
}

#################

buildEMmodel <- function(LR, G, by){
	'
	Called by EMnormalize(object)
	Model the distribution as a gaussian mixture.
	'
	model <- Mclust(LR[seq(1, length(LR), by = by)], G = G)
	nG <- model$G
	p <- model$parameters$pro
	m <- model$parameters$mean
	s <- model$parameters$variance$sigmasq
	if(length(s)<length(m)) s <- rep(s, length(m))
	p <- p[order(m)]
	s <- s[order(m)]
	m <- m[order(m)]
	return(list(nG = nG, m = m, p = p, s =s))
}

#################

mergePeaks <- function(nG, m, s, p, MergeVal){
	'
	Called by EMcentr(object)
	Given the initial EM parameters, recompute the gaussian mixture by merging the groups for which the between-group distance < MergeVal
	'
	Raw <- c(1, 1)
	while(length(Raw)!=0){
		Mdist <- matrix(0, nG, nG)
		for(i in 1:nG)
			for(j in 1:nG){
				Mdist[i, j] <- abs(m[i] - m[j])
				}
			diag(Mdist) <- NA
			Raw <- ceiling(which(Mdist<MergeVal)/nG)
			cat('Merging', Raw, "\n")
			if(length(Raw)!=0){
				C1 <- Raw[1]
				C2 <- Raw[2]
				m[C1] <- ((p[C1])*m[C1] + (p[C2])*m[C2])/(p[C1] + p[C2])
				s[C1] <- ((p[C1])*s[C1] + (p[C1])*s[C2])/(p[C1] + p[C2])
				p[C1] <- p[C1] + p[C2]
				m <- m[-C2]; s <- s[-C2]; p <- p[-C2]
				nG <- length(m)
				cat("means:", m, "\nVar:", s, "\nprops:", p, "\n\n")
				}
			}
		return(list(nG = nG, m = m, s = s, p = p))
		}


#################

computeDensities <- function(n, m, p, s){
	'
	Called in EMnormalize(object)
	Simulates the mixture model according to the returned EM paramaters.
	'
	dList = list()
	peaks <- kurt <- c()
	for(i in 1:length(m)){
		tmp <- rnorm(n*p[i], m[i], sqrt(s[i]))
		tmpD <- density(tmp, na.rm = T)
		tmpD$y = tmpD$y *p[i]
		dList[[i]] <- tmpD
		kurt <- c(kurt, kurtosis(tmp, type = 2))
		peaks <- c(peaks, max(tmpD$y))
		}
	return(list(dList = dList, peaks = peaks, kurt = kurt))
}

#################

chooseBestPeak <- function(peaks, m, peakThresh){
	'
	Called in EMnormalize(object)
	Estimates what peak as to be used as the centralization value.
	'
	best <- which(peaks>=max(peaks)*peakThresh & m<=0)
	if (length(best) != 0)
		cat('\nLeft peak at', m[best], 'has been chosen.')
	else{
		# Centered peak next
		best <- which(peaks>=max(peaks)*peakThresh)
		best <- best[which.min(abs(m)[best])]
		if (length(best) != 0)
			cat('\nCentral peak at', m[best], 'has been chosen.')
		else{
			# Right peak last
			best <- which(peaks>=max(peaks)*peakThresh & m>=m[which.max(peaks)])
			if (length(best) != 0)
					cat('\nRight peak at', m[best], 'has been chosen.\n')
				}
			}
	bestPeak = best[which.max(m[best])]
	return(bestPeak)
}

#################

plotEMmodel <- function(LR, dList, m, bestPeak, cut, Title){
	'
	Called in EMnormalize(object)
	Visualization of the mixture model
	'
	dLR <- density(LR)
	currentPlot = xyplot(dLR$y~dLR$x, type = "n", main = Title,
						xlab = "Log2R", ylab = "Density",
						#xlim = range(-0.6, 0.6), ylim = range(0, max(dLR$y)*1.5),
						xlim = cut, ylim = range(0, max(dLR$y)*1.5),
						panel = function(x, y){
											lpolygon(dLR$x, dLR$y, col = 'grey90')
											n = length(LR)
											nG = length(dList)
											for (i in 1:nG){
												tmp <- dList[[i]]
												llines(tmp$x, tmp$y, lwd = 1, col = rgb(i/nG, 0.2, (nG-i)/nG, 0.75))
												lpolygon(tmp$x, tmp$y, col = rgb(i/nG, 0.2, (nG-i)/nG, 0.25))
												ltext(	x = mean(tmp$x), y = min(max(tmp$y*1.5, na.rm = TRUE), max(dLR$y)*1.25), labels = round(m[i], 3),
														cex = ifelse(i == bestPeak, 1.5, 1.25), font =  ifelse(i == bestPeak, 2, 1))
												}
									}
					)
	return(currentPlot)
}

# To Do

#################
# QCSegm

#################
MedSegm <- function(seg.cna.obj, cnSet){
	'
	Called by SegmentCGH()
	'
	seg.start <- seg.cna.obj$output$loc.start
	seg.end <- seg.cna.obj$output$loc.end			
	seg.len <- seg.cna.obj$output$num.mark
	cum.seg.len <- cumsum(seg.len)
	s <- cnSet$ChrStart[cum.seg.len]
	e <- cnSet$ChrEnd[cum.seg.len]
	if(!is.null(e)) seg.end <- seg.end + (e-s)


	seg.med <- c()
	for(i in 1:length(seg.start)){
		index <- which(cnSet$genomicPos>= seg.start[i] & cnSet$genomicPos<= seg.end[i])
		tmpLR = cnSet$Log2Ratio[index]
		tmpMed <- tukey.biweight(tmpLR[!is.na(tmpLR)])
		seg.med <- c(seg.med, tmpMed)
		}
	seg.cna.obj$output <- cbind.data.frame(seg.cna.obj$output, seg.med = seg.med)
	return(seg.cna.obj)
}

#################
AddSegments <- function(seg.cna.obj, cnSet, cutGainLoss = 10, use.medians = TRUE){
	'
	Called by SegmentCGH()
	'
	seg.start <- seg.cna.obj$output$loc.start
	seg.end <- seg.cna.obj$output$loc.end			# !!!! verifier seg-end
	seg.mean <- seg.cna.obj$output$seg.mean
	if(use.medians) seg.mean <- seg.cna.obj$output$seg.med
	seg.len <- seg.cna.obj$output$num.mark
	cum.seg.len <- cumsum(seg.len)
	s <- cnSet$ChrStart[cum.seg.len]
	e <- cnSet$ChrEnd[cum.seg.len]
	if(!is.null(e)) seg.end <- seg.end + (e-s)

	# Proportion d'aberrations
	cutoff <- log2(1 + cutGainLoss/100)
	seg.delta <- seg.end - seg.start
	Prop <- sum(seg.delta[which(abs(seg.mean)>=cutoff)])/sum(seg.delta)*100
	Prop <- round(Prop, 2)

	seg.val <- rep(NA, nrow(cnSet))
	for(i in 1:length(seg.mean)){
		index <- which(cnSet$genomicPos>= seg.start[i] & cnSet$genomicPos<=seg.end[i])
		seg.val[index] <- seg.mean[i]
		}
	cnSet <- cbind.data.frame(cnSet, Segm = seg.val)
	return(cnSet)
	}

#################
dLRsd <- function(LR){
	'
	Not used yet
	'
	n <- length(LR)
	V1 <- LR[-1]
	V2 <- LR[-n]
	dLR <- V2-V1
	q1 <- quantile(dLR, 0.25, na.rm = TRUE)
	q3 <- quantile(dLR, 0.75, na.rm = TRUE)
	s <- sd(dLR[which(dLR > q1 & dLR < q3)], na.rm = TRUE)/sqrt(2)
	return(s)
	}

#################
# Gene list functions
#################

# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi

# esearch: return pubmed Ids
eSearch <- function (term, n){
	'
	Not used yet
	'
  # term = keywords. Same form as in a PubMed query
  # n = max number of pubmed Ids (papers) to return
  srch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
  srch.mode <- paste('db=pubmed&retmax=', n, '&retmode=xml&term=', sep = "")
  doc <-xmlTreeParse(paste(srch.stem,srch.mode,term,sep=""), isURL = TRUE, useInternalNodes = TRUE)
  sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
}

# gSearch: return gene Ids
gSearch <- function (geneSymb, database){  
  '
	Called by geneRequest.v7()
  '
	require(XML)
  # ! This function can return more than one Id ! 
  gsrch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
  gsrch.mode <- paste0("db=", database, "&retmode=xml","&term=")
  URL <- paste0(gsrch.stem, gsrch.mode, geneSymb)
  doc <- xmlTreeParse(URL, isURL = TRUE, useInternalNodes = TRUE)
  sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
}

getItem <- function(doc, item){
  '
	Called by gSummary()
  '
	expr = paste0('<Item Name=\"', item, '\"')
  if(any(grepl(expr, doc))){
    String = doc[grep(expr, doc)]
    r = regexec('>(.*?)<', String)
    return(unlist(regmatches(String, r))[2])
  }
  return(NA)
}

gSummary <- function(id, database, itemList){
  '
	Called by geneRequest.v7()
  '
	require(XML)
	sum.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
  	sum.mode <- paste0("db=", database, "&id=")
  	urlFile = url(paste0(sum.stem, sum.mode, id), 'r')
  	doc = readLines(urlFile)
  	close(urlFile)
  	geneInfo = c()
  	foreach(item = iter(itemList)) %do% {
    	geneInfo = cbind(geneInfo, getItem(doc, item))
    	colnames(geneInfo)[length(geneInfo)] = item
  		}
  return(as.data.frame(geneInfo))
}


geneRequest.v7 <- function(geneId, database, verbose = TRUE){
	itemList = c(	'Orgname', 'NomenclatureStatus', 
					'NomenclatureSymbol', 'OtherAliases', 'Description',
             		'Chromosome', 'MapLocation', 'ChrStart', 'ChrStop')
  	output = data.frame()
  	notFound = cbind.data.frame(NA, t(rep(NA, length(itemList))))
    colnames(notFound) = c(itemList, 'entrezgene')
    if(is.character(geneId)){
    	geneId = toupper(geneId)
	    id = gSearch(paste0(geneId, '[symbol]+homo+sapiens[Organism]'), database)
	    id = unlist(id)
	    }
    if(length(id) == 0){
      if(verbose) cat('\n', geneId, '\t*** not found ***')
      output = rbind.data.frame(output, notFound)
    	}
    else{
      # should have NomenclatureStatus = 'Official'
      Official = FALSE
      k = 1
      while (!Official & k <= length(id)){
        tmp = gSummary(paste0(id[k], '%5BUID%5D'), database, itemList)
        tmp = cbind.data.frame(tmp, entrezgene = id[k])
        Official <- tmp$NomenclatureStatus == 'Official'
        k = k + 1
      	}
      if(!Official) tmp = notFound
      output = rbind.data.frame(output, tmp)
      if(verbose){
      	if (is.character(geneId)) cat('\n', geneId, 'found. entrezgene:', as.character(output$entrezgene))
      	else cat('\n', geneId, 'found as', as.character(output$NomenclatureSymbol))
      	}
    }
 
	# Add genomic position.
  	arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
	hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
	cumLen = cumsum(as.numeric(hg19$length))
	cumLen = c(0, cumLen[-length(cumLen)])
	output = cbind.data.frame(output[,-ncol(output)], genomicStart = rep(NA, nrow(output)), genomicStop = rep(NA, nrow(output)), entrezgene = output[,ncol(output)])
	chr = as.numeric(as.character(output$Chromosome))
	output$genomicStart = as.numeric(as.character(output$ChrStart))+ cumLen[chr]
	output$genomicStop = as.numeric(as.character(output$ChrStop)) + cumLen[chr]
  	#cat('\n')
  	#rownames(output) = seq(1, nrow(output))
  	return(output)
}



####################
#	HTML Table functions
####################

keggLink <- function(id){
	return(paste0('http://www.genome.jp/dbget-bin/www_bget?hsa:', id))
}
ncbiLink <- function(id){
	return(paste0('http://www.ncbi.nlm.nih.gov/gene?term=', id, '%5BUID%5D'))
}
sangerLink <- function(id){
	return(paste0('http://www.sanger.ac.uk/perl/genetics/CGP/cosmic?action=gene&ln=', id))
}
pharmGkbLink <- function(id){
	return(paste0('http://www.pharmgkb.org/gene/', id))
}
keggToDrugLink <- function(id){
	return(paste0('http://www.genome.jp/kegg-bin/get_htext?htext=br08303_target.keg&query=', id))
}
CTDbaseLink <- function(id){
	# Use entrezgene as id
	return(paste0('http://ctdbase.org/detail.go?type=gene&acc=', id))
}
clinTrialsLink <- function(id){
	return(paste0('http://clinicaltrials.gov/ct2/results?term=', id))
}
htmlLink <- function(tag, link){
	return(paste0("<font size='4'><a href=", link, " target=_blank >", as.character(tag), "</a></font>"))
}

html_css <- function(filename){
		write.table("<style type=\"text/css\">", file = filename, quote=F, append=T, col.names=F, row.names=F)
		write.table("table{width:100%;height:10px;background-color:#FBFBEF;}", file = filename, quote=F, append=T, col.names=F, row.names=F)
		write.table("tr.firstline{background-color:#FFBD9D;}", file = filename, quote=F, append=T, col.names=F, row.names=F)					# fond entêtes de sous-tables
		write.table("caption.captiondataframe{font-style:italic;font-size:15pt;}", file = filename, quote=F, append=T, col.names=F, row.names=F)	# fond entêtes de sous-tables
		write.table("a:link{text-decoration:none;color:blue;}", file = filename, quote=F, append=T, col.names=F, row.names=F)					# couleur du lien
		write.table("a:visited{text-decoration:none;color:#8A008A;}", file = filename, quote=F, append=T, col.names=F, row.names=F)		# couleur du lien après activation
		write.table("a:hover{text-decoration:underline;color:red;}", file = filename, quote=F, append=T, col.names=F, row.names=F)			# couleur du lien au passage de souris
		write.table("h2{background-color:#FFA366;text-align:center;}", file = filename, quote=F, append=T, col.names=F, row.names=F)	# fond entête principal
		write.table("span{font-weight:bold;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
		write.table("#Norm{color:#A1A0A0;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
		write.table("#Gain{color:#3075ED;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
		write.table("#Ampli{color:#1100A6;}",file=filename, quote=F, append=T,col.names=F, row.names=F)	
		write.table("#Loss{color:#E50810;}",file=filename, quote=F, append=T,col.names=F, row.names=F)												# 
		write.table("</style> ",file = filename, quote=F, append=T,col.names=F, row.names=F)
	}
	
textStyle <- function(x, align = 'center'){
	#return(paste0('<font size="4", align="center">', x, "</font>"))
	return(paste0("<p style=font-size:18;text-align:", align, ">", x, "</p>"))
}

logStyle<- function(LR, thresh){
	# Lower = red
	styleL = function(x){paste0("<p style=color:#E62E00;font-size:18;font-weight:bold;text-align:center>", round(x, 3), "</p>")}
	# Normal = grey
	styleN = function(x){paste0("<p style=color:grey;font-size:18;text-align:center>", round(x, 3), "</p>")}
	# Higher = blue
    styleH = function(x){paste0("<p style=color:#0000FF;font-size:18;font-weight:bold;text-align:center>", round(x, 3), "</p>")}
    LR  = ifelse(LR<(-thresh), styleL(LR), ifelse(LR>thresh, styleH(LR), styleN(LR)))
	return(LR)
}

formatTable <- function(Table){
	Table$Chr = as.numeric(as.character(Table$Chr))
	Table$chrStart = as.numeric(as.character(Table$chrStart))
	Table$chrStop = as.numeric(as.character(Table$chrStop))
	Table = Table[order(Table$Chr, Table$chrStart),]
	isSymbol = grep('Symbol', colnames(Table))
	desc = grep('Description', colnames(Table))
	chr = grep('Chr', colnames(Table))
	mapLoc = grep('Cytoband', colnames(Table))
	cStart = grep('genomicStart', colnames(Table))
	cStop = grep('genomicStop', colnames(Table))
	toNcbi <- htmlLink(Table[,isSymbol], ncbiLink(Table$entrezgene))
	toKegg <- htmlLink(Table[,isSymbol], keggLink(Table$entrezgene))
	Table = data.frame('Symbol' = toNcbi,
												'Description' = textStyle(Table[,desc], align = 'left'),
												'Chr' = textStyle(Table[,chr]),
												'MapLoc' = textStyle(as.character(Table[,mapLoc])),
												'ChrStart (Kb)' = textStyle(round(Table[,cStart]/1e3)),
												'ChrStop (Kb)' = textStyle(round(Table[,cStop]/1e3)),
												'entrezgene' = textStyle(Table$entrezgene),
												'Log2Ratio' = logStyle(Table$Log2Ratio, 0.1),
												check.names = FALSE)
	#colnames(Table)[5:6] = c('ChrStart (Kb)', 'ChrEnd (Kb)')
	return(Table)	
}

buildHtml <- function(geneTable, filePath, fileName, GeneOfInt = TRUE){
	Table <- formatTable(geneTable)
	if(GeneOfInt){
		TitleName <- paste0("<h1 align=center>Genes of Interest</h1><h2 align=center>Sample ", fileName, '</h2>')
		output = HTMLInitFile(filePath, filename = fileName, Title = 'Genes Of Interest')
		HTML("<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />", file = output)
		HTML(as.title(TitleName), file=output)
		HTML(Table, file = output, row.names = F, innerBorder = 1, caption = NULL, captionalign = "top")
		html_css(output)
		HTMLEndFile(output)
		}	
	else{
		TitleName <- paste0("Something else", fileName, '</h2>')
		}
}


FullHtml <- function(object, gain = log2(2.25/2), loss = log2(1.75/2),
								Select = c("Both", "Gain", "Loss"), CHR = 1:23,
								sangercensus = TRUE, pharmgkb = TRUE, keggtodrug = TRUE, ctdbase = TRUE, clintrials = TRUE, Restrict = TRUE,
								filePath){
		# Load information tables
		# arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
		cat('Sanger:', sangercensus, 'pharmgKB:', pharmgkb, 'keggToDrug:',keggtodrug,
			'ctdbase:', ctdbase, 'clintrails:', clintrials, '\n')
			
		# Replace with myGeneDB2
		# geneDB <- readRDS(paste0(arrayInfoPath, 'myGeneDB2.rds'))
		#if(sangercensus) census.db <- read.table(paste0(arrayInfoPath, "CensusTable_Annot_FC.txt"), header = T, sep = "\t")
		#if(pharmgkb) PGKB.db <- read.csv(paste0(arrayInfoPath, "relationships.tsv"), header = T, sep = "\t")
		#if(keggtodrug) KeggToDrug.db <- read.csv(paste0(arrayInfoPath, "Kegg_Drugs_Table.txt"), header = T, sep = "\t")
		#if(clintrials) clinTrials.db <- read.csv(paste0(arrayInfoPath, "ClinicalTrials_Drugs_Table.txt"), header = T, sep = "\t")
		#if(ctdbase){
		#	CTDb <- read.csv(paste0(arrayInfoPath, "CTD_chem_gene_ixns_SansDuplic.txt"), header = T, sep = "\t")
		#	CTDbList <- CTDb$GeneSymbol
		#	}	

		# Define what to consider
		segTable <- getSegTable(object)
		selectTable <- segTable[which(segTable$chrom %in% CHR),]
		SegmentsValues <- selectTable$seg.med
		Select <- match.arg(Select)
		switch(Select,
					Both = (whichSeg = which(SegmentsValues >= gain | SegmentsValues <= loss)),
					Gain = (whichSeg = which(SegmentsValues >= gain)),
					Loss = (whichSeg = which(SegmentsValues <= loss))
					)
		cat(Select, '\n')
		gainSubTitl <- paste0("gained (>", round(gain, 3), ")")
		lossSubTitl <- paste0("lost (<", round(loss, 3), ")")
		switch(Select,
					Gain = (selecType = gainSubTitl),
					Loss = (selecType = lossSubTitl),
					Both = (selecType = paste(gainSubTitl, "or", lossSubTitl))
					)
		Explore <- paste("Chr", CHR[1], "to", CHR[length(CHR)], sep = "")

		# Initialize html
		fileName <- paste(getInfo(object, 'sampleId'), getInfo(object, 'barCode'), getInfo(object, 'analyseDate'), getInfo(object, 'platform'), Explore, Select, "SupplTable", sep = "_")
		titleName <- paste("Supplementary tables: ", getInfo(object, 'sampleId'), " / ", getInfo(object, 'barCode'), "<br>",
						getInfo(object, 'platform'), " / ", getInfo(object, 'analyseDate'), " / ", paste("Chr", CHR[1], "to", CHR[length(CHR)]), ": ", paste("segments", selecType))
		output <- HTMLInitFile(filePath, filename = fileName, Title = titleName)
		HTML(as.title(titleName), file = output)
		HTML("<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />", file = output)

		if(length(whichSeg)<=1) Append = FALSE
		k = 1
		for(segment in whichSeg){	#whichSeg
			cat('segment', k, 'of', length(whichSeg), '\t')
			subTable <- selectTable[segment,]	#[segment,]
			chrom <- subTable$chrom
			startAt = subTable$loc.start
			stopAt = subTable$loc.end
			segValue = subTable$seg.med
			cat(', Segment value:', segValue,'\n')
			tmpTable <- .createSubTable(subTable, sangercensus = sangercensus, pharmgkb = pharmgkb,
														keggtodrug = keggtodrug, ctdbase = ctdbase, clintrials = clintrials)
			if(!is.null(tmpTable)){
				if(Restrict){
					cat('Restric...\t')
					ofInterest <- which(tmpTable$SangerCancer != '-' | tmpTable$PgKB != '-' | tmpTable$KeggToDrugs != '-' |
												tmpTable$CTDbase != '-' | tmpTable$ClinicalTrials != '-')
					#ofInterest <- which(!is.na(tmpTable$SangerCancer) | !is.na(tmpTable$PgKB) | !is.na(tmpTable$KeggToDrugs) |
					#							!is.na(tmpTable$CTDbase) | !is.na(tmpTable$ClinicalTrials))
					tmpTable <- tmpTable[ofInterest,]
					cat(length(ofInterest), 'genes\n')
					}
				if(nrow(tmpTable)>0){
					for(i in 4:ncol(tmpTable))
						tmpTable[,i] <- as.data.frame(paste0('<div style="text-align:center;">', tmpTable[,i], '</div>'))
					HTML(tmpTable, file = output, row.names = F, innerBorder = 1,
							caption = paste('<b >Segment on Chr', chrom, ": ", round(startAt/1e6, 3), "-",
											round(stopAt/1e6, 3), "(Mb), Log2R = ", signif(segValue, 3), ", nb genes = ",
											nrow(tmpTable),"</b>", sep = ""),
							captionalign = "top")
					}
				}
			k = k + 1
			cat('\n')
			}
		html_css(output)
		HTMLEndFile(output)
}

.createSubTable <- function(subTable, sangercensus, pharmgkb, keggtodrug, ctdbase, clintrials){
			chr = subTable$chrom
			startAt = subTable$loc.start
			stopAt = subTable$loc.end
			segValue = subTable$seg.med
			# Get the list of genes included in segment
			cat('In geneDB:\t')
			tmpTable <- geneDB[which(geneDB$genomicStart >= startAt & geneDB$genomicStop <= stopAt |
													geneDB$genomicStop >= startAt & geneDB$genomicStart <= stopAt),
											-c(1, 5, 10:12)]
			cat(nrow(tmpTable), '\n')
			if(nrow(tmpTable) > 0){
				Symbols = as.character(tmpTable$Symbol)
				entrezgene = as.character(tmpTable$entrezgene)
				# Add first html entregene & html KeggPathway.
				tmpTable$Symbol <- htmlLink(Symbols, ncbiLink(entrezgene))
				tmpTable <- cbind.data.frame(Symbol = htmlLink(Symbols, ncbiLink(entrezgene)), Kegg = htmlLink(Symbols, keggLink(entrezgene)), tmpTable[,-1])

				if(sangercensus){
					censusIds <- filtrSangerCensus(Symbols, census.db)
					tmpTable <- cbind.data.frame(tmpTable, SangerCancer = ifelse(!is.na(censusIds), htmlLink(censusIds, sangerLink(censusIds)), "-"))
					}
					if(pharmgkb){
						pgkbIds <- filtrPharmGkb(Symbols, PGKB.db)
						tmpTable <- cbind.data.frame(tmpTable, PgKB = ifelse(!is.na(pgkbIds), htmlLink(pgkbIds, pharmGkbLink(pgkbIds)), "-"))
						}
					if(keggtodrug){
						keggIds <- filtrKeggToDrug(Symbols, KeggToDrug.db)				
						tmpTable <- cbind.data.frame(tmpTable, KeggToDrugs = ifelse(!is.na(keggIds), htmlLink(keggIds, keggToDrugLink(keggIds)), "-"))
						}
					if(ctdbase){ # should return entrezgenes but not symbols
						ctdIds <- filtrCTDbase(Symbols, CTDbList)
						tmpTable <- cbind.data.frame(tmpTable, CTDbase = ifelse(!is.na(ctdIds), htmlLink(Symbols, CTDbaseLink(entrezgene)), "-"))
						}
					if(clintrials){
						clinTrialsIds <- filtrClinTrials(Symbols, clinTrials.db)
						tmpTable <- cbind.data.frame(tmpTable, ClinicalTrials = ifelse(!is.na(clinTrialsIds), htmlLink(clinTrialsIds, clinTrialsLink(clinTrialsIds)), "-"))
						}
					#cat('Return tmpTable from create\n')
					return(tmpTable)
					}
				return(NULL)
}

arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
census.db <- read.table(paste0(arrayInfoPath, "CensusTable_Annot_FC.txt"), header = T, sep = "\t")
PGKB.db <- read.csv(paste0(arrayInfoPath, "relationships.tsv"), header = T, sep = "\t")
PGKB.db <- PGKB.db[grepl('Drug', PGKB.db$Entity2_id),]
#PGKB.db <- PGKB.db[order(PGKB.db$Entity1_name),]
KeggToDrug.db <- read.csv(paste0(arrayInfoPath, "Kegg_Drugs_Table.txt"), header = T, sep = "\t")
KeggToDrug.db <- KeggToDrug.db[KeggToDrug.db$nDrug > 0,]
CTDb <- read.csv(paste0(arrayInfoPath, "CTD_chem_gene_ixns_SansDuplic.txt"), header = T, sep = "\t")
CTDbList <- CTDb$GeneSymbol
clinTrials.db <- read.csv(paste0(arrayInfoPath, "ClinicalTrials_Drugs_Table.txt"), header = T, sep = "\t")
clinTrials.db <- clinTrials.db[clinTrials.db$CTfound > 0,]

filtrSangerCensus <- function(Symbols, census.db){
	cat('Searching in Sanger Cancer...\t')
	#censusIds <- ifelse(as.character(Symbols) %in% census.db$Symb, Symbols, NA)
	censusIds <- Symbols
	censusIds[!is.element(Symbols, census.db$Symb)] <- NA
	cat(sum(!is.na(censusIds)), '\n')
	return(as.character(censusIds))
}

# filtrPharmGkb <- function(Symbols, PGKB.db){
	# cat('Searching in PharmgKB...\t')
	# # Search for PGKB annotations
	# pgkbIds <- rep("-", length(Symbols))
	# pgkbGenes <- ifelse(as.character(Symbols) %in% PGKB.db$Entity1_name, as.character(Symbols), NA)
	# pgkbIds <- rep(NA, length(pgkbGenes))
	# if(any(!is.na(pgkbGenes))){
		# for(pg in 1:length(pgkbGenes))
			# if(!is.na(pgkbGenes[pg])){
				# tmp <- as.character(unique(PGKB.db$Entity1_id[which(PGKB.db$Entity1_name == pgkbGenes[pg])]))
				# pgkbIds[pg] <- substr(tmp, 6, 50)
				# }
		# }
	# cat(sum(!is.na(pgkbIds)), '\n')
	# return(as.character(pgkbIds))	
# }

filtrPharmGkb <- function(Symbols, PGKB.db){
	cat('Searching in PharmgKB...\t')
	# Search for PGKB annotations, retunrs PharmgKB ids
	pgkbIds <- lapply(Symbols, function(gene){
		if(is.element(gene, PGKB.db$Entity1_name)){
			entId <- as.character(unique(PGKB.db$Entity1_id[PGKB.db$Entity1_name == gene]))
			gsub('Gene:', '', entId)
		}
		else NA
	})
	pgkbIds <- do.call(c, pgkbIds)
	cat(sum(!is.na(pgkbIds)), '\n')
	return(as.character(pgkbIds))	
}

# filtrKeggToDrug <- function(Symbols, KeggToDrug.db){
	# cat('Searching in Kegg-Drugs...\t')
	# KeggIds <- rep(NA, length(Symbols))
	# is.Kegg <- ifelse(as.character(Symbols) %in% KeggToDrug.db$Symb, Symbols, NA)
	# KeggIds <- apply(as.data.frame(is.Kegg), 1,
						# function(x){
							# if(as.character(x) %in% Symbols){
								# kegg.index <- which(KeggToDrug.db$Symb == x)
								# if(KeggToDrug.db$nDrug[kegg.index]!=0) return (KeggToDrug.db$GeneId[kegg.index])
								# else return(NA)
								# }
							# else return (NA)
						# })
	# cat(sum(!is.na(KeggIds)), '\n')
	# return(as.character(KeggIds))
# }

filtrKeggToDrug <- function(Symbols, KeggToDrug.db){
	cat('Searching in Kegg-Drugs...\t')
	#KeggIds <- rep(NA, length(Symbols))
	#is.Kegg <- ifelse(as.character(Symbols) %in% KeggToDrug.db$Symb, Symbols, NA)
	KeggIds <- lapply(Symbols, function(gene){
							if(is.element(gene, KeggToDrug.db$Symb))
								KeggToDrug.db$GeneId[KeggToDrug.db$Symb == gene]
							else NA
						})
	KeggIds <- do.call(c, KeggIds)
	cat(sum(!is.na(KeggIds)), '\n')
	return(as.character(KeggIds))
}

filtrCTDbase <- function(Symbols, CTDbList){
	cat('Searching in CTD...\t')
	# ctdIds <- rep(NA, length(Symbols))
	# ctdIds <- ifelse(as.character(Symbols) %in% CTDbList, Symbols, NA)
	ctdIds <- Symbols
	ctdIds[!is.element(Symbols, CTDbList)] <- NA
	cat(sum(!is.na(ctdIds)), '\n')
	return(as.character(ctdIds))
}

# filtrClinTrials <- function(Symbols, clinTrials.db){
	# cat('Searching in clinical Trials...\t')
	# ctIds <- c()
	# foreach(gene = iter(Symbols)) %do%{
		# iselem <- is.element(gene, clinTrials.db$Symb)
		# if(iselem) found = clinTrials.db$CTfound[which(clinTrials.db$Symb == gene)]
		# else found = 0
		# if(iselem & found != 0)
			# ctIds <- c(ctIds, gene)
		# else ctIds <- c(ctIds, NA)
	# }
	# cat(sum(!is.na(ctIds)), '\n')
	# return(ctIds)
# }

# filtrClinTrials <- function(Symbols, clinTrials.db){
	# cat('Searching in clinical Trials...\t')
	# ctIds <- c()
	# foreach(gene = iter(Symbols)) %do%{
		# if(is.element(gene, clinTrials.db$Symb))
			# ctIds <- c(ctIds, gene)
		# else ctIds <- c(ctIds, NA)
	# }
	# cat(sum(!is.na(ctIds)), '\n')
	# return(ctIds)
# }

filtrClinTrials <- function(Symbols, clinTrials.db){
	cat('Searching in clinical Trials...\t')
	ctIds <- Symbols
	ctIds[!is.element(Symbols, clinTrials.db$Symb)] <- NA
	cat(sum(!is.na(ctIds)), '\n')
	return(ctIds)
}
