# gSearch: return gene Ids
gSearch <- function (geneSymb, database){  
  '
	Called by geneRequest.v7()
  '
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

gSummary <- function(id, database){
  '
	Called by geneRequest.v7()
  '
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

itemList = c(	'Orgname', 'NomenclatureStatus', 
					'NomenclatureSymbol', 'OtherAliases', 'Description',
             		'Chromosome', 'MapLocation', 'ChrStart', 'ChrStop')

geneRequest.v7 <- function(geneList, database, verbose = TRUE){
  output = data.frame()
  foreach(gene = iter(geneList)) %do% {
    notFound = cbind.data.frame(gene, t(rep(NA, length(itemList))))
    colnames(notFound) = c(itemList, 'entrezgene')
    id = gSearch(paste0(gene, '[symbol]+homo+sapiens[Organism]'), database)
    id = unlist(id)
    if(length(id) == 0){
      if(verbose) cat('\n', gene, '\t*** not found ***')
      output = rbind.data.frame(output, notFound)
    }
    else{
      # should have NomenclatureStatus = 'Official'
      Official = FALSE
      k = 1
      while (!Official & k <= length(id)){
        tmp = gSummary(paste0(id[k], '%5BUID%5D'), database)
        tmp = cbind.data.frame(tmp, entrezgene = id[k])
        Official <- tmp$NomenclatureStatus == 'Official'
        k = k + 1
      	}
      if(!Official) tmp = notFound
      output = rbind.data.frame(output, tmp)
      if(verbose) cat('\n', gene, 'found')
    }
  }
  	arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
	hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
	cumLen = cumsum(as.numeric(hg19$length))
	cumLen = c(0, cumLen[-length(cumLen)])
	output = cbind.data.frame(output[,-ncol(output)], genomicStart = rep(NA, nrow(output)), genomicStop = rep(NA, nrow(output)), entrezgene = output[,ncol(output)])
	foreach(i = 1:nrow(output)) %do%{
		chr = as.numeric(output$Chromosome[i])
		output$genomicStart[i] = as.numeric(as.character(output$ChrStart[i]))+ cumLen[chr]
		output$genomicStop[i] = as.numeric(as.character(output$ChrStop[i])) + cumLen[chr]
		}
  cat('\n')
  rownames(output) = seq(1, nrow(output))
  return(output)
}


geneRequest.v7 <- function(geneId, database, verbose = TRUE){
  	output = data.frame()
  	notFound = cbind.data.frame(geneId, t(rep(NA, length(itemList))))
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
        tmp = gSummary(paste0(id[k], '%5BUID%5D'), database)
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
  	cat('\n')
  #rownames(output) = seq(1, nrow(output))
  	return(output)
}

gene = 'EGFR'
database = 'gene'

geneRequest.v7('egfr', 'gene')
