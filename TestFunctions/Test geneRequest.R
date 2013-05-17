

# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi
require(XML)
require(foreach)
require(iterators)


# esearch: return pubmed Ids
eSearch <- function (term, n){
  # term = keywords. Same form as in a PubMed query
  # n = max number of pubmed Ids (papers) to return
  srch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
  srch.mode <- paste0('db=pubmed&retmax=', n, '&retmode=xml&term=')
  doc <-xmlTreeParse(paste(srch.stem,srch.mode,term,sep=""), isURL = TRUE, useInternalNodes = TRUE)
  sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
}

itemList = c('NomenclatureSymbol', #'NomenclatureName',
				'Description',
				'NomenclatureStatus', 'OtherAliases', 'Chromosome',
				'MapLocation', 'ChrStart', 'ChrStop', 'Orgname')

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

geneRequest.v7 <- function(geneList, database, verbose = TRUE){
	output = data.frame()
	foreach(gene = iter(geneList)) %do% {
		notFound = cbind.data.frame(gene, t(rep(NA, length(itemList))))
		colnames(notFound) = c(itemList, 'geneId')
		id = gSearch(paste0(gene, '[symbol]+homo+sapiens[Organism]'), database)
		id = unlist(id)
			if(length(id) == 0){
				cat('\n', gene, '\t*** not found ***')
				output = rbind.data.frame(output, notFound)
						}
				else{
					# should have NomenclatureStatus = 'Official'
					Official = FALSE
					k = 1
					while (!Official & k <= length(id)){
						tmp = gSummary(paste0(id[k], '%5BUID%5D'), database)
						tmp = cbind.data.frame(tmp, geneId = id[k])
						Official <- tmp$NomenclatureStatus == 'Official'
						k = k + 1
						}
					if(!Official) tmp = notFound
					output = rbind.data.frame(output, tmp)
					cat('\n', gene, 'found')
					}
	}
	cat('\n\n')
	rownames(output) = seq(1, nrow(output))
	return(output)
}

myGenes = c('EGFR', 'ERCC1', 'BRCA1', 'FGFR1', 'A_gene_Called_Fred', 'BRCA2', 'FGFR2', 'SageBionetworks', 'MYC', 'PTEN', 'RPTOR', 'FGF4')
geneRequest.v7(myGenes, 'gene')

