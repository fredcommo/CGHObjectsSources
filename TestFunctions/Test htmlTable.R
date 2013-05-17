keggLink <- function(id){
	return(paste0('http://www.genome.jp/dbget-bin/www_bget?hsa:', id))
}
ncbiLink <- function(id){
	return(paste0('http://www.ncbi.nlm.nih.gov/gene?term=', id, '%5BUID%5D'))
}
htmlLink <- function(tag, link){
	return(paste0("<font size='4'><a href=", link, " target=_blank >", as.character(tag), "</a></font>"))
}

html_css <- function(filename){
		write.table("<style type=\"text/css\">", file = filename, quote=F, append=T, col.names=F, row.names=F)
		write.table("table{width:100%;height:10px}", file = filename, quote=F, append=T, col.names=F, row.names=F)
		write.table("tr.firstline{background-color:#FFBD9D;}", file = filename, quote=F, append=T, col.names=F, row.names=F)				# fond entêtes de sous-tables
		write.table("a:link{text-decoration:none;color:blue;}", file = filename, quote=F, append=T, col.names=F, row.names=F)				# couleur du lien
		write.table("a:visited{text-decoration:none;color:#8A008A;}", file = filename, quote=F, append=T, col.names=F, row.names=F)			# couleur du lien après activation
		write.table("a:hover{text-decoration:underline;color:red;}", file = filename, quote=F, append=T, col.names=F, row.names=F)			# couleur du lien au passage de souris
		write.table("h2{background-color:#FFA366;text-align:center;}", file = filename, quote=F, append=T, col.names=F, row.names=F)		# fond entête principal
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
	# Lower
	styleL = function(x){paste0("<p style=color:#0000FF;font-size:18;font-weight:bold;text-align:center>", round(x, 3), "</p>")}
	# Normal 
    styleN = function(x){paste0("<p style=color:grey;font-size:18;text-align:center>", round(x, 3), "</p>")}
	# Higher
    styleH = function(x){paste0("<p style=color:#E62E00;font-size:18;font-weight:bold;text-align:center>", round(x, 3), "</p>")}
    LR  = ifelse(LR<(-thresh), styleL(LR), ifelse(LR>thresh, styleH(LR), styleN(LR)))
	return(LR)
}

buildHtml <- function(geneTable, filePath, fileName){ #, cssFile
	toNcbi <- htmlLink(geneTable$Symbol, ncbiLink(geneTable$entrezgene))
	toKegg <- htmlLink(geneTable$Symbol, keggLink(geneTable$entrezgene))
#	TitleName <- paste0("Genes of Interest<br>Sample: ", fileName, "<br>")
	TitleName <- paste0("<h1 align=center>Genes of Interest</h1><h2 align=center>Sample ", fileName, '</h2>')
	output = HTMLInitFile(filePath, filename = fileName, Title = 'Gene Of Interest')
	Table = cbind.data.frame(Symbol = toNcbi,
												Description = textStyle(geneTable$Description, align = 'left'),
												Chr = textStyle(geneTable$Chr),
												MapLocation = textStyle(geneTable$Cytoband),
												ChrStart = textStyle(round(geneTable$chrStart/1e3)),
												ChrStop = textStyle(round(geneTable$chrStop/1e3)),
												entrezgene = textStyle(geneTable$entrezgene),
												Log2Ratio = logStyle(geneTable$Log2Ratio, 0.1))
	colnames(Table)[5:6] = c('ChrStart (Kb)', 'ChrEnd (Kb)')
	HTML("<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />", file = output)
	HTML(as.title(TitleName), file=output)
	HTML(Table, file = output, row.names = F, innerBorder = 1, caption = NULL, captionalign = "top")#, CSSFile = cssFile)
	html_css(output)
	HTMLEndFile(output)		
}

	# create and save the table
geneList = c("CCND1", "ALK", "MDM2", "FRS2", "MET", "RPTOR", 'EGFR')
geneTable <- geneOfInt(object5, geneList)
for(Col in c('Chromosome', 'ChrStart', 'ChrStop', 'genomicStart', 'genomicStop'))
	geneTable[,colnames(geneTable) == Col] <- as.numeric(as.character(geneTable[,colnames(geneTable) == Col]))
	
setwd('/Users/fredcommo/Documents/Projet Safir')
filePath = paste0('./Output/', getInfo(object5, 'sampleId'))
fileName = paste0(getInfo(object5, 'sampleId'))
buildHtml(geneTable, filePath, fileName)


