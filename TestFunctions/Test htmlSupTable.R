EditSupTab.v7.6 <- function(object, ReqFunc = Request, use.medians = TRUE, 
							Select = c("Gain", "Loss", "Both"), ThreshGain = 1, ThreshLoss = -1, CHR = 1:23, Restrict = FALSE,
							sangercensus = T, pharmgkb = T, keggtodrug = T, ctdbase = T, clintrials = T, Root = "E"){

require(tcltk)
require(R2HTML)

	# Check for system
	system <- Sys.info()["sysname"]
	workingDir <- getwd()

	##############################
	## Liste des liens partiels ##
	##############################
	linkModel <- NULL
	linkModel$Entrez = "http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids="					# Link to EntrezGene
	linkModel$Kegg = "http://www.genome.jp/dbget-bin/www_bget?hsa:"																# Link to Kegg database
	if(sangercensus) linkModel$SangerCensus = "http://www.sanger.ac.uk/perl/genetics/CGP/cosmic?action=gene&ln="				# Link to Sanger Census Cancer
	if(pharmgkb) linkModel$PharmgKB = "http://www.pharmgkb.org/gene/"															# Link to Pharmacogenomics Knowledge base
	if(keggtodrug) linkModel$KeggToDrugs = "http://www.genome.jp/kegg-bin/get_htext?htext=br08303_target.keg&query="			# Link to Kegg database 'genes to drugs'
	if(ctdbase) linkModel$CTDbase = "http://ctdbase.org/detail.go?type=gene&acc="												# Link to CTdatabase (genes/drugs information, direct or indirect interactions)
	if(clintrials) linkModel$ClinicalTrials = "http://clinicaltrials.gov/ct2/results?term="										# Link to ClinicalTrials.gov (official site for clinical trials registry)
	
	linkModel <- as.list(linkModel)								
	
	linkNames <- names(linkModel)

################################################################





################################################################


	####################
	## Bloc principal ##
	####################

	# Build path regarding system used (quite different for linux and windows)
	if(system=="Linux"){
		setwd('/Projet Safir/Data SAFIR/Arrays Infos')
		FileDir <- paste(Root, "/Projet Safir/Data SAFIR/Safir Output/Safir SupplGeneList", sep = "")
	}
	if(system=="Windows"){
		setwd(paste(Root, ":/Projet Safir/Data SAFIR/Arrays Infos", sep = ""))
		FileDir <- paste(Root, ":/Projet Safir/Data SAFIR/Safir Output/Safir SupplGeneList", sep = "")
	}
	if(system=="Darwin"){
		setwd("/Users/fredcommo/Documents/Projet Safir/Arrays Infos")
		FileDir <- "/Users/fredcommo/Documents/Projet Safir/Arrays Infos/supTables"
	}

	
	# load annotation tables
	# full <- read.table("022060_4x180K_hg19_20110628_GenomePos&GC_GNames_FC.txt", header = T, sep = "\t")
	human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	FullDB = getBM(attributes=c('hgnc_symbol', 'description', 'chromosome_name', 'band', 'start_position','end_position', 'entrezgene'), mart = human)
	FullDB$description =  gsub(' \\[.*\\]', '', FullDB$description)
	FullDB = as.data.frame(FullDB)

	if(sangercensus) census <- read.table("CensusTable_Annot_FC.txt", header = T, sep = "\t")
	if(pharmgkb) PGKB.db <- read.csv("relationships.tsv", header = T, sep = "\t")
	if(keggtodrug) KeggToDrug.db <- read.csv("Kegg_Drugs_Table.txt", header = T, sep = "\t")
	if(ctdbase){
		CTDb <- read.csv("CTD_chem_gene_ixns_SansDuplic.txt", header = T, sep = "\t")
		CTDbList <- CTDb$GeneSymbol
	}	
	if(clintrials) CT.db <- read.csv("ClinicalTrials_Drugs_Table.txt", header = T, sep = "\t")

	
	setwd(FileDir)

	# Select the sTab rows containing the Chr to explore
	sTab = getSegTable(object)
	sTab <- sTab[which(sTab$chromosome %in% CHR),]

	# Select what segment values: means or medians
	segValues <- sTab$seg.mean
	if(use.medians) segValues <- sTab$seg.med

	# What imbalances to consider: Gain, Loss or Both, and define both values and table title
	Select <- match.arg(Select)
	switch(Select,  Gain = (whichseg = which(Values >= ThreshGain)),
					Loss = (whichseg = which(Values <= ThreshLoss)),
					Both = (whichseg = which(Values >= ThreshGain | Values <= ThreshLoss))
					)
	gainSubTitl <- paste("gained (>", round(ThreshGain, 3), ")", sep = "")
	lossSubTitl <- paste("lost (<", round(ThreshLoss, 3), ")", sep = "")
	switch(Select,  Gain = (selecType = gainSubTitl),
					Loss = (selecType = lossSubTitl),
					Both = (selecType = paste(gainSubTitl, "or", lossSubTitl))
					)
					
	# Define File name and main table title
	Explore <- paste("Chr", CHR[1], "to", CHR[length(CHR)], sep = "")
	FileName <- paste(getInfo(object, 'sampleId'), getInfo(object, 'barCode'), getInfo(object, 'analysisDate'), getInfo(object, 'platform'), Explore, Select, "SupplTable", sep = "_")
	TitleName <- paste("Supplementary tables: ", getInfo(object, 'sampleId'), " / ", getInfo(object, 'barCode'), "<br>",
						getInfo(object, 'platform'), " / ", getInfo(object, 'analysisDate'), " / ", paste("Chr", CHR[1], "to", CHR[length(CHR)]), ": ", paste("segments", selecType))

	# Define progressBar
	nseg = length(whichseg)
	PB <- tkProgressBar(title = "The R'atWork BaBar", min = 0, max = nseg, width = 500)
	k = 0

	## Initialisation du html		
	html_init(init = T, filedir = FileDir, filename = FileName, titlename = TitleName)	# changer titlename = tabName


	if(length(whichseg)>1) Append <- T
	for(i in whichseg){
	
		Sys.sleep(0.1)
		# launch & increment the pBar
		k = k + 1
		chr <- sTab$ChrNum[i]
		setTkProgressBar(PB, k, label = paste("Segment #", i, " on Chr#", chr,": ", k," of ", nseg, " in progress...", sep = ""))
	
		cat(paste("Segment#", i, " on Chr", chr,": ", k," of ", nseg, sep = ""), "\n")
		start <- sTab$loc.start[i]
		end <- sTab$loc.end[i]
		seg.value <- round(Values[i], 3)

		# What genes are included in a given segment
		index.full = which(FullDB$genomicPos >= start & FullDB$genomicPos <= end)
		Symbols = as.character(unique(FullDB$hgnc_symbol[index.full]))
		geneIds = as.character(unique(FullDB$entrezgene[index.full]))
		if(any(Symbols == ""))
			Symbols <- Symbols[-which(Symbols == "")]
		if(length(Symbols)>0){
			cat("Symbols\n", Symbols, "\n")
			refseq <- c()
			for(j in 1:length(Symbols)) refseq <- c(refseq, as.character(full$RefSeq[which(full$Symbol == Symbols[j])][1]))

			# Search for PGKB annotations
			pgkb.Ids <- rep("-", length(Symbols))
			if(pharmgkb){
				pgkb.genes <- ifelse(as.character(Symbols) %in% PGKB.db$Entity1_name, as.character(Symbols), "-")
				cat("PgKB\n", ifelse(pgkb.genes!="-", pgkb.genes, "---"), "\n")
				pgkb.Ids <- rep("-", length(pgkb.genes))
				if(any(pgkb.genes != "-")){
					for(pg in 1:length(pgkb.genes))
						if(pgkb.genes[pg] != "-"){
							tmp <- as.character(unique(PGKB.db$Entity1_id[which(PGKB.db$Entity1_name == pgkb.genes[pg])]))
							pgkb.Ids[pg] <- substr(tmp, 6, 50)
						}
				}
			}
			cat("PkGB ", pgkb.Ids, "\n")
			
			# Search for Census Cancer annotations
			is.census <- rep(NA, length(Symbols))
			if(sangercensus) is.census <- ifelse(as.character(Symbols) %in% census$Symb, Symbols, NA)
			cat("Census ", is.census, "\n")
			
			# search for KeggToDrugs
			Kegg.GI <- rep(NA, length(Symbols))
			if(keggtodrug){
				is.Kegg <- ifelse(as.character(Symbols) %in% KeggToDrug.db$Symb, Symbols, NA)
				Kegg.GI <- apply(as.data.frame(is.Kegg), 1, function(x){
																		if(as.character(x) %in% Symbols){
																			kegg.index <- which(KeggToDrug.db$Symb == x)
																			if(KeggToDrug.db$nDrug[kegg.index]!=0) return (KeggToDrug.db$GeneId[kegg.index])
																			else return(NA)
																			}
																		else return (NA)
																		})
			}
			cat("KeggD ", Kegg.GI, "\n")
			
			# Search for CTDbase
			is.ctd <- rep(NA, length(Symbols))
			if(ctdbase) is.ctd <- ifelse(as.character(Symbols) %in% CTDbList, Symbols, NA)
			cat("is.ctd ", is.ctd, "\n")
			
			# search for ClinicalTrials
			is.CT <- rep(NA, length(Symbols))
			if(clintrials){
				is.CT <- ifelse(as.character(Symbols) %in% CT.db$Symb, Symbols, NA)
				CT.Symb <- apply(as.data.frame(is.CT), 1, function(x){
																		if(as.character(x) %in% Symbols){
																			CT.index <- which(CT.db$Symb == x)
																			if(CT.db$CTfound[CT.index]!=0) return (x)
																			else return(NA)
																			}
																		else return (NA)
																	})
			}
			cat("ClinTrials ", CT.Symb, "\n")

			# Search for restrictions: if Restrict = TRUE
			if(Restrict){
				of.interest <- which(!is.na(is.census) | pgkb.Ids != "-" | !is.na(Kegg.GI) | !is.na(is.ctd) | !is.na(CT.Symb))
				Symbols <- as.character(Symbols[of.interest])
				refseq <- as.character(refseq[of.interest])
				is.census <- as.character(is.census[of.interest])
				pgkb.Ids <- as.character(pgkb.Ids[of.interest])
				Kegg.GI <- Kegg.GI[of.interest]
				is.ctd <- as.character(is.ctd[of.interest])
				CT.Symb <- CT.Symb[of.interest]
				cat(length(of.interest), "\n")
			}

			# Building output table
			if(length(Symbols)>0){
				# Search for NCBI annotations
				tmp.request <- ReqFunc(as.character(Symbols), hg19.info, verbose = FALSE)
				tmp.request <- cbind.data.frame(tmp.request, RefSeq = refseq)

				tmpTable = getBM(attributes=c('hgnc_symbol', 'description', 'chromosome_name', 'band', 'start_position','end_position', 'entrezgene'),
              					filters = 'hgnc_symbol', values = Symbols, mart = human)
				tmpTable$description =  gsub(' \\[.*\\]', '', output$description)
				if(sangercensus) tmp.request <- cbind.data.frame(tmpTable, Cancer.Census.List = ifelse(!is.na(is.census), as.character(is.census), "-"))
				if(pharmgkb) tmp.request <- cbind.data.frame(tmpTable, PgKB = pgkb.Ids)
				if(keggtodrug) tmp.request <- cbind.data.frame(tmpTable, KeggToDrugs = ifelse(!is.na(Kegg.GI), as.character(Kegg.GI), "-"))
				if(ctdbase) tmp.request <- cbind.data.frame(tmpTable, CTDbase = ifelse(!is.na(is.ctd), as.character(tmp.request$GeneId), "-"))
				if(clintrials) tmp.request <- cbind.data.frame(tmpTable, ClinicalTrials = ifelse(!is.na(CT.Symb), as.character(CT.Symb), "-"))
												

				ord <- order(tmpTable$Chr.start)
				tmpTable <- tmpTable[ord,]
				
				
				# Build a list of Ids to complete http links
				GeneId <- as.list(tmp.request[,c(12, 2, 14:ncol(tmp.request))])					# Symbols is used for CliniclaTrials.gov
				print(GeneId)			

				# Edition html
				######################################
				##	Impression des data dans html	##
				######################################
				
				## Building complete links
				linkListe <- linkModel
				linkListe$Entrez <- http_link(linkListe$Entrez, tmpTable$entrezgene)
				linkListe$Kegg <- http_link(linkListe$Kegg, tmpTable$entrezgene)
				if(sangercensus) linkListe$SangerCensus <- http_link(linkListe$SangerCensus, tmpTable$Cancer.Census.List)
				if(pharmgkb) linkListe$PharmgKB <- http_link(linkListe$PharmgKB, tmpTable$PgKB)
				if(keggtodrug) linkListe$KeggToDrugs <- http_link(linkListe$KeggToDrugs, tmpTable$KeggToDrugs)
				if(ctdbase) linkListe$CTDbase <- http_link(linkListe$CTDbase, tmpTable$CTDbase)
				if(clintrials) linkListe$ClinicalTrials <- http_link(linkListe$ClinicalTrials, tmpTable$ClinicalTrials)
				
				
				## Unlist linkList to build subtables of links
				# linkListe <- sapply(linkListe, function(x){paste("<a href=", x, " target=_blank >", tmp.request$Symb, "</a>", sep="")})
				linkList2 <- c()
				for(L in 1:length(linkListe)){
					tmplist <- sapply(as.data.frame(linkListe[[L]]), function(x){paste("<a href=", x, " target=_blank >", GeneId[[L]], "</a>", sep="")})
					linkList2 <- cbind(linkList2, tmplist)
				}
			
				linkList2 <- as.data.frame(linkList2)
				colnames(linkList2) <- linkNames

				# Liste des noms dans la table
				#  names(tmp.request)
				# [1] "Query"			"Symb"			"Name"				"Org"				"Chr"               
				# [6] "Cytoband"		"Chr.start"		"Chr.end"			"Genom.start"		"Genom.end"         
				# [11] "RangeGB.Id"		"GeneId"		"RefSeq"			"PgKB"				"Cancer.Census.List"
				# [16] "KeggToDrugs"	"CTDbase"		"ClinicalTrials"    

				# tab_name <- c(linkNames, colnames(tmp.request[,c(2:3, 6, 12:13)]))	##LA VARIABLE linkNames était remplacée par link_name
				
				# bug de formatage si il n'y a qu'un gÃšne sur le fragment...
				if(ncol(linkList2)>1){
					linkList2 <- cbind.data.frame(tmp.request[,c(2, 3, 6, 12, 13)], linkList2)			# Information first(see above for #col concordance), then links.
				}
				else{
						linkList2 <- cbind.data.frame(tmp.request[,c(2, 3, 6, 12, 13)], t(linkList2))		# Information first(see above for #col concordance), then links.
				}
					
					html_tab(filename = FileName, data = linkList2, caption = paste("<b>Segment on Chr", chr, ": ", round(start/1e6, 3), "-", round(end/1e6, 3), "(Mb), Log2R = ", seg.value, ", nb genes = ", length(Symbols),"</b>", sep = ""))
			}
		}
	}
	close(PB)
	rm(census, full, PGKB.db, KeggToDrug.db, CT.db, CTDb, CTDbList)
	cat("Supplementary table saved: \n\t", paste(getwd(), FileName, sep = "/"), "\n")
	
	## 
	html_css(filename = FileName, css_liste = T)
	html_init(end = T, filedir = FileDir, filename = FileName)
	setwd(workingDir)
}

