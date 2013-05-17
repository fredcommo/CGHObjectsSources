

object <- object5
gain = log2(2.25/2)
loss = log2(0.8/2)
whichseg = which(Segments >= gain | Segments <= loss)
selecType = paste(gainSubTitl, "or", lossSubTitl)
CHR = 1:23
sangercensus = TRUE; pharmgkb = TRUE; keggtodrug = TRUE; ctdbase = TRUE; clintrials = TRUE
filePath
fileName
cssFile

		# Load information tables
		arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
		geneDB <- readRDS(paste0(arrayInfoPath, 'myGeneDB.rds'))
		if(sangercensus) census <- read.table(paste0(arrayInfoPath, "CensusTable_Annot_FC.txt"), header = T, sep = "\t")
		if(pharmgkb) PGKB.db <- read.csv(paste0(arrayInfoPath, "relationships.tsv"), header = T, sep = "\t")
		if(keggtodrug) KeggToDrug.db <- read.csv(paste0(arrayInfoPath, "Kegg_Drugs_Table.txt"), header = T, sep = "\t")
		if(clintrials) clinTrials.db <- read.csv(paste0(arrayInfoPath, "ClinicalTrials_Drugs_Table.txt"), header = T, sep = "\t")
		if(ctdbase){
			CTDb <- read.csv(paste0(arrayInfoPath, "CTD_chem_gene_ixns_SansDuplic.txt"), header = T, sep = "\t")
			CTDbList <- CTDb$GeneSymbol
			}	

		# Define what to consider
		segTable <- getSegTable(object)
		selectTable <- segTable[which(segTable$chrom %in% CHR),]
		Segments <- selectTable$seg.med
		Select <- match.arg(Select)
		switch(Select,
					Gain = (whichseg = which(Segments >= gain)),
					Loss = (whichseg = which(Segments <= loss)),
					Both = (whichseg = which(Segments >= gain | Segments <= loss))
					)
		gainSubTitl <- paste("gained (>", round(gain, 3), ")", sep = "")
		lossSubTitl <- paste("lost (<", round(loss, 3), ")", sep = "")
		switch(Select,
					Gain = (selecType = gainSubTitl),
					Loss = (selecType = lossSubTitl),
					Both = (selecType = paste(gainSubTitl, "or", lossSubTitl))
					)
		Explore <- paste("Chr", CHR[1], "to", CHR[length(CHR)], sep = "")

		FileName <- paste(getInfo(object, 'sampleId'), getInfo(object, 'barCode'), getInfo(object, 'analyseDate'), getInfo(object, 'platform'), Explore, Select, "SupplTable", sep = "_")
		TitleName <- paste("Supplementary tables: ", getInfo(object, 'sampleId'), " / ", getInfo(object, 'barCode'), "<br>",
						getInfo(object, 'platform'), " / ", getInfo(object, 'analyseDate'), " / ", paste("Chr", CHR[1], "to", CHR[length(CHR)]), ": ", paste("segments", selecType))
		fullTable <- HTMLInitFile(filePath, filename = fileName, Title = TitleName)

		if(length(Segments)>1) Append = TRUE
		foreach(Segment = iter(Segments)) %do% {
			
			# Get the list of genes included in Segment
			tmpGeneDB <- geneDB[which(geneDB$genomicStart >= selectTable$loc.start[Segment] & geneDB$genomicStop <= selectTable$loc.end[Segment]),]
			Symbols = as.character(tmpGeneDB$NomenclatureSymbol)
			entregene = as.character(tmpGeneDB$entrezgene)
			
			# Find those which are in the Drug DBs
			censusIds <- filtrSangerCensus(Symbols, census)
			pgkbIds <- filtrPharmGkb(Symbols, PGKB.db)
			KeggIds <- filtrKeggToDrug(Symbols, KeggToDrug.db)
			ctdIds <- filtrCTDbase(Symbols, CTDbList)
			clinTrialsIds <- filtrClinTrials(Symbols, clinTrials.db)

			# Build a table
			
			# Add first html entregene & html KeggPathway.
			
			if(sangercensus)
				tmpGeneDB <- cbind.data.frame(tmpGeneDB, SangerCancer = ifelse(!is.na(censusIds), htmlLink(censusIds, sangerLink(censusIds)), "-"))
			if(pharmgkb)
				tmpGeneDB <- cbind.data.frame(tmpGeneDB, PgKB = htmlLink(pgkbIds, pharmGkbLink(pgkbIds)))
			if(keggtodrug)
				tmpGeneDB <- cbind.data.frame(tmpGeneDB, KeggToDrugs = ifelse(!is.na(KeggIds), htmlLink(KeggIds, keggToDrugLink(KeggIds)), "-"))
			if(ctdbase)
				tmpGeneDB <- cbind.data.frame(tmpGeneDB, CTDbase = ifelse(!is.na(ctdIds), htmlLink(Symbols, CTDbaseLink(Symbols)), "-"))
			if(clintrials)
				tmpGeneDB <- cbind.data.frame(tmpGeneDB, ClinicalTrials = ifelse(!is.na(clinTrialsIds), htmlLink(clinTrialsIds, clinTrialsLink(clinTrialsIds)), "-"))
			# Restrict the table if requested
			if(Restrict){
				ofInterest <- which(!is.na(censusIds) | pgkbIds != "-" | !is.na(KeggIds) | !is.na(ctdIds) | !is.na(clinTrialsIds))
				tmpGeneDB <- tmpGeneDB[ofInterest,]
				}
			
			# Add the current tmpTable+\n to the main output table
			}
		
		# Save the final table using html function and css
}

filtrPharmGkb <- function(Symbols, PGKB.db){
	# Search for PGKB annotations
	pgkbIds <- rep("-", length(Symbols))
	pgkbGenes <- ifelse(as.character(Symbols) %in% PGKB.db$Entity1_name, as.character(Symbols), "-")
	#cat("PgKB\n", ifelse(pgkb.genes!="-", pgkb.genes, "---"), "\n")
	pgkbIds <- rep("-", length(pgkbGenes))
	if(any(pgkbGenes != "-")){
		for(pg in 1:length(pgkbGenes))
			if(pgkbGenes[pg] != "-"){
				tmp <- as.character(unique(PGKB.db$Entity1_id[which(PGKB.db$Entity1_name == pgkbGenes[pg])]))
				pgkbIds[pg] <- substr(tmp, 6, 50)
				}
		}
	return(as.character(pgkbIds))	
}

filtrSangerCensus <- function(Symbols, census){
	censusIds <- ifelse(as.character(Symbols) %in% census$Symb, Symbols, NA)
	return(as.character(censusIds))
}

filtrKeggToDrug <- function(Symbols, KeggToDrug.db){
	# search for KeggToDrugs
	KeggIds <- rep(NA, length(Symbols))
	is.Kegg <- ifelse(as.character(Symbols) %in% KeggToDrug.db$Symb, Symbols, NA)
	KeggIds <- apply(as.data.frame(is.Kegg), 1,
						function(x){
							if(as.character(x) %in% Symbols){
								kegg.index <- which(KeggToDrug.db$Symb == x)
								if(KeggToDrug.db$nDrug[kegg.index]!=0) return (KeggToDrug.db$GeneId[kegg.index])
								else return(NA)
								}
							else return (NA)
						})
	return(as.character(KeggIds))
}

filtrCTDbase <- function(Symbols, CTDbList){
	# Search for CTDbase
	ctdIds <- rep(NA, length(Symbols))
	ctdIds <- ifelse(as.character(Symbols) %in% CTDbList, Symbols, NA)
	return(as.character(ctdIds))
}

filtrClinTrials <- function(Symbols, clinTrials.db){
	# search for ClinicalTrials
	is.CT <- rep(NA, length(Symbols))
	is.CT <- ifelse(as.character(Symbols) %in% clinTrials.db$Symb, Symbols, NA)
	clinTrialsIds <- apply(as.data.frame(is.CT), 1,
							function(x){
								if(as.character(x) %in% Symbols){
									CT.index <- which(clinTrials.db$Symb == x)
									if(clinTrials.db$CTfound[CT.index]!=0) return (x)
									else return(NA)
									}
								else return (NA)
								})
	return(as.character(clinTrialsIds))
}
