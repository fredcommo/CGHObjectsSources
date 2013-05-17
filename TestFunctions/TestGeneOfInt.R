# Test geneRequest using myGeneDB2

arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'		
geneDB <- readRDS(paste0(arrayInfoPath, 'myGeneDB_2013_Mar_26.rds'))
geneOfInt <- function(object, geneList, DB = geneDB){
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

geneList <- c('FGFR1', 'EGFR', 'ERBB2')
geneOfInt(object5, geneList)

		# foreach(gene = iter(geneList)) %do% {
			# tmp = geneDB[which(geneDB$Symbol == gene),]
			# start = tmp$genomicStart
			# stop = tmp$genomicStop
			# containStart = which(segTable$loc.start<=start & segTable$loc.end>=start)
			# containStop = which(segTable$loc.start<=stop &  segTable$loc.end>=stop)
			# if(is.na(start) | is.na(stop))
				# tmp = cbind.data.frame(tmp, Log2Ratio = NA)
			# else if(containStart == containStop)
				# tmp = cbind.data.frame(tmp, Log2Ratio = segTable$seg.med[containStart])
			# else{
				# Log2Ratio =  segTable$seg.med[union(containStart, containStop)]
				# tmp = cbind.data.frame(tmp, Log2Ratio)
				# }
			# output = rbind(output, tmp)
			# }
