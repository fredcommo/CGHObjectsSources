myBiomaRt <- function(){
	require(biomaRt)
	require(foreach)
	require(iterators)
	source('/Users/fredcommo/Documents/Projet Safir/Safir R Sources/CGHObjects/AllHelperFunctions.R')
	human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	cat('human built\n')
	output = getBM(attributes=c('hgnc_symbol', 'description', 'chromosome_name', 'band', 'start_position','end_position', 'entrezgene'),
			, mart = human)
	output$description =  gsub(' \\[.*\\]', '', output$description)
	output = as.data.frame(output)
	# Add genomic positions
	arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
	hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
	cumLen = cumsum(as.numeric(hg19$length))
	cumLen = c(0, cumLen[-length(cumLen)])
	output$chromosome_name <- gsub('X', 23, output$chromosome_name)										
	output$chromosome_name <- gsub('Y', 24, output$chromosome_name)										
	output$chromosome_name <- as.numeric(as.character(output$chromosome_name))
	output <- output[order(output$chromosome_name, output$start_position),]
	chrTab <- table(output$chromosome_name, useNA = 'ifany')
	gStart = gEnd = c()
	foreach(chr = 1:24) %do%{
		cat('Compute Chr:', chr, '\n')
		index = which(output$chromosome_name == chr)
		gStart = c(gStart, output$start_position[index] + rep(cumLen[chr], length(index)))
		gEnd = c(gEnd, output$end_position[index] + rep(cumLen[chr], length(index)))
		}
	index = which(is.na(output$chromosome_name))
	gStart = c(gStart, rep(NA, length(index)))
	gEnd = c(gEnd, rep(NA, length(index)))
	output = cbind.data.frame(output[,-ncol(output)], genomicStart = gStart, genomicEnd = gEnd, entrezgene = output[,ncol(output)])

	cat('Run GeneRequest\n')
	myGeneDB <- c()
	foreach(i = 1:nrow(output)) %do%{
		gene = output$hgnc_symbol[i]
		if(gene == '') gene = output$entrezgene[i]
		if(!is.na(gene)){
			tmp <- geneRequest.v7(gene, 'gene', verbose = TRUE)
			myGeneDB <- rbind(myGeneDB, tmp)
			}
		}		
	saveRDS(output, file = paste0(arrayInfoPath,'myBiomaRt.rds'))
	cat('myBiomaRt saved\n')
	
	saveRDS(myGeneDB, file = paste0(arrayInfoPath,'myGeneDB.rds'))
	cat('myGeneDB saved\n')
}
myBiomaRt()

