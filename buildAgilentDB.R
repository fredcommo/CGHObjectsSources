# Build Agilent_022060_4x180K_hg19_20110628

setwd('/Users/fredcommo/Documents/Projet Safir/Arrays Infos')
tabSeq <- read.csv('Agilent_022060_4x180_hg19_20130116_FC.txt', header = TRUE, sep = '\t', stringsAsFactor = FALSE)
tabSeq = tabSeq[order(tabSeq$TargetID),]
hg19 <- read.csv('human.chrom.info.hg19.FC.txt', header = TRUE, sep = '\t')

# Rebuild Chr:
		isChr = grep('^chr[^Mrandom]*$', tabSeq$TargetID)
		tabSeq <- tabSeq[isChr, ]
		# cnSet <- cnSet[order(cnSet$ProbeName), ]
		systNames = sapply(tabSeq$TargetID, function(x){unlist(strsplit(x, ':'))})
		cat('\n\tStep3...')
		systNames = gsub('X', '23', systNames); systNames = gsub('Y', '24', systNames)
		cat('\n\tStep4...')
		chrNum = systNames[seq(1, length(systNames), by = 2)]
		cat('\n\tStep5...')
		chrNum = as.numeric(gsub('chr', '', chrNum))
		cat('\n\tStep6...')
		positions = systNames[seq(2, length(systNames), by = 2)]
		cat('\n\tStep7...')
		positions = sapply(positions, function(x){unlist(strsplit(x, '-'))})
		start = as.numeric(positions[seq(1, length(systNames), by = 2)])
		end = as.numeric(positions[seq(2, length(systNames), by = 2)])
		
# Check positions
# Check overlaps
check = c()
for (i in 1:24){
	first = min(start[chrNum==i], na.rm = TRUE)
	last = max(start[chrNum==i], na.rm = TRUE)
	delta = last - first
	L = hg19$length[i]
	cat(i, 'start at', first, '\tend at', last, '\tDelta =', delta, '\tLenght in table:', L, '\n')
	check = c(check, ifelse(delta<L, 'shorter', 'longer'))
}
# Chromosome coverage is lower than chromosome length: check is expected as shorter for all.
check

myTab <- cbind.data.frame(tabSeq[,1:4], ChrNum = chrNum, ChrStart = start, ChrEnd = end, tabSeq[,-c(1:4)])
myTab = myTab[order(myTab$ChrNum, myTab$ChrStart),]

cumLen = cumsum(as.numeric(hg19$length))
cumLen = c(0, cumLen[-length(cumLen)])
chrTable = table(myTab$ChrNum, useNA = 'ifany')
genomic = c()
for(i in 1:length(chrTable)){
	n = chrTable[i]
	genomic = c(genomic, rep(cumLen[i], n))
	cat('Chr', i, n,'\n')
}

genomic = ifelse(!is.na(myTab$ChrStart), genomic + myTab$ChrStart, NA)

# Check overlaps
check = c()
for (i in 1:23){
	last = max(genomic[myTab$ChrNum==i], na.rm = TRUE)
	first = min(genomic[myTab$ChrNum==i+1], na.rm = TRUE)
	#cat(i, 'start at', first, '\tend at', last, '\tlen', hg19$length[i], '\n')
	cat('Last pos', i, 'at', last, '\tFirst pos', i+1, 'at', first, ifelse(last<first, "Ok", "Error"),'\n')
	check = c(check, last<first)
}
# If no overlap, chech returns True for all.
check

myTab = cbind.data.frame(myTab[,1:7 ], genomicPos = genomic, myTab[,-c(1:7)])
head(myTab)


# Save AgilentDB
#save(Agilent, file = 'Agilent_022060_4x180K_hg19.RData')
#mydb = load('Agilent_022060_4x180K_hg19.RData')

saveRDS(myTab, file = 'Agilent_022060_4x180_hg19_20130116_FC.RData')

# mydb = readRDS('Agilent_022060_4x180K_hg19_dataFrame.RData')



# Define Class
setClass('Agilent_022060_4x180K_hg19', representation(probeName = 'character',
																				Sequence = 'character',
																				GCpercent = 'numeric',
																				TargetID = 'character',
																				ChrNum = 'numeric',
																				ChrStart = 'numeric',
																				ChrEnd = 'numeric',
																				genomicPos = 'numeric',
																				Symbol = 'vector',
																				GeneName = 'character',
																				Accessions = 'character',
																				Description = 'character'
																				)
			)

# Constructor
Agilent_022060_4x180K_hg19 = function(probeName, Sequence, GCpercent, TargetID,
															ChrNum, ChrStart, ChrEnd, genomicPos,
															Symbol, GeneName, Accessions, Description)
		{
		new('Agilent_022060_4x180K_hg19', 	probeName = probeName,
																Sequence = Sequence,
																GCpercent = GCpercent,
																TargetID = TargetID,
																ChrNum = ChrNum,
																ChrStart = ChrStart,
																ChrEnd = ChrEnd,
																genomicPos = genomicPos,
																Symbol = Symbol,
																GeneName = GeneName,
																Accessions = Accessions,
																Description = Description)
		}


# Build R object
Agilent = Agilent_022060_4x180K_hg19(as.character(myTab$ProbeID),
															myTab$Sequence,
															myTab$GCpercent,
															myTab$TargetID,
															myTab$ChrNum,
															as.numeric(myTab$ChrStart),
															as.numeric(myTab$ChrEnd),
															as.numeric(myTab$genomicPos),
															myTab$GeneSymbol,
															myTab$GeneName,
															as.character(myTab$Accessions),
															as.character(myTab$Description))


save(Agilent, file = 'Agilent_022060_4x180K_hg19.RData')

# Show method
setMethod('show', signature = 'Agilent_022060_4x180K_hg19',
					function(object){
						slotnames = slotNames(object) 
						n = length(slotnames)
						cat('Agilent_022060_4x180K_hg19 DB annotated with genomic positions and probes GC content.\n')
						cat('Available items by calling getAnnot(object, \'item\'):\n\n')
						for(i in 1:n){
							cat('\t', slotnames[i], '\n')
							}
							cat('\n')
							}
				)

# getAnnot method
setGeneric('getAnnot', function(object,...) standardGeneric('getAnnot'))
setMethod('getAnnot', 'Agilent_022060_4x180K_hg19',
	function(object, arg){
		output = as.matrix(slot(object, arg))
		names(output) = slot(object, 'probeName')
			return(output)
		}
)
