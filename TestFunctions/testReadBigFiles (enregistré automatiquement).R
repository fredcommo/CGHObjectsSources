#system.time(test <- read.csv('S355_BC31_ADN0118_s1h1_SNP6.CN5.CNCHP.txt'))

getLines <- function(x){
	return(unlist(strsplit(x, '\\t')))
}

getTagValue <- function(x){
	return(unlist(strsplit(x, '='))[2])
}

setwd('/Users/fredcommo/Documents/Projet Safir/RawData_RData')
testFile = 'S357_US82900153_252206024159_S01_CGH_107_Sep09_1_2.txt.bz2'
system.time(test <- read.csv(testFile, skip = 9, header = TRUE, sep = '\t', stringsAsFactors = FALSE))

systNames = as.character(test$SystematicName)
isChr = grep('^chr[^Mrandom]*$', test$SystematicName)
test <- test[isChr, ]

chrNames = sapply(test$SystematicName, function(x){unlist(strsplit(x, ':'))})
chrNames = gsub('X', '23', chrNames); chrNames = gsub('Y', '24', chrNames)
chrNum = chrNames[seq(1, length(chrNames), by = 2)]
chrNum = as.numeric(gsub('chr', '', chrNum))
positions = chrNames[seq(2, length(chrNames), by = 2)]
start = sapply(positions, function(x){unlist(strsplit(x, '-'))[1]})
plot(chrNum, start)

start = positions[seq(1, length(positions), by = 2)]


setwd('/Users/fredcommo/Documents/Projet Safir/RawData_RData')
testFile = 'S357_US82900153_252206024159_S01_CGH_107_Sep09_1_2.txt.bz2'
system.time(cnSet <- read.csv(testFile, skip = 9, header = TRUE, sep = '\t', stringsAsFactors = FALSE))

		keepCol <- which(as.character(colnames(cnSet)) %in% c(	"ProbeName", "SystematicName",
																								"gMedianSignal", "rMedianSignal",
																								"gIsSaturated", "rIsSaturated",
																								"gIsFeatNonUnifOL", "rIsFeatNonUnifOL",
																								"gIsWellAboveBG", "rIsWellAboveBG")
																								)
															
		# Select QC columns and rows containing probes with Ids of type 'A_xxx'
		cnSet <- cnSet[grep('^A', cnSet$ProbeName), keepCol]
		cnSet <- cnSet[order(cnSet$ProbeName), ]
	
		# use grep methods
		cat('\n\tFiltering probes...')
		isChr = grep('^chr[^Mrandom]*$', cnSet$SystematicName)
		cnSet <- cnSet[isChr, keepCol]
		cat('\n\tDefining chr nums...')
		chrNames = sapply(cnSet$SystematicName, function(x){unlist(strsplit(x, ':'))})
		chrNames = gsub('X', '23', chrNames); chrNames = gsub('Y', '24', chrNames)
		chrNum = chrNames[seq(1, length(chrNames), by = 2)]
		chrNum = as.numeric(gsub('chr', '', chrNum))
		cat('\n\tDefining start pos...')
		positions = chrNames[seq(2, length(chrNames), by = 2)]
		start = sapply(positions, function(x){unlist(strsplit(x, '-'))[1]})
