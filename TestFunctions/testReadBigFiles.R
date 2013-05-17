#system.time(test <- read.csv('S355_BC31_ADN0118_s1h1_SNP6.CN5.CNCHP.txt'))

getLines <- function(x){
	return(unlist(strsplit(x, '\\t')))
}

getTagValue <- function(x){
	return(unlist(strsplit(x, '='))[2])
}

setwd('/Users/fredcommo/Documents/Projet Safir/RawData_RData')
testFile = 'S357_US82900153_252206024159_S01_CGH_107_Sep09_1_2.txt.bz2'
system.time(test <- read.csv(testFile, skip = 9, header = TRUE, sep = '\t'))

systName = as.character(test$SystematicName)
isChr = grep('^chr', systName)
chrNames = sapply(systName[isChr], function(x){unlist(strsplit(x, ':'))[1]})
positions = sapply(systName[isChr], function(x){unlist(strsplit(x, ':'))[2]})
chrNames = sapply(systName[isChr], function(x){strsplit(x, ':')})
chrNames = gsub('X', '23', chrNames)
testAnnot = gsub('Y', '24', chrNames)
testAnnot[grep('random', testAnnot)] = 'chrNA'
testAnnot[grep('M', testAnnot)] = 'chrNA'
ChrNum = as.numeric(substr(testAnnot, 4, 5))
table(chrNum)

system.time(cnSet <- read.csv(testFile, skip = max(isInfo), header = TRUE, sep = '\t'))
