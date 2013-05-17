
# Test Affy workflow
scriptPath = "/Users/fredcommo/Documents/Projet Safir/Safir R Sources/CGHObjects/"
setwd(scriptPath)
source('SourceCodeObj.R')

# Analysis workflow
object1 <- buildCGHObj.v01()
object2 <- adjustSignal(object1)
object3 <- EMnormalize(object2)
object4 <- SegmentCGH(object3, UndoSD = 2)
object5 <- createProfile(object4)
getProfile(object5, ylim = range(-1.5, 1.5))					# the genomic profile


# Save plots in a directory
setwd('/Users/fredcommo/Documents/Projet Safir')
path = paste0('./Output/', getInfo(object5, 'sampleId'))
dir.create(path)

pdf(file = paste0(path, '/', getInfo(object4, 'sampleId'), '_Density.pdf'), width = 5, height = 5)
getDensity(object5)
dev.off()

pdf(file = paste0(path, '/', getInfo(object4, 'sampleId'), '_plot.pdf'), width = 10, height = 7.5)
getProfile(object5, ylim = range(-1.5, 1.5))
dev.off()

# Save .Rdata
cghProfile <- object5
save(cghProfile, file = paste0(path, '/', getInfo(object5, 'sampleId'), '.RData'))

# Save in synapse: it works using .rds or .rbin
fileName <- getInfo(object5, 'sampleId')
assign(fileName, list(gProfile = getProfile(object5), segTable = getSegTable(object5)))
save(list = fileName, file = file.path(tempdir(), paste0(fileName, '.rds')))

file <- File(file.path(tempdir(), paste0(fileName, '.rds')), parentId = 'syn1864121')
file <- synStore(file)

# Generate an save the results table for a given list of genes on interest
	# create a list of symbols
geneList = c("CCND1", "ALK", "MDM2", "FRS2", "MET", "RPTOR", "ESR1", "PGR", "FGFR1", "FGFR2",
					"MYC", "FGF4", "FGF9", "EGFR", "ERBB2", "TOP2A", "IGF1", "IGF1R", "BRCA1", "BRCA2",
					"NOTCH4", "VEGFA", "PTEN", "PIK3CB", "PAK1")
	# retrieve the annotations
geneTable <- geneOfInt(object5, geneList)

	# create and save the HTML tables
setwd('/Users/fredcommo/Documents/Projet Safir')
path = paste0('./Output/', getInfo(object5, 'sampleId'))
buildHtml(geneTable, path, paste0(getInfo(object5, 'sampleId')))
FullHtml(object5, filePath = path)

require(shiny)
shinyPath = '/Users/fredcommo/Documents/Projet Safir/Safir R Sources/CGHObjects/'
runApp(paste0(shinyPath, 'shinyAppCGH/'))



# Accessors
getInfo(object5)															# retreive the information
getInfo(object5, 'sampleId')											# retreive a specific information
getParam(object5)														# retreive the analysis parameters
head(getSegTable(object5))											# the segmentation table	
getDensity(object5)													# Display the density plot
getProfile(object5, ylim = range(-1.5, 1.5))					# the genomic profile
tagMyGene(object5, c('PIK3CB', 'MET', 'FGF4', 'EGFR','ERCC1'))							# Add a specified gene
zoomIn(object5, 1)														# Zoom in a specified chromosome



############################
projectId <- 'syn1864121'

# note, we recommend using "File" rather than "Data"
myPath <- '/Users/fredcommo/Documents/Projet Safir/RawData_RData/'
fileId <- 'Example_Bx_xxx_s1h1_xxx_GenomeWideSNP6_CN5_CNCHP.txt'
file <- Data(list(name = fileId, parentId=projectId))
file <- addFile(file, fileId)
file <- synStore(file)
sourceData <- getEntity(propertyValue(file, 'id'))
sourceData <- getEntity('syn1864341')

myPath <- '/Users/fredcommo/Documents/Projet Safir/Output/Example/'
fileId <- 'Example.RData'
setwd(myPath)
file <- Data(list(name = fileId, parentId=projectId))
file <- addFile(file, fileId)
file <- synStore(file)
rData <- getEntity(propertyValue(file, 'id'))

# Create provenance
used(rData)<-list(list(entity = sourceData, wasExecuted = FALSE))
rData <- synStore(rData)

# Push profile to wiki
folderWikiUri <- sprintf("/entity/%s/wiki", projectId)
    # Start a wiki, and initialize the markdown as empty.
folderWiki <- list()
folderWiki$attachmentFileHandleIds <- list()
folderWiki$markdown <- c()
folderWiki <- synRestPOST(folderWikiUri, folderWiki)










# Check overlaps in genomic positions
cnSet = getCNset(object)
for(chr in 1:23){
	g1 = cnSet$genomiPos[cnSet$ChrNum == chr]
	g2 = cnSet$genomiPos[cnSet$ChrNum == chr+1]
	cat('\n', chr, '/', chr+1, '\toverlap:', any(g1>g2, na.rm = TRUE))
}

cnSet = getCNset(object)
snp = getSNPset(object)

left = 1e6
right = 1e7
snpIndex = which(snp$genomicPos > left & snp$genomicPos < right)
cnIndex = which(cnSet$genomicPos > left & cnSet$genomicPos < right)
length(snpIndex); length(cnIndex)
plot(cnSet$genomicPos[cnIndex], cnSet$Log2Ratio[cnIndex]-2, pch = 19, cex = 0.2, col = 'grey75', ylim = range(-5, 5))
points(snp$genomicPos[snpIndex], snp$Log2Ratio[snpIndex]+2, pch = 19, cex = 0.2, col = 'red3')


cnSet = getCNset(object4)
l = cnSet$Log2Ratio
g = cnSet$genomicPos
s = cnSet$Segm
Title = getInfo(object4, 'sampleId')
subTitle = paste(getInfo(object4, 'platform'), getInfo(object4, 'gridGenomicBuild'), getInfo(object4, 'analyseDate'), sep = ' - ')
xRange = range(g, na.rm = TRUE); yRange = range(-1.5, 1.5)

currentPlot = xyplot(yRange~xRange, type = 'n', xlab = 'Genomic Position', ylab = 'Log2Ratio', main = Title, sub = subTitle,
panel = function(){panel.xyplot(xRange, yRange, type = 'n')
			locateChr(max(yRange))
			locateProbes(g, l)
			locateSegments(g, s, 0.1)
			}
			)
			
currentPlot



locateGenes(Table)

segTab = getSegTable(object)
ratio = segTab$seg.med/abs(segTab$loc.end - segTab$loc.start)
score = log10(abs(ratio))
par(mfrow = c(2, 1))
xyplot(yRange ~ xRange, type = 'n', xlab = 'Genomic Position', ylab = 'Log2Ratio', main = Title, sub = subTitle)
locateChr(max(yRange))
locateProbes(g, l)
locateSegments(g, s, 0.1)
#locateGenes(Table)

plot(xRange, range(ratio, na.rm = TRUE), type = 'n', xlab = 'Genomic Position', ylab = 'Log2R/Length', main = Title, sub = subTitle)
segments(x0 = segTab$loc.start, y0 = ratio, x1 = segTab$loc.end, lwd = 3)

plot(xRange, range(score, na.rm = TRUE), type = 'n', xlab = 'Genomic Position', ylab = 'Log2Ratio', main = Title, sub = subTitle)
segments(x0 = segTab$loc.start, y0 = score, x1 = segTab$loc.end, lwd = 3)

mu = mean(score, na.rm = T)
pl = pnorm(score[score<mu], mean = mu, sd = sd(score,na.rm=T), lower.tail = T)
pR = pnorm(score[score>=mu], mean = mu, sd = sd(score,na.rm=T), lower.tail = F)
P = c(pl, pR)
plot(score, P)
plot(xRange, range(score, na.rm = TRUE), type = 'n', xlab = 'Genomic Position', ylab = 'Log2Ratio', main = Title, sub = subTitle)
segments(x0 = segTab$loc.start, y0 = score, x1 = segTab$loc.end, lwd = 3, col = ifelse(P<0.05, 'red','black'))

# Testing adjustSignal
init = log2(cnv$rMedianSignal/cnv$gMedianSignal)
object <- adjustSignal(object)
post <- getCNset(object)$Log2Ratio

# Step by step
f = try(getFile(), silent = TRUE); f
object = getAnnot(f)
object

object <- readInfo(object)
object

object <- readCN(object)
object

object <- supressFlags(object)
object

object <- supressDuplic(object)	
object

object <- preset(object)
object


cnSet1 <- getCNset(object)
cnSet2 <- CyAdjust(cnSet1, 1, TRUE) 
cnSet3 <- GCadjust(cnSet2)
LR = log2(cnSet$rMedianSignal/cnSet$gMedianSignal)

Samp = sample(1:nrow(cnSet), 5000)
par(mfrow = c(2, 2))
plot(LR[Samp], cnSet2$Log2Ratio[Samp], cex = 0.1)
plot(cnSet2$Log2Ratio[Samp], cnSet3$Log2Ratio[Samp], cex = 0.1)
plot(LR[Samp], cnSet3$Log2Ratio[Samp], cex = 0.1)
par(mfrow = c(1, 1))
summary(lm(cnSet2$Log2Ratio[Samp]~LR[Samp]))
summary(lm(cnSet3$Log2Ratio[Samp]~cnSet2$Log2Ratio[Samp]))
summary(lm(cnSet3$Log2Ratio[Samp]~LR[Samp]))


