# Test in batch from synapse
require(synapseClient)
synapseLogin('frederic.commo@sagebase.org', 'Se@ttle7')

source('/Users/fredcommo/Documents/Projet Safir/Safir R Sources/CGHObjects/fromSynapse/fromSynapse.R')

listFiles <- synapseQuery('select id, name from entity where entity.parentId == "syn1715847"')
listFiles
synId <- listFiles$entity.id[1]
object1 <- buildAgilent(synId)
object2 <- adjustSignal(object1, Fract = 0.15)
object3 <- EMnormalize(object2, cut = c(-0.4, 0.4), MergePeaks = TRUE)
object4 <- SegmentCGH(object3, UndoSD = 1)
object5 <- createProfile(object4)

getDensity(object5)
getProfile(object5, ylim = range(-1.5, 1.5))

# Save:
	# - densityPlot,
	# - profile
	# - segmentation table
	# - full table
	# - object

savePath <- "/Users/fredcommo/Documents/Projet Safir/Safir R Sources/CGHObjects/fromSynapse/RData/"
save(object5, file = paste0(savePath, getInfo(object5, 'sampleId'), '.RData'))

source('/Users/fredcommo/Documents/Projet Safir/Safir R Sources/CGHObjects/fromSynapse/fromSynapse.R')
setwd("/Users/fredcommo/Documents/Projet Safir/Safir R Sources/CGHObjects/fromSynapse/RData/")
load(file='syn1725770.RData')