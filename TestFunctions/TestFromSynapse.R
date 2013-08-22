# Test in batch from synapse
# source('/Users/fredcommo/Documents/Projet Safir/CGHObjectsSources/fromSynapse/fromSynapse_v2.R')

# Download R (class def & codes) from synapse
require(synapseClient)
workflow <- synGet('syn2128342')
source(workflow@filePath)

listFiles <- synapseQuery("select id, name from entity where entity.parentId == 'syn2025161'")

lapply(listFiles$entity.id[1:nrow(listFiles)], function(synId){
  synId <- listFiles$entity.id[1]
  cat('\nRunning on:', synId, '\n')
  object1 <- buildAgilent(synId)
  object2 <- adjustSignal(object1)
  object3 <- EMnormalize(object2, cut = c(-0.5, 0.5), MergePeaks = TRUE)
  object4 <- SegmentCGH(object3)
  object5 <- createProfile(object4)
  Profile <- pushToSynapse(object5)
  }
)

# Check the genomic profile before saving: adjust the parameters when necessary.
# getDensity(object5)
# getProfile(object5, ylim = range(-1.5, 1.5))


# Use this to create a unique shinyApps page
newShiny <- buildShinyData(object5, builtFrom = Profile)
onWeb(newShiny)

#Or check the shiny App localy
require(shiny)
scriptPath = "/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/shinyAppsCGH"
runApp(file.path(scriptPath, paste0('shiny_', propertyValue(shiny, 'id'))))

# Multi-samples
require(shiny)
setwd("/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/shinyAppsCGH")
runApp('./shinyCGHselectSample/')

# Multi-samples with gain/loss sliders
require(shiny)
setwd("/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/shinyAppsCGH")
runApp('./shinyCGHwithSliders/')

# Multi-samples with gain/loss sliders
require(shiny)
setwd("/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/shinyAppsCGH")
runApp('./shinyCGHwithSliders_3/')

# # Save:
# 	# - densityPlot,
# 	# - profile
# 	# - segmentation table
# 	# - full table
# 	# - object
# 
# savePath <- "/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/fromSynapse/RData/"
# save(object5, file = paste0(savePath, getInfo(object5, 'synapseId'), '.RData'))
# 
# source('/Users/fredcommo/Documents/Projet_Safir/Safir R Sources/CGHObjects/fromSynapse/fromSynapse.R')
# setwd("/Users/fredcommo/Documents/Projet_Safir/Safir R Sources/CGHObjects/fromSynapse/RData/")
# load(file='syn1725770.RData')
