# push CGH workflow on synapse

parentId <- 'syn2117484'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/')
files <- list.files()
Rfiles <- files[grepl('All', files)]
lapply(Rfiles, function(f){
  file <- File(file.path('.', f), parentId = parentId)
  file <- synStore(file)
})

parentId <- 'syn2117484'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources')
file <- File('./AllAccessors.R', parentId = parentId)
file <- synStore(file)

parentId <- 'syn2117484'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources')
file <- File('./AllClasses.R', parentId = parentId)
file <- synStore(file)

parentId <- 'syn2117484'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources')
file <- File('./AllGenerics.R', parentId = parentId)
file <- synStore(file)

parentId <- 'syn2117484'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources')
file <- File('./AllMethods.R', parentId = parentId)
file <- synStore(file)

parentId <- 'syn2117484'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources')
file <- File('./AllHelperFunctions.R', parentId = parentId)
file <- synStore(file)

parentId <- 'syn2117484'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources')
file <- File('./AllRpackages.R', parentId = parentId)
file <- synStore(file)

parentId <- 'syn2117484'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources')
file <- File('./AllShowMethods.R', parentId = parentId)
file <- synStore(file)

parentId <- 'syn2117484'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources')
file <- File('./loadWorkflow.R', parentId = parentId)
file <- synStore(file)

parentId <- 'syn2117484'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/fromSynapse/')
file <- File('./fromSynapse_v2.R', parentId = parentId)
file <- synStore(file)

# Shiny templates for individual visualization
parentId <- 'syn2117484'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/shinyAppsCGH/shinyCGHtemplates/')
file <- File('./plotCGHtemplate_server.R', parentId = parentId)
file <- synStore(file)

parentId <- 'syn2117484'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/shinyAppsCGH/shinyCGHtemplates/')
file <- File('./plotCGHtemplate_ui.R', parentId = parentId)
file <- synStore(file)

# shiny templates for multi-visualization on the main results page: includes samples selection
parentId <- 'syn2145860'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/shinyAppsCGH/shinyCGHselectSample/')
file <- File('./plotCGHselectSample_server.R', parentId = parentId)
file <- synStore(file)

parentId <- 'syn2145860'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/shinyAppsCGH/shinyCGHselectSample/')
file <- File('./plotCGHselectSample_ui.R', parentId = parentId)
file <- synStore(file)

# shiny templates for multi-visualization & sliders on the main results page: includes samples selection
parentId <- 'syn2145860'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/shinyAppsCGH/shinyCGHwithSliders/')
file <- File('./plotCGHwithSliders_server.R', parentId = parentId)
file <- synStore(file)

parentId <- 'syn2145860'
setwd('/Users/fredcommo/Documents/Projet_Safir/CGHObjectsSources/shinyAppsCGH/shinyCGHwithSliders/')
file <- File('./plotCGHwithSliders_ui.R', parentId = parentId)
file <- synStore(file)
