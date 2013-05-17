# Create a project
metProject2 <- Project(list(name = 'metPoject2', description = 'Correlations between MET exprs. and ETs family genes'))

# Then create an entity
metProject2 <- creatEntity(metProject2)

# Get the synapse Id
parentId = propertyValue(metProject2, "id")

# Create folders attached to the main entity 
# Doesn't work !!!!

metPlots<-Folder(list(name = "Correlation plots", parentId=parentId, description = "Correlation_heatmaps"))
metTables<-Folder(list(name = "Correlation tables", parentId=parentId, description = "Correlation_tables"))

# Create new entities (use the web client)
metPlots <- createEntity(metPlots)
metTables <- createEntity(metTables)

# Get the entity Ids
plotsId = propertyValue(metPlots, 'id')
tablesId = propertyValue(metTables, 'id')

# Push all the plots in the directory as data entities
setwd("/gluster/home/jcommo/MetProject/Plots")
listFiles = list.files()
for (l in 1:length(listFiles)){
  tmpName = listFiles[l]
  myData <- Data(list(name = tmpName, parentId = plotsId))
  myData <- addFile(myData, tmpName)
  myData <- storeEntity(myData)
}

# Push all the tables in the directory as data entities
setwd("/gluster/home/jcommo/MetProject/Tables")
listFiles = list.files()
for (l in 1:length(listFiles)){
  tmpName = listFiles[l]
  myData <- Data(list(name = tmpName, parentId = tablesId))
  myData <- addFile(myData, tmpName)
  myData <- storeEntity(myData)
}


########################################################
# Create a zip file
# Define this zip file as a data entity
# Push the data entity in the parent entity
########################################################

# Create a project
metProject2 <- Project(list(name = 'metPoject2', description = 'Correlations between MET exprs. and ETs family genes on 172 data sets'))
# Then create an entity
metProject2 <- createEntity(metProject2)
# Get the synapse Id
parentId = propertyValue(metProject2, "id")

# Or
parentId = 'syn1621698'

# Push annotation tables
myData <- Data(list(name = 'ETsFamilyList', parentId = parentId))
myData <- addFile(myData, paste(getwd(), '/ETsFamily.txt', sep = ''))
myData <- storeEntity(myData)

myData <- Data(list(name = 'dataSetsList', parentId = parentId))
myData <- addFile(myData, paste(getwd(), '/SynapseQuery.txt', sep = ''))
myData <- storeEntity(myData)

myData <- Data(list(name = 'METCorrelationResults', parentId = parentId))
myData <- addFile(myData, paste(getwd(), '/METvsETSfamily_correl_table.xls', sep = ''))
myData <- storeEntity(myData)

myData <- Data(list(name = 'HGFCorrelationResults', parentId = parentId))
myData <- addFile(myData, paste(getwd(), '/HGFvsETSfamily_correl_table.xls', sep = ''))
myData <- storeEntity(myData)


# Push plots
setwd("/gluster/home/jcommo/MetProject2/Plots")
listFiles = list.files()
zip(paste(getwd(),'/Plots', sep = ""), listFiles)
myData <- Data(list(name = 'Plots.zip', parentId = parentId))
myData <- addFile(myData, 'Plots.zip')
myData <- storeEntity(myData)

# Push result-tables
setwd("/gluster/home/jcommo/MetProject2/Tables")
listFiles = list.files()
zip(paste(getwd(),'/Tables', sep = ""), listFiles)
myData <- Data(list(name = 'Tables.zip', parentId = parentId))
myData <- addFile(myData, 'Tables.zip')
myData <- storeEntity(myData)

# Push Rcode
setwd("/gluster/home/jcommo/MetProject2")
myCode <- Code(list(name = "Rcode", parentId = parentId))
myCode <- addFile(myCode, paste(getwd(), '/Rcode_METproject_Brig.R', sep = ''))
myCode <- storeEntity(myCode)

