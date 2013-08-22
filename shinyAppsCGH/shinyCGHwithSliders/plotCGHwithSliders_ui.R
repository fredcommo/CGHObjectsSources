require(synapseClient)
query <- synapseQuery("select id, name from entity where parentId == 'syn2116884'")
sampleIds <- gsub('_(.*)', '', query$entity.name[grep('Profile', query$entity.name)])
synIds <- query$entity.id[grep('Profile', query$entity.name)]
idList <- paste(sampleIds, synIds, sep = '-')
##idList <- c("syn2121688", "syn2138168", "syn2144145")

# Define UI for application that plots random distributions 
shinyUI(
  bootstrapPage(	 	
  #  	pageWithSidebar(
    
    headerPanel(textOutput('Id')),

    mainPanel(
      h4('Sample list'),
      selectInput(inputId = "Samp", label = "", choices = idList, selected = idList[1]),
      h4('Predefined gene list'),
      selectInput(inputId = "PredefSymb", label = "",
                  choices = c("CCND1", "ALK", "MDM2", "FRS2", "MET", "RPTOR",
                              "ESR1", "PGR", "FGFR1", "FGFR2", "MYC", "FGF4",
                              "FGF9", "EGFR", "ERBB2", "TOP2A", "IGF1", "IGF1R",
                              "BRCA1", "BRCA2", "NOTCH4", "VEGFA", "PTEN", "PIK3CB",
                              "PAK1", "NONE"), selected = 'NONE'),
      h4('Other symbol'), textInput("OtherSymb", '', ''),
      sliderInput("gain", "Gain threshold:", min=0, max=1, value=log2(1.25), step = .01),
      sliderInput("loss", "Loss threshold:", min=-1, max=0, value=log2(.85), step = .01),
      h4('Gene table'), tableOutput('geneTable'),
      plotOutput("Profile", width = "900px", height = '500px')
    )
  )
)
