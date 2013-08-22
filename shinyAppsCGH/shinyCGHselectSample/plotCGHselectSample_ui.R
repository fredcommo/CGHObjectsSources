require(synapseClient)
query <- synapseQuery("select id, name from entity where parentId == 'syn2116884'")
idList <- query$entity.id[grep('Profile', query$entity.name)]
##idList <- c("syn2121688", "syn2138168", "syn2144145")

# Define UI for application that plots random distributions 
shinyUI(
  bootstrapPage(	 	
    #	pageWithSidebar(
    
    # Display the sampleId
    #	headerPanel(getInfo(cghProfile, 'sampleId')),
#    h4('Choose a sample'),
    headerPanel(textOutput('Id')),

    h4('Choose samples'),
    selectInput(inputId = "Samp", label = "",
                choices = idList, selected = idList[1]),
    
    mainPanel(
      #			tabsetPanel(
      # tabPanel("Plot",
      h4('Predefined list'),
      selectInput(inputId = "PredefSymb", label = "",
                  choices = c("CCND1", "ALK", "MDM2", "FRS2", "MET", "RPTOR",
                              "ESR1", "PGR", "FGFR1", "FGFR2", "MYC", "FGF4",
                              "FGF9", "EGFR", "ERBB2", "TOP2A", "IGF1", "IGF1R",
                              "BRCA1", "BRCA2", "NOTCH4", "VEGFA", "PTEN", "PIK3CB",
                              "PAK1", "NONE"), selected = 'NONE'),
      h4('Other symbol'),
      textInput("OtherSymb", '', ''),
      submitButton("Update view"),
      h4('Gene table'), tableOutput('geneTable'),
      h4('Genomic Profile'), plotOutput("Profile", width = "900px", height = '500px')
      #					 ),
      #				 tabPanel("Summary", tableOutput('summaryTable'))
      #			)
    )
  )
)
