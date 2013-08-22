# Define UI for application that plots random distributions 
shinyUI(
  bootstrapPage(	 	
    #	pageWithSidebar(
    
    # Display the sampleId
    #	headerPanel(getInfo(cghProfile, 'sampleId')),
    headerPanel(textOutput('Id')),
    
    mainPanel(
      #			tabsetPanel(
      # tabPanel("Plot",
      conditionalPanel(condition = "input.userId == ''",
                       tableOutput('Access')),
      h4('Predefined list'),
      selectInput(inputId = "PredefSymb",
                  label = "",
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
