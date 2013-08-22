require(synapseClient)
query <- synapseQuery("select id, name from entity where parentId == 'syn2166893'")
sampleIds <- sapply(query$entity.name, function(Name) unlist(strsplit(Name, '_'))[2])
synIds <- query$entity.id
idList <- paste(sampleIds, synIds, sep = '-')
##idList <- c("syn2121688", "syn2138168", "syn2144145")

# Define UI for application that plots random distributions 
shinyUI(
  bootstrapPage(	 	
  #  	pageWithSidebar(
    p('ShinyCGHApp - F. Commo'),
    br(),
    h1(textOutput('Id')),
#    h4(textOutput('dlrs')),
    h4(textOutput('status')),

    mainPanel(
      h4('Sample list'),
      selectInput(inputId = "Samp", label = "", choices = idList, selected = idList[1]),
      radioButtons(inputId = "modtype", label = "Choose a model",
                   choices = list("Non centered profile" = "ncp",
                                  "Left-Centered profile" = "lcp",
                                  "Zero-Centered profile" = "zcp",
                                  "Centromeric centered profile" = "ccp")
                   ),
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
