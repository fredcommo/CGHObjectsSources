# Define UI for application that plots random distributions 
 shinyUI(
	bootstrapPage(	 	
#	pageWithSidebar(
	
		# Display the sampleId
  	#	headerPanel(getInfo(cghProfile, 'sampleId')),
	    headerPanel('Example'),
	  
		# Predefined lsit of genes
#      sidebarPanel(
#  		textInput("userId", 'Enter password', 'foo'),
      textInput("OtherSymb", '', ''), # Search for:
  		selectInput(inputId = "PredefSymb",
  		            label = "display genes",
  		            choices = c("CCND1", "ALK", "MDM2", "FRS2", "MET", "RPTOR", "ESR1", "PGR", "FGFR1", "FGFR2",
  		                        "MYC", "FGF4", "FGF9", "EGFR", "ERBB2", "TOP2A", "IGF1", "IGF1R", "BRCA1", "BRCA2",
  		                        "NOTCH4", "VEGFA", "PTEN", "PIK3CB", "PAK1", "NONE"), selected = 'NONE')
#     )
      ,
  		
  		#submitButton("Update View"),
  # Show a plot of the generated profile

		mainPanel(
#			tabsetPanel(
				# tabPanel("Plot",
 		  conditionalPanel(condition = "input.userId == ''",
 		                   tableOutput('Access')),
 		  h4('Gene table'),
 		  conditionalPanel(condition = "input.PredefSymb != 'NONE' || input.OtherSymb != ''",
 		                   tableOutput('geneTable')),
 		  h4('Genomic Profile'),
			plotOutput("Profile", width = "600px", height = '300px')
					# ),
				# tabPanel("Summary", tableOutput('summaryTable'))
#			)
		)
	)
)

