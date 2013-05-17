# Define UI for application that plots random distributions 
 shinyUI(
	#bootstrapPage(	 	
	pageWithSidebar(
	
		# Display the sampleId	
  		headerPanel(getInfo(object5, 'sampleId')),

		# Predefined lsit of genes	 		
  		selectInput(inputId = "PredefSymb",
  						label = "display genes",
  						choices = c("CCND1", "ALK", "MDM2", "FRS2", "MET", "RPTOR", "ESR1", "PGR", "FGFR1", "FGFR2",
										"MYC", "FGF4", "FGF9", "EGFR", "ERBB2", "TOP2A", "IGF1", "IGF1R", "BRCA1", "BRCA2",
										"NOTCH4", "VEGFA", "PTEN", "PIK3CB", "PAK1", "NONE"), selected = 'NONE'),
  
  		textInput("OtherSymb", 'Search for:', ''),
  		#submitButton("Update View"),
  # Show a plot of the generated profile

		mainPanel(
#			tabsetPanel(
				# tabPanel("Plot",
					h4('Genomic Profile'),
					plotOutput("Profile", width = "1600px", height = '600px'),
					h4('Gene table'),
					conditionalPanel(condition = "input.PredefSymb != 'NONE' || input.OtherSymb != ''",
											tableOutput('geneTable'))
					# ),
				# tabPanel("Summary", tableOutput('summaryTable'))
#			)
		)
	)
)

