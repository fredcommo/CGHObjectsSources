# library(shiny)

require(synapseClient)
# Load class def & accessors
cat('Loading functions...')
# synGet() doesn't work with codeEntity
eHelper <- synGet('syn1889772') # codeEntity: syn1894643, fileEntity: syn1889772
source(eHelper@filePath)

# Load annot tables
eAnnot <- synGet('syn1877556')
HG19 <- read.csv(eAnnot@filePath, header = TRUE, sep = '\t')
eDB <- synGet('syn1877601')
geneDB <- readRDS(eDB@filePath)
rm(eHelper, eAnnot, eDB)

# Load data
cat('\n\nLoading data...')
cgh  <- synGet('syn1871353') # RDS= syn1871353, RData = syn1864342
load(cgh@filePath)
cat('\n\n')

#
shinyServer(function(input, output) {

  # Function that generates a plot of the genomic profile and provide the information about specified genes.
  # It's updated each time a new gene symbol is called.
 	# createPlot is the reactive function called by renderPlot() to generate the main plot.
  # createTable is the reactive function called by renderTable()
     
  createPlot <- reactive({
    gene <- toupper(input$PredefSymb)
  	if(input$OtherSymb != '') gene <- toupper(input$OtherSymb)
  	if(gene != 'NONE' & gene != '')
      gPlot <- .tagMyGene(Example, geneDB, tag = gene)
    else
      gPlot <- .plotProfile(Example)
    print(gPlot)
  	})
  output$Profile <- renderPlot({createPlot()})
  
  createTable <- reactive({
		gene <- toupper(input$PredefSymb)
  	if(input$OtherSymb != '') gene = toupper(input$OtherSymb)
		if(gene == 'NONE' | gene == '')
      data.frame(Symbol = 'Choose a gene in the "Predefined list" or enter a valid symbol in "Other symbol",\n
                 then click on "Update view"')
		else{
      if (any(gene %in% geneDB$Symbol)){
        tmp <- .getValue(Example, geneDB, gene)
			  data.frame(Symbol = tmp$Symbol, entrezgene = tmp$entrezgene,
                   Description = tmp$Description, Chr = tmp$Chr, Cytoband = tmp$Cytoband,
                   Log2Ratio = round(tmp$LR, 3))
        }
      else data.frame(Symbol = paste0('Too bad!!, "', gene, '" is not a valid symbol. Try again!'))
  	  }
		})
  output$geneTable <- renderTable({createTable()})
	  }
  )


#setwd('/Users/fredcommo/Documents/Stats/Docs R/shinyApps/')
#runApp('./shinyAppCGH/')