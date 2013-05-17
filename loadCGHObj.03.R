

# The function activates an interactive selection of file to load
# Collect microarray informations and Cy5/Cy3 intensity values,
# and suppress flags and duplicated probes.
# Root : A letter indicating the hard drive partition to use. May run with other indications (~/home)

##########################################################################################
###########################

# Class are defined in AllClasses.R
# Accessors 	are defined in AllAccessors.R
# Generic methids are defined in AllGenerics.R
# Methods are in AllMethods .R

source('SourceCodeObj.R')
source("AllClasses.R")
cat('All classes... Done.\n')
source('AllAccessors.R')
cat('All accessors... Done.\n')
source('AllShowMethods.R')
cat('All show methods... Done.\n')
source('AllGenerics.R')
cat('All generics... Done.\n')
source('AllMethods.R')
cat ('All methods... Done.\n')
source('AllHelperFunctions.R')
cat ('Helper functions... Done.\n')

# validity method : To do !

# Common functions : To do !
	# setGeneric('Norm', function(object) standardGeneric('Norm'))
	# setGeneric('Segment', function(object) standardGeneric('Segment'))
	# setGeneric('Plot', function(object) standardGeneric('Plot'))
	# setGeneric('Save', function(object) standardGeneric('Save'))
	# setGeneric('EditSupTab', function(object) standardGeneric('EditSubTab'))



###########################
###########################

# loadCGHObj.03 <- function(){
	# f = try(getFile(), silent = TRUE)
	# if(class(f) == 'try-error') stop ('No selected file!\n')
	
	# cat("\nGetting Array Information...\n")
	# object = getAnnot(f)
	# object <- readInfo(object)
	
	# cat('\nReading', f, '...\n')
	# object <- readCN(object)
	
	# cat('\nPreprocessing...')	
	# object <- supressFlags(object)
	# object <- supressDuplic(object)	
	# object <- preset(object)
	# cat('\tDone.')
	
	# cat("\n\n\t*** Everything is fine. Lucky you :-) ***\n\n")
	# return(object)
# }

