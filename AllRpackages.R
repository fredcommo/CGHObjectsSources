	# Load R packages
	cat('Loading R packages...')
	source("http://bioconductor.org/biocLite.R")
	biocList <- c('limma', 'DNAcopy', 'preprocessCore', 'affy', 'annotate')
  cranList <- c('robustbase', 'mclust', 'XML', 'R2HTML', 'plyr', 'e1071',
                'lattice', 'foreach', 'iterators')
  currentPack <- installed.packages()
  for(pack in biocList)
    if(!pack %in% currentPack){
      cat('Installing required package:', pack, '\n')
      biocLite(pack)
    }
	for(pack in cranList)
	  if(!pack %in% currentPack){
	    cat('Installing required package:', pack, '\n')
	    install.packages(pack)
	  }
  
  require(limma, quietly = TRUE)
	require(DNAcopy, quietly = TRUE)
	require(preprocessCore, quietly = TRUE)
#	require(AMORE, quietly = TRUE)
	require(affy, quietly = TRUE)
#	require(robustbase, quietly = TRUE)
	require(mclust, quietly = TRUE)
#	require(tcltk, quietly = TRUE)
	require(annotate, quietly = TRUE)
	require(XML, quietly = TRUE)
	require(R2HTML, quietly = TRUE)
	require(plyr, quietly = TRUE)
	require(e1071, quietly = TRUE)
	require(lattice, quietly = TRUE)
	require(foreach, quietly = TRUE)
	require(iterators, quietly = TRUE)
	cat('\tDone\n')
	