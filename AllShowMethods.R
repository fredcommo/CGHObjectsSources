#showMethod
setMethod('show', signature = 'cghObj',
					function(object){
						d <- dim(object@cnSet)
						infoTab = getInfo(object)
						cat('\nInstance of class', class(object), '\n\n')
						cat('CNSet with', d[1], 'probes and', d[2], 'columns\n\n')
						cat('Array information:\n\n')
						for(i in 1:nrow(infoTab)){
							item = rownames(infoTab)[i]
							extraTabul = ifelse(nchar(item)<10, '\t\t\t\t:\t', '\t\t:\t')
							value = as.character(infoTab[i,1])
							cat('\t', item, extraTabul, value, '\n')
							}
							cat('\n')
							cat('Use getInfo(object) to get array information\n')
							cat('Use getCNset(object) to get the CGH matrix\n')
							cat('Use getSNPset(object) to get the SNP matrix (when available)\n')
							cat('Use getParam(object) to get segmentation parameters\n')
							cat('Use getSegTable(object) to get segmentation table\n')
							cat('Use getProfile(object) to display the genomic profile\n\n')
							}
				)
