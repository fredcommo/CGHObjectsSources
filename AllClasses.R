# Class cghObj with 2 children AgilentObj, AffyObj.

setClass('cghObj', representation(	info = 'character',									# a vector of char containing.... specify what!
													cnSet = 'data.frame',							# a data.frame containing the CN probes values
													snpSet = 'data.frame', 						# a data.frame containing the SNIP probes values, when available.
													param = 'vector',								# a vector of mixted values (boolean/numeric) containing the parameters used by CBS algorithm.
													segTable = 'data.frame',						# a data.frame containing the segment values. Availabe after the segmenation step.
													probesDensity = 'ANY',						# a xyplot (of class 'trellis') representing the probe density as a gaussian mixture.
													gProfile = 'ANY')									# a xyplot (of class 'trellis') representing the genomic profile.
						)
cghObj = function(info, cnSet = data.frame(), snpSet = data.frame(), param = NA, segTable = data.frame(),
							probesDensity = xyplot(c(0, 1)~c(0, 1), type  ='n'), gProfile = xyplot(c(0, 1)~c(0, 1), type  ='n')){
	new('cghObj', info = info, cnSet = cnSet, snpSet = snpSet, param = param, segTable = segTable,
						probesDensity = probesDensity, gProfile = gProfile)
}

setClass('AffyObj', contains = 'cghObj')
AffyObj <- function(info, cnSet = data.frame(), snpSet = data.frame(), param = NA, segTable = data.frame(),
							probesDensity = xyplot(c(0, 1)~c(0, 1), type  ='n'), gProfile = xyplot(c(0, 1)~c(0, 1), type  ='n')){
		new('AffyObj', cghObj(info = info, cnSet = cnSet, snpSet = snpSet, param = param, segTable = segTable,
										probesDensity = probesDensity, gProfile = gProfile))
		}

setClass('AgilentObj', contains = 'cghObj')
AgilentObj <- function(info, cnSet = data.frame(), snpSet = data.frame(), param = NA, segTable = data.frame(),
									probesDensity = xyplot(c(0, 1)~c(0, 1), type  ='n'), gProfile = xyplot(c(0, 1)~c(0, 1), type  ='n')){
		new('AgilentObj', cghObj(info = info, cnSet = cnSet, snpSet = snpSet, param = param, segTable = segTable,
											probesDensity = probesDensity, gProfile = gProfile))
		}

# To do : build a Nimblegen Class