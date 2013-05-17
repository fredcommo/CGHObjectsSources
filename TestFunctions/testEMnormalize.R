# Test EMnormalise.

# cghObj = function(info, cnSet = data.frame(), snpSet = data.frame(), param = NA){
	# new('cghObj', info = info, cnSet = cnSet, snpSet = snpSet, param = param)

N = 170334
p = c(0.2288778, 0.07877624, 0.15586947, 0.1271292, 0.1934367, 0.1844995)
p = c(p, 1 - sum(p))
m = c(-0.1547146, -0.07450492, -0.01949387,  0.001256451, 0.05373657, 0.1250104, 0.3049194)
s = c(0.001303411, 0.002754407, 0.05224733,  0.0005247088, 0.001008751, 0.00204692, 0.001393142)

LR = c()
for (i in 1:length(m)){
	LR = c(LR, rnorm(N*p[i], m[i], 4*sqrt(s[i])))
}
if(length(LR)<N) LR = c(LR, rep(NA, N - length(LR)))

testData = data.frame(ProbeName = paste("A", seq(1,N), sep = "_"),
								ChrNum = rep(1, N),
								ChrStart = rep(1, N),
								genomicPos = rep(1, N),
								Log2Ratio = LR)
testObject = cghObj(info = c(platform = 'Affymetrix'), cnSet = testData)

scriptPath = "/Users/fredcommo/Documents/Projet Safir/Safir R Sources/CGHObjects/"
setwd(scriptPath)
source('SourceCodeObj.R')
testRes <- EMnormalize(testObject)
