# .Iqr <- function(x, n){
	# xprim = x[-(n+1)]
	# q = quantile(xprim, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
	# iqr = q[3] - q[1]
	# if(x[n+1] < q[2]-1.5*iqr | x[n+1] > q[2]+1.5*iqr)
		# x[n+1] <- q[2]
	# return(x)
# }

# .supprOutliers <- function(x, n = 5){
	# k = 0
	# p = length(x)
	# xnew <- xprim <- c(x[1:n], x, x[(p-n+1):p])
	# for(i in (n+1):p){
		# tmp <- xprim[(i-n):(i+n)]
		# xnew[i] <- .Iqr(tmp, n)[n+1]
		# }
	# xnew <- xnew[(n+1):(p+n)]
	# return(xnew)
# }

# .supprOutliers2 <- function(x, n = 5){
	# k = 0
	# p = length(x)
	# xnew <- xprim <- c(x[1:n], x, x[(p-n+1):p])
	# for(i in (n+1):p){
		# tmp <- xprim[(i-n):(i+n)]
		# if(.Iqr(tmp, n)){
			# xnew[i] <- median(tmp[-i])
			# k = k + 1
			# }
		# }
	# xnew <- xnew[(n+1):(p+n)]
	# cat(k, 'replacements\n')
	# return(xnew)
# }

###########################
.Iqr <- function(x, n){
	xprim = x[-(n+1)]
	q = quantile(xprim, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
	iqr = q[3] - q[1]
	return(x[n+1] < q[2]-1.5*iqr | x[n+1] > q[2]+1.5*iqr)
}

.supprOutliers3 <- function(x, n = 5){
	k = 0
	p = length(x)
	xnew <- c(x[1:n], x, x[(p-n+1):p])
	for(i in (n+1):p){
		tmp <- xnew[(i-n):(i+n)]
		if(.Iqr(tmp, n)){
			xnew[i] <- median(tmp[-i])
			k = k + 1
			}
		}
	xnew <- xnew[(n+1):(p+n)]
	cat(k, 'replacements\n')
	return(xnew)
}



filterOutliers <- function(x, n = 5) {
	skip = n+1
	Q = runquantile(x[-skip], 2*n+1, c(0.25, 0.5, 0.75))
	IQR = (Q[,3]-Q[,1])*1.5
	if(is.na(x[skip]) | ((x < Q[,2] - IQR) | (x > Q[,2] + IQR))){
		x[skip]=Q[2]
		}
	return(x)
	}

s = 0.5
n1 = 1e3; n2 = 5e3; n3 = 1.25e4; n4 = 5e3; n5 = 1.25e4; n6 = 1.5e3; n7 = 5e3; n8 = 5e3; n9 = 1e3; n10 = 4e3
y <- c(rnorm(n1, 0, s),
			rnorm(n2, 0.8, s),
			rnorm(n3, 0, s),
			rnorm(n4, 0.9, s),
			rnorm(n5, 0, s),
			rnorm(n6, -0.6, s),
			rnorm(n7, 0, s),
			rnorm(n8, -0.8, s),
			rnorm(n9, 0.75, s),
			rnorm(n10, 0, s))
tags = c(rep('A', n1) , rep('B', n2), rep('C', n3), rep('D', n4), rep('E', n5), rep('F', n6), rep('G', n7), rep('G', n8), rep('H', n9), rep('I', n10))
replic <- rep(c('1', '2'), each = length(y))
tags = factor(tags)
x <- seq(1, length(y))

Q <- runquantile(y, 11, probs = c(0.25, 0.5, 0.75))
IQR <- (Q[,3] - Q[,1])*1.5
idx <- which(y<Q[,2]-IQR | y>Q[,2]+IQR)
length(idx)
y4 <- y
y4[idx] <- Q[,2][idx]

system.time(y2 <- .supprOutliers2(y))
system.time(y3 <- .supprOutliers3(y))
system.time(out <- filterOutliers(x))

par(mfrow = c(1, 3))
plot(y, y2, cex = 0.1, xlab = 'Orignal values', ylab  = 'filtered skipping values')
legend('topleft', legend = '#filtered: 7885')
plot(y, y3, cex = 0.1, xlab = 'Orignal values', ylab  = 'filtered without skipping values')
legend('topleft', legend = '#filtered: 15661')
plot(y, y4, cex = 0.1, xlab = 'Orignal values', ylab  = 'filtered using runquantile')
legend('topleft', legend = '#filtered: 4369')
par(mfrow = c(1, 1))

par(mfrow = c(1, 4))
plot(x, y, cex = 0.1, main = 'Original profile')
plot(x, y2, ylim = range(y), cex = 0.1, main = 'Filtered by skipping values')
plot(x, y3, ylim = range(y), cex = 0.1, main = 'Filtered without skipping values')
plot(x, y4, ylim = range(y), cex = 0.1, main = 'Filtered using runquantile')
par(mfrow = c(1, 1))


Y <- c(y, y2, y3, y4)
Tags <- factor(c(tags, tags, tags, tags)); levels(Tags) = levels(tags)
replic <- rep(c('1', '2', '3', '4'), each = length(y))
boxplot(Y ~ replic:Tags, col = c('steelblue1', 'indianred1', 'red3', 'seagreen'), outpch = 8, outcol = c('steelblue1', 'indianred1', 'red3', 'seagreen'), outcex = 0.5)
legend('bottomleft', legend = c('Original', 'skipping value', 'without skipping values', 'using runquantile'), fill = c('steelblue1', 'indianred1', 'red3', 'seagreen'), ncol = 2)

S1 = S2 = S3 = S4 = c()
for(tag in levels(tags)){
	S1 <- c(S1, rep(sd(y[tags == tag]), length(which(tags == tag))))
	S2 <- c(S2, rep(sd(y2[tags == tag]), length(which(tags == tag))))
	S3 <- c(S3, rep(sd(y3[tags == tag]), length(which(tags == tag))))
	S4 <- c(S4, rep(sd(y4[tags == tag]), length(which(tags == tag))))
}

par(mfrow = c(1, 4))
plot(x, y, cex = 0.1, main = 'Original profile'); lines(x, log(S1)*2, col = 'steelblue4', lwd = 3)
plot(x, y2, ylim = range(y), cex = 0.1, main = 'Filtered by skipping values'); lines(x, log(S2)*2, col = 'indianred1', lwd = 3)
plot(x, y3, ylim = range(y), cex = 0.1, main = 'Filtered without skipping values'); lines(x, log(S3)*2, col = 'red3', lwd = 3)
plot(x, y4, ylim = range(y), cex = 0.1, main = 'Filtered using runquantile'); lines(x, log(S4)*2, col = 'seagreen', lwd = 3)
par(mfrow = c(1, 1))

plot(S1, ylim = range(c(S1, S2, S3)), type = 'l'); lines(S2, col = 'indianred1'); lines(S3, col = 'seagreen')
par(mfrow = c(2, 2))
plot(x, y, cex = 0.1)
plot(x, newy, ylim = range(y), cex = 0.1)
plot(d1 <- density(y), ylim = range(d1$y)*1.5, col = 'blue', lwd = 3); lines(density(newy), col = 'red3', lwd = 3)
plot(y, newy, asp = 1)
par(mfrow = c(1, 1))

cnSet <- cnSet <- getCNset(object2)
lr <- cnSet$Log2Ratio
