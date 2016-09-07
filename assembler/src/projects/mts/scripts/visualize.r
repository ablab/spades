x <- read.csv("canopy/points_matrix.csv", header=FALSE)
mx <- as.matrix(x)
vecs <- mx[,-1]
d <- 1 - abs(cor(t(vecs)))
# d <- as.matrix(dist(vecs, diag=TRUE, upper=TRUE))
fit <- cmdscale(d,eig=TRUE, k=2) 
df <- data.frame(fit$points)
plot(df)
sdf <- split(df, mx[,1])
len <- length(sdf)
cols <- rainbow(len)
for (i in 1:len) {
	points(sdf[[i]], col=cols[i])
}