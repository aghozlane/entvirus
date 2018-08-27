library(Matrix)
library(biclust)
library(blockcluster)

m2 = as.matrix(read.csv(file = "CVA13_contigs_bmge.tsv", sep="\t", header = FALSE))
m = as.matrix(read.csv(file = "CVA13_contigs_bmge_binary.tsv", sep="\t", header = FALSE))
m = as.matrix(forceSymmetric(m,uplo="U"))
plot.new()
par(mar=rep(0, 4))
plot.window(xlim=c(0, ncol(m)), ylim=c(0, nrow(m)), asp=1)
o <- cbind(c(row(m)), c(col(m))) - 1
rect(o[, 1], o[, 2], o[, 1] + 1, o[, 2] + 1, col=t(m)[, ncol(m):1])
res1 <- biclust(m, method=BCCC(), delta=1.5,  alpha=1, number=10)
res1 <- biclust(m, method=BCBimax(), minr=4, minc=4, number=100)
heatmapBC(m, res1)            
drawHeatmap2(m, bicResult=res1,1)

# clustering hierarchique
d = dist(m, method = "binary")
hc = hclust(d, method="ward")
plot(hc)
cluster.means = aggregate(m,by=list(cutree(hc, k = 6)), mean)
rect.hclust(hc, k = 3, border = "red")

# blockcluster
out=cocluster(m,"binary",nbcocluster = c(2,3))
summary(out)
plot(out)