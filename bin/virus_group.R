library(plyr)
library(gplots)
#resume = read.table("result_summary.csv", sep="\t", header=T)
#matrix = read.table("matrix_LISTE_MAD_3D.txt", sep="\t", skip=1, header=F)
#matrix = read.table("matrix_LISTE_MAD_5UTR.txt", sep="\t", skip=1, header=F)
matrix = read.table("matrix_LISTE_MAD_VP1.txt", sep="\t", skip=1, header=F)
# Select column of interest
#interest = resume[, c("X3D_sequences","Serotype_3D")]
# Remove empty lines
#interest = interest[!apply(interest == "", 1, all),]
#interest$X3D_sequences = gsub('|', '_', interest$X3D_sequences,fixed = TRUE)

header=matrix[,1]
matrix=matrix[,2:dim(matrix)[2]]
colnames(matrix)=header
rownames(matrix)=header
species=unique(sub("-.+", "",header))

res = matrix(NA, length(species), length(species))
colnames(res) = species
rownames(res) = species
for(s in species){
	for(s2 in species){
		res[s,s2] = mean(as.matrix(matrix[grepl(s, header),grepl(s2, header)]))
	}	
}
#pdf("enterovirus_distance_3D.pdf")
#pdf("enterovirus_distance_5UTR.pdf")
pdf("enterovirus_distance_VP1.pdf")
heatmap.2(res,tracecol=NA)
dev.off()
