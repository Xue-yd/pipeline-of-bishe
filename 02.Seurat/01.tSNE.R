library(Seurat)

load(file = "outfile/human_PCA.RData")
human <- readRDS("outfile/human_PCA.rds")

# calculate the nearest neighbor graph based on pca
human <- FindNeighbors(human, dims = 1:18)
#The resolution parameter determines the number of clusters obtained from downstream clustering analysis
#suggestion: 3K cells~resolution:0.4--1.2
human <- FindClusters(human, resolution = 1.3)
#Idents() for cells clusters
head(Idents(human), 8)
#for the cell numbers of each cluster
table(human@active.ident) 

#Nonlinear dimensionalization,t-SNE
human <- RunTSNE(human, dims = 1:18)
pdf("out/t-SNE_human.pdf")
tsneplot<-TSNEPlot(human,label = TRUE, pt.size = 0.4)
tsneplot
dev.off()
#save data
saveRDS(human, file = "outfile/human_tSNE.rds")
save(human,file="outfile/human_tSNE.RData")
