library(Seurat)

human <- readRDS("outfile/human_Markers.rds")
cluster.averages <- AverageExpression(human,return.seurat = T)
head(cluster.averages[["RNA"]][, 1:5])

pdf("out/human_cluster_average_heatmap.pdf")
DoHeatmap(cluster.averages, features = unlist(TopFeatures(human[["pca"]],balanced = TRUE)), size = 3, draw.lines = FALSE,label = F)
dev.off()

write.table(cluster.averages[["RNA"]][, 1:20],"outfile/cluster_average_human.txt")
