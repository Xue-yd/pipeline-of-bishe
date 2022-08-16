library(pheatmap)

human_clusters <- read.table("average_genes_human.txt")
monkey_clusters <- read.table("average_genes_monkey.txt")

Pearson_mat <- matrix(nrow = 14,ncol = 14)

rownames(Pearson_mat) <- c("H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12","H13","H14")
colnames(Pearson_mat) <- c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","M13","M14")

for(h in c(1,2,3,4,5,6,7,8,9,10,11,12,13,14))
{
  Xh <-human_clusters[,c(h)]
  for (m in c(1,2,3,4,5,6,7,8,9,10,11,12,13,14))
  {  Ym <- monkey_clusters[,c(m)]
     cc <- cor(Xh, Ym, method = "pearson")
     Pearson_mat[h,m] = cc
  }
}  

pdf("pheatmap.pdf")
d <- pheatmap(Pearson_mat, cluster_rows = T,
              cluster_cols = T,show_rownames = T,show_colnames = T)
d
dev.off()
