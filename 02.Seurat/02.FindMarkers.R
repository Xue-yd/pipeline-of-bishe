library(Seurat)
library(dplyr)

load(file = "outfile/human_tSNE.RData")
human <- readRDS("outfile/human_tSNE.rds")

human.markers <- FindAllMarkers(human, only.pos = TRUE, 
                                 min.pct = 0.25, logfc.threshold = 0.25) 

human_markers_top20 <- human.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
human_markers_top10 <- human.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

write.table(human_markers_top20,"outfile/human_top20_markers")
write.table(human_markers_top10,"outfile/human_top10_markers")

saveRDS(human, file = "outfile/human_Markers.rds")
save(human,file="outfile/human_Markers.RData")
