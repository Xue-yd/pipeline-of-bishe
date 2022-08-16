library(Seurat)
#read data & create Seurat object
human <- read.table(file = "human_orthologs_DGE.txt",header = T)
human <- CreateSeuratObject(human,project = "Human_orthologs",assay = "RNA",
                            min.cells = 3,min.features = 200,names.delim = "-")
#calculate the distribution frequency of mitochondrial genes & Draw a violin map
human[["percent.mt"]] <- PercentageFeatureSet(human,pattern = "^MT-")
pdf("out/violin_map.pdf")
VlnPlot(object = human, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        idents = NULL,assay = NULL,pt.size = 0,ncol = 3)
dev.off()

human$log10GenesPerUMI <- log10(human$nFeature_RNA)/log10(human$nCount_RNA)

#filter the data according to the number of genes and the percent of mitochondria
human <- subset(human, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 
                & percent.mt < 10 &log10GenesPerUMI > 0.80)

human

#draw a scatter map of data
pdf("out/Feature_scatter.pdf")
plot1 <- FeatureScatter(human,feature1 = "nCount_RNA",
                        feature2 = "percent.mt")+ NoLegend()
plot2 <- FeatureScatter(human, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")+ NoLegend()
CombinePlots(plots = list(plot1, plot2))
dev.off()

#keep the genes(expressed in >10 cells)
#counts <- GetAssayData(object = human, slot = "counts")
#nonzero <- counts > 0
#keep_genes <- Matrix::rowSums(nonzero) >= 10
#counts <- counts[keep_genes,]
#human <- CreateSeuratObject(counts,meta.data = human@meta.data)

#normalization
human <- NormalizeData(human, normalization.method = "LogNormalize",scale.factor = 10000)

#identify the highly variable genes
human <- FindVariableFeatures(human,selection.method = "vst")
top10 <- head(VariableFeatures(human), 10) 
pdf("out/Variable_feature.pdf")
plot3 <- VariableFeaturePlot(human)
plot3
dev.off()
pdf("out/Variable_feature_top10.pdf")
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE,xnudge=0,ynudge=0)
plot4
dev.off()

#CombinePlots(plots = list(plot3, plot4))
#dev.off()

#scaling data, just highly variable genes
all.genes <- rownames(human)
#long time,several hours
human <- ScaleData(human,features = all.genes,vars.to.regress = "percent.mt")

#Linear dimensionalization PCA, defalut: highly variable genes
human<- RunPCA(human, features = VariableFeatures(human))

#some PCA methods
pdf("out/pca1.pdf")
VizDimLoadings(object = human, dims = 1:6, reduction = "pca")
dev.off()


pdf("out/pca_dotplot.pdf")
DimPlot(human, reduction = "pca")
dev.off()

pdf("out/pca_heatmap.pdf")
DimHeatmap(human, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()

#Identify the available dimensions of the dataset, and after PCA is complete, 
#need to determine which genes represented by the main component 
#can go downstream for analysis, where JackStraw can be used for resamping analysis. 
#JackStrawPlot visualization to see which pcs can be analyzed downstream.

human<- JackStraw(human, num.replicate = 100)  #long time
human<- ScoreJackStraw(human, dims = 1:20)
#above the line is useful pcs
pdf("out/pca_pvalue.pdf")
JackStrawPlot(human, dims = 1:20)
dev.off()

#Ranked based on the Variance rate of each pc
pdf("out/pca_SD.pdf")
ElbowPlot(object = human)
dev.off()

saveRDS(human, file = "outfile/human_PCA.rds")
save(human,file="outfile/human_PCA.RData")
