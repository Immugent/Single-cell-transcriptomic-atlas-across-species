# --------- Prog Info ---------
# R version 3.6.3
# Purpose: Cell clustering.

# 00. Prepare-------------------
print("library...")
library(Seurat)
library(ggplot2)
library(cowplot)
setwd("scRNA")
cell <- read.table("01.Rscript/cells.txt",header=TRUE)

# 01. Load the data-------------------
print("Loading...")
# Human
Human.data <- read.table("01.Data/Human.txt",header=TRUE)
colnames(Human.data) <- paste(colnames(Human.data),"Human",sep="_")
Human.data <- Human.data[,rownames(cell[which(cell$orig.ident == "Human"),])]
Human <- CreateSeuratObject(counts = Human.data, project = "Human", min.cells = 3, min.features = 200)
Human$species <- "Human"
# Monkey
Monkey.data <- read.table("01.Data/Macaca.txt",header=TRUE)
colnames(Monkey.data) <- paste(colnames(Monkey.data),"Monkey",sep="_")
Monkey.data <- Monkey.data[,rownames(cell[which(cell$orig.ident == "Monkey"),])]
Monkey <- CreateSeuratObject(counts = Monkey.data, project = "Monkey", min.cells = 3, min.features = 200)
Monkey$species <- "Monkey"
# Pig
Pig.data <- read.table("01.Data/Sus.txt",header=TRUE)
colnames(Pig.data) <- paste(colnames(Pig.data),"Pig",sep="_")
Pig.data <- Pig.data[,rownames(cell[which(cell$orig.ident == "Pig"),])]
Pig <- CreateSeuratObject(counts = Pig.data, project = "Pig", min.cells = 3, min.features = 200)
Pig$species <- "Pig"
# Rat
Rat.data <- read.table("01.Data/Rat.txt",header=TRUE)
colnames(Rat.data) <- paste(colnames(Rat.data),"Rat",sep="_")
Rat.data <- Rat.data[,rownames(cell[which(cell$orig.ident == "Rat"),])]
Rat <- CreateSeuratObject(counts = Rat.data, project = "Rat", min.cells = 3, min.features = 200)
Rat$species <- "Rat"
# Mouse
Mouse.data <- read.table("01.Data/Mus.txt",header=TRUE)
colnames(Mouse.data) <- paste(colnames(Mouse.data),"Mouse",sep="_")
Mouse.data <- Mouse.data[,rownames(cell[which(cell$orig.ident == "Mouse"),])]
Mouse <- CreateSeuratObject(counts = Mouse.data, project = "Mouse", min.cells = 3, min.features = 200)
Mouse$species <- "Mouse"
# Frog
Frog.data <- read.table("01.Data/Xenopus.txt",header=TRUE)
colnames(Frog.data) <- paste(colnames(Frog.data),"Frog",sep="_")
Frog.data <- Frog.data[,rownames(cell[which(cell$orig.ident == "Frog"),])]
Frog <- CreateSeuratObject(counts = Frog.data, project = "Frog", min.cells = 3, min.features = 200)
Frog$species <- "Frog"
# Fish
Fish.data <- read.table("01.Data/Danio.txt",header=TRUE)
colnames(Fish.data) <- paste(colnames(Fish.data),"Fish",sep="_")
Fish.data <- Fish.data[,rownames(cell[which(cell$orig.ident == "Fish"),])]
Fish <- CreateSeuratObject(counts = Fish.data, project = "Fish", min.cells = 3, min.features = 200)
Fish$species <- "Fish"

save.image("02.Results/RData/data.RData")

# 02. Normalized-------------------
print("Normalized...")
sample.list <- list(Human,Monkey,Pig,Rat,Mouse,Frog,Fish)
for (i in 1:length(sample.list)){
	sample.list[[i]] <- NormalizeData(sample.list[[i]], verbose = FALSE)
	sample.list[[i]] <- FindVariableFeatures(sample.list[[i]], selection.method = "vst", nfeatures = 4000, verbose = FALSE)
}

# 03. Integration-------------------
print("Integration...")
sample.anchors <- FindIntegrationAnchors(object.list = sample.list, dims = 1:40,anchor.features = 4000)
sample.integrated <- IntegrateData(anchorset = sample.anchors, dims = 1:40)

# 04. analysis (Run UMAP)-------------------
print("TSNE/UMAP...")
DefaultAssay(sample.integrated) <- "integrated"
sample.integrated <- ScaleData(sample.integrated, verbose = FALSE)
sample.integrated <- RunPCA(sample.integrated, npcs = 40, verbose = FALSE)
####
rm(sample.list)
####
sample.integrated <- FindNeighbors(sample.integrated, reduction = "pca", dims = 1:30)
sample.integrated <- FindClusters(sample.integrated, resolution = 0.2)
sample.integrated <- RunTSNE(sample.integrated, reduction = "pca", dims = 1:30)
save.image("02.Results/RData/cross.res0.2.RData")

tsne <- as.data.frame(sample.integrated@reductions$tsne@cell.embeddings)
write.table(tsne,"02.Results/TSNE/00.TSNE-res0.2.txt",row.names=TRUE,sep="\t",quote=TRUE)


# 05. Plot TSNE-------------------
sample.integrated$species <- factor(sample.integrated$species ,levels = c("Fish","Fish","Fish","Rat","Pig","Pig","Human"))
 #TSNE
t1 <- DimPlot(sample.integrated, reduction = "tsne", group.by = "species",pt.size =0.01) 
t2 <- DimPlot(sample.integrated, reduction = "tsne", group.by = "species",pt.size =0.01)  + NoLegend()
t3 <- DimPlot(sample.integrated, reduction = "tsne", label = TRUE, repel = TRUE,pt.size =0.01)
t4 <- DimPlot(sample.integrated, reduction = "tsne", label = TRUE, repel = TRUE,pt.size =0.01) + NoLegend()
pdf("02.Results/TSNE/00.FigS1B.pdf",width = 6, height = 5)
t1
t2
dev.off()
pdf("02.Results/TSNE/01.TSNE-bycluster.pdf",width = 6, height = 5)
t3
t4
dev.off()
t5 <- DimPlot(object = sample.integrated, reduction = "tsne", split.by = "species",label = TRUE,label.size = 5,ncol=4,pt.size =0.01)
pdf("02.Results/TSNE/02.TSNE.splitbyspecies.pdf",width = 15, height = 8)
t5
dev.off()


# 06. Change the name of the cell type-------------------
load("02.Results/RData/cross.res0.2.RData")
sample.integrated$species <- factor(sample.integrated$species ,levels = c("Fish","Fish","Fish","Rat","Pig","Pig","Human"))
newnames <- c("C1","C1","C2","C3","C4","C1","C1","C2","C5","C1","C1","C6","C7","C8")
names(newnames) <- levels(sample.integrated)
sample.integrated <- RenameIdents(sample.integrated,newnames)
n_t1 <- DimPlot(sample.integrated, reduction = "tsne", label = TRUE, repel = TRUE,pt.size =0.01)
n_t2 <- DimPlot(sample.integrated, reduction = "tsne", label = TRUE, repel = TRUE,pt.size =0.01) + NoLegend()
pdf("02.Results/TSNE/01.Fig1B.pdf",width = 6, height = 5)
n_t1
n_t2
dev.off()

### marker heatmap
markers <- read.table("02.Results/TSNE/Final-paper/03.Marker-celltype-0419.txt",header=TRUE,sep="\t")
write.table(markers,"02.Results/markers.txt",sep="\t",row.names=TRUE)

top30 <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
genes <- as.vector(top30$gene)
h1 <- DoHeatmap(sample.integrated, features = genes) + theme(axis.text.y = element_text(size = 0)) #+ NoLegend()
pdf("02.Results/Heatmap/02.Fig1C.pdf",width = 20, height = 10)
h1
dev.off()

newnames <- c("T cell","T cell","B cell","Mc-1","Ne","T cell","T cell","B cell","NK","T cell","T cell","Monocyte","PlasmaB","Mc-2")
names(newnames) <- levels(sample.integrated)
sample.integrated <- RenameIdents(sample.integrated,newnames)
n_t1 <- DimPlot(sample.integrated, reduction = "tsne", label = TRUE, repel = TRUE,pt.size =0.01)
n_t2 <- DimPlot(sample.integrated, reduction = "tsne", label = TRUE, repel = TRUE,pt.size =0.01) + NoLegend()
pdf("02.Results/TSNE/01.TSNE-bycelltype.pdf",width = 6, height = 5)
n_t1
n_t2
dev.off()
n_s1 <- DimPlot(object = sample.integrated, reduction = "tsne", split.by = "species",label = TRUE,label.size = 5,ncol=4,pt.size =0.8)
pdf("02.Results/TSNE/03.Fig1E.pdf",width = 15, height = 8)
n_s1
dev.off()

n_s2 <- DimPlot(object = sample.integrated, reduction = "tsne", split.by = "species",label = FALSE,label.size = 5,ncol=4,pt.size =0.8)
n_s3 <- DimPlot(object = sample.integrated, reduction = "tsne", split.by = "species",label = TRUE,label.size = 5,ncol=4,pt.size =0.8) + NoLegend()
pdf("02.Results/TSNE/02.TSNE.split-0.8_newnames-nolengend-0417.pdf",width = 15, height = 8)
n_s2
n_s3
dev.off()

# 07. Marker plot-------------------
markers_P <- c("CD3E","CD19","C1QC","C1QC","ADAM8","ITGA2","CD68","IGHG3")
n <- length(markers_P)
for (i in 1:n){
  marker <- markers_P[i]
  f <- FeaturePlot(object = sample.integrated,features = marker,cols = c("#C0C0C0", "#E41C12"),pt.size = 0.01,sort.cell=TRUE,reduction="tsne",slot = "data")
  pdf(paste("02.Results/Markers/",marker,".pdf",sep=""),width = 6, height = 5)
  print(f)
  dev.off()
}
