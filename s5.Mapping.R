#Analysis the spleen data from GSE159929

rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(Seurat)
library(ggplot2)
library(openxlsx)
library(monocle)
library(tidyverse)   
library(patchwork)

# Load data 
raw <- read.csv("rds_data/GSM4850589_Spleen_Counts.csv", header=T, row.names = 1)
dim(raw)
raw[1:4,1:4]
# Load metadata 
#metadata <- read.csv("GSE132465_annotation.csv",  header=T)
spleen <- CreateSeuratObject(counts = raw, min.cells = 3, min.features = 200)
spleen
median(spleen$nFeature_RNA)
pbmc <- spleen

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc[["percent.rb"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[LS]")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(pbmc))
pbmc[["percent.HB"]]<-PercentageFeatureSet(pbmc, features=HB.genes) 

minGene=300
maxGene=8000
maxUMI=22000
pctMT=10
pctHB=1

pbmc <- subset(pbmc, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & 
                 nCount_RNA < maxUMI & percent.mt < pctMT & percent.HB < pctHB)
table(pbmc$orig.ident)

plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(pbmc,  pt.size = 0.01,cols = c("#265297","#DD322E"),
                       features = plot.featrures[i]) + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
violin

pbmc <- NormalizeData(pbmc)  # 解决每个细胞测序深度不同的问题
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(pbmc))
s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(pbmc))
pbmc <- CellCycleScoring(pbmc, g2m.features=g2m_genes, s.features=s_genes)

pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000) %>% ScaleData(vars.to.regress=c("percent.mt","percent.rb","S.Score","G2M.Score"))

pbmc <- RunPCA(pbmc, verbose = F)
ElbowPlot(pbmc, ndims = 50)
pc.num=1:30
pbmc <- pbmc %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
pbmc <- pbmc %>% FindNeighbors(dims = pc.num) %>% FindClusters(resolution=0.3) %>% FindClusters(resolution=0.5)
saveRDS(pbmc,"rds_data/one_spleen.rds")

#pbmc <- readRDS("rds_data/one_spleen.rds")  #read rds data directly
#pbmc

reference1 <- readRDS("rds_data/reference.rds")#直接读取注释文件
reference1 <- subset(reference1, subset = primary_type %in% c('T','B','Mc','PlasmaB'))
reference1 <- subset(reference1, subset = secondary_type != c('UndefT'))
table(reference1$secondary_type)

reference <- readRDS("rds_data/integrate_classified.rds") #read integrated data as reference
table(reference$species,reference$secondary_type)
colnames(reference@meta.data)
reference1 <- subset(reference, subset = species =='Human')
head(reference1@meta.data)
table(reference1$secondary_type)

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(reference1))
s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(pbmc))
reference1 <- CellCycleScoring(reference1, g2m.features=g2m_genes, s.features=s_genes)
reference1 <- FindVariableFeatures(reference1, nfeatures = 2000) %>% ScaleData(vars.to.regress=c("percent.mt","percent.rb","S.Score","G2M.Score"))

###Map two datasets
###### Look for anchors
anchors <- FindTransferAnchors(reference = reference1,
                               query = pbmc,
                               normalization.method = "LogNormalize",
                               reduction = "pcaproject",
                               reference.reduction = "pca",
                               dims = 1:30)

###### Anchor conversion
pbmc <- MapQuery(
  anchorset = anchors,
  query = pbmc,
  reference = reference1,
  refdata = list(celltype.l1 = "primary_type",  
                 celltype.l2 = "secondary_type"  
  ),       
  reference.reduction = "pca")
colnames(pbmc@meta.data)
table(pbmc$predicted.celltype.l2)
pbmc
pbmc <- subset(pbmc, subset = predicted.celltype.l1 %in% c('T','B','Mc','PlasmaB'))
pbmc <- subset(pbmc, subset = predicted.celltype.l2 != c('UndefT'))
table(pbmc$predicted.celltype.l2)

#dir.create("Mapping")
p1 <- DimPlot(reference1,group.by = "primary_type", label = T, label.size = 5) + NoLegend()
p2 <- DimPlot(reference1,group.by = "secondary_type", label = T, label.size = 5) + NoLegend()
p3 <- DimPlot(pbmc, group.by = "predicted.celltype.l1", label = T, label.size = 5) + NoLegend()
p4 <- DimPlot(pbmc, group.by = "predicted.celltype.l2", label = T, label.size = 5) + NoLegend()
p <- (p1|p2)/(p3|p4)
p
ggsave("Mapping/Celltype_filter.pdf", p, width = 12, height = 10)

table(pbmc$predicted.celltype.l1)
table(reference1$primary_type)
round(table(reference1$primary_type)/sum(table(reference1$primary_type)),2)

###Calculate relevance
pbmc <- readRDS("rds_data/one_spleen.rds")  

pbmc <- readRDS("rds_data/0h_spleen.rds")  
#table(pbmc$Donor,pbmc$predicted.celltype.l2)
table(pbmc$Donor,pbmc$predicted.celltype.l1)
pbmc <- subset(pbmc, subset = Donor == '356C')
pbmc

library(ape)
library(dplyr)
library(psych)
library(ComplexHeatmap)
library(edgeR)
library(limma)
library(grid)
library(circlize)
library("GetoptLong")
library(corrplot)

head(pbmc@meta.data)
head(reference1@meta.data)

df1<-as.data.frame(AverageExpression(pbmc,features = rownames(pbmc),group.by = "predicted.celltype.l1")[["RNA"]])
colnames(df1) <- paste(colnames(df1),"Query",sep = "_")

#secondary_type #primary_type
df2<-as.data.frame(AverageExpression(reference1,features = rownames(pbmc),group.by = "primary_type")[["RNA"]])
colnames(df2) <- paste(colnames(df2),"Ref",sep = "_")
df1$gene <- rownames(df1)
df2$gene <- rownames(df2)
head(df1)
head(df2)
df <- merge(df1,df2[,c(-4,-5,-6)],by = "gene")#[,c(-3,-4)]
df <- merge(df1[,c(-6,-7,-8,-11)],df2[,c(-4,-7,-6,-8,-9,-10,-13,-14,-15,-16)],by = "gene")

rownames(df) <- df$gene
df <- df[,-1]
dim(df)

head(df)
#head(df)

#df.s <- scale(df[genes$V1,c(1:3,7:15,4:6)],center = T,scale = T) 
df.s <- scale(df,center = T,scale = T)

df.t <- as.matrix(t(df.s))
#the distance measure #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski" 
d <- dist(df.t,method = "euclidean") 
#d <- as.dist(1-cor(t(df.t)))  
hc <- hclust(d,method = "ward.D2") 
py <- as.phylo(hc)
row_dend = as.dendrogram(hc)
plot(py,type="phylogram",use.edge.length = 0,
     cex = 0.8,edge.width = 2,font = 3,label.offset=0.01,
     adj = 0.1)

###Calculate the correlation coefficient
df.pearson <- cor(df.s,method="spearman") #pearson, spearman
df.pearson
n <- dim(df)[2]-1
colnames(df.pearson) <- gsub("X","",colnames(df.pearson))
colnames(df.pearson) <- gsub("_Mean","",colnames(df.pearson))
rownames(df.pearson) <- gsub("X","",rownames(df.pearson))
rownames(df.pearson) <- gsub("_Mean","",rownames(df.pearson))

# heatmap
f1 <- colorRamp2(breaks = c(0.7,0.85,0.95,1), 
                 c("#114A85", "#FFFFFF","#B61F2E","#8D1114"))
p <- Heatmap(df.pearson, name ="Correlation", # coefficient
             column_title = " ",row_title = " ",
             # rect_gp = gpar(col = "white"),
             col = f1,
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(df.pearson[i, j] > -1)
                 grid.text(sprintf("%.2f", df.pearson[i, j]),
                           x, y, gp = gpar(fontsize = 10))
             }, 
             cluster_rows = row_dend,cluster_columns = row_dend,
             # row_names_gp = gpar(fontsize = 10,font=4),
             show_column_names = TRUE,show_row_names = TRUE,row_names_side = "right",
             #        row_dend_width = unit(4, "mm"),column_dend_height = unit(2,"mm"),
             gap = unit(1,"mm"),row_dend_side = "left",row_title_side = "left",
             column_dend_side = "top",column_title_rot = 0,column_title_side = "top")
p

pdf("0.all_Correlation_byRep.pdf",width = 5.1,height = 3.8)
p
dev.off()

library(ape)
library(dplyr)
library(psych)
library(ComplexHeatmap)
library(edgeR)
library(limma)
library(grid)
library(circlize)
library("GetoptLong")
library(corrplot)
library(factoextra)

df.s <- scale(df,center = T,scale = T)
# 03. 层次聚类
df.t <- as.matrix(t(df.s))
#d <- dist(df.t,method = "euclidean") 
#?get_dist
d <- get_dist(df.t,method = "spearman") 
#d <- as.dist(1-cor(t(df.t)))  
hc <- hclust(d,method = "ward.D2")
py <- as.phylo(hc)
row_dend = as.dendrogram(hc)
plot(py,type="phylogram",use.edge.length = 0,
     cex = 0.8,edge.width = 2,font = 3,label.offset=0.01,
     adj = 0.1)

df.pearson <- cor(df.s,method="spearman")
df.pearson
n <- dim(df)[2]-1
colnames(df.pearson) <- gsub("X","",colnames(df.pearson))
colnames(df.pearson) <- gsub("_Mean","",colnames(df.pearson))
rownames(df.pearson) <- gsub("X","",rownames(df.pearson))
rownames(df.pearson) <- gsub("_Mean","",rownames(df.pearson))

f1 <- colorRamp2(breaks = c(0.7,0.8,0.9,1), 
                 c("#114A85", "#FFFFFF","#B61F2E","#8D1114"))#set color

##heatmap
p <- Heatmap(df.pearson, name ="Correlation", 
             column_title = " ",row_title = " ",
             # rect_gp = gpar(col = "white"),
             col = f1,
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(df.pearson[i, j] > -1)
                 grid.text(sprintf("%.2f", df.pearson[i, j]),
                           x, y, gp = gpar(fontsize = 10))
             }, 
             cluster_rows = row_dend,cluster_columns = row_dend,
             show_column_names = TRUE,show_row_names = TRUE,row_names_side = "right",
             gap = unit(1,"mm"),row_dend_side = "left",row_title_side = "left",
             column_dend_side = "top",column_title_rot = 0,column_title_side = "top")



p 
