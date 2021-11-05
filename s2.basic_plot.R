# --------- Prog Info ---------
# R version 3.6.3
# Purpose: Basic plot (TSNE, Doheatmap, Fraction) for a given cluster.

# 0. Prepare -------------------
library(optparse)
# construct options
option_list <- list(make_option(opt_str = "--cluster",action = "store",dest = "cluster",help = "Cluster need basic plots"),
					make_option(opt_str = "--outdir",action = "store",dest = "outdir",metavar = "OUTDIR"),
					make_option(opt_str = "--species",action = "store",dest = "species",type = "numeric",default = 127,help = "flag of species which to plot top30 heatmap, default all (127)"))
parse <- OptionParser(option_list = option_list) # instance
args <- parse_args(object = parse) # parse cmd
if (is.null(args$cluster) | is.null(args$outdir)){
	stop("--cluster/--outdir are required and expected one argument, please check!")
}
if (args$species < 1 | args$species > 127){
	stop("--species out of range, should be [1,127]")
}
library(ggplot2)
library(ggpubr)
library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)

# 1. Function Definition ------------------
parse.species <- function(num){
# Parse species flag
# Args:
#	num: species flag, a number in [1,127]
# Returns:
#	species name character
	species <- c()
	bites <- list(Fish = 2^0,Frog = 2^1,Mouse = 2^2,Rat = 2^3,Pig = 2^4,Monkey = 2^5,Human = 2^6)
	for (species_name in names(bites)){
		if (bitwAnd(bites[[species_name]],num)){
			species <- c(species,species_name)
		}
	}
	return(species)
}

which.type <- function(obj,name){
# Determin which cluster types the name belongs to
# Args:
#	obj: Seurat object
#	name: the name to determine
# Returns:
#	the current and next cluster types (primary_type, secondary_type, tertiary_type)
	if(name %in% obj$tertiary_type){
		write(paste("WARNING: No subcluster for:",name),stdout())
		return(list(this = "tertiary_type",following = NULL))
	}else if(name %in% obj$secondary_type){
		return(list(this = "secondary_type",following = "tertiary_type"))
	}else if(name %in% obj$primary_type){
		return(list(this = "primary_type",following = "secondary_type"))
	}else{
		stop(paste("Wrong cluster name:",name))
	}
}

do.heatmap <- function(obj,spe,prefix,cutoff){
# Top 30 Heatmap for given spe
# Args:
#	obj: Seurat object
#	spe: species
#	prefix: prefix of outfile
#	cutoff: a list of p, min.pct, fold change cutoff
	print(unlist(cutoff))
	subdata <- subset(obj,subset = species %in% spe)
	degs <- FindAllMarkers(subdata,slot = "data",only.pos = T,logfc.threshold = log(cutoff$fc),min.pct = cutoff$min.pct,return.thresh = cutoff$p)
	degs <- degs %>% 
		group_by(cluster) %>%
		top_n(30,wt = avg_logFC)
	write.table(x = degs,file = paste0(prefix,"_markers_top30.txt"),sep = "\t",col.names = T,row.names = F,quote = F)
	subdata <- ScaleData(subdata,features = degs$gene)
	p1 <- DoHeatmap(subdata,features = degs$gene) +
		guides(color = FALSE,fill = guide_colorbar(override.aes = list(size = 5))) +  # no color legend and set fiil lengend key size
		theme(legend.title = element_text(size = 15),legend.text = element_text(size = 15))
	p2 <- p1 + theme(axis.text.y = element_blank()) # no gene labeled
	ggsave(filename = paste0(prefix,"_heatmap.pdf"),plot = p1,device = "pdf",width = 20,height = 20)
	ggsave(filename = paste0(prefix,"_heatmap_no_genes.pdf"),plot = p2,device = "pdf",width = 20,height = 10)
}

# 2. Set global parameters or from cmd --------------------
rds_file <- "scRNA/02.Results/RData/integrate_classified.rds"
cluster <- args$cluster
species <- parse.species(args$species)
outdir <- args$outdir
## construct outdir
if (!dir.exists(outdir)){
	dir.create(outdir,recursive = T) # mkdir -p
}
setwd(outdir)
## log 
write("================= 1. Parameters List ===================",stdout())
write(paste("Raw RDS:",rds_file),stdout())
write(paste("Cluster:",cluster),stdout())
write(paste("Species:",paste(species,collapse = ",")),stdout())
write(paste("Outdir:",outdir),stdout())

# 3. Pre-process scRNA data --------------
write("================= 2. Pre-Process ===================",stdout())
scrna <- readRDS(rds_file)
DefaultAssay(scrna) <- "integrated"
types <- which.type(scrna,cluster)
Idents(scrna) <- scrna@meta.data[[types$this]]
scrna <- subset(scrna,idents = cluster)
scrna <- RunTSNE(scrna,reduction = "pca",dims = 1:30)
if (is.null(types$following)){
	Idents(scrna) <- scrna@meta.data[[types$this]]
} else {
	Idents(scrna) <- factor(scrna@meta.data[[types$following]])
}
DefaultAssay(scrna) <- "RNA"
print(addmargins(table(Idents(scrna),scrna$species))) # log

# 4. TSNE Plot -----------------
write("================= 3. TSNE Plot ===================",stdout())
if (!dir.exists("TSNE")){
	dir.create("TSNE",recursive = T)
}
p_integrate <- DimPlot(scrna,reduction = "tsne",group.by = "ident",label = F,pt.size = 0.05)
p_split <- DimPlot(scrna,reduction = "tsne",group.by = "ident",split.by = "species",pt.size = 0.05,ncol = 7)
ggsave(filename = "TSNE/integrate.pdf",plot = p_integrate,width = 5,height = 4,device = "pdf")
ggsave(filename = "TSNE/split.pdf",plot = p_split,width = 15,height = 3,device = "pdf")
ggsave(filename = "TSNE/integrate_nolegend.pdf",plot = p_integrate + NoLegend(),width = 5,height = 4,device = "pdf")
ggsave(filename = "TSNE/split_nolegend.pdf",plot = p_split + NoLegend(),width = 15,height = 3,device = "pdf")

# 5. Fraction Calculate and Visualization -----------------
write("================= 4. Fraction Calculate & Visualization ===================",stdout())
if (!dir.exists("Fraction")){
	dir.create("Fraction",recursive = T)
}
df_frac <- as.data.frame(table(scrna$species,Idents(scrna)))
colnames(df_frac) <- c("Species","Type","Frequency")
df_frac$Species <- factor(df_frac$Species,levels = levels(scrna$species))
write.table(x = df_frac,file = "Fraction/frac_df.txt",sep = "\t",col.names = T,row.names = F,quote = F)
## pie chart
p <- list()
for (current_species in levels(df_frac$Species)){
    data <- df_frac[df_frac$Species == current_species,]
    p[[current_species]] <- ggpie(data,x = "Frequency",fill = "Type",color = "Type",lab.pos = "in",lab.font = c(0,"bold","black")) + 
        labs(fill = "",title = current_species) + 
        guides(color = FALSE) +
        theme(plot.title = element_text(size = 15,hjust = 0.5)
              ,plot.margin = margin(-1,-5,-1,-5,unit = "pt")
              ,legend.title = element_text(size = 15),legend.text = element_text(size = 15)
              ,legend.box.margin = margin(0,0,18,0,unit = "pt"))
    if (current_species == "Human"){
        # the last one(Human) owns legend
        p[[current_species]] <- p[[current_species]] + theme(legend.position = "right")
    } else {
        p[[current_species]] <- p[[current_species]] + theme(legend.position = "none")
    }
}
ggsave(filename = "Fraction/piechart.pdf",plot = wrap_plots(p,ncol = 7),device = "pdf",width = 12,height = 3)
## construct out data
df_sum <- aggregate(df_frac$Frequency,by = list(df_frac$Species),sum)
colnames(df_sum) <- c("Species","Num")
df_frac[["label"]] <- ""
for (i in 1:nrow(df_frac)){
    df_frac$label[i] <- round(100*df_frac$Frequency[i]/df_sum$Num[df_sum$Species==df_frac$Species[i]],2)
    df_frac$label[i] <- paste0(df_frac$label[i],"%(",df_frac$Frequency[i],")")
}
df <- dcast(df_frac,Type~Species)
write.table(x = df,file = "Fraction/fraction_table.csv",sep = ",",col.names = T,row.names = F,quote = F)

# 6. Top 30 Heatmap -------------------
write("================= 5. Top 30 Heatmap ===================",stdout())
if (!dir.exists("Heatmap")){
	dir.create("Heatmap",recursive = T)
}
threshold <- list(p = 0.05,fc = exp(0.25),min.pct = 0.1)
if (is.null(types$following)){
	write(paste("No subcluster for:",types$following,"skip heatmap!"),stdout())
} else {
	write("Integrate:",stdout())
	do.heatmap(obj = scrna,spe = levels(scrna$species),prefix = "Heatmap/integrate",cutoff = threshold)
	for (current_species in species){
		write(paste0(current_species,":"),stdout())
		do.heatmap(obj = scrna,spe = current_species,prefix = paste0("Heatmap/",current_species),cutoff = threshold)
	}
}
