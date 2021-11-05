# --------- Prog Info ---------
# R version 3.6.3
# Purpose: Differential Gene Analysis (DEA) and scatter plot for given groups

# 0. Prepare -------------------
library(optparse)
# construct options
option_list <- list(make_option(opt_str = c("-t","--treat"),action = "store",dest = "treat",metavar = "TREAT"),
					make_option(opt_str = c("-c","--control"),action = "store",dest = "control",metavar = "CTRL"),
					make_option(opt_str = c("-o","--outdir"),action = "store",dest = "outdir",metavar = "OUTDIR"),
					make_option(opt_str = "--species",action = "store",dest = "species",type = "numeric",default = 127,help = "flag of species, default all (127)"),
					make_option(opt_str = "--speciesType",action = "store",dest = "speciesType",default = NULL,help = "Type name of given species. If not NULL, given species will see as a whole."))
parse <- OptionParser(option_list = option_list) # instance
args <- parse_args(object = parse) # parse cmd
if (is.null(args$treat) | is.null(args$control) | is.null(args$outdir)){
	stop("-t/-c/-o are required and expected one argument, please check!")
}
if (args$species < 1 | args$species > 127){
	stop("--species out of range, should be [1,127]")
}
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)

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
#	obj: Seurat objective
#	name: the name to determine
# Returns:
#	the cluster types (primary_type, secondary_type, tertiary_type), if name exist in several types, return the low-level type
	if(name %in% obj$tertiary_type){
		return("tertiary_type")
	}else if(name %in% obj$secondary_type){
		return("secondary_type")
	}else if(name %in% obj$primary_type){
		return("primary_type")
	}else{
		stop(paste("Wrong Treat or Control name:",name))
	}
}

scatter.plot <- function(data,x,y,labels = NULL){
# Scatter plot for DEGs, color by change types (UP,DOWN,NOT)
# Args: 
#	data: data frame contains following columns: Gene, Change and 2 columns with loged expression (x,y)
#	x: a column name of data indicating values of x axis
#	y: a column name of data indicating values of y axis
#	labels: genes to label in graph, defalut (NULL) for top 20 avg_logfc (each 10 for UP/DOWN)
# Returns:
#	a ggplot2 object
	p <- ggplot(data,aes(.data[[x]],.data[[y]],label = Gene))+
			geom_point(aes(color = Change)) +
			geom_text_repel(data = subset(data,Gene %in% labels),fontface = "italic",size = 5,box.padding = unit(0.45, "lines")) + 
			scale_color_manual(breaks = c("NOT","UP","DOWN"),values = c("grey","red","blue")) + # assign color of points
			guides(color = guide_legend(override.aes = list(size = 4))) + # set legend key size
			labs(x = paste("Log-scaled Expression in",x),y = paste("Log-scaled Expression in",y)) +
			theme_classic() + 
			theme(axis.text = element_text(size = 20),axis.title = element_text(size = 20),
				legend.text = element_text(size = 15),legend.title = element_text(size = 15))
	return(p)
}

# 2. Set global parameters or from cmd --------------------
rds_file <- "a1.data/integrate_classified.rds"
treat <- args$treat
control <- args$control
outdir <- args$outdir
given_species <- parse.species(args$species)
speciesType <- args$speciesType
## construct outdir
if (!dir.exists(outdir)){
	dir.create(outdir,recursive = T) # mkdir -p
}
setwd(outdir)
## log 
write("================= 1. Parameters List ===================",stdout())
write(paste("Raw RDS:",rds_file),stdout())
write(paste("Comparison:",treat,"vs",control),stdout())
write(paste("Species:",paste(given_species,collapse = ",")),stdout())
write(paste("speciesType:",speciesType),stdout())
write(paste("Outdir:",outdir),stdout())

# 3. Pre-process scRNA data --------------
scrna <- readRDS(rds_file)
## grouping cell by Treat/Control
scrna$tmp <- "Passenger" # Passenger will not use while analysis
scrna$tmp[scrna[[which.type(scrna,treat)]] == treat] <- treat
if(grepl(pattern = "^Non-",x = control)){
	control_name <- gsub(pattern = "^Non-",replacement = "",x = control)
	scrna$tmp[scrna[[which.type(scrna,control_name)]] != control_name] <- control
}else if(grepl(pattern = "-Others$",x = control)){
	control_name <- gsub(pattern = "-Others$",replacement = "",x = control)
	scrna$tmp[scrna[[which.type(scrna,control_name)]] == control_name & scrna[[which.type(scrna,treat)]] != treat] <- control
}
Idents(scrna) <- scrna$tmp
DefaultAssay(scrna) <- "RNA"
write("================= 2. Pre-Process ===================",stdout())
print(addmargins(table(Idents(scrna),scrna$species))) # log

# 4. DEG analysis: output all genes ------------------------------------------
write("================= 3. DEA START ===================",stdout())
if (!dir.exists("DEGs_Table")){
	dir.create("DEGs_Table",recursive = T)
}
if (is.null(speciesType)){
	write("DEA inner Species one by one!",stdout())
	for (current_species in given_species){
		write(paste("\n",current_species),stdout())
		subdata <- subset(scrna,subset = species == current_species)
		print(table(subdata$species,Idents(subdata)))
		deg <- FindMarkers(subdata,ident.1 = treat,ident.2 = control,logfc.threshold = 0,min.pct = 0)
		deg <- data.frame(gene = rownames(deg),deg)
		write.table(x = deg,file = paste0("DEGs_Table/",current_species,"_nofilter.txt"),quote = F,row.names = F,col.names = T,sep = '\t')
	}
} else {
	write("DEA, Species as a Whole!",stdout())
	subdata <- subset(scrna,subset = species %in% given_species)
	print(table(Idents(subdata),subdata$species))
	deg <- FindMarkers(subdata,ident.1 = treat,ident.2 = control,logfc.threshold = 0,min.pct = 0)
	deg <- data.frame(gene = rownames(deg),deg)
	write.table(x = deg,file = paste0("DEGs_Table/",speciesType,"_nofilter.txt"),quote = F,row.names = F,col.names = T,sep = '\t')
}

# 5. Scatter plot for DEGs ----------------------------
write("================= 4. DEGs Scatter Plot START ===================",stdout())
if (!dir.exists("DEGs_Scatter")){
	dir.create("DEGs_Scatter",recursive = T)
}
cutoff <- list(p = 0.05,fc = exp(0.25),min.pct = 0.1)
write("Cutoff List:",stdout())
print(unlist(cutoff))

if (is.null(speciesType)){
	for (current_species in given_species){
		# extract average expression
		subdata <- subset(scrna,subset = species == current_species)
		exprs <- AverageExpression(subdata,assays = "RNA",slot = "data")[["RNA"]]
		exprs <- log1p(exprs) # when the expression were not log scaled
		exprs[["gene"]] <- rownames(exprs)
		write.table(x = exprs,file = paste0("DEGs_Scatter/",current_species,"_loged_exprs.txt"),col.names = T,row.names = F,sep = "\t",quote = F)
		# get DEGs
		degs <- read.table(file = paste0("DEGs_Table/",current_species,"_nofilter.txt"),sep = "\t",header = T,stringsAsFactors = F)
		degs[["Change"]] <- "NOT"
		degs$Change[degs$avg_logFC > log(cutoff$fc) & degs$p_val < cutoff$p & (degs$pct.1 > cutoff$min.pct | degs$pct.2 > cutoff$min.pct)] <- "UP"
		degs$Change[degs$avg_logFC < -log(cutoff$fc) & degs$p_val < cutoff$p & (degs$pct.1 > cutoff$min.pct | degs$pct.2 > cutoff$min.pct)] <- "DOWN"
		write.table(x = degs[degs$Change != "NOT",-grep("Change",colnames(degs))],
					file = paste0("DEGs_Table/",current_species,"_filtered.txt"),
					col.names = T,row.names = F,sep = "\t",quote = F)
		# construct data
		df <- merge(degs,exprs,by = "gene")
		df <- df[,c("gene","p_val","p_val_adj",treat,control,"Change","avg_logFC")]
		colnames(df) <- c("Gene","p_val","p_val_adj",treat,control,"Change","avg_logFC")
		# top 10 DEGs(each 10 of up/down)
		tmp <- df[df$Change!="NOT",]
		label_genes <- group_by(tmp,Change) %>%
			top_n(10,wt = abs(avg_logFC)) %>%
			pull(Gene)
		p <- scatter.plot(data = df,x = treat,y = control,labels = label_genes)
		ggsave(filename = paste0("DEGs_Scatter/",current_species,".pdf"),plot = p,device = "pdf",width = 8,height = 7)
	}
} else {
	subdata <- subset(scrna,subset = species %in% given_species)
	exprs <- AverageExpression(subdata,assays = "RNA",slot = "data")[["RNA"]]
	exprs <- log1p(exprs) # when the expression were not log scaled
	exprs[["gene"]] <- rownames(exprs)
	write.table(x = exprs,file = paste0("DEGs_Scatter/",speciesType,"_loged_exprs.txt"),col.names = T,row.names = F,sep = "\t",quote = F)
    # get DEGs
    degs <- read.table(file = paste0("DEGs_Table/",speciesType,"_nofilter.txt"),sep = "\t",header = T,stringsAsFactors = F)
	degs[["Change"]] <- "NOT"
	degs$Change[degs$avg_logFC > log(cutoff$fc) & degs$p_val < cutoff$p & (degs$pct.1 > cutoff$min.pct | degs$pct.2 > cutoff$min.pct)] <- "UP"
	degs$Change[degs$avg_logFC < -log(cutoff$fc) & degs$p_val < cutoff$p & (degs$pct.1 > cutoff$min.pct | degs$pct.2 > cutoff$min.pct)] <- "DOWN"
	write.table(x = degs[degs$Change != "NOT",-grep("Change",colnames(degs))],
				file = paste0("DEGs_Table/",speciesType,"_filtered.txt"),
				col.names = T,row.names = F,sep = "\t",quote = F)
	# construct data
	df <- merge(degs,exprs,by = "gene")
	df <- df[,c("gene","p_val","p_val_adj",treat,control,"Change","avg_logFC")]
	colnames(df) <- c("Gene","p_val","p_val_adj",treat,control,"Change","avg_logFC")
	# top 10 DEGs(each 10 of up/down)
	tmp <- df[df$Change!="NOT",]
	label_genes <- group_by(tmp,Change) %>%
		top_n(10,wt = abs(avg_logFC)) %>%
		pull(Gene)
	p <- scatter.plot(data = df,x = treat,y = control,labels = label_genes)
	ggsave(filename = paste0("DEGs_Scatter/",speciesType,".pdf"),plot = p,device = "pdf",width = 8,height = 7)
}
