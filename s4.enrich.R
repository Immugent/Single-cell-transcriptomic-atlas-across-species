# --------- Prog Info ---------
# R version 3.6.3
# Author: Songming Liu
# E-Maile: smliu-bio@qq.com
# Date: 2020/11/19
# Purpose: Gene Enrichment (ORA) for DEGs

# 0. Prepare -------------------
library(optparse)
# construct options
option_list <- list(make_option(opt_str = c("-i","--indir"),action = "store",dest = "indir",metavar = "INDIR",help = "Dir contains DEGs files"),
					make_option(opt_str = c("-s","--suffix"),action = "store",dest = "suffix",help = "Suffix to recognize DEGs files [.txt, _filtered.txt]"),
					make_option(opt_str = c("-o","--outdir"),action = "store",dest = "outdir",metavar = "OUTDIR"),
					make_option(opt_str = "--seven",action = "store_true",dest = "seven",default = FALSE,help = "Whether to enrich for Species Specific DEGs (will be divided to 7 types)"),
					make_option(opt_str = "--single",action = "store_true",dest = "single",default = FALSE,help = "Whether to visualize one by one, default FALSE"),
					make_option(opt_str = "--C7",action = "store_true",dest = "c7",default = FALSE,help = "Whether to do MSigDB C7 enrichment, defalut FALSE. only works when set --single"))
parse <- OptionParser(option_list = option_list) # instance
args <- parse_args(object = parse) # parse cmd
if (is.null(args$indir) | is.null(args$suffix) | is.null(args$outdir)){
	stop("-i/-s/-o are required and expected one argument, please check!")
}
if (!args$single & args$c7){
	write("WARNING: --C7 was set but ignored since not set --single",stdout())
	args$c7 <- FALSE
}
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(ggplot2)
library(msigdbr)

# 1. Function Definition ------------------
check.dir <- function(dir){
# Check dir and create if not exists
# Args:
#	dir: dircetory
	if (!dir.exists(dir)){
		dir.create(dir,recursive = T) # mkdir -p
	}
}

check.empty <- function(file,header){
# Check file is empty or not
# Args:
#	file: file to check
#	header: T/F
# Returns:
#	NULL when error occurs while reading, data frame while successfully read
	df <- try(read.table(file = file,header = header,stringsAsFactors = F,sep = "\t"),silent = T)
	if ("try-error" %in% class(df)){
		simpleWarning(df)
		return(NULL)
	} else {
		return(df)
	}
}

parse.seven <- function(data,outfile){
# Parse Union Results from LogfcHeatmap
# Args:
#	data: data frame of results
#	outfile: name of outfile
# Returns:
#	a list of 7 types
	ls <- vector(mode = "list",length = 7)
	names(ls) <- paste0("Type",1:7)
	for (i in 1:nrow(data)){
		switch (as.character(data[i,"Type1"]),
			"127" = {ls$Type1 <- c(ls$Type1,rownames(data)[i])},
			"126" = {ls$Type2 <- c(ls$Type2,rownames(data)[i])},
			"124" = {ls$Type3 <- c(ls$Type3,rownames(data)[i])},
			"120" = {ls$Type4 <- c(ls$Type4,rownames(data)[i])},
			"112" = {ls$Type5 <- c(ls$Type5,rownames(data)[i])},
			"96" = {ls$Type6 <- c(ls$Type6,rownames(data)[i])},
			"64" = {ls$Type7 <- c(ls$Type7,rownames(data)[i])}
		)
	}
	df <- as.data.frame(do.call(cbind,lapply(ls,`length<-`,max(lengths(ls))))) # turn list to data.frame
	write.table(file = outfile,x = df,sep = "\t",quote = F,col.names = T,row.names = F,na = "")
	return(ls)
}

find.enrich <- function(symbol,pcut = 1,qcut = 1,out_prefix,ont = "ALL"){
# GO/KEGG enrich for given genes
# Args:
#	symbol: gene symbols
#	pcut: cutoff of pvalue, p.adj
#	qcut: cutoff of qvalue
#	out_prefix: prefix of outfile
#	ont: GO ontology to enrich
	# get ENTREZID, some may fail
	write(paste("\t1. Before converting gene number(SYMBOL):",length(symbol)),stdout())
	geneid <- bitr(symbol,fromType = 'SYMBOL',toType = "ENTREZID",OrgDb = org.Hs.eg.db)
	write(paste("\t2. After converting gene number(ENTREZID):",length(unique(geneid$ENTREZID))),stdout())
	write(paste("\t\tThe following gene(s) failed converting:",paste(unique(setdiff(symbol,geneid$SYMBOL)),collapse = ", ")),stdout())
	if(length(unique(geneid$ENTREZID)) < 3){
		write("\tToo few eligible genes (<3), skip GO/KEGG enrichment!",stdout())
		return(NULL) # stop function
	}
	# GO enrich
	go <- enrichGO(gene = geneid$ENTREZID,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',ont = ont,pvalueCutoff = pcut,qvalueCutoff = qcut,readable = TRUE)
	if (is.null(go)){
		write("\tWARNING: GO finished but NO RESULT!",stdout())
	}else{
		write.table(go@result,file = paste0(out_prefix,"_go.txt"),,sep = '\t',quote = F,col.names = T,row.names = F)
	}
	# KEGG enrich
	kegg <- enrichKEGG(gene = geneid$ENTREZID,organism = "hsa",keyType = "kegg",pvalueCutoff = pcut,qvalueCutoff = qcut)
	if (is.null(kegg)){
		write("\tWARNING: KEGG finished but NO RESULT!",stdout())
	}else{
		kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")# set readable for gene
		write.table(kegg@result,file = paste0(out_prefix,"_kegg.txt"),,sep = '\t',quote = F,col.names = T,row.names = F)
	}
	write("\t3. Enrichment successfully finished!",stdout())
}

vlz.enrich <- function(data,out_prefix,title = "",show.num,col = 'pvalue',single = F){
# Visualize GO/KEGG results
# Args: 
#	data: Enrichment results for visualization
#	out_prefix: prefix of outfile
#	title: title of grpah
#	show.num: number of pathway with minimum pvalue to show in graph
#	col: column to color by
#	single: whether to plot one by one
	plot_params <- list(y_len = 40,base_size = 22,width = 11,height = 7)
    data$GeneRatio <- sapply(as.character(data$GeneRatio), function(x){return(eval(parse(text = x)))})
    if (!single){
		data <- data %>% group_by(Group)
		plot_params <- list(y_len = 50,base_size = 25,width = 16,height = 12)
	}
	pathway_id <- data %>%
		top_n(show.num,wt = -log10(pvalue)) %>%
		top_n(show.num,wt = GeneRatio) %>%
		pull(ID)
    df <- as.data.frame(data[data$ID %in% pathway_id,],stringsAsFactors = F)
    # select data using for plotting
    if (col == 'p.adjust'){
        df[['color']] <- -log10(df$p.adjust)
        colorname <- "-log10(p.adjust)"
    } else {
        df[['color']] <- -log10(df$pvalue)
        colorname <- "-log10(pvalue)"
    }
    # max length of y lab
    string_cut <- function(x,len){
        if(nchar(x) > len){
            return(paste(str_sub(x,1,len),"..."))
        }else {
            return(x)
        }
    }
	# plot
	df$Description <- unlist(sapply(df$Description,string_cut,plot_params$y_len)) # max length
	if (single){
		p <- ggplot(data = df,aes(GeneRatio,Description))
	} else {
		p <- ggplot(data = df,aes(Group,Description))
	}
	p <- p + geom_point(aes(size = Count,color = color)) +
		labs(title = title,x = '',y = '',color = colorname) +
		scale_color_gradient(low = "blue",high = "red") +
		scale_size(range = c(3,10)) + # size range of "Count"
		guides(size = guide_legend(reverse = T,nrow = 5,order = 1),color = guide_colorbar(order = 2)) +
		theme_bw(base_size = plot_params$base_size) +
		theme(axis.text.x.bottom = element_text(angle = 25,hjust = 0.9,vjust = 0.9))
    ggsave(filename = paste0(out_prefix,".pdf"),plot = p,device = "pdf",width = plot_params$width,height = plot_params$height)
    ggsave(filename = paste0(out_prefix,".png"),plot = p,device = "png",width = plot_params$width,height = plot_params$height)
}

find.c7 <- function(symbol,outfile){
# C7 Enrichment
# Args:
#	symbol: gene symbol
#	outfile: outfile name
# Returns:
#	NULL when too less gene
	m_t2g <- msigdbr(species = "Homo sapiens",category = "C7") %>% 
		select(gs_name,gene_symbol) # use gene symbol
	if (length(symbol) < 3){
		write("\tToo few eligible genes (<3), skip C7 enrichment!",stdout())
		return(NULL)
	}
	res <- enricher(gene = symbol,TERM2GENE = m_t2g,pvalueCutoff = 1,qvalueCutoff = 1)
	if(is.null(res)){
		write("\tWARNING: MSigDB C7 finished but NO RESULT!",stdout())
	}else{
		write.table(x = res@result,file = outfile,col.names = T,row.names = F,quote = F,sep = "\t")
		write("\tMSigDB C7 successfully finished!",stdout())
	}
}

# 2. Set global parameters or from cmd --------------------
indir <- args$indir
suffix <- args$suffix
outdir <- args$outdir
seven <- args$seven
single <- args$single
c7 <- args$c7
## construct outdir
check.dir(outdir)
setwd(outdir)
## log 
write("================= 1. Parameters List ===================",stdout())
write(paste("Indir:",indir),stdout())
write(paste0("Filenames: *",suffix),stdout())
write(paste("Outdir:",outdir),stdout())
write(paste("Seven Types Or Not:",seven),stdout())
write(paste("Single Visualization:",single),stdout())
write(paste("C7 Enrichment:",c7),stdout())

# 3. Read DEGs ------------------
write("=================== 2. Reading DEGs ======================",stdout())
files <- list.files(path = indir,pattern = paste0(suffix,"$"),full.names = T)
degs <- list(up = list(),down = list())
for (infile in files){
	prefix <- gsub(pattern = paste0(".+/|logfc_|",suffix,"$"),replacement = "",x = infile)
	if (suffix == "_filtered.txt"){
	# for DEA results
		df <- check.empty(file = infile,header = T)
		if (is.null(df)){
			next
		}
		degs$up[[prefix]] <- df$gene[df$avg_logFC > 0]
		degs$down[[prefix]] <- df$gene[df$avg_logFC < 0]
	} else if (seven){
	# for specific results and divided into 7 types
		df <- check.empty(file = infile,header = T)
		if (is.null(df)){
			next
		}
		degs[[prefix]] <- parse.seven(data = df,outfile = paste0(prefix,".txt"))
	} else if (single){
	# for Shared DEGs
		df <- check.empty(file = infile,header = F)
		if (is.null(df)){
			next
		}
		degs[[prefix]][["Shared"]] <- df[[1]]
	} else {
	# for Specific DEGs
		df <- check.empty(file = infile,header = T)
		if (is.null(df)){
			next
		}
		degs[[prefix]] <- lapply(X = as.list(df),FUN = function(x){x[x != "" & !is.na(x)]})
	}
}
write("Up-Regulatory Genes:",stdout())
print(lengths(degs$up))
write("Down-Regulatory Genes:",stdout())
print(lengths(degs$down))

# 4. Coustruct Dir ---------------
write("=================== 3. Constructing Directory =====================",stdout())
check.dir(dir = "All_Results")
check.dir(dir = "Visualization")
if (c7){
	check.dir(dir = "C7_Results")
}

# 5. Enrichment ------------------
write("=================== 4. Enrichment ======================",stdout())
for (reg_type in names(degs)){
	for (group in names(degs[[reg_type]])){
		write(paste(reg_type,group,"START:"),stdout())
		genes <- degs[[reg_type]][[group]]
		prefix <- paste0("All_Results/",reg_type,"_",group)
		if (length(genes) == 0){
			write("\tNo genes, skip!",stdout())
			next
		}
		find.enrich(symbol = genes,pcut = 1,qcut = 1,out_prefix = prefix,ont = "ALL")
		if (c7){
			write(paste(reg_type,group,"C7 START:"),stdout())
			find.c7(symbol = genes,outfile = paste0("C7_Results/",reg_type,"_",group,"_C7.txt"))
		}
	} 
}

# 6. Visualization (GO/KEGG) ---------------------------------
write("=================== 5. Visualization (GO/KEGG) ======================",stdout())
files <- list.files(path = "All_Results/",pattern = ".txt$",full.names = T)
## read file
df <- list(up = list(go = rbind(),kegg = rbind()),down = list(go = rbind(),kegg = rbind()))
for (resfile in files){
	prefix <- gsub(pattern = paste0(".+/|.txt$"),replacement = "",x = resfile)
	info <- unlist(strsplit(x = prefix,split = "_")) # prefix: up_Human_go
	reg_type <- info[1]
	group <- info[2]
	enrich_type <- info[3]
	tmp <- read.table(file = resfile,header = T,sep = "\t",stringsAsFactors = F,quote = "")
	tmp[["Group"]] <- group
	df[[reg_type]][[enrich_type]] <- rbind(df[[reg_type]][[enrich_type]],tmp)
}
## visualization
for (reg_type in names(df)){
	for (enrich_type in names(df[[reg_type]])){
		tmp <- df[[reg_type]][[enrich_type]]
		if (single){
		# when single
			for (group in unique(tmp$Group)){
				prefix <- paste0("Visualization/",reg_type,"_",group,"_",enrich_type)
				if (enrich_type == "go"){
					vlz.enrich(data = tmp[tmp$Group == group & tmp$ONTOLOGY == "BP", ],out_prefix = paste0(prefix,"_BP"),show.num = 10,col = 'pvalue',single = T)
					vlz.enrich(data = tmp[tmp$Group == group & tmp$ONTOLOGY == "MF", ],out_prefix = paste0(prefix,"_MF"),show.num = 10,col = 'pvalue',single = T)
					vlz.enrich(data = tmp[tmp$Group == group & tmp$ONTOLOGY == "CC", ],out_prefix = paste0(prefix,"_CC"),show.num = 10,col = 'pvalue',single = T)
				} else {
					vlz.enrich(data = tmp[tmp$Group == group, ],out_prefix = prefix,show.num = 10,col = 'pvalue',single = T)
				}
				write(paste(reg_type,group,enrich_type,"done (single)"),stdout())
			}
		} else {
		# when not single
			if (seven){
				tmp$Group <- factor(tmp$Group,levels = paste0("Type",1:7))
			} else {
				tmp$Group <- factor(tmp$Group,levels = c("Fish","Frog","Mouse","Rat","Pig","Monkey","Human"))
			}
			prefix <- paste0("Visualization/",reg_type,"_",enrich_type)
			if (enrich_type == "go"){
				vlz.enrich(data = tmp[tmp$ONTOLOGY == "BP", ],out_prefix = paste0(prefix,"_BP"),show.num = 5,col = 'pvalue',single = F)
				vlz.enrich(data = tmp[tmp$ONTOLOGY == "MF", ],out_prefix = paste0(prefix,"_MF"),show.num = 5,col = 'pvalue',single = F)
				vlz.enrich(data = tmp[tmp$ONTOLOGY == "CC", ],out_prefix = paste0(prefix,"_CC"),show.num = 5,col = 'pvalue',single = F)
			} else {
				vlz.enrich(data = tmp,out_prefix = prefix,show.num = 5,col = 'pvalue',single = F)
			}
			write(paste(reg_type,enrich_type,"done (multi)"),stdout())
		}
	}
}

