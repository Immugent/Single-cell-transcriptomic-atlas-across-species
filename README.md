# Single-cell-transcriptomic-atlas-across-species
This code is used for annotating and analyzing single-cell sequencing data derived from the spleen across 7 species.
## Introduction: 
Both innate and adaptive immune system undergo evolution from low to high vertebrates.Due to the limitation of conventional approaches in identifying broader spectrum of immune cells and
molecules from various vertebrates, it remains unclear how immune molecules evolve among vertebrates.
## Objectives: 
Here, we utilized carry out comparative transcriptome analysis in various immune cellsacross seven vertebrate species.
## Methods: 
Single-cell RNA sequencing (scRNA-seq).
## Results: 
We uncovered both conserved and species-specific profiling of gene expression in innate andadaptive immunity. Macrophages exhibited highly-diversified genes and developed sophisticated molecular
signaling networks along with evolution, indicating effective and versatile functions in higher species.In contrast, B cells conservatively evolved with less differentially-expressed genes in analyzed
species. Interestingly, T cells represented a dominant immune cell populations in all species and uniqueT cell populations were identified in zebrafish and pig. We also revealed compensatory TCR cascade components utilized by different species. Inter-species comparison of core gene programs demonstratedmouse species has the highest similarity in immune transcriptomes to human.
#  Conclusions: 
Therefore, our comparative study reveals gene transcription characteristics across multiplevertebrate species during the evolution of immune system, providing insights for species-specific immunity
as well as the translation of animal studies to human physiology and disease.

# Scripts

- *s1.Integrated-scRNA.R*: Perform analysis of single-cell transcriptome data for all species in the paper, including cell clustering, cell type annotations, and a heatmap of the top30 key genes for each cell typ.  It generates the following figs: 1b-e, S1b, S2 and S3a
- *s2.basic_plot.R*: The basic plot (TSNE, Doheatmap, Fraction) for a given cluster. it generates figs. 2b, 2h, 3b, 4f, S1c,  
- *s3.DEA.R*:  R script for finding genes that differ between two cell types (cell type is the input parameter). It generates data for a heatmap of different genes
- *s4.enrich.R*: Gene Enrichment for DEGs, Including go, KEGG and C7 enrichment analysis. 
- *s5.Mapping.R*: Mapping other human spleen scRNAseq datasets to our data, and calculate similarity between them.


# Plots

The scripts and data in this repository can generate almost all figures in the main text and supplement, as follows:

- Main text: fig.1 (b-e); fig.2 (a-d,f); fig. 3(a-g); fig. 4(a,b,d)
- Supplementary Information: figs. S2, S6-S8

# Data
All sequencing data generated in this paper have been deposited at Gene Expression Omnibus, under GSE186158 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186158). Other human spleen datasets used in benchmarking were taken from GSE159929 and PRJEB31843.

The following useful data structures are included in the repository:

- *cells.txt*: The ID of the filtered cells

# REFERENCE 
Please cite this paper if you use our data or code.

Jiao A, Zhang C, Wang X, Sun L, Liu H, Su Y, Lei L, Li W, Ding R, Ding C, Dou M, Tian P, Sun C, Yang X, Zhang L, Zhang B. Single-cell sequencing reveals the evolution of immune molecules across multiple vertebrate species. J Adv Res. 2024 Jan;55:73-87. doi: 10.1016/j.jare.2023.02.017. Epub 2023 Mar 4. PMID: 36871615; PMCID: PMC10770119.

# Contact me
Please feel free to contact us at the following email addresses: cg_zhang2021@163.com.
