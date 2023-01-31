# Single-cell-transcriptomic-atlas-across-species
This repository contains materials for the paper [**Single-cell sequencing reveals the evolution of immune molecules across multiple vertebrate species**](), by Anjun Jiao, Cangang Zhang et al.

## Scripts

- *s1.Integrated-scRNA.R*: Perform analysis of single-cell transcriptome data for all species in the paper, including cell clustering, cell type annotations, and a heatmap of the top30 key genes for each cell typ.  It generates the following figs: 1b-e, S1b, S2 and S3a
- *s2.basic_plot.R*: The basic plot (TSNE, Doheatmap, Fraction) for a given cluster. it generates figs. 2b, 2h, 3b, 4f, S1c,  
- *s3.DEA.R*:  R script for finding genes that differ between two cell types (cell type is the input parameter). It generates data for a heatmap of different genes
- *s4.enrich.R*: Gene Enrichment for DEGs, Including go, KEGG and C7 enrichment analysis. 



## Plots

The scripts and data in this repository can generate almost all figures in the main text and supplement, as follows:

- Main text: fig.1 (b-e); fig.2 (a-d,f); fig. 3(a-g); fig. 4(a,b,d)
- Supplementary Information: figs. S2, S6-S8

## Data

The following useful data structures are included in the repository:

- *cells.txt*: The ID of the filtered cells
