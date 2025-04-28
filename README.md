# quick_gsea

This script calculates pathway enrichment scores based on a list of differentially expressed genes. It works in two stesp:

1. Raw gene counts are normalized and a list of differentially expressed genes is produced using the DESeq2 algorithm
2. All differentially expressed genes are ranked and fed into the fGSEA algorithm which calculates and plots enriched pathways

!!! It is very important to set the correct reference for your experiment. Since this is a comparative (conrast) study, there will be two states in which to compare (i.e. reference vs. treatment). Results are always computed in treatment, meaning set the your control as the reference.

Example:

We are looking for enriched patwhays in treated explants. In this case, we have two classes 'control' and 'treated'. To see enriched genes and pathways in 'treated', set argument ref <- 'control' to ensure it works properly. These names are arbitrary and can reflect any kind of experiment, this is set by the user in the file 'colData.csv'. 

```
How to use:

Simply open quick_gsea.R with Rstudio and run the code line by line.


rnaseq.counts   <- read.csv(file.choose()) # Opens a dialogue to choose your RNA-seq expression data (counts)
experiment.data <- read.csv(file.choose()) # Opens a dialogue to choose your mapping table (colData)
pathways        <- file.choose()           # Opens a dialogue to choose which reference gene set collection you wish to compute enrichments for
outfolder       <- "./"                    # Folder path in which to save your results


 rnaseq.counts file:

This should be a csv file that contains your gene symbols in the first column, and rows as expression data (raw count data). The first column has to be named 'gene_id'! All the genes have to be unique (quick_gsea.R does not allow duplicated genes).

experiment.data:

This should be a simple two column csv file, where first column are your sample names (these must be the same names as sample names in rnaseq.counts.csv file) and the second column is your experiment condition (i.e. control vs. treatment).

```

There are two example files for rnaseq.counts and experiment.data that can be used as example file formats. A 2024 collection of MSigDB molecular sets can be found under the 'msigdb' folder, howver its best to always check and download the latest molecular sets from https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp.
