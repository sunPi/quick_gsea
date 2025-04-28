# quick_gsea

This script calculates pathway enrichment scores based on a list of differentially expressed genes. It works in two stesp:

1. Raw gene counts are normalized and a list of differentially expressed genes is produced using the DESeq2 algorithm
2. All differentially expressed genes are ranked and fed into the fGSEA algorithm which calculates and plots enriched pathways

!!! It is very important to set the correct reference for your experiment. Since this is a comparative (conrast) study, there will be two states in which to compare (i.e. reference vs. treatment). Results are always computed in treatment, meaning set the your control as the reference.

Example:

We are looking for enriched patwhays in treated explants. In this case, we have two classes 'control' and 'treated'. To see enriched genes and pathways in 'treated', set argument ref <- 'control' to ensure it works properly. These names are arbitrary and can reflect any kind of experiment, this is set by the user in the file 'colData.csv'. 
