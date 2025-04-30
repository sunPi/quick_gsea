# Load Functions
handleRequirements <- function(pkgs) {
  ipkgs <- sapply(pkgs, function(...) require(..., character.only = TRUE))
  if (any(!ipkgs)) {
    BiocManager::install(pkgs[!ipkgs])
    install.packages(pkgs[!ipkgs])
  }
  else {
    message("\n\nCool! your machine has everything is needed.\n\n")
  }
  print("Loading required packages...")
  pacman::p_load(pkgs, install = TRUE, character.only = TRUE)
  return(pacman::p_loaded())
}
prepareCounts      <- function(cts, verbose){
  if(verbose){
    print("Preparing RNA-seq counts for DEG and GSEA...")
  }
  
  # # Remove duplicated values
  # dup.marker <- all(duplicated(cts[,1]))
  # 
  # if(dup.marker){
  #   if(verbose){
  #     print("Found duplicates, removing...")
  #   }
  
  # Move gene names to row names
  if(verbose){
    print("Moving first column to row names...")
  }
  cts <- column_to_rownames(cts, var = "gene_id")
  cts <- cts[,order(colnames(cts))]
  
  # } else{
  #   if(verbose){
  #     print("Found no duplicates!")
  #   }
  # }
  # Round values just in case
  # if(verbose){
  #   print("Rounding the integers...")
  # }
  # cts <- round(cts)
  
  return(cts)
}
getColData         <- function(exp.data, cts, reference){
  
  cts <- cts[rowSums(cts == 0) == 0, ] # Remove rows with 0 values
  if(length(exp.data[[1]]) != length(colnames(cts))){
    print("Detecting a difference in sample size between count data and experiment design data, conforming...")
    exp.data <- dplyr::filter(exp.data, exp.data[[1]] %in% colnames(cts))
  }
  
  
  
  coldata <- exp.data[order(exp.data[[1]]), ]
  names(coldata)[2] <- "condition"
  # coldata$condition <- relevel(factor(coldata$condition), ref = reference)
  
  return(coldata)
  
}
RankGenes          <- function(rnaseq.counts, experiment.data, reference, verbose){
  
  # Prepare RNA-seq expression data
  cts <- prepareCounts(rnaseq.counts, verbose)
  
  # Create a 'coldata' object from the experiment design
  coldata <- getColData(exp.data = experiment.data, cts = cts, reference = reference)
  
  # Ensure the order of samples in coldata matches the counts data
  all(colnames(cts) == coldata[[1]])
  
  
  
  # Create a DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
  
  rowsum.threshold <- 150 # user chosen
  # fdr.threshold <- 0.1 # user chosen
  rs <- rowSums(counts(dds))
  dds <- dds[ rs > rowsum.threshold ,]
  
  if(verbose){
    print("Running DESeq2...")
  }
  
  # Step 2: Relevel the condition (for reference)
  dds$condition <- relevel(dds$condition, ref = reference)
  
  # Perform the DESeq analysis
  dds <- DESeq(dds)
  
  # Extract the results
  # Contrast is always built such as - "name of the column in coldata" "test" "reference"
  
  test <- levels(dds$condition)[levels(dds$condition) != reference]
  ctr <- c("condition", test, reference)
  
  if(verbose){
    print(paste0("Extracting list of DEG in ", test, "..."))
  }
  
  gene.list <- as.data.frame(results(dds))
  # names(gene.list)[1] <- "gene_symbol"
  
  if(all(grepl(pattern = "^ENSG*", rownames(gene.list)))){
    gene.list$gene_symbol <- mapIds(org.Hs.eg.db, keys = row.names(gene.list), keytype = "ENSEMBL", column = "SYMBOL")
    
  } else{
    gene.list$gene_symbol <- rownames(gene.list)
  }
  # gene.list <- na.omit(gene.list)
  
  # Filter out the significant genes based on p.value and order/sort them
  sig <- dplyr::filter(gene.list, gene.list$padj < 0.05)
  sig <- sig[order(sig$padj),]
  sig <- sig[order(sig$log2FoldChange, decreasing = T), ]
  sig <- sig %>% as.data.frame() %>% arrange(desc(log2FoldChange), desc(padj))
  top.10.genes <- head(sig[1:10, ])
  bot.10.genes <- tail(sig[1:10, ])
  
  if(nrow(sig) > 0){
    any.sig <- T
    print(paste0("DESeq2 analysis found ", nrow(sig), " significant genes!"))
    
  } else{
    any.sig <- F
    print("DESeq2 found 0 significant genes!")
  }
  
  p <- EnhancedVolcano(gene.list,
                       lab = gene.list$gene_symbol,
                       x = 'log2FoldChange',
                       y = 'pvalue') + 
    theme_prism() +
    labs(title = NULL, subtitle = NULL)
  
  
  return(list(gene.list   = gene.list,
              significant = sig,
              top.10      = top.10.genes,
              bot.10      = bot.10.genes,
              design      = list(ref = reference,
                                 test = test),
              any.sig     = any.sig,
              volcano.p   = p)
  )
}
RunfGSEA           <- function(gene.list, pathways, design, verbose){
  res <- gene.list %>% 
    dplyr::select(gene_symbol, stat) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(gene_symbol) %>% 
    summarize(stat=mean(stat))
  
  ranks <- deframe(res)
  head(ranks, 20)
  
  # Load the pathways into a named list
  pathways <- gmtPathways(pathways)
  # names(pathways) <- gsub("HALLMARK_", "", names(pathways))
  
  fgseaRes <- fgsea(pathways=pathways, stats=ranks)
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy <- fgseaResTidy[fgseaResTidy$padj < 0.1,]
  fgseaResTidy$Direction <- recode(factor(fgseaResTidy$NES > 0), "TRUE" = "Upregulated", "FALSE" = "Downregulated")
  
  # Convert leading edge list into vectors for convenience
  fgseaResTidy <- fgseaResTidy %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ", "))) %>%
    mutate(leadingEdge = as.character(leadingEdge)) # Force character
  
  # Create custom pallete in accordance with our lab's theme
  custom_colors <- c("Upregulated" = "#077f97", "Downregulated" = "#400257")
  
  # title <- design$test
  
  # Create a sideways barplot of pathway enrichments
  p <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=Direction)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score", title = design$test) + # title=title
    theme(text = element_text(size = 8), 
          title = element_text(size = 18), 
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = custom_colors) +
    theme_prism()
  
  return(list(enrichments = fgseaResTidy,
              plot        = p)) 
}

# Load Libraries
pkgs <- c("DESeq2", "tidyverse", "dplyr", "fgsea", "ggprism", 
          "EnhancedVolcano", "org.Hs.eg.db", "AnnotationDbi", "stringr", "writexl")
suppressMessages(handleRequirements(pkgs))

#----- Main -----
# Load Arguments
rnaseq.counts   <- read.csv(file.choose())
experiment.data <- read.csv(file.choose())
pathways        <- file.choose()
outfolder       <- "./" # Specify path where to save the results - "./" means it will use current directory as output directory

# Load additional arguments 
p.value         <- 0.05 
verbose         <- T
serialise       <- T
reference       <- "..." # Set your reference (e.g. if you want results IN responders, set reference = Non-Responders)
dname           <- pathways %>% basename()
outfolder       <- here(outfolder, dname)
dir.create(outfolder, F, T)

# Run differential gene expression
print("Ranking Genes...")
DEG             <- RankGenes(rnaseq.counts, experiment.data, reference, verbose) # DESeq2 is used to generate diff. exp. genes

# Run gene set enrichment analysis
print("Running fGSEA...")
fGSEA           <- RunfGSEA(DEG$gene.list, pathways, DEG$design, verbose) # fGSEA algorithm computes enriched pathways

if(serialise){ # Run this block of code to save all your results
  print(paste("Serialising results to...", outfolder))
  write_xlsx(list(significant = DEG$significant,
                  top10       = DEG$top.10,
                  bot10       = DEG$bot.10),
                      here(outfolder,paste0(DEG$design$test, '_DEG.xlsx'))
  )

  png(here(outfolder,paste0(DEG$design$test, '_volcano.png')),
      width = 2400, height = 1800, res = 300)
  print(DEG$volcano.p)
  dev.off()

 write_xlsx(fGSEA$enrichments,
            here(outfolder,paste0(DEG$design$test, '_pathways.xlsx'))
  )

  png(here(outfolder,paste0(DEG$design$test, '_patwhays.png')),
      width = 3200, height = 1800, res = 300)
  print(fGSEA$plot)
  dev.off()
}


