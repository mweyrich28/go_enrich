# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("fgsea")
library(fgsea)
library(data.table)
library(msigdbr)

genes <- fread("genes.tsv")
#      id      fc signif
#           <char>   <num> <lgcl>
# 1: DNAJC25-GNG10 -1.3420  FALSE
# 2:      IGKV2-28 -2.3961  FALSE
# 3:      IGHV3-64  0.6092  FALSE
# 4:     IGKV2D-30 -0.2356  FALSE
# 5:  LOC100509620  0.6891  FALSE
# 6:       PPIAL4E -3.2035  FALSE


genes <- genes[order(genes$fc, decreasing = TRUE),]
geneList <- genes$fc
names(geneList) <- genes$id

gene_sets <- msigdbr(species = "Homo sapiens", category = "C5")
gene_sets  <- as.data.table(gene_sets)
pathways <- split(gene_sets$gene_symbol, gene_sets$gs_name)
 

fgsea_results <- fgsea(pathways = pathways, stats = geneList, minSize = 50, maxSize = 500)

go_id_map <- unique(gene_sets[, .(gs_name, gs_exact_source)])

fgsea_results <- as.data.table(fgsea_results)

fgsea_results <- merge(fgsea_results, 
                      go_id_map, 
                      by.x = "pathway", 
                      by.y = "gs_name", 
                      all.x = TRUE)

filtered_results <- fgsea_results[
  grepl("^GOBP_", pathway) & 
  padj < 0.05
][order(padj)]

head(filtered_results)
filtered_results <- filtered_results[, .(gs_exact_source, size, padj, ES)]
filtered_results <- filtered_results[size >= 50 & size <= 500]
write.table(filtered_results, "gsea_results.tsv", row.names = FALSE, sep="\t")


