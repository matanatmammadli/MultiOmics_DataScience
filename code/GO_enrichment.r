library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(gprofiler2)

data <- read.csv("/Users/matanatmammadli/Desktop/MultiOmic_Datascience-main/all_gene_annotation.csv", sep = "\t", header = TRUE)

genes_to_test <- data$ensembl_gene_id

##goResults <- gprofiler(query=genes_to_test, 
##                       organism = 'hsapiens', 
##                       src_filter = 'GO', 
##                       hier_filtering = 'moderate')
  
GO_results <- enrichGO(gene = genes_to_test,
                       OrgDb = "org.Hs.eg.db", 
                       keyType = "ENSEMBL", 
                       ont = "MF",
                       ##pvalueCutoff = 1,
                       ##qvalueCutoff = 1,
                       readable = TRUE
                       )
GO_results



fit <- plot(barplot(GO_results, showCategory = 15))
fit


GO_results2 <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "MF")
##as.data.frame(GO_results2)

fit2 <- plot(barplot(GO_results2, showCategory = 15))
fit2


GO_results3 <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "CC")
##as.data.frame(GO_results3)

fit3 <- plot(barplot(GO_results3, showCategory = 15))
fit3



## Biological theme comparison

library(clusterProfiler)
library(GOSemSim)
##BiocManager::install("enrichplot")

ensembl_gene_id <- data$ensembl_gene_id
ensembl_gene_id

gene_biotype <- data$gene_biotype
gene_biotype

hgnc_symbol <- data$hgnc_symbol
hgnc_symbol

xx <- compareCluster(ensembl_gene_id ~ gene_biotype, 
                     data=data, 
                     fun="enrichGO",
                     OrgDb = "org.Hs.eg.db", 
                     keyType = "ENSEMBL", 
                     ont = "MF")

clusterProfiler::dotplot(xx)

xx <- enrichplot::pairwise_termsim(xx)                     
p1 <- enrichplot::emapplot(xx)
p2 <- enrichplot::emapplot(xx, pie.params = list(legend_n = 2)) 
p3 <- enrichplot::emapplot(xx, pie.params = list(pie = "count"))
p4 <- enrichplot::emapplot(xx, pie.params = list(pie = "count"), 
               cex.params = list(category_node = 1.5), 
               layout.params = list(layout = "kk"))
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])


##xx2 <- compareCluster(ensembl_gene_id ~ gene_biotype + hgnc_symbol, 
##                     data=data, 
##                     fun="enrichGO",
##                     OrgDb = "org.Hs.eg.db", 
##                     keyType = "ENSEMBL", 
##                     ont = "MF")

##head(as.data.frame(xx2))

##clusterProfiler::dotplot(xx2)
##clusterProfiler::dotplot(xx2, x=~gene_biotype) + 
##    ggplot2::facet_grid(~hgnc_symbol)

##xx2 <- pairwise_termsim(xx2)                     
##p1 <- emapplot(xx2)
##p2 <- emapplot(xx2, pie.params = list(legend_n = 2)) 
##p3 <- emapplot(xx2, pie.params = list(pie = "count"))
##p4 <- emapplot(xx2, pie.params = list(pie = "count"), 
##               cex.params = list(category_node = 1.5), 
##               layout.params = list(layout = "kk"))
##cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

