##### Cluster Pathway Enrichment Script v 1.4 #####
## Serdar Turkarslan                            ##
## Institute for Systems Biology, 2020          ##
## Last update: 08/20/2020                      ##
## Input set of genes with gene symbol either   ##
## as vector or data frame with Ensembl/Symbol  ##
## Enrichments are performed against different  ##
## MSigSdb signatures                           ##
#################################################

library(pathfindR)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)


cluster_enrichment <- function(my.genes){
 

  ## Input gene list could be data frame or vector of genes
  if(is.vector(my.genes)){
    genes.input <- my.genes
  }
  if(is.data.frame(my.genes)){
    genes.input <- my.genes$SYMBOL
  }

  ## Hallmarks Enrichment (Uses MSIGDB Hallmark signatures)
  h_df = msigdbr(category = "H")
  h_t2g = h_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  edo.h <- enricher(gene = genes.input, TERM2GENE = h_t2g)
  if(!is.null(edo.h)) {
     if(nrow(edo.h) > 0) {
        edo.h <- edo.h %>%
        mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))
     }
  }
  
  
  ## Immunologic Signatures (MSIGDb)
  c7_df = msigdbr(category = "C7")
  c7_t2g = c7_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  edo.c7 <- enricher(gene = genes.input, TERM2GENE = c7_t2g)
  if(!is.null(edo.c7)) {
     if(nrow(edo.c7) > 0) {
        edo.c7 <- edo.c7 %>%
          mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))
     }
  }
  
  ## Oncogenic Signatures (MSIGDb)
  c6_df = msigdbr(category = "C6")
  c6_t2g = c6_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  edo.c6 <- enricher(gene = genes.input, TERM2GENE = c6_t2g)
  if(!is.null(edo.c6)) {
     if(nrow(edo.c6) > 0) {
        edo.c6 <- edo.c6 %>%
          mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))
     }
  }

  ## Curated datasets (MSIGDb)
  c2_df = msigdbr(category = "C2")
  c2_t2g = c2_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  edo.c2 <- enricher(gene = genes.input, TERM2GENE = c2_t2g)
  if(!is.null(edo.c2)) {
     if(nrow(edo.c2) > 0) {
        edo.c2 <- edo.c2 %>%
          mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))
     }
  }  

  ## KEGG  enrichments (MSIGDb)
  kegg_df = msigdbr(category = "C2", subcategory = "CP:KEGG")
  kegg_t2g = kegg_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  edo.kegg <- enricher(gene = genes.input, TERM2GENE = kegg_t2g)
  if(!is.null(edo.kegg)) {
     if(nrow(edo.kegg)) {
         edo.kegg <- edo.kegg %>%
         mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))
     }
  }
  

  ## Curated datasets (MSIGDb)
  cm_df = msigdbr(category = "C4", subcategory = "CM")
  cm_t2g = cm_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  edo.cm <- enricher(gene = genes.input, TERM2GENE = cm_t2g)
  if(!is.null(edo.cm)) {
     if(nrow(edo.cm) > 0) {
        edo.cm <- edo.cm %>%
          mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))
     }
  }

  ## GO BP (MSIGDb)
  gobp_df = msigdbr(category = "C5", subcategory = "BP")
  gobp_t2g = gobp_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  edo.gobp <- enricher(gene = genes.input, TERM2GENE = gobp_t2g)
  if(!is.null(edo.gobp)) {
    if(nrow(edo.gobp) > 0) {
        edo.gobp <- edo.gobp %>%
          mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))
    }
  }
  
  
### Dotplots for each enrichments
  if(!is.null(edo.h)) {
     if(nrow(edo.h) > 0) {
        p1 <- enrichplot::dotplot(edo.h,font.size = 8)
     } else {
        p1 <- NULL
     } 
  } else {
    p1 <- NULL
  }

  p2 <- NULL
  if(!is.null(edo.c7)) {  
     if(nrow(edo.c7) > 0) {
        p2 <- enrichplot::dotplot(edo.c7, font.size = 8)
     }
  }
  
  
  p3 <- NULL
  if(!is.null(edo.c6)) {
     if(nrow(edo.c6) > 0) {
       p3 <- enrichplot::dotplot(edo.c6,font.size = 8)
     } 
  }
  
  p4 <- NULL
  if(!is.null(edo.c2)) {
     if(nrow(edo.c2) > 0) {
       p4 <- enrichplot::dotplot(edo.c2,font.size = 8)
     }
  }
  
  p5 <- NULL
  if(!is.null(edo.cm)) {
     if(nrow(edo.cm) > 0) {
       p5 <- enrichplot::dotplot(edo.cm,font.size = 8)
     }
  }
  
  p6 <- NULL
  if(!is.null(edo.gobp)) {
     if(nrow(edo.gobp) > 0) {
       p6 <- enrichplot::dotplot(edo.gobp,font.size = 8)
     } 
  }
  
  

## if no enrichment handle empty results
  q2 =""
  if(!is.null(edo.h)) {
     if(length(data.frame(edo.h)$ID) <= 2){
        q2 = ""
     } else{
        q2 <- try(cnetplot(edo.h, font.size=6))
     }
  }

  q3 = ""
  if(!is.null(edo.kegg)) {
     if(length(data.frame(edo.kegg)$ID) <= 2){
       q3 = ""
     } else{
       q3 <- try(cnetplot(edo.kegg, font.size=6))
     }
  }

 
    ## Enrichment tables for the above comparisons.
  d1 <- as.data.frame(edo.h)

  d2 <- as.data.frame(edo.c7)

  d3 <- as.data.frame(edo.c6)

  d4 <- as.data.frame(edo.c2)

  d5 <- as.data.frame(edo.cm)

  d6 <- as.data.frame(edo.gobp)
  
  d7 <- as.data.frame(edo.kegg)

  ## Return result list
  return(list(hallmark=p1, immuno=p2, onco=p3, curated=p4, cancer.modules=p5,
              GO.BP=p6, cnetplot=q2, kegg.plot=q3, 
              hallmarkd=d1, immunod=d2, oncod=d3, curatec2d=d4, curatecmd=d5,
              god=d6, keggd=d7))

}
