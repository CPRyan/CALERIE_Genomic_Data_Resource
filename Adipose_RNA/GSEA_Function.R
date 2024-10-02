
# From
# https://bioinformaticsbreakdown.com/how-to-gsea/
# Also found in my zotero

GSEA = function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO, 
                        stats = gene_list,
                        minSize=15,
                        maxSize=600,
                        nperm=10000) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) 
  print(dim(fgRes))
  
  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(15,600))
  
  ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  print(dim(rbind(ups,downs)))
  ## Define up / down pathways which are significant in both tests
  keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
  keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
  
  fgRes = fgRes[ !is.na(match(fgRes$pathway, 
                              c( keepups$pathway, keepdowns$pathway))), ] %>% 
    arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  
  
  upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
  downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
  colos = c(upcols, downcols)
  names(colos) = 1:length(colos)
  filtRes$Index = as.factor(1:nrow(filtRes))
  
  g = ggplot(filtRes %>% mutate(pathway = sub("^GOBP_", "", pathway))
, aes(reorder(pathway, NES), NES)) +
    geom_col( aes(fill = Index )) +
    scale_fill_manual(values = colos ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes") + 
    theme_minimal()+
    theme(legend.position = 'none')
  
  
  output = list("Results" = fgRes, "Plot" = g)
  return(output)
}