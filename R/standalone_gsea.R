# GSEA using DEGs from clusters ---------------------------
#' GSEA for clusters.
#'
#' A function to perform GSEA on differentially expressed gene lists from clusters and plot as a heatmap.
#'
#' Expects as input a data.frame with column names consistent with output of Seurat::FindAllMarkers(). Essentially, a dataframe where each row corresponds to a gene, and columns consists of the following: "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", and "gene". Minimally requires "cluster", "gene", and "avg_logFC". Pvalue and/or adjusted pvalue can be included if weighting by those values.
#'
#' Pathways must be a named list of character vectors; by default uses mouse orthologs of Msigdb Hallmark sets, accessible via tamlabscpipeline::hallmark.
#'
#' See FGSEA package docs for details. Also check out the Broad's Msigdb to view a large database of pathways and the MsigdbR package to get orthologs of the Msigdb pathways for your species of interest.
#'
#' Output plots can be difficult to read if many pathways are used. Recommended 50-75 or so pathways as the maximum, many more will generate very ugly heatmaps.
#' Utilizes special characters, so pdf printing may be problematic. Alternative PDF devices should be able to handle this though, one can use Quartz if on Mac, or CairoPDF if on Windows.
#' @param degdf A dataframe in the format of the output of Seurat::FindAllMarkers.
#' @param pathways A named "list of lists" of genes. The format of this object is important for functionality of this function and can be previewed in the tamlabscpipeline::hallmark object. The names of each list element provides the Y axis; the genes within each list element are the target genes used for GSEA. If NULL, defaults to the Msigdb's Hallmark pathways for mouse.
#' @param nperm An integer denoting how many permutations FGSEA will run. Default = 10000.
#' @param weightmethod A string, one of either "pvalue", "padj", or "foldchange". Default = "pvalue".
#' @param onlypos T/F. Whether to filter each cluster to only include upregulated genes. Default = F.
#' @param filter_nonsig_pathways T/F. Whether to remove non-significant rows from resulting heatmap. Useful if running very large numbers of pathways, or to generate finalized plots after exploratory analysis. Default = F.
#' @return Returns a ggplot object, and prints a heatmap to the standard out.
#' @examples
#' \dontrun{
#' sobjmarkers <- Seurat::FindAllMarkers(sobj)
#' gsea.clusters(sobjmarkers)
#' }
gsea.clusters <- function(degdf,
                          pathways=NULL,
                          nperm=NULL,
                          weightmethod=NULL,
                          onlypos=NULL,
                          filter_nonsig_pathways=NULL){

  set.seed(500)

  #if(is.null( clustercolname )) {clustercolname <- 'cluster'}
  clustercolname <- 'cluster'
  if(is.null( pathways )) {pathways <- tamlabscpipeline::hallmark}
  if(is.null( nperm )) {nperm <- 10000}
  if(is.null( weightmethod )) {weightmethod <- 'pvalue'}
  if(is.null( onlypos )) {onlypos <- F}

  #weightmethod <- 'pvalue'
  if(is.null( filter_nonsig_pathways )) {filter_nonsig_pathways <- F}

  res <- data.frame()
  for(i in unique(degdf[,clustercolname]) ){

    tmp <- degdf[degdf[,clustercolname] == i, ]


    if(onlypos == T){
      tmp[tmp$avg_logFC >0,]
    }


    if(weightmethod == 'pvalue') {
      scores <- log(tmp$p_val)

      #fix underflow
      if(any(scores == -Inf)){

        numuf = length(scores[scores==-Inf])
        adder = rev(1:numuf)
        for(uf in rev(1:numuf)){
          scores[uf] <- scores[numuf + 1] + (adder[uf] * -1)
        }

      }

      scores <- (-1 * scores) * sign(tmp$avg_logFC)
      names(scores) <- tmp$gene
      rm(tmp)
      scores <- sort(scores, decreasing = T)
    }

    if(weightmethod == 'padj') {
      scores <- log(tmp$p_val_adj)

      #fix underflow
      if( any(scores == -Inf) ){

        numuf = length(scores[scores==-Inf])
        adder = rev(1:numuf)
        for(uf in rev(1:numuf)){
          scores[uf] <- scores[numuf + 1] + (adder[uf] * -1)
        }

      }

      scores <- (-1 * scores) * sign(tmp$avg_logFC)
      names(scores) <- tmp$gene
      rm(tmp)
      scores <- sort(scores, decreasing = T)
    }

    if(weightmethod == 'foldchange'){
      scores <- 10^tmp$avg_logFC

      names(scores) <- tmp$gene
      rm(tmp)
      scores <- sort(scores, decreasing = T)
    }


    suppressWarnings( fgseaRes <- fgsea(pathways=pathways, stats=scores, nperm=nperm) )
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))


    rm(fgseaRes)


    fgseaResTidy$cluster <- rep(paste0('cluster', i), length(fgseaResTidy$pathway))

    res <- rbind(res, fgseaResTidy)
    rm(fgseaResTidy)
    rm(scores)
    message('\tCompleted ', i)

  }


  #fix NaN errors... origin unclear, may be related to low nperm?
  res$NES[is.na(res$NES)] <- 0


  #to properly order clusters above 10 (so stupid)
  numorder <- stringr::str_sort(res$cluster, numeric = T)
  res$cluster <- factor(res$cluster, levels = unique(numorder))


  #to remove nonsignificant rows...
  if(filter_nonsig_pathways == T){
    res2 <-  data.frame()
    for(set in unique(res$pathway)){
      res_tmp <- res[res$pathway == set,]
      if(any(res_tmp$padj <= 0.25)){
        res2 <- rbind(res2, res_tmp)
      } else {NULL}
    }


    if( nrow(res2) > 0 ) {
      res <- res2
      rm(res2) } else { message('No significant results for ', j, ' detected.') }

  }

  #plot
  gseaout <- ggplot(res, aes( cluster, forcats::fct_rev(pathway) ) )+
    geom_tile(aes(fill = NES), color = "white")+
    scale_fill_gradient2(low = "steelblue", high = "red",
                         midpoint = 0, mid = "white")+
    theme_minimal()+
    theme(legend.title = element_text(size = 10, face="bold"),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
    )+
    labs(#subtitle = j,
      caption = paste0(sprintf("\u25cb"), ' = FDR < 0.25\n',
                       sprintf("\u25cf"), ' = FDR < 0.05'))

  if( nrow( res[res[,3] <= 0.25 & res[,3] >= 0.05 ,] ) > 0 ){

    res25 <- res[res[,3] <= 0.25 & res[,3] >= 0.05,]
    gseaout <- gseaout + geom_text(data = res25,
                                   aes(x = cluster, y = forcats::fct_rev(pathway),
                                       label = sprintf("\u25cb")),
                                   family = "Arial Unicode MS")

  }

  if( nrow( res[res[,3] < 0.05,] ) > 0 ){
    res05 <- res[res[,3] < 0.05,]
    gseaout <- gseaout + geom_text(data = res05,
                                   aes(x = cluster, y = forcats::fct_rev(pathway),
                                       label = sprintf("\u25cf")),
                                   family = "Arial Unicode MS")



  }


  print(gseaout)

}


# GSEA using DEG across conditions ---------------------------
#' GSEA for conditions.
#'
#' Perform and plot GSEA heatmap for comparison across condtions.
#'
#' Expects input as list of dataframes. Each dataframe corresponds to a particular condtion; expected dataframe format matches output of Seurat::FindMarkers(). Multiple conditions (ie >= 2) can be run.
#'
#' Pathways must be a named list of character vectors; by default uses mouse orthologs of Msigdb Hallmark sets, accessible via tamlabscpipeline::hallmark.
#'
#' See FGSEA package docs for details. Also check out the Broad's Msigdb to view a large database of pathways and the MsigdbR package to get orthologs of the Msigdb pathways for your species of interest.
#'
#' Output plots can be difficult to read if many pathways are used. Recommended 50-75 or so pathways as the maximum, many more will generate very ugly heatmaps.
#' Utilizes special characters, so pdf printing may be problematic. Alternative PDF devices should be able to handle this though, one can use Quartz if on Mac, or CairoPDF if on Windows.
#' @param inputlist A list of dataframes. Must be named. Names will be used as x axis. Order of names retained in axis.
#' @param pathways A named "list of lists" of genes. The format of this object is important for functionality of this function and can be previewed in the tamlabscpipeline::hallmark object. The names of each list element provides the Y axis; the genes within each list element are the target genes used for GSEA. If NULL, defaults to the Msigdb's Hallmark pathways for mouse.
#' @param nperm An integer denoting how many permutations FGSEA will run. Default = 10000.
#' @param weightmethod A string, one of either "pvalue", "padj", or "foldchange". Default = "pvalue".
#' @param onlypos T/F. Whether to filter each cluster to only include upregulated genes. Default = F.
#' @param filter_nonsig_pathways T/F. Whether to remove non-significant rows from resulting heatmap. Useful if running very large numbers of pathways, or to generate finalized plots after exploratory analysis. Default = F.
#' @return Returns a ggplot object, and prints a heatmap to the standard out.
#' @examples
#' \dontrun{
#' inputlist <- list(treatment = treatmentdf, control = controldf)
#' gsea.clusters(sobjmarkers)
#' }
gsea.conditions <- function(inputlist,
                            pathways=NULL,
                            nperm=NULL,
                            weightmethod=NULL,
                            onlypos=NULL,
                            filter_nonsig_pathways=NULL){

  set.seed(500)


  if(is.null( pathways )) {pathways <- tamlabscpipeline::hallmark}
  if(is.null( nperm )) {nperm <- 10000}
  if(is.null( weightmethod )) {weightmethod <- 'pvalue'}
  if(is.null( onlypos )) {onlypos <- F}
  if(is.null( filter_nonsig_pathways )) {filter_nonsig_pathways <- F}




  res <- data.frame()
  for(i in names(inputlist) ){

    tmp <- inputlist[[i]]

    if(onlypos == T){
      tmp[tmp$avg_logFC >0,]
    }

    if(weightmethod == 'pvalue') {
      scores <- log(tmp$p_val)

      #fix underflow
      if(any(scores == -Inf)){

        numuf = length(scores[scores==-Inf])
        adder = rev(1:numuf)
        for(uf in rev(1:numuf)){
          scores[uf] <- scores[numuf + 1] + (adder[uf] * -1)
        }

      }

      scores <- (-1 * scores) * sign(tmp$avg_logFC)
      names(scores) <- rownames(tmp)
      rm(tmp)
      scores <- sort(scores, decreasing = T)
    }

    if(weightmethod == 'padj') {
      scores <- log(tmp$p_val_adj)

      #fix underflow
      if(any(scores == -Inf)){

        numuf = length(scores[scores==-Inf])
        adder = rev(1:numuf)
        for(uf in rev(1:numuf)){
          scores[uf] <- scores[numuf + 1] + (adder[uf] * -1)
        }

      }

      scores <- (-1 * scores) * sign(tmp$avg_logFC)
      names(scores) <- rownames(tmp)
      rm(tmp)
      scores <- sort(scores, decreasing = T)
    }

    if(weightmethod == 'foldchange'){
      scores <- 10^tmp$avg_logFC
      names(scores) <- rownames(tmp)
      rm(tmp)
      scores <- sort(scores, decreasing = T)
    }


    suppressWarnings( fgseaRes <- fgsea(pathways=pathways, stats=scores, nperm=nperm) )
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))


    rm(fgseaRes)


    fgseaResTidy$cluster <- rep(i, length(fgseaResTidy$pathway))

    res <- rbind(res, fgseaResTidy)
    rm(fgseaResTidy)
    rm(scores)
    message('Completed ', i)

  }



  res$NES[is.na(res$NES)] <- 0


  #to properly order clusters above 10 (so stupid)
  #numorder <- stringr::str_sort(res$cluster, numeric = T)
  #res$cluster <- factor(res$cluster, levels = unique(numorder))

  res$cluster <- factor(res$cluster, levels = names(inputlist) )

  #to remove nonsignificant rows...
  if(filter_nonsig_pathways == T){
    res2 <-  data.frame()
    for(set in unique(res$pathway)){
      res_tmp <- res[res$pathway == set,]
      if(any(res_tmp$padj <= 0.25)){
        res2 <- rbind(res2, res_tmp)
      } else {NULL}
    }


    if( nrow(res2) > 0 ) {
      res <- res2
      rm(res2) } else { message('No significant results for ', j, ' detected.') }

  }

  #plot
  gseaout <- ggplot(res, aes( cluster, forcats::fct_rev(pathway) ) )+
    geom_tile(aes(fill = NES), color = "white")+
    scale_fill_gradient2(low = "steelblue", high = "red",
                         midpoint = 0, mid = "white")+
    theme_minimal()+
    theme(legend.title = element_text(size = 10, face="bold"),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
    )+
    labs(
      caption = paste0(sprintf("\u25cb"), ' = FDR < 0.25\n',
                       sprintf("\u25cf"), ' = FDR < 0.05'))

  if( nrow( res[res[,3] <= 0.25 & res[,3] >= 0.05 ,] ) > 0 ){

    res25 <- res[res[,3] <= 0.25 & res[,3] >= 0.05,]
    gseaout <- gseaout + geom_text(data = res25,
                                   aes(x = cluster, y = forcats::fct_rev(pathway),
                                       label = sprintf("\u25cb")),
                                   family = "Arial Unicode MS")

  }

  if( nrow( res[res[,3] < 0.05,] ) > 0 ){
    res05 <- res[res[,3] < 0.05,]
    gseaout <- gseaout + geom_text(data = res05,
                                   aes(x = cluster, y = forcats::fct_rev(pathway),
                                       label = sprintf("\u25cf")),
                                   family = "Arial Unicode MS")



  }


  return(gseaout)
  #rm(res, res_tmp)


}
