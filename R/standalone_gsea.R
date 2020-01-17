# GSEA using DEGs from clusters ---------------------------
#' GSEA for clusters.
#'
#' A function to perform GSEA on differentially expressed gene lists from clusters. Relies for input on the output format from deg.acrossclusters(), and gene pathways in the form of named lists. Also relies on the FGSEA package by Alexey Sergushichev.
#' unfortunately, currently this function only works for Mac users due to PDF encoding limitations special symbols (ie, points denoting statistical signifcance on heatmaps), but windows portability will come very soon.
#' Output plots can be difficult to read if many pathways are used. Recommended 50-75 or so pathways as the maximum, many more will generate very ugly heatmaps.
#' @param inputfolder A string containing the path to a directory of .rds files storing dataframes outputted by Seurat::FindMarkers(). This is the format of the output of tamlabscpipeline::deg.acrossclusters().
#' @param pathways A named "list of lists" of genes. The format of this object is important for functionality of this function and can be previewed in the tamlabscpipeline::hallmark object. The names of each list element provides the Y axis; the genes within each list element are the target genes used for GSEA. If NULL, defaults to the Msigdb's Hallmark pathways for mouse.
#' @param nperm An integer denoting how many permutations FGSEA will run. Default = 10000.
#' @param makepdf T/F. Whether to print to an output PDF or not. Default = F.
#' @param pdfname A string denoting the name of putput pdf. No effect if makepdf == F. Default is 'clusterGSEA_(date).pdf'
#' @param filter_nonsig_pathways T/F. Whether to remove non-significant rows from resulting heatmap. Useful if running very large numbers of pathways, or to generate finalized plots after exploratory analysis. Default = F.
#' @return Prints a heatmap to the standard out, or to a pdf.
#' @examples
#' \dontrun{
#' gseapipline.clusters(inputfolder = 'deg.acrossclusters.outputdir/')
#' }
gsea.clusters <- function(inputfolder,
                                  pathways=NULL,
                                  nperm=NULL,
                                  #weightmethod=NULL,
                                  makepdf=NULL,
                                  pdfname=NULL,
                                  filter_nonsig_pathways=NULL){

  set.seed(500)

  if(is.null( pathways )) {pathways <- tamlabscpipeline::hallmark}
  if(is.null( nperm )) {nperm <- 10000}
  #if(is.null( weightmethod )) {weightmethod <- 'pvalue'}
  weightmethod <- 'pvalue'
  if(is.null( makepdf )) { makepdf <- F}
  if(is.null( filter_nonsig_pathways )) {filter_nonsig_pathways <- F}
  if(is.null( pdfname )) { pdfname <- paste0('clusterGSEA_', Sys.Date(), '.pdf') }


  if( makepdf == T) {quartz(type = 'pdf', file = pdfname, width = 8)}


  #indir <- "clusters_mast/"
  files <- stringr::str_sort(list.files(inputfolder), numeric = T)

  res <- data.frame()
  for(i in files){
    #message('Running ', i)
    tmp <- readRDS( paste0(inputfolder, '/',  i) )

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

    filename <- strsplit(x = basename(i), split = "\\.")[[1]][1]
    cluster <- strsplit(x = filename, split = "\\_")[[1]][2]

    rm(fgseaRes)


    fgseaResTidy$cluster <- rep(cluster, length(fgseaResTidy$pathway))

    res <- rbind(res, fgseaResTidy)
    rm(fgseaResTidy)
    rm(scores)
    message('\tCompleted ', cluster)

  }



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
  #rm(res, res_tmp)

  if( makepdf == T) {dev.off()}

}


# GSEA using DEG across conditions ---------------------------




gseapipeline.conditions <- function(inputfolder,
                                    pathways,
                                    nperm=NULL,
                                    weightmethod=NULL,
                                    clusterlabs=NULL,
                                    makepdf=NULL,
                                    pdfname=NULL,
                                    filter_nonsig_pathways=NULL){

  set.seed(500)

  if(is.null( nperm )) {nperm <- 10000}
  if(is.null( weightmethod )) {weightmethod <- 'pvalue'}
  if(is.null( makepdf )) { makepdf <- F}
  if(is.null( filter_nonsig_pathways )) {filter_nonsig_pathways <- F}

  if( makepdf == T) {quartz(type = 'pdf', file = pdfname, width = 8)}

  for(j in stringr::str_sort(list.files(inputfolder), numeric = T) ){
    indir <- paste0(inputfolder, '/', j, '/')

    #indir <- "clusters_mast/"
    files <- list.files(indir)

    message('\nRunning FGSEA on ',j, ':')

    res <- data.frame()
    for(i in files){
      tmp <- readRDS(paste0(indir, '/',  i))

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

      filename <- strsplit(x = basename(i), split = "\\.")[[1]][1]
      cluster <- strsplit(x = filename, split = "\\_")[[1]][2]

      rm(fgseaRes)


      fgseaResTidy$cluster <- rep(cluster, length(fgseaResTidy$pathway))

      res <- rbind(res, fgseaResTidy)
      rm(fgseaResTidy)
      rm(scores)
      message('\tCompleted ', cluster)

    }

    if(!is.null( clusterlabs )) {
      res$cluster <- rep(clusterlabs, each = length(pathways))}


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
      labs(subtitle = j,
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
    #rm(res, res_tmp)

  }

  if( makepdf == T) {dev.off()}

}
