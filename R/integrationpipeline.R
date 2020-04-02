# Integration ---------------------------





seuratpipeline.integration <- function(objectdir,
                                       keepwholeintmat=NULL,
                                       project=NULL,
                                       conditions=NULL,
                                       PCAgenelist=NULL,
                                       jackstraw=NULL,
                                       dims=NULL,
                                       res=NULL,
                                       plotout.dimplot=NULL
){

  set.seed(500)

  if(is.null( project )) {project <- "IntegratedSeuratObject"}

  #if(is.null( PCAgenelist )) {VariableFeatures(sobjint)} #must be set later?
  if(is.null( jackstraw )) {jackstraw <- TRUE}
  if(is.null( keepwholeintmat )) {keepwholeintmat <- FALSE}
  if(is.null( dims )) {dims <- 1:30}
  if(is.null( res )) {res <- c(0.5, 1.5, 1)}

  if(is.null( plotout.dimplot )) {plotout.dimplot <- T}


  sobjlist <- list()
  for(sobjfile in list.files(objectdir) ){
    sobj <-  readRDS(paste0(objectdir, '/', sobjfile))
    sobjlist[[sobjfile]] <- sobj

    message('Reading in ', sobjfile, '\n')

    rm(sobj, sobjfile)
  }

  message('\nSelecting Features for integration\n')
  #new sct integration steps: select var features
  sobj.features <- SelectIntegrationFeatures(object.list = sobjlist, nfeatures = 3000)
  sobjlist <- PrepSCTIntegration(object.list = sobjlist, anchor.features = sobj.features,
                                 verbose = TRUE)

  message('\nInitiating CCA pipeline to detect integration anchors\n')
  #CCA: find pairwise anchors
  sobjanchors <- FindIntegrationAnchors(object.list = sobjlist,
                                        normalization.method = 'SCT',
                                        anchor.features = sobj.features)

  message('\nIntegrating Datasets, constructing corrected matrix\n')

  if( keepwholeintmat == T) {

    genes <- Reduce(intersect, lapply(sobjlist, rownames))

    sobjint <- IntegrateData(anchorset = sobjanchors, preserve.order = T,
                             normalization.method = "SCT",
                             features.to.integrate = genes)
    rm(genes)

  } else{
    #integrate using anchors, generated corrected matrix
    sobjint <- IntegrateData(anchorset = sobjanchors, preserve.order = T,
                             normalization.method = "SCT"
    )

  }

  rm(sobjanchors, sobjlist, sobj.features)

  #conditions column
  if(is.null( conditions )) {conditions <- unique(sobjint$orig.ident)}
  conditionvector <- c()
  for(condindex in 1:length(conditions)){
    conditionvector <-  c(conditionvector,
                          rep(conditions[condindex], table(sobjint$orig.ident)[condindex] )  )
  }

  sobjint$condition <- conditionvector
  rm(conditionvector, condindex)

  sobjint@project.name <- project

  #PCA
  if(is.null( PCAgenelist )) {PCAgenelist <- VariableFeatures(sobjint)}

  sobjint <- RunPCA(object = sobjint, features = PCAgenelist, verbose = F)
  rm(PCAgenelist)

  #jackstraw
  if(jackstraw == T){
    message('Initiating Jackstraw Procedure\n')
    sobjint <- JackStraw(sobjint, dims = 50, verbose = T)
    sobjint <- ScoreJackStraw(sobjint, dims = 1:50)

    dims <- sobjint@reductions$pca@jackstraw$overall.p.values[sobjint@reductions$pca@jackstraw$overall.p.values[,2] < 0.05,1]

    message('\nJackstraw Procedure completed.\nRetaining PCs: ', list(dims), '\n')
  }

  message('Initiating clustering of integrated dataset for resolution: ', list(res))

  sobjint <- FindNeighbors(object = sobjint, dims = dims, verbose = F)
  sobjint <- FindClusters(object = sobjint, resolution = res, verbose = T)

  message('\nRunning tSNE / UMAP:\n')
  sobjint <- RunTSNE(sobjint, dims = dims)
  sobjint <- RunUMAP(sobjint, dims = dims, verbose = F)

  if(plotout.dimplot == T){
    print(  DimPlot(sobjint, label = F, group.by = 'condition') )
  }

  message('Clustering procedure complete.')
  DefaultAssay(sobjint) <- 'SCT'

  return(sobjint)
}





# DEG across integrated conditions ---------------------------






deg.acrossconditions <- function(object,
                                 assay=NULL,
                                 test=NULL,
                                 latent.vars=NULL,
                                 outdir=NULL){

  set.seed(500)


  if(is.null( assay )) {assay <- 'SCT'}
  if(is.null( latent.vars )) {latent.vars <- NULL}
  if(is.null( test )) {test <- "MAST"}

  if(is.null( outdir )) {outdir <- "conditions_mast/"}
  if(!dir.exists(outdir)){dir.create(outdir)}


  DefaultAssay(object) <- assay

  #fix latent vars for mast...
  #cannot have charactervector / categorical variable as latent var
  #so convert to numeric.... not sure if makes sense...
  if( !is.null( latent.vars ) & test=='MAST' ) {
    for(latent in latent.vars){

      if( !is.numeric(object@meta.data[,latent]) ){
        object@meta.data[,latent] <- unclass( factor(object$orig.ident.tag) )
      }

    }
  }

  #create condition column in object@meta.data
  if( !('condition' %in% names(object@meta.data)) ) {
    object$condotions <- object$orig.ident }




  for(cluster in levels(object@active.ident) ){
    clusterdir <- paste0(outdir, '/cluster', cluster)


    tmp_sobj <- subset(object, idents = cluster)
    Idents(tmp_sobj) <- tmp_sobj@meta.data$condition


    if( length( unique(Idents(tmp_sobj)) ) == length( unique(object$condition ) ) ) {

      #test to see if clusters are done; pickup where left off
      if(dir.exists(clusterdir)){
        if( length(list.files(clusterdir)) == length( unique(Idents(tmp_sobj)) ) ) {
          message('Cluster ', cluster, ' already completed')
        } else {

          message('\nCluster ', cluster, ' partially completed:')

          for(cond in levels(tmp_sobj@active.ident)){

            c1 <- levels(tmp_sobj@active.ident)[levels(tmp_sobj@active.ident) == cond]
            c2 <- levels(tmp_sobj@active.ident)[levels(tmp_sobj@active.ident) != cond]

            if( !file.exists(paste0(clusterdir, "/markers_", c1, ".rds")) ){

              message('\tInitating DEG test ', test, ' for ', c1, ' vs other conditions:')

              tmp <- FindMarkers(object = tmp_sobj, test.use = test,
                                 ident.1 = c1, ident.2 = c2,
                                 logfc.threshold = -Inf, min.pct = -Inf,
                                 min.diff.pct = -Inf, min.cells.group = 1,
                                 latent.vars = latent.vars,
                                 assay = 'SCT')

              saveRDS(tmp, file = paste0(clusterdir, "/markers_", c1, ".rds"))
            } else{ message('\t', 'Cluster ', cluster, ' condition "', c1, '" already completed') }

          }



        }
      }

      #make cluster dir and run DEG
      if(!dir.exists(clusterdir)){
        dir.create(clusterdir)

        message('\tInitating DEG test ', test, ' for cluster ', cluster,':\n')

        for(cond in levels(tmp_sobj@active.ident)){

          c1 <- levels(tmp_sobj@active.ident)[levels(tmp_sobj@active.ident) == cond]
          c2 <- levels(tmp_sobj@active.ident)[levels(tmp_sobj@active.ident) != cond]

          message('\tDEG test ', test, ' for ', c1, ' vs other conditions:')

          tmp <- FindMarkers(object = tmp_sobj, test.use = test,
                             ident.1 = c1, ident.2 = c2,
                             logfc.threshold = -Inf, min.pct = -Inf,
                             min.diff.pct = -Inf, min.cells.group = 1,
                             latent.vars = latent.vars,
                             assay = 'SCT')

          saveRDS(tmp, file = paste0(clusterdir, "/markers_", c1, ".rds"))
        }

        message('\nCompleted DEG test for cluster ', cluster,'.\n')

      }


    } else{ message("\nSkipping cluster ", cluster, ' due to lack of cells in one condition\n') }
  }

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

  #loop thru inputdir
  for(j in stringr::str_sort(list.files(inputfolder), numeric = T) ){
    indir <- paste0(inputfolder, '/', j, '/')

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
      cluster <- strsplit(x = filename, split = "markers_")[[1]][2]

      rm(fgseaRes)


      fgseaResTidy$cluster <- rep(cluster, length(fgseaResTidy$pathway))

      res <- rbind(res, fgseaResTidy)
      rm(fgseaResTidy)
      rm(scores)
      message('\tCompleted ', cluster)

    }

    if(!is.null( clusterlabs )) {
      res$cluster <- rep(clusterlabs, each = length(pathways))}

    #fix broken tibble issue...
    res <- as.data.frame(res)

    #prevent NA...
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
