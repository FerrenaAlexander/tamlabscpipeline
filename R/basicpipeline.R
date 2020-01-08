# Process single dataset ---------------------------
#'Basic SC Pipeline for QC and clustering
#'
#' Basic SC pipeline for a single input dataset.
#' Perform QC and clustering of an input dataset.
#'
#' @param data a string containing a filepath with format connoted by the format parameter.
#' @param format a string, either 'dir' for cellranger dir output; 'h5' for cellranger h5 output, or 'kallisto' for the kallisto|bustools pipeline output
#' @param transcript_gene_file a string containing a filepath for the "transcript_gene" conversion file. Only used for Kallisto|bustools workflow.
#' @param baselinefilter.mad T/F; whether to perform "global" QC, ie without pre-clustering. Default = False.
#' @param baseline.mito.filter T/F; whether to perform global maximum mitochondrial content filtration using median absolute deviation threshold; will only work if baselinefilter.mad is set to True. Default is True.
#' @param madmax.dist.percentmito.baseline a numeric, or a string reading 'predict'. If numeric is provided, will use this as median absolute deviation threshold for global mito cutoff. If set to 'predict', will attempt to learn cutoff from data. Default = 'predict'
#' @param baseline.libsize.filter T/F; whether to perform global minimal lib size filtration using median absolute deviation threshold; will only work if baselinefilter.mad is set to True. Default is True.
#' @param madmax.dist.nCount_RNA.baseline a numeric, or a string reading 'predict'. If numeric is provided, will use this as median absolute deviation threshold for global libsize cutoff. If set to 'predict', will attempt to learn cutoff from data. Default = 'predict'
#' @param removemitomaxclust T/F ; whether to identify and remove abnormally high mitochondrial content clusters after first-pass clustering; if no baseline mito filtration is used, there will almost certainly be a mito-driven cluster. Identification is via Grubbs' test for outliers, based on Lukasz Komsta's implementation in the outliers package. Default is True.
#' @param iterativefilter T/F ; whether to perform iterative filtering on first-pass filtering. Default is True.
#' @param iterativefilter.libsize either one of two strings ('twosided' or 'lefttail') or False. Twosided will attempt to learn median abs. dev. cutoffs for both min and max to catch debris and, ostensibly, doublets. Lefttail will attempt to learn median abs. dev. threshold for min to catch debris; doublets may not accurately be captured by max cutoffs as this is more an artifact of sequencing than cell suspension. False will skip. Default is 'lefttail'.
#' @param iterativefilter.mito T/F; whether to learn right-tail median abs. dev. thresholds and filter maximal mitochondrial content from each cluster. May incorrectly remove mito-okay cells while missing true mito-hi cells. Default = F.
#' @param cellcycleregression either one of two strings ('total', 'difference'), or False. If false, still calculates cell cycle score but does not attempt correction. Default is False. See here: https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html
#' @param PCAgenelist a character vector of gene names to use for PCA. If null, defaults to highly variable genes called by SCT. default is NULL.
#' @param jackstraw T/F; whether to perform jackstraw to score significant PCs for use in clustering / dimreduction. May be incompatible with SCT. Default is false.
#' @param dims an integer range. Controls graph construction prior to clustering and dimensionality reduction for visualization. Connotes which PC dimensions to use in clustering. Defaults to 1:30.
#' @param res a numeric, vector of numerics, or range of numerics. Controls assignment of cells to clusters. Connotes the "resolution" paramter used as correction in Louvain clustering. Defaults to c(0.5, 1.0, 1.5).
#' @return will return a bunch of plots related to QC and an output in the form of a Seurat object to the standard out.
#' @examples
#' \dontrun{
#' pdf('qcplots.pdf')
#' sobj <- seuratpipeline('datafilepath.h5', format=h5)
#' dev.off()
#' }
seuratpipeline <- function(data,
                           format,
                           transcript_gene_file,
                           project=NULL,

                           baselinefilter.mad=NULL,
                           baseline.mito.filter=NULL,
                           madmax.dist.percentmito.baseline=NULL,
                           baseline.libsize.filter=NULL,
                           madmax.dist.nCount_RNA.baseline=NULL,

                           removemitomaxclust=NULL,
                           iterativefilter=NULL,
                           iterativefilter.libsize=NULL,
                           iterativefilter.libsize.twosided=NULL,
                           iterativefilter.libsize.lefttail=NULL,
                           iterativefilter.mito=NULL,


                           cellcycleregression=NULL,
                           PCAgenelist=NULL,
                           jackstraw=NULL,
                           dims=NULL,
                           res=NULL){


  set.seed(500)

  if(is.null( project )) {project <- "SeuratProject"}

  # baseline (global) filtration
  if(is.null( baselinefilter.mad )) {baselinefilter.mad <- F}
  if(is.null( baseline.mito.filter )) {baseline.mito.filter <- T}
  if(is.null( madmax.dist.percentmito.baseline )) {madmax.dist.percentmito.baseline <- 'predict'}

  if(is.null( baseline.libsize.filter )) {baseline.libsize.filter <- T}
  if(is.null( madmax.dist.nCount_RNA.baseline )) {madmax.dist.nCount_RNA.baseline <- 'predict'}

  # baseline (global) filtration
  if(is.null( iterativefilter )) {iterativefilter <- T}
  if(is.null( removemitomaxclust )) {removemitomaxclust <- T}
  if(is.null( iterativefilter.libsize )) {iterativefilter.libsize <- 'lefttail'}
  if(is.null( iterativefilter.mito )) {iterativefilter.mito <- F}

  if(is.null( cellcycleregression )) {cellcycleregression <- F}

  #if(is.null( PCAgenelist )) {VariableFeatures(tmp)} #must be set later!
  if(is.null( jackstraw )) {jackstraw <- F}
  if(is.null( dims )) {dims <- 1:30}
  if(is.null( res )) {res <- c(0.5, 1.5, 1)}


  #read in
  if(format == 'dir'){
    tmp <- CreateSeuratObject(Read10X(data),
                              project = project,
                              min.features = 200, min.cells = 3)
  }

  if(format == 'h5'){
    tmp <- CreateSeuratObject(Read10X_h5(data),
                              project = project,
                              min.features = 200, min.cells = 3)
  }

  if(format == 'kallisto'){
    #read in data from kallisto | bustools
    res_mat <- BUSpaRse::read_count_output(dir = data, name = 'gene', tcc = F)

    #initial cell vs debris classification
    bc_rank <- barcodeRanks(res_mat)

    #keep cells above inflection point
    tot_counts <- Matrix::colSums(res_mat)
    res_mat <- res_mat[, tot_counts > metadata(bc_rank)$inflection]
    dim(res_mat)

    rm(bc_rank, tot_counts)

    #rename genes from ensebl IDs to MGI gene symbols
    tr2g <- read.table(transcript_gene_file, sep = '\t',
                       stringsAsFactors = F,
                       col.names = c('transcript', 'gene', 'gene.symbol'))

    mg <- rownames(res_mat)

    tr2g <- tr2g[tr2g[,2] %in% mg, 2:3]
    tr2g <-  tr2g[match(mg, tr2g[,1]),]

    tr2g[,2] <- make.unique(tr2g[,2])

    rownames(res_mat) <- tr2g[,2]

    rm(mg, tr2g)

    tmp <- CreateSeuratObject(res_mat, project = project,
                              min.cells = 3, min.features = 200)
  }

  message('Reading in ', project ,'; num cells = ', ncol(tmp)  , '\n')


  #mito content, add to metadata
  mito.features <- grep(pattern = "^mt-", x = rownames(x = tmp), value = TRUE)
  percent.mito <- Matrix::colSums(
    x = GetAssayData(object = tmp, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = tmp, slot = 'counts')
    )
  tmp[['percent.mito']] <- percent.mito
  rm(mito.features, percent.mito)

  #whether to perform baseline filtering at all.
  if(baselinefilter.mad == T) {

    #find mito cutoff high
    if(baseline.mito.filter == T){

      message('Initiating global baseline mito hi filtration')

      x <- tmp$percent.mito
      m         <- median(x)
      abs.dev   <- abs(x - m)

      right.mad <- median(abs.dev[x>=m])

      x.mad <- rep(right.mad, length(x))

      mad.distance <- abs(x - m) / x.mad
      mad.distance[x==m] <- 0

      #select bad cells (outside of MAD cutoff)
      if(madmax.dist.percentmito.baseline == 'predict'){

        message('Predicting baseline mito MAD threshold...')

        #predict mad cutoff. for right tail, x > m
        y=sort(mad.distance[x>m])
        ecpout <- ecp::e.divisive(diff(as.matrix(y)), k = 2, min.size = 2)
        #this step is non-trivial... may have to play around with k value and which estimate to use.


        ecpp = ecpout$estimates[2]
        madmax.dist.percentmito.baseline = as.numeric(y[ecpout$estimates[2]])
        ecpplot <- data.frame(x = 1:length(y), y = y) #should have called this ecpdf...

        ecpplot1 <- ggplot(data=ecpplot, aes(x = x, y = y)) +
          geom_point(shape=1, size = 0.1) +
          geom_hline(yintercept = madmax.dist.percentmito.baseline, col = 'red', linetype = 'dashed') +
          geom_vline(xintercept = ecpp, col = 'red', linetype = 'dashed')+
          geom_point(mapping = aes(x = ecpp, y = madmax.dist.percentmito.baseline),col = 'red') +
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')


        ecpplot2 <- ggplot(data=ecpplot, aes(x = 0, y = y)) +
          geom_violin(fill = 'steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = madmax.dist.percentmito.baseline, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')+
          scale_x_discrete(limits=0)
        #scaling x this way so the line is at the same level in the cowplot

        title <- ggdraw() +
          draw_label(
            paste0('ECP Prediction = ', ecpp),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(ecpplot1, ecpplot2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(y, ecpout, ecpplot, ecpplot1, ecpplot2, ecpp, plotrow, title)
      } #end predict block

      bad <- x[mad.distance > madmax.dist.percentmito.baseline]

      #for plotting / message output: get mito upper bound
      mitohi <- ifelse( test = length(bad[bad > median(x)]) == 0,
                        yes = Inf,
                        no = min(bad[bad > median(x)])
      )


      print(
        ggplot(tmp@meta.data, aes(x = 0, y = percent.mito))+
          geom_violin(fill='steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = mitohi, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          labs(title = paste0("percent mito cutoff"),
               subtitle = paste0(length(rownames(tmp@meta.data)), " cells total",
                                 "; cutoff for percent mito is ", as.character(round(mitohi, digits = 3))),
               caption = paste0('Mito max MAD = ', round(madmax.dist.percentmito.baseline, 3) ,
                                '\npercent.mito cutoff = ', round(mitohi, digits = 3),
                                '\nNum cells presubset = ', ncol(tmp),
                                '\nNum cells remaining = ', length(tmp$percent.mito[tmp@meta.data$percent.mito < mitohi]), '\n'))+
          xlab('Cells')
      )


      filteredcells <- names(x[!names(x) %in% names(bad)] )


      message('Mito max MAD = ', round(madmax.dist.percentmito.baseline, 3) ,
              '\npercent.mito cutoff = ', round(mitohi, digits = 3),
              '\nNum cells presubset = ', ncol(tmp),
              '\nNum cells remaining = ', length(filteredcells), '\n')



      #save cutoff params...
      tmp@commands$baselinemitofilter <- list('MadThresh' = madmax.dist.percentmito.baseline,
                                              'PercMitoCutoff' = mitohi,
                                              'CellsPreSub' = ncol(tmp),
                                              'CellsPostSub' = length(filteredcells)
      )

      tmp <- subset(x = tmp,
                    cells = filteredcells)

      rm(abs.dev, bad, right.mad, m, x, x.mad, mad.distance,filteredcells,
         percent.mito, madmax.dist.percentmito.baseline, mitohi, mito.features)


    } #end mito baseline block



    #find libsize cutoff low.
    if(baseline.libsize.filter == T){

      message('Initiating global baseline lib size low filtration')

      x <- tmp$nCount_RNA
      m         <- median(x)
      abs.dev   <- abs(x - m)

      left.mad <- median(abs.dev[x<=m])

      x.mad <- rep(left.mad, length(x))

      mad.distance <- abs(x - m) / x.mad
      mad.distance[x==m] <- 0

      #select bad cells (outside of MAD cutoff)
      if(madmax.dist.nCount_RNA.baseline == 'predict'){

        message('Predicting baseline lib size threshold...')

        #predict mad cutoff. for left tail, x < m...
        y=sort(mad.distance[x<m])
        ecpout <- ecp::e.divisive(diff(as.matrix(y)), k = 1, min.size = 2)
        #this step is non-trivial... may have to play around with k value and which estimate to use.


        ecpp = ecpout$estimates[2]
        madmax.dist.nCount_RNA.baseline = as.numeric(y[ecpout$estimates[2]])
        ecpplot <- data.frame(x = 1:length(y), y = y) #should have called this ecpdf...

        ecpplot1 <- ggplot(data=ecpplot, aes(x = x, y = y)) +
          geom_point(shape=1, size = 0.1) +
          geom_hline(yintercept = madmax.dist.nCount_RNA.baseline, col = 'red', linetype = 'dashed') +
          geom_vline(xintercept = ecpp, col = 'red', linetype = 'dashed')+
          geom_point(mapping = aes(x = ecpp, y = madmax.dist.nCount_RNA.baseline),col = 'red') +
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')


        ecpplot2 <- ggplot(data=ecpplot, aes(x = 0, y = y)) +
          geom_violin(fill = 'steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = madmax.dist.nCount_RNA.baseline, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')+
          scale_x_discrete(limits=0)
        #scaling x this way so the line is at the same level in the cowplot

        title <- ggdraw() +
          draw_label(
            paste0('ECP Prediction = ', ecpp),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(ecpplot1, ecpplot2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(y, ecpout, ecpplot, ecpplot1, ecpplot2, ecpp, plotrow, title)
      }

      bad <- x[mad.distance > madmax.dist.nCount_RNA.baseline]


      ncountslo <- ifelse( test = length(bad[bad < median(x)]) == 0,
                           yes = 0,
                           no = max(bad[bad < median(x)])
      )

      #plot cutoffs for library size / UMIs
      print(
        ggplot(tmp@meta.data) +
          geom_histogram(aes(nCount_RNA),
                         color="black", fill = "steelblue",
                         binwidth = 0.05)+
          #geom_vline(xintercept = ncountshi, linetype = "dashed", colour = "red")+
          geom_vline(xintercept = ncountslo, linetype = "dashed", colour = "red")+
          geom_vline(xintercept = median(tmp@meta.data$nCount_RNA),
                     linetype = "dotted", colour = "red", size = 1.2)+
          scale_x_log10()+
          labs(title = paste0("Library Size per Cell"),
               x = "Library size (UMI, aka 'nCount'), log10 scale",
               y = "Number of cells" ,
               subtitle = paste0(length(rownames(tmp@meta.data)), " cells; median lib size = ",
                                 median(tmp@meta.data$nCount_RNA)),
               caption = paste0("cells remaining = ",
                                length(x) - length(bad)
               )
          )+
          theme_linedraw()
      )

      print(
        ggplot(tmp@meta.data, aes(x = 0, y = nCount_RNA))+
          geom_violin(fill='steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = ncountslo, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          scale_y_log10()+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          labs(title = paste0("LibSize Cutoff"),
               subtitle = paste0(length(rownames(tmp@meta.data)), " cells total",
                                 "; cutoff for libsize low is ", as.character(round(ncountslo, digits = 3))),
               caption = paste0('Mito max MAD = ', round(madmax.dist.nCount_RNA.baseline, 3) ,
                                '\npercent.mito cutoff = ', round(ncountslo, digits = 3),
                                '\nNum cells presubset = ', ncol(tmp),
                                '\nNum cells remaining = ', length(tmp$percent.mito[tmp@meta.data$nCount_RNA > ncountslo]), '\n'))+
          xlab('Cells')
      )



      filteredcells <- names(x[!names(x) %in% names(bad)] )


      message('Libsize global max MAD = ', round(madmax.dist.nCount_RNA.baseline, 3) ,
              '\nLibSize Low Cutoff = ', ncountslo,
              '\nNum cells presubset = ', ncol(tmp),
              '\nNum cells remaining = ', length(filteredcells), '\n')

      #save cutoff params...
      tmp@commands$baselinelibsizefilter <- list('MadThresh' = madmax.dist.nCount_RNA.baseline,
                                                 'PercMitoCutoff' = ncountslo,
                                                 'CellsPreSub' = ncol(tmp),
                                                 'CellsPostSub' = length(filteredcells)
      )

      tmp <- subset(x = tmp,
                    cells = filteredcells)

      rm(abs.dev, bad, left.mad, m, x, x.mad, mad.distance, y, ncountslo, madmax.dist.nCount_RNA.baseline, filteredcells)


    }

  }





  message('Initiating first-pass clustering\n')

  #normalize and cluster, first pass, for lib size filtration
  suppressWarnings(tmp <- SCTransform(tmp, verbose = T))

  tmp <- RunPCA(object = tmp, verbose = F)

  tmp <- FindNeighbors(object = tmp, dims = 1:30, verbose = F)
  tmp <- FindClusters(object = tmp, resolution = 1.0, verbose = T)

  #record params for all clusters...
  md <- tmp@meta.data
  avgs <- aggregate(data = md, cbind(nCount_RNA , nFeature_RNA , percent.mito) ~ seurat_clusters, FUN=mean)
  rm(md)

  print(VlnPlot(tmp, c('percent.mito', 'nCount_RNA', 'nFeature_RNA'), ncol=1, pt.size = 0.1))

  ### remove the mito cluster and recluster ###
  if(removemitomaxclust == T){
    message('\nRemoving mito cluster(s)')

    avgstmp <- avgs
    gt <- outliers::grubbs.test(avgstmp$percent.mito)

    #grubbs test while loop for right-tail outliers, w/ bonferroni correction
    iter=1
    while(
      str_split(gt$alternative, ' ')[[1]][1] == 'highest' &
      gt$p.value < (0.05 / nrow(avgstmp))
    ) {
      message('\tgrubbs outlier test iteration: ',iter)

      avgstmp <- avgstmp[-which.max(avgstmp$percent.mito),]

      gt <- grubbs(avgstmp$percent.mito)

      iter=iter+1
    }

    mitohiclusts <- as.vector(avgs$seurat_clusters[!(avgs$seurat_clusters %in% avgstmp$seurat_clusters)])

    mitoclustcells <- WhichCells(tmp, idents = mitohiclusts)
    tmp <- tmp[,!(colnames(tmp) %in% mitoclustcells)]

    message('Reclustering without mito clusters\n')
    suppressWarnings(tmp <- SCTransform(tmp, verbose = T))

    tmp <- RunPCA(object = tmp, verbose = F)

    tmp <- FindNeighbors(object = tmp, dims = 1:30, verbose = F)
    tmp <- FindClusters(object = tmp, resolution = 1.0, verbose = T)
  }

  #record params for all clusters...
  #md <- tmp@meta.data
  #avgs <- aggregate(data = md, cbind(nCount_RNA , nFeature_RNA , percent.mito) ~ seurat_clusters, FUN=mean)
  #rm(md)

  print(VlnPlot(tmp, c('percent.mito', 'nCount_RNA', 'nFeature_RNA'), ncol=1, pt.size = 0.1))











  message('Initiating lib-size QC of first-pass clustering\n')


  ########################### filtration loop ####################################


  if(iterativefilter == T) {



    #loop through clusters and get cells with proper library size, based on MAD
    filteredcells <- c()
    params <- list()
    for(clust in levels(tmp$seurat_clusters)){        #print(clust) }

      message('Filtering cluster ', clust)
      tmpclust <- subset(tmp, idents = clust)


      if(ncol(tmpclust) < 100) {
        message('  -Under 100 cells; calculating params, but will not subset')
      }


      ### LIBSIZE FILTER ###

      #two sided, left sided, or False;
      #if false, else statement will just add all cells to filtered cells obj


      if(iterativefilter.libsize ==  'twosided'){

        #calculate median absolute deviation (MAD) of library size (nCount_RNA)
        x <- tmpclust$nCount_RNA
        m <- median(x)

        #get absolute deviation from the median
        abs.dev   <- abs(x - m)

        #get left and right side of absolute deviation
        left.mad  <- median(abs.dev[x<=m])
        right.mad <- median(abs.dev[x>=m])
        two.sided.mad <- c(left.mad, right.mad)

        #for each cell, assign to tail of median
        x.mad <- rep(two.sided.mad[1], length(x))
        x.mad[x > m] <- two.sided.mad[2]

        #for each cell, get distance from the MAD for the cell's assigned tail
        mad.distance <- abs(x - m) / x.mad
        mad.distance[x==m] <- 0

        #predict mad cutoff. for left tail, x < m...
        y=sort(mad.distance[x<m])
        ecpout <- ecp::e.divisive(diff(as.matrix(y)), k = 1, min.size = 2)
        #this step is non-trivial... may have to play around with k value and which estimate to use.


        ecpp = ecpout$estimates[2]
        madmax.dist.nCount_RNA = as.numeric(y[ecpout$estimates[2]])
        ecpplot <- data.frame(x = 1:length(y), y = y) #should have called this ecpdf...

        ecpplot1 <- ggplot(data=ecpplot, aes(x = x, y = y)) +
          geom_point(shape=1, size = 0.1) +
          geom_hline(yintercept = madmax.dist.nCount_RNA, col = 'red', linetype = 'dashed') +
          geom_vline(xintercept = ecpp, col = 'red', linetype = 'dashed')+
          geom_point(mapping = aes(x = ecpp, y = madmax.dist.nCount_RNA),col = 'red') +
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')


        ecpplot2 <- ggplot(data=ecpplot, aes(x = 0, y = y)) +
          geom_violin(fill = 'steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = madmax.dist.nCount_RNA, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')+
          scale_x_discrete(limits=0)
        #scaling x this way so the line is at the same level in the cowplot

        title <- ggdraw() +
          draw_label(
            paste0('Cluster: ', clust, '\nLibsize MAD threshold prediction = ', ecpp),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(ecpplot1, ecpplot2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(y, ecpout, ecpplot, ecpplot1, ecpplot2, ecpp, plotrow, title)


        #select bad cells (outside of MAD cutoff)
        bad <- x[mad.distance > madmax.dist.nCount_RNA]


        #for plotting: get hi and lo bounds for histogram based on hi and lo cells
        ncountshi <- ifelse( test = length(bad[bad > median(x)]) == 0,
                             yes = Inf,
                             no = min(bad[bad > median(x)])
        )

        ncountslo <- ifelse( test = length(bad[bad < median(x)]) == 0,
                             yes = 0,
                             no = max(bad[bad < median(x)])
        )

        #plot cutoffs for library size / UMIs
        g1 <- ggplot(tmpclust@meta.data) +
          geom_histogram(aes(nCount_RNA),
                         color="black", fill = "steelblue",
                         binwidth = 0.05)+
          geom_vline(xintercept = ncountshi, linetype = "dashed", colour = "red")+
          geom_vline(xintercept = ncountslo, linetype = "dashed", colour = "red")+
          geom_vline(xintercept = median(tmpclust@meta.data$nCount_RNA),
                     linetype = "dotted", colour = "red", size = 1.2)+
          scale_x_log10()+
          labs(title = paste0("Library Size per Cell for Cluster", clust),
               x = "Library size (UMI, aka 'nCount'), log10 scale",
               y = "Number of cells" ,
               subtitle = paste0(length(rownames(tmpclust@meta.data)), " cells; median lib size = ",
                                 median(tmpclust@meta.data$nCount_RNA)),
               caption = paste0("cells remaining = ",
                                length(x) - length(bad)
               )
          )+
          theme_linedraw()




        g2 <-  ggplot(tmpclust@meta.data, aes(x = 0, y = nCount_RNA))+
          geom_violin(fill='steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = ncountslo, col = 'red', linetype = 'dashed')+
          geom_hline(yintercept = ncountshi, col = 'red', linetype = 'dashed')+
          geom_hline(yintercept = median(tmpclust@meta.data$nCount_RNA),
                     linetype = "dotted", colour = "red", size = 1.2)+
          theme_linedraw()+
          scale_y_log10()+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          labs(title = paste0("Library Size per Cell for Cluster", clust),
               subtitle = paste0(length(rownames(tmpclust@meta.data)), " cells total",
                                 "; cutoff for percent mito is ", as.character(round(ncountslo, digits = 3))),
               caption = paste0('libsize max MAD = ', round(madmax.dist.nCount_RNA, 3) ,
                                '\nNum cells presubset = ', ncol(tmpclust),
                                '\nNum cells remaining = ',
                                length(x) - length(bad), '\n'))+
          xlab('Cells') + ylab('LibSize, Log10 scale')



        title <- ggdraw() +
          draw_label(
            paste0('Cluster: ', clust, '\nLibsize cutoffs'),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(g1, g2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(g1,g2,title, plotrow)

        #add good cells to iterator list
        filteredclust <- names(x[!names(x) %in% names(bad)] )

        #save params
        #params[[clust]] <- data.frame(libsizemad = madmax.dist.nCount_RNA,
        #                              libsizelo = ncountslo,
        #                              libsizehi = ncountshi)

        #remove objects
        rm(abs.dev, bad, left.mad, right.mad, m, x, x.mad,
           two.sided.mad, mad.distance,
           ncountshi, ncountslo,
           madmax.dist.nCount_RNA)

      } #end two sided iterative lib size filt


      if(iterativefilter.libsize == 'lefttail'){

        #calculate median absolute deviation (MAD) of library size (nCount_RNA)
        x <- tmpclust$nCount_RNA
        m <- median(x)

        #get absolute deviation from the median
        abs.dev   <- abs(x - m)

        left.mad  <- median(abs.dev[x<=m])

        x.mad <- rep(left.mad, length(x))

        #for each cell, get distance from the MAD for the cell's assigned tail
        mad.distance <- abs(x - m) / x.mad
        mad.distance[x==m] <- 0

        #predict mad cutoff. for left tail, x < m...
        y=sort(mad.distance[x<m])
        ecpout <- ecp::e.divisive(diff(as.matrix(y)), k = 1, min.size = 2)
        #this step is non-trivial... may have to play around with k value and which estimate to use.


        ecpp = ecpout$estimates[2]
        madmax.dist.nCount_RNA = as.numeric(y[ecpout$estimates[2]])
        ecpplot <- data.frame(x = 1:length(y), y = y) #should have called this ecpdf...

        ecpplot1 <- ggplot(data=ecpplot, aes(x = x, y = y)) +
          geom_point(shape=1, size = 0.1) +
          geom_hline(yintercept = madmax.dist.nCount_RNA, col = 'red', linetype = 'dashed') +
          geom_vline(xintercept = ecpp, col = 'red', linetype = 'dashed')+
          geom_point(mapping = aes(x = ecpp, y = madmax.dist.nCount_RNA),col = 'red') +
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')


        ecpplot2 <- ggplot(data=ecpplot, aes(x = 0, y = y)) +
          geom_violin(fill = 'steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = madmax.dist.nCount_RNA, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')+
          scale_x_discrete(limits=0)
        #scaling x this way so the line is at the same level in the cowplot

        title <- ggdraw() +
          draw_label(
            paste0('Cluster: ', clust, '\nLibsize MAD threshold prediction = ', ecpp),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(ecpplot1, ecpplot2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(y, ecpout, ecpplot, ecpplot1, ecpplot2, ecpp, plotrow, title)


        #select bad cells (outside of MAD cutoff)
        bad <- x[mad.distance > madmax.dist.nCount_RNA]


        #for plotting: get lo bounds for histogram based on lo cells

        ncountslo <- ifelse( test = length(bad[bad < median(x)]) == 0,
                             yes = 0,
                             no = max(bad[bad < median(x)])
        )

        #plot cutoffs for library size / UMIs
        g1 <- ggplot(tmpclust@meta.data) +
          geom_histogram(aes(nCount_RNA),
                         color="black", fill = "steelblue",
                         binwidth = 0.05)+
          geom_vline(xintercept = ncountslo, linetype = "dashed", colour = "red")+
          geom_vline(xintercept = median(tmpclust@meta.data$nCount_RNA),
                     linetype = "dotted", colour = "red", size = 1.2)+
          scale_x_log10()+
          labs(title = paste0("Library Size per Cell for Cluster", clust),
               x = "Library size (UMI, aka 'nCount'), log10 scale",
               y = "Number of cells" ,
               subtitle = paste0(length(rownames(tmpclust@meta.data)), " cells; median lib size = ",
                                 median(tmpclust@meta.data$nCount_RNA)),
               caption = paste0("cells remaining = ",
                                length(x) - length(bad)
               )
          )+
          theme_linedraw()




        g2 <-  ggplot(tmpclust@meta.data, aes(x = 0, y = nCount_RNA))+
          geom_violin(fill='steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = ncountslo, col = 'red', linetype = 'dashed')+
          geom_hline(yintercept = median(tmpclust@meta.data$nCount_RNA),
                     linetype = "dotted", colour = "red", size = 1.2)+
          theme_linedraw()+
          scale_y_log10()+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          labs(title = paste0("Library Size per Cell for Cluster", clust),
               subtitle = paste0(length(rownames(tmpclust@meta.data)), " cells total",
                                 "; cutoff for percent mito is ", as.character(round(ncountslo, digits = 3))),
               caption = paste0('libsize max MAD = ', round(madmax.dist.nCount_RNA, 3) ,
                                '\nNum cells presubset = ', ncol(tmpclust),
                                '\nNum cells remaining = ',
                                length(x) - length(bad), '\n'))+
          xlab('Cells') + ylab('LibSize, Log10 scale')



        title <- ggdraw() +
          draw_label(
            paste0('Cluster: ', clust, '\nLibsize cutoffs'),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(g1, g2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(g1,g2,title, plotrow)

        #add good cells to iterator list
        filteredclust <- names(x[!names(x) %in% names(bad)] )

        #save params
        # params[[clust]] <- data.frame(libsizemad = madmax.dist.nCount_RNA,
        #                                libsizelo = ncountslo)

        #remove objects
        rm(abs.dev, bad, left.mad, m, x, x.mad,
           mad.distance,  ncountslo,
           madmax.dist.nCount_RNA)

      } #end left iterative lib size filt

      if(iterativefilter.libsize ==  F){
        filteredclust <- colnames(tmp)
        }

      #end iterative lib size filt








      ### MITO FILTER ###

      if(iterativefilter.mito == T){

        ### find mito cutoff high ###
        x <- tmpclust$percent.mito
        m         <- median(x)
        abs.dev   <- abs(x - m)

        right.mad <- median(abs.dev[x>=m])

        x.mad <- rep(right.mad, length(x))

        mad.distance <- abs(x - m) / x.mad
        mad.distance[x==m] <- 0

        #predict mad cutoff. for right tail, x > m
        y=sort(mad.distance[x>m])
        ecpout <- ecp::e.divisive(diff(as.matrix(y)), k = 3, min.size = 2)
        #this step is non-trivial... may have to play around with k value and which estimate to use.


        ecpp = ecpout$estimates[2]
        madmax.dist.percentmito = as.numeric(y[ecpout$estimates[2]])
        ecpplot <- data.frame(x = 1:length(y), y = y) #should have called this ecpdf...

        ecpplot1 <- ggplot(data=ecpplot, aes(x = x, y = y)) +
          geom_point(shape=1, size = 0.1) +
          geom_hline(yintercept = madmax.dist.percentmito, col = 'red', linetype = 'dashed') +
          geom_vline(xintercept = ecpp, col = 'red', linetype = 'dashed')+
          geom_point(mapping = aes(x = ecpp, y = madmax.dist.percentmito),col = 'red') +
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')


        ecpplot2 <- ggplot(data=ecpplot, aes(x = 0, y = y)) +
          geom_violin(fill = 'steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = madmax.dist.percentmito, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          xlab(label = 'Cells') + ylab(label = 'Med. Abs. Dev.')+
          scale_x_discrete(limits=0)
        #scaling x this way so the line is at the same level in the cowplot

        title <- ggdraw() +
          draw_label(
            paste0('Cluster: ', clust, '\nMito MAD threshold prediction = ', ecpp),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(ecpplot1, ecpplot2)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(y, ecpout, ecpplot, ecpplot1, ecpplot2, ecpp, plotrow, title)


        #select bad cells (outside of MAD cutoff)
        bad <- x[mad.distance > madmax.dist.percentmito]



        #for plotting / message output: get mito upper bound
        mitohi <- ifelse( test = length(bad[bad > median(x)]) == 0,
                          yes = Inf,
                          no = min(bad[bad > median(x)])
        )



        m1 <- ggplot(tmpclust@meta.data, aes(x = 0, y = percent.mito))+
          geom_violin(fill='steelblue')+
          geom_jitter(height = 0, width = 0.25, size = 0.1)+
          geom_hline(yintercept = mitohi, col = 'red', linetype = 'dashed')+
          theme_linedraw()+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          labs(title = paste0("percent mito cutoff for cluster ", clust),
               subtitle = paste0(length(rownames(tmpclust@meta.data)), " cells total",
                                 "; cutoff for percent mito is ", as.character(round(mitohi, digits = 3))),
               caption = paste0('Mito max MAD = ', round(madmax.dist.percentmito, 3) ,
                                '\npercent.mito cutoff = ', round(mitohi, digits = 3),
                                '\nNum cells presubset = ', ncol(tmpclust),
                                '\nNum cells remaining = ', length(tmpclust$percent.mito[tmpclust@meta.data$percent.mito < mitohi]), '\n'))+
          xlab('Cells')

        title <- ggdraw() +
          draw_label(
            paste0('Cluster: ', clust, '\nMito Hi Cutoff'),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )

        plotrow <- cowplot::plot_grid(m1)

        print(
          plot_grid(
            title, plotrow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1)
          )
        )

        rm(m1, title, plotrow)

        #save mito filtered cells
        mitofilteredcells <- names(x[!names(x) %in% names(bad)] )

        params[[clust]]$mitomad <- madmax.dist.percentmito
        params[[clust]]$mitohi  <- mitohi

        #remove objects
        rm(abs.dev, bad, right.mad, m, x, x.mad, mad.distance,mitohi, madmax.dist.percentmito)

      } else{mitofilteredcells <- colnames(tmp)} #end iterative mito filt



      #get the good cell barcodes for this cluster
      filteredclust <- intersect(filteredclust, mitofilteredcells)

      #save cell numbers pre and post
      params[[clust]]$prefilt <- ncol(tmpclust)
      params[[clust]]$postfilt <- length(filteredclust)

      #if under 100, do not actually subset.
      if(ncol(tmpclust) > 100) {

        #barcodes of this cluster to global good barcode list
        filteredcells <- c(filteredcells, filteredclust)

      } else{

        filteredcells <- c(filteredcells, colnames(tmpclust))
      }

      rm(filteredclust, mitofilteredcells, tmpclust)


    } #close for clust loop

  } #close iterative filter conditional loop


  #avgs <- cbind(avgs, rbindlist(params))

  ### second pass normalization and clustering, after removing bad cells ###
  #read in
  if(format == 'dir'){
    tmp <- CreateSeuratObject(Read10X(data),
                              project = project,
                              min.features = 200, min.cells = 3)
  }

  if(format == 'h5'){
    tmp <- CreateSeuratObject(Read10X_h5(data),
                              project = project,
                              min.features = 200, min.cells = 3)
  }

  if(format == 'kallisto'){
    tmp <- CreateSeuratObject(res_mat, project = project,
                              min.cells = 3, min.features = 200)
  }

  #subset
  tmp <- tmp[,filteredcells]

  #mito content, add to metadata
  mito.features <- grep(pattern = "^mt-", x = rownames(x = tmp), value = TRUE)
  percent.mito <- Matrix::colSums(
    x = GetAssayData(object = tmp, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = tmp, slot = 'counts')
    )
  tmp[['percent.mito']] <- percent.mito

  rm(percent.mito, mito.features)

  #cluster
  message('\nInitiating second-pass clustering:\n')

  suppressWarnings(tmp <- SCTransform(tmp, verbose = T))

  tmp <- RunPCA(object = tmp, verbose = F)
  tmp <- FindNeighbors(object = tmp, dims = 1:30, verbose = F)
  tmp <- FindClusters(object = tmp, resolution = 1.0, verbose = T)


  message('Initiating DF Parameter Sweep:\n')

  #param sweep...
  sweepres <- paramSweep_v3(seu = tmp, PCs = 1:30, sct = T)

  message('DF Parameter Sweep Completed')

  sweepstats <- summarizeSweep(sweepres)
  bcmvn <- find.pK(sweepstats)

  maxscorepk <- bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric),2]
  maxscorepk <- as.numeric( levels(maxscorepk)[maxscorepk] )





  #homotypic doublet modelling

  ### using 10x table, use linear regression --> important for predicting homotypic / total doublet number
  tamlabscpipeline::dratedf
  #dratedf <- read.delim('/Users/ferrenaa/Documents/tam/scripts/doublets/doubletrate.txt', header = T)
  dratedf[,1] <- as.numeric(sub("%","",dratedf[,1]))/100

  names(dratedf) <- c('MultipletRate', 'CellsLoaded_100%Viability', 'CellsRecovered')

  dbmodel <- lm(MultipletRate ~ CellsRecovered, data = dratedf)

  predicteddoubletrate <- as.numeric((dbmodel$coefficients['CellsRecovered'] * ncol(tmp)) + dbmodel$coefficients[1])

  homotypicprop <- modelHomotypic(tmp$seurat_clusters)
  nexppoi <- round(predicteddoubletrate * length(rownames(tmp@meta.data)))
  nexppoiadj <- round(nexppoi * (1 - homotypicprop))

  #classify doublets
  message('Initiating third-pass clustering for doublet estimation:\n')

  tmp <- suppressWarnings(doubletFinder_v3(seu = tmp, PCs = 1:30, pN = 0.25, pK = maxscorepk, nExp = nexppoi, sct = T))

  ddf <- data.frame(cells = rownames(tmp@meta.data),
                    orig.ident = tmp$orig.ident,
                    classification = tmp@meta.data[,ncol(tmp@meta.data)])

  print( table(ddf$classification ))

  ddf <- ddf[ddf$classification == 'Singlet',]

  ### 4th pass clustering ###

  message('\nInitiating final clustering:\n')

  #read in
  if(format == 'dir'){
    tmp <- CreateSeuratObject(Read10X(data),
                              project = project,
                              min.features = 200, min.cells = 3)
  }

  if(format == 'h5'){
    tmp <- CreateSeuratObject(Read10X_h5(data),
                              project = project,
                              min.features = 200, min.cells = 3)
  }

  if(format == 'kallisto'){
    tmp <- CreateSeuratObject(res_mat, project = project,
                              min.cells = 3, min.features = 200)
  }

  tmp <- tmp[,ddf$cells]

  #mito content, add to metadata
  mito.features <- grep(pattern = "^mt-", x = rownames(x = tmp), value = TRUE)
  percent.mito <- Matrix::colSums(
    x = GetAssayData(object = tmp, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = tmp, slot = 'counts')
    )
  tmp[['percent.mito']] <- percent.mito
  rm(mito.features, percent.mito)

  #cell cycle scoring

  if(cellcycleregression == 'difference' | cellcycleregression == 'total'){
    message('Initating cell cycle regression procedure\n')

    #get mouse homologs
    #mart <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
      #              dataset ='hsapiens_gene_ensembl',
       #             host = 'www.ensembl.org')

    #s.genes <- getBM(attributes = c("mmusculus_homolog_associated_gene_name"),
     #                filters = 'hgnc_symbol',
      #               #values = genes$ensemble.genes,
       #              value = cc.genes[[1]],
        #             mart = mart)[,1]

    #g2m.genes <- getBM(attributes = c("mmusculus_homolog_associated_gene_name"),
     #                  filters = 'hgnc_symbol',
      #                 #values = genes$ensemble.genes,
       #                value = cc.genes[[2]],
        #               mart = mart)[,1]

    s.genes <- tamlabscpipeline::s.genes
    g2m.genes <- tamlabscpipeline::g2m.genes

    tmp <- NormalizeData(tmp, verbose = T)
    tmp <- CellCycleScoring(tmp, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

    rm(s.genes, g2m.genes)
  }

  if(cellcycleregression == 'difference'){
    message('Regressing out S - G2M differences and normalizing\n')

    tmp$CC.Difference <- tmp$S.Score - tmp$G2M.Score
    suppressWarnings(tmp <- SCTransform(tmp, verbose = T, return.only.var.genes = F,
                                        vars.to.regress = c('percent.mito', 'CC.Difference') ) )
  }


  if(cellcycleregression == 'total'){
    message('Regressing out total cell cycle score and normalizing\n')

    suppressWarnings(tmp <- SCTransform(tmp, verbose = T, return.only.var.genes = F,
                                        vars.to.regress = c('percent.mito', "S.Score", "G2M.Score") ) )
  }


  if(cellcycleregression == F){
    message('Calculating cell cycle score, but normalizing without cell cycle regression\n')

    #still calculate cell cycle score!!!
    s.genes <- tamlabscpipeline::s.genes
    g2m.genes <- tamlabscpipeline::g2m.genes

    tmp <- NormalizeData(tmp, verbose = T)
    tmp <- CellCycleScoring(tmp, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

    rm(s.genes, g2m.genes)


    suppressWarnings(tmp <- SCTransform(tmp, verbose = T, return.only.var.genes = F,
                                        vars.to.regress = 'percent.mito') )
  }


  #PCA
  if(is.null( PCAgenelist )) {PCAgenelist <- VariableFeatures(tmp)}

  tmp <- RunPCA(object = tmp, features = PCAgenelist, verbose = F)

  #jackstraw
  if(jackstraw == T){
    message('Initiating Jackstraw Procedure\n')
    tmp <- JackStraw(tmp, dims = 50, verbose = T)
    tmp <- ScoreJackStraw(tmp, dims = 1:50)
    dims <- tmp@reductions$pca@jackstraw$overall.p.values[tmp@reductions$pca@jackstraw$overall.p.values[,2] < 0.05,1]
    message('\nJackstraw Procedure completed.\nRetaining PCs: ', list(dims), '\n')
  }

  message('Initiating final clustering for resolution: ', list(res))



  tmp <- FindNeighbors(object = tmp, dims = dims, verbose = F)
  tmp <- FindClusters(object = tmp, resolution = res, verbose = T)

  message('\nRunning tSNE / UMAP:\n')
  tmp <- RunTSNE(tmp, dims = dims)
  tmp <- RunUMAP(tmp, dims = dims, verbose = F)

  print(  DimPlot(tmp, label = T) )

  message('Clustering procedure complete.')
  return(tmp)

}





# DEGs from each clusters, output compatible with GSEA fxn ---------------
#' Differential expression testing for clusters in a single sample
#'
#' A function to run differentially expressed gene (DEG) testing and save to disk the output in a format compatible with gseapipeline.clusters() function.
#' This function takes a while. To help deal with the risk of interruptions, this function is written in a way that tries to pick up where previously left off if stopped part-way through.
#' @param sobj A seurat object, usually one on which Seurat::FindClusters() has been run.
#' @param idents A string referring to a categorical column of the seurat metadata, such as the output of Seurat::FindClusters(). Defaults to 'seurat_project'.
#' @param test A string connoting which DEG test to perform. See ?Seurat::FindMarkers() for details on options. Defaults to "MAST".
#' @param latent.vars A string or character vector referring to which columns to use a "latent variables", which some of the tests will attempt to regress out the effect of. Some tests (such as MAST) can only perform this for quantitative, continuous variables. Defaults to ('nFeature_RNA', 'nCount_RNA', 'percent.mito').
#' @param outdir a string connoting which directory to save results in. Will create directory if it does not exist. Will not overwrite. Defaults to the following: "ClusterDEG_(Seurat object project name)_(date)_(test)"
#' @return Saves output files to outdir. Will not return anything else besides printing progress reports to the standard out.
#' @examples
#' \dontrun{
#' pdf('qcplots.pdf')
#' deg.acrossclusters(sobj = sobj)
#' dev.off()
#' }
deg.acrossclusters <- function(sobj,
                               idents=NULL,
                               test=NULL,
                               latent.vars=NULL,
                               outdir=NULL){

  if(is.null( idents )) {idents <- "seurat_clusters"}
  if(is.null( test )) {test <- "MAST"}

  if(is.null( latent.vars )) {latent.vars <- c('nFeature_RNA', 'nCount_RNA', 'percent.mito') }

  if(is.null( outdir )) {outdir <- paste0('clusterDEG_', sobj@project.name, '_', Sys.Date(), '_', test) }
  if(!dir.exists(outdir)){dir.create(outdir)}

  #check if input ident is already the active ident; fix if not; will reset later if not.
  if( !(identical(as.vector(sobj@meta.data[,idents]), as.vector(sobj@active.ident)) ) ){
    ditest = T
    di <- sobj@active.ident
    sobj <- SetIdent(sobj, value = sobj@meta.data[,idents])
  }

  #set default assay to SCT, reset to previous later
  da <- DefaultAssay(sobj)
  DefaultAssay(sobj) <- 'SCT'

  for(cluster in levels(sobj@active.ident) ){
    message('\nDEG for cluster ', cluster, ':')

    #test if file exists
    if( !file.exists( paste0(outdir, "/markers_cluster", cluster, ".rds") ) ){

      tmp <- FindMarkers(object = sobj, test.use = test, ident.1 = cluster,
                         assay = 'SCT', slot = 'data', verbose = T,
                         logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf,
                         latent.vars = latent.vars)

      saveRDS(tmp, file = paste0(outdir, "/markers_cluster", i, ".rds"))

    }else{ message('\t', 'Cluster ', cluster, ' already completed ', sprintf('\u2714')) }

  }
  #return to defaults
  DefaultAssay(sobj) <- da

  if( ditest == T ){

    sobj <- SetIdent(sobj, value = di)
  }

}




# GSEA using DEGs from clusters ---------------------------
#' A function for GSEA using the output of differential gene expression testing
gseapipeline.clusters <- function(inputfolder,
                                  pathways,
                                  nperm=NULL,
                                  weightmethod=NULL,
                                  makepdf=NULL,
                                  pdfname=NULL,
                                  filter_nonsig_pathways=NULL){

  set.seed(500)

  if(is.null( nperm )) {nperm <- 10000}
  if(is.null( weightmethod )) {weightmethod <- 'pvalue'}
  if(is.null( makepdf )) { makepdf <- F}
  if(is.null( filter_nonsig_pathways )) {filter_nonsig_pathways <- F}

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
