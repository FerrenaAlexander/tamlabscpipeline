htodemux <- function(sobjdir,
                     hashouts,
                     pos.quant.sweep=NULL,
                     plotout.pos.quant.sweep=NULL,
                     positive.quantile=NULL,
                     plotout.htoheatmap=NULL,
                     predict.opt.posquant=NULL,
                     subset.demux=NULL,
                     move.predemux=NULL,
                     hto.predemux=NULL,
                     return.classification=NULL) {

  if(is.null( pos.quant.sweep )) {pos.quant.sweep <- T}
  if(is.null( plotout.pos.quant.sweep )) {plotout.pos.quant.sweep <- F}
  if(is.null( predict.opt.posquant )) {predict.opt.posquant <- T}
  if(is.null( plotout.htoheatmap )) {plotout.htoheatmap <- F}
  if(is.null( subset.demux )) {subset.demux <- F}
  if(is.null( move.predemux )) {move.predemux <- F}
  if(is.null( hto.predemux )) {hto.predemux <- 'hto.predemux'}
  if(is.null( return.classification )) {return.classification <- F}

  if(return.classification == T){
    retclass = list()
  }

  for(sobjfile in list.files(sobjdir) ){

    sobjfile <- paste0(sobjdir, '/', sobjfile)

    sobjbasename <- tools::file_path_sans_ext(basename(sobjfile))

    message('\nInitiating ', sobjbasename, '\n')

    htodir <- paste0(hashouts, '/', sobjbasename, '/umi_count/')

    sobj <- readRDS(sobjfile)
    htos <- Read10X(htodir, gene.column = 1)

    if( 'unmapped' %in% rownames(htos) ){
      message('\nRemoving unmapped from HTO matrix\n')
      htos <- htos[rownames(htos) != 'unmapped',]
    }

    #deal with cell barcoedes having number... ugh
    if( any(grepl('-', x = colnames(sobj))) | any( grepl('-', x = colnames(htos))) ){
      cellbcnumed <- T
      addon <- unique(str_sub(colnames(sobj), start=-2))
      if(length(addon) > 1) {stop("more than one cell labels attached to barcodes of cell names / colanems of sample ", sobjbasename)}

      colnames(htos) <- paste0(colnames(htos), addon)
    }


    #not ideal; simply throw out hto cell barcodes not found in filtered UMI matrix...
    jointbcs <- intersect(colnames(sobj), colnames(htos))
    htos <- htos[,jointbcs]
    sobj <- sobj[,jointbcs]


    ### add in HTO data ###
    sobj[['HTO']] <- CreateAssayObject(counts = htos)

    #fix for low counts, suggested by Eleni Mimitou NYGC
    sobj@assays$HTO@counts <- (as.matrix(sobj@assays$HTO@counts) + 1)

    #normalize htos using CLR method
    sobj <- NormalizeData(sobj, assay = 'HTO', normalization.method = 'CLR')


    if(pos.quant.sweep == T){
      ### positive quantile param sweep

      #use 10x cell user guide table... for 10x v2, need to update for v3 chem?
      dratedf <- read.delim('/Users/ferrenaa/Documents/tam/scripts/doublets/doubletrate.txt', header = T)
      dratedf[,1] <- as.numeric(sub("%","",dratedf[,1]))/100

      names(dratedf) <- c('MultipletRate', 'CellsLoaded_100%Viability', 'CellsRecovered')

      #predict number of doublets used num cells recovered
      dbmodel <- lm(MultipletRate ~ CellsRecovered, data = dratedf)

      predicteddoubletrate <- as.numeric((dbmodel$coefficients['CellsRecovered'] * ncol(sobj)) + dbmodel$coefficients[1])
      predicteddoublets <- round(predicteddoubletrate * length(rownames(sobj@meta.data)))


      message('\nInitiating Parameter Sweep for ', sobjbasename, ' :')
      pdb <- list()
      for(param in seq(0.95, 0.999, by = 0.001)) {
        message('\t', 'setting pos.quant param = ', param)
        suppressMessages( sobj <- HTODemux(sobj, assay = 'HTO', positive.quantile = param) )
        pdb[[as.character(param)]] <- as.data.frame( t(as.matrix(table(sobj$HTO_classification.global)))  )
      }

      message('Parameter sweep completed for ', sobjbasename, '\n')

      pdb <- data.table::rbindlist(pdb)
      pdb$positive.quantile <- seq(0.95, 0.999, by = 0.001)

      pl <- tail(pdb[pdb$Doublet >= predicteddoublets,], 1)
      pl <- rbind(pl, head(pdb[pdb$Doublet <= predicteddoublets,], 1))


      if(predict.opt.posquant == T) {

        tryCatch({
          if( any(duplicated(pl[,1])) ){
            predictedposquan = mean(pl$positive.quantile)
          } else {
            f <- approxfun(pl$Doublet, pl$positive.quantile)
            predictedposquan <- f(predicteddoublets)
          }

          message('\nPos.quant prediction succesful, predicted pos.quant = ', predictedposquan, '\n')

        }, error = function(e) {message('\nWarning: Predicted Pos-quant failed.\nSetting to 0.95.\n')
          predictedposquan = 0.95
        }
        ) #end tryCatch block

      } #end predict.opt.posquant block


      if(plotout.pos.quant.sweep==T){


        pdb3 <- data.frame(Cells = c(pdb$Singlet, pdb$Negative, pdb$Doublet),
                           Classification = c(rep('Singlet', 50), rep('Negative', 50), rep('Doublet', 50)),
                           positive.quantile = rep(pdb$positive.quantile, 3)
        )

        if(predict.opt.posquant == F) {

          g2 <- ggplot(data = pdb3, aes(x = positive.quantile, y = Cells, col = Classification))+
            geom_point() +
            geom_line() +
            scale_x_continuous(breaks = pretty(pdb$positive.quantile, n = 50) )+
            #scale_y_log10()+
            #geom_point(aes(x = predictedposquan, y = predicteddoublets), col = 'red', size = 1.5) +
            geom_hline(yintercept = predicteddoublets, col = 'red', size = 0.25, linetype = 'dashed')+
            geom_vline(xintercept= positive.quantile, col = 'red', size = 0.25, linetype = 'dashed')+
            theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))+
            scale_color_brewer(palette = 'Set2', direction = 1)+
            labs(title = sobjbasename,
                 subtitle = paste0('positive.quantile manually set to ', positive.quantile),
                 caption = paste0(predicteddoublets, ' predicted doublets (horizontal dashed line).') )

        }



        if(predict.opt.posquant == T) {

          g2 <- ggplot(data = pdb3, aes(x = positive.quantile, y = Cells, col = Classification))+
            geom_point() +
            geom_line() +
            scale_x_continuous(breaks = pretty(pdb$positive.quantile, n = 50) )+
            #scale_y_log10()+
            #geom_point(aes(x = predictedposquan, y = predicteddoublets), col = 'red', size = 1.5) +
            geom_hline(yintercept = predicteddoublets, col = 'red', size = 0.25, linetype = 'dashed')+
            #geom_vline(xintercept= predictedposquan, col = 'red', size = 0.25, linetype = 'dashed')+
            theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))+
            scale_color_brewer(palette = 'Set2', direction = 1)+
            labs(title = sobjbasename,
                 caption = paste0(predicteddoublets, ' predicted doublets (horizontal dashed line).') )


          if( is.null(predictedposquan) ){
            predictedposquan = 0.95
            g2 <- g2 + #geom_point(aes(x = predictedposquan, y = predicteddoublets), col = 'red', size = 1.5)+
              geom_vline(xintercept= 0.95, col = 'red', size = 0.25, linetype = 'dashed')+
              labs(title = sobjbasename,
                   subtitle = paste0("Optimal pos.quant prediction failed. Set to 0.95."))

          } else{ #else block for catch failure of prediction

            g2 <- g2 + geom_point(aes(x = predictedposquan, y = predicteddoublets), col = 'red', size = 1.5)+
              geom_vline(xintercept= predictedposquan, col = 'red', size = 0.25, linetype = 'dashed')+
              labs(title = sobjbasename,
                   subtitle = paste0("Predicted positive.quantile parameter = ", round(predictedposquan, digits = 5)),
                   caption = paste0(predicteddoublets, ' predicted doublets (horizontal dashed line).\n',
                                    predictedposquan, ' =  predicted pos.quant (vertical dashed line)'))

          }
        }

        print(g2)
      }

      #set pos quant; fix if failed prediciton
      if(predict.opt.posquant == T) {
        if( is.null(predictedposquan) ){ predictedposquan = 0.95}
        positive.quantile = predictedposquan}


    } else { if(is.null( positive.quantile )) { stop('Parameter Sweep set to false, but no positive.quantile paramter provided!') } }


    sobj <- HTODemux(sobj, assay = 'HTO', positive.quantile = positive.quantile)



    if( plotout.htoheatmap == T) { print(HTOHeatmap(sobj, assay = 'HTO', raster = F) + ggtitle(sobjbasename) ) }



    if(subset.demux == T){

      DefaultAssay(sobj) <- 'SCT'

      #remove doublets / negatives
      singletcells <- rownames( sobj@meta.data[sobj@meta.data$HTO_classification.global == 'Singlet',] )
      sobj <- subset(sobj, cells = singletcells)

      #subset good tags based on counts. hash with more than 50 hits retained...
      goodtags <- names( table(sobj$hash.ID)[table(sobj$hash.ID) > 50] )

      #subset and store each isolated tag
      for(hashid in goodtags) {
        hashidcells <- rownames( sobj@meta.data[sobj@meta.data$hash.ID == hashid,] )

        tmp <- subset(sobj, cells = hashidcells)

        vst <- tmp[["SCT"]]@misc$vst.out
        vst$cell_attr <- vst$cell_attr[Cells(tmp),]
        vst$cells_step1 <- intersect(vst$cells_step1, Cells(tmp))
        tmp[["SCT"]]@misc$vst.out <- vst


        tagshort <- str_split(hashid, pattern = '-')[[1]][1]

        tmp$orig.ident.tag <- paste0(tmp$orig.ident, '.', tagshort)


        if(move.predemux == T) {

          if( !file.exists(paste0(hto.predemux, '/', basename(sobjfile))) ) {

            suppressWarnings( dir.create(hto.predemux) )
            file.copy(from = sobjfile, to = paste0(hto.predemux, '/', basename(sobjfile)) ) }

          unlink(sobjfile)
        }

        tmpfile <- saveRDS(tmp, file = paste0(sobjdir, '/', sobjbasename, '_', tagshort, '.rds') )

      }

    }#subset demux

    if(return.classification == T){
      rcx <- data.frame(cells = colnames(sobj),
                        hash.id = sobj$hash.ID,
                        hash.global = sobj$HTO_classification.global)

      retclass[[sobjfile]] <- rcx

    }

  } #end forsobjfile in list.files(outdir)

  if(return.classification == T){
    return(retclass)

  }

} #end function
