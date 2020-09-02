# tamlabscpipeline
A basic package for scRNAseq data analysis.

Update (September 2 2020) - this package will soon be updated to support new versions of R and dependency packages. This is in the works and will hopefully be done in the coming weeks.


This pipeline was designed to automatically perform various QC steps such as mitochondrial content cutoffs and doublet detection.

A big motivation for this package was to develop a QC approach that both accurately removes low quality cells in the context of very high levels of heterogeneity and cell-type specific differences in mito-content and even library size. Thus, the main philosophy of the QC approach here is to first cluster, then perform analyze quality on each cluster, with the aim of more accurately applying QC thresholds than simple global cutoffs.

The following describes the default pipeline for individual sample QC and processing, as implemented in the `seuratpipeline()` function:
The path to Cellranger output in the H5 or directory formats is specified. Kallisto output is also supported but has some extra requirements (see documentation).


* Mito cutoff and library size outlier cutoffs are applied after an initial round of clustering. First, "mito clusters" (driven entirely by mito content) are identified using the Grubb's test for outliers and removed, as implemented in the "outliers" package. 
* Next, the data is reclustered. In each remaining cluster, the distribution of mitochondrial content is assessed. Cutoffs for high mitochondrial content are made based on median absolute deviation. The upper threshold for median absolute deviation for each cluster is decided based on "changepoint analysis" (via the package ecm), based on the fact that most cells have low mito content, but often a few cells in each cluster have high amounts of mito content.
* A similar analysis is performed for library size, but this analysis is performed in a two-tail manner and does not include grubbs test or whole cluster removal.
* For both mito and lib-size, if each cluster is under 100 cells, the QC metrics are calculated, but cells are not removed.
Next, the data is reclustered after exlcuding the cells called as outliers in these ways. Then, doublet detection (via DoubletFinder) is applied.
* After this, a final clustering is undertaken, and a "clean and processed" Seurat object is returned.

Other functions included in this package include convenience functions for other common scRNAseq analysis methods including:

- GSEA implemented as a wrapper around the FGSEA package. Also includes Mouse orthologs for the MSIGDB Hallmarks set, as provided by the msigdbr package.
- integration, wrapper around Seurat
- differential expression wrapper around Seurat::FindMarkers() function
- HTO demultiplexing, a wrapper around Seurat's HTODemux() package that also includes a function for parameter sweep of the "positive.quantile" parameter for optimal classification




Installation instructions:
```
install.packages("devtools") #if devtools not already installed
devtools::install_github('apf2139/tamlabscpipeline', build_vignettes = T)
```

Use `browseVignettes("tamlabscpipeline")` and select the HTML option for a discussion of functions and analysis tips.


Quickstart:

```
seurat_obj <- tamlabscpipeline::seuratpipeline(data = "path_to_cellranger_output.h5", 
                                               format = "h5")
```

Currently, single-sample pipeline is up and running. Functions are mostly set for integration (multi-sample comparison) pipeline, but documentation / vignettes are not ready yet. Upcoming updates in terms of code include HTO Demux; smarter doublet calling if Hashing was used; parallelization; better Windows compatibility (almost everything currently works); and better compatibility for human data (almost everything already currently works). (AF, 2020.01.08)

Dependencies are being tweaked / tested. but some required ones along with their versions include:

* Seurat vers 3.1.2
* sctransform vers 0.2.1
* ecp vers 3.1.2
* tidyverse vers 1.3.0
* cowplot vers 1.0.0
* DoubletFinder vers 2.0.2
* MAST (suggested) vers 1.12.0
* fgsea vers 1.12.0



```
install.packages("tidyverse")
install.packages('Seurat')
devtools::install_github(repo = 'ChristophH/sctransform')
install.packages("ecp")
install.packages("cowplot")
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')

#BioconductorPkgs
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MAST")
BiocManager::install("fgsea")
```


A big thank you to to the developers of these packages, especially the Satija lab (Seurat); Gottardo lab (MAST); Chris McGinnis and the Gartner lab (DoubletFinder); and Alexey Sergushichev (fgsea). Also, thanks to Peer lab at MSKCC; Zheng lab at Einstein (Dr.s Deyou Zheng and Yang Liu); and all the members of the Tammela lab.


Please enjoy this nice dolphin.

<img src="vignettes/embed/dolphins.jpg" width="300" height="250">
