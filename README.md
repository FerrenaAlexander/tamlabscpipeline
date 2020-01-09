# tamlabscpipeline
A basic package for scRNAseq data analysis.


Installation instructions:
```
install.packages("devtools") #if devtools not already installed
devtools::install_github('apf2139/tamlabscpipeline', build_vignettes = T)
```

Use `browseVignettes("tamlabscpipeline")` and select the HTML option for a quickstart guide, discussion of functions and analysis tips.


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
