# Date: Dec 25 2020 
# Packages update manually

BiocManager::install("ade4", lib = "/mnt/sas0/AD/jli2274/R/x86_64-redhat-linux-gnu-library/3.6")

BiocManager::install(c("ade4","annotate", "AnnotationDbi", "ape", "assertthat",
                       "backports", "BH", "Biobase", "BiocGenerics", "BiocParallel", "bit", "bit64", "blob", "boot", "checkmate", "class", "cli", "cluster",
                       "coda", "codetools", "colorspace", "data.table", "DBI", "DelayedArray", "deldir", "DESeq2", "digest", "doParallel", "dplyr", "e1071",
                       "evaluate", "expm", "fit.models", "foreach", "formatR", "Formula", "genefilter", "geneplotter", "generics", "GenomeInfoDb",
                       "GenomeInfoDbData", "GenomicRanges", "ggplot2", "glue", "GO.db", "gtable", "hierfstat", "highr", "Hmisc", "htmlTable", "htmltools",
                       "htmlwidgets", "igraph", "impute", "IRanges", "iterators", "jsonlite", "KernSmooth", "knitr", "labeling", "lambda.r", "lattice",
                       "latticeExtra", "lazyeval", "locfit", "magrittr", "markdown", "MASS", "Matrix", "matrixStats", "mgcv", "mime", "mvtnorm", "nlme", "nnet",
                       "pheatmap", "pillar", "pkgconfig", "pls", "plyr", "preprocessCore", "prettyunits", "R6", "raster", "Rcpp", "RcppArmadillo", "RCurl",
                       "reshape2", "rlang", "robust", "robustbase", "rrBLUP", "rrcov", "RSQLite", "rstudioapi", "S4Vectors", "scales", "segmented", "seqinr",
                       "sf", "shiny", "sp", "spatial", "spData", "spdep", "stringi", "stringr", "SummarizedExperiment", "survival", "tibble", "units", "vctrs",
                       "vegan", "WGCNA", "withr", "XML", "xtable", "XVector", "yaml", "zlibbioc"), lib = "/mnt/sas0/AD/jli2274/R/x86_64-redhat-linux-gnu-library/3.6")
BiocManager::install("MAST", lib = "/mnt/sas0/AD/jli2274/R/x86_64-redhat-linux-gnu-library/3.6")

BiocManager::install(c("nloptr"), lib = "/mnt/sas0/AD/jli2274/R/x86_64-redhat-linux-gnu-library/3.6")
