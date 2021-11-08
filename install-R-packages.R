#!/usr/bin/env Rscript

.libPaths("/opt/Rlibs/")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")

BiocManager::install("SWATH2stats", ask=FALSE)
BiocManager::install("PECA", ask=FALSE)
install.packages("tidyr",repos = "http://cran.us.r-project.org")
install.packages("argparse",repos = "http://cran.us.r-project.org")
install.packages("corrplot",repos = "http://cran.us.r-project.org")

