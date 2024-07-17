load(file = "Data/LGGdata.rda")
clinical_data <- read.csv("Data/clinical_data.csv")
clinical_drug <- read.csv("Data/clinical_drug.csv")
protein_data <- read.csv("Data/protein_data.csv")

#remove na columns
protein_data <- protein_data[ , colSums(is.na(protein_data)) == 0]#plot histogram for protein data
protein_data$bcr_patient_barcode<- NULL
library(ggplot2)
library(gridExtra)
library(patchwork)
library(reshape2)

ggplot(melt(protein_data[, 301:400]), aes(x = value)) + 
    facet_wrap(~ variable, scales = "free", ncol = 10) + 
    geom_histogram(binwidth = .5)
ggsave('hist_4.pdf')


library("Rcpp")
library("pracma")
library("tidyverse")
library("ggpubr")
library("uwot")
library("irlba")
library("mclust")
library("mcclust")
library("lpSolve")

Sys.setenv(CC = "/opt/homebrew/opt/llvm/bin/clang")
Sys.setenv(CXX = "/opt/homebrew/opt/llvm/clang++")
Sys.setenv(CXXFLAGS = "-fopenmp -I/usr/local/opt/libomp/include")
Sys.setenv(LDFLAGS = "-L/usr/local/opt/libomp/lib -lomp")

sourceCpp("R/test.cpp")
test_openmp()
sourceCpp("R/DL_linear_split_merge_package.cpp",verbose = TRUE)
Rcpp::evalCpp("1 + 13 - 1")
