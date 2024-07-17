rm(list=ls())

# Install necessary libraries
if (!"BiocManager" %in% rownames(installed.packages()))
  install.packages("BiocManager")
if (!"TCGAWorkflow" %in% rownames(installed.packages()))
  BiocManager::install("TCGAWorkflow")
if (!"TCGAWorkflowData" %in% rownames(installed.packages()))
  BiocManager::install("TCGAWorkflowData")
if (!"TCGAbiolinks" %in% rownames(installed.packages()))
  BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(TCGAWorkflow)
library(TCGAWorkflowData)
library(DT)
library(tidyverse)

# Download and prepare clinical data
query_clinical <- GDCquery(
  project = "TCGA-LGG",
  data.format = "bcr xml",
  data.category = "Clinical"
)

GDCdownload(query_clinical)
clinical_data <- GDCprepare_clinic(
  query = query_clinical, 
  clinical.info = "patient"
  )

datatable(
  data = clinical_data, 
  options = list(scrollX = TRUE, keys = TRUE), 
  rownames = FALSE
)

# Extract treatment information
clinical_drug <- GDCprepare_clinic(
  query = query_clinical, 
  clinical.info = "drug"
)

clinical_drug |>
  datatable(
    options = list(scrollX = TRUE, keys = TRUE), 
    rownames = FALSE
  ) 

# Note: There is a need to determine how to handle multiple treatment cycles:
# Should a composite therapy be considered as advanced if it features 
# radiotherapy or molecularly targeted (even at a later moment), 
# or should only the first treatment be taken into account?
# eg TCGA-DH-5140

# Download and prepare protein expression data
query_protein <- GDCquery(
  project = "TCGA-LGG", 
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification", 
  sample.type = "Primary Tumor"
)

GDCdownload(query_protein)
protein_data <- GDCprepare(query_protein)

datatable(
  data = protein_data, 
  options = list(scrollX = TRUE, keys = TRUE), 
  rownames = FALSE
)

# Protein data name is a bit different from clinical data bcr_patient_barcode: 
# I think simply removing portion and analyte is safe here: 
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
protein_datanew <- protein_data[,-c(1:4)] %>%
  rename_with(~ str_sub(., 1, -5), -1)
## > length(colnames(protein_datanew))
## [1] 430
## > length(unique(colnames(protein_datanew)))
## [1] 430

# Transpose the tibble keeping the names
data_matrix <- t(as.matrix(protein_datanew[,-1]))
colnames(data_matrix) <- protein_datanew$peptide_target 
protein_data_new <- data.frame(data_matrix) %>%
  rownames_to_column(var = "bcr_patient_barcode")

# Save datasets to CSV files
write.csv(clinical_data, "clinical_data.csv", row.names = FALSE)
write.csv(clinical_drug, "clinical_drug.csv", row.names = FALSE)
write.csv(protein_data_new, "protein_data.csv", row.names = FALSE)

session_info <- sessionInfo()
writeLines(capture.output(session_info), "session_info.txt")
