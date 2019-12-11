#### Data download: nuclear-encoded mitochondrial gene mutations ####

# load packages
# https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
library(TCGAbiolinks)
library(maftools)
library(tidyverse)

#### assess data availability ####

sink("data_info.txt", append = TRUE)
sessionInfo()
print("LUSC, lung (NSCLC/non-small cell lung cancer, SCC/squamous-cell carcinoma is the histologic type)")
TCGAbiolinks:::getProjectSummary("TCGA-LUSC")
print("PAAD, pancreatic")
TCGAbiolinks:::getProjectSummary("TCGA-PAAD")
print("BRCA, breast")
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")
print("PRAD, prostate")
TCGAbiolinks:::getProjectSummary("TCGA-PRAD")
print("READ, rectum adenocarcinoma")
TCGAbiolinks:::getProjectSummary("TCGA-READ")
print("COAD, colon adenocarcinoma")
TCGAbiolinks:::getProjectSummary("TCGA-COAD")
print("OV, ovary")
TCGAbiolinks:::getProjectSummary("TCGA-OV")
print("KIRC, kidney")
TCGAbiolinks:::getProjectSummary("TCGA-KIRC")
print("HNSC, head and neck")
TCGAbiolinks:::getProjectSummary("TCGA-HNSC")
sink()

#### obtain variant calling + clinical data ####

# function to download and format variant + clinical data
create_maf <- function(cancer){
  # download nucleotide variant data
  query_maf <- GDCquery_Maf(paste(cancer), 
                            directory = "GDCdata", 
                            pipelines = "mutect2", save.csv = TRUE)
  # use mutect2: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#pipeline-descriptions 
  
  # extract barcodes from maf query
  barcodes <- sort(query_maf$Tumor_Sample_Barcode) %>%
    unique
  #str_trunc(barcodes, 12, side = "right", ellipsis = "")
  
  # download clinical data 
  clinical <- GDCquery_clinic(project = paste("TCGA", cancer, sep = "-"), 
                              type = "clinical")
  
  # extract barcodes from clinical data
  sub_id <- sort(clinical$submitter_id) %>%
    unique()
  clinical <- rename(clinical, Tumor_Sample_Barcode = submitter_id)
  
  # create maf object with clinical data
  maf <- read.maf(query_maf, clinicalData = clinical, useAll = FALSE)
  # defaults to only somatic mutations
  
  # save maf data to file
  write.mafSummary(maf, basename = paste(cancer))
  # sample summary
  return(getSampleSummary(maf))
}

#### canned visualizations ####

# maf summary
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
# oncoplot
oncoplot(maf = maf, top = 10, removeNonMutated = TRUE)
# transitions and transversions
maf_titv <- titv(maf = maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = maf_titv)
