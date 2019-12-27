#### Data download: nuclear-encoded mitochondrial gene mutations ####

# load packages
# https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
library(TCGAbiolinks)
library(maftools)
library(tidyverse)

# create directories for output
dir.create("data_summary")
dir.create("maf")

#### assess data availability ####

sink("data_summary/data_info.txt", append = TRUE)
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

# define data download and assessment function

create_maf <- function(cancer){
  # download nucleotide variant data (if data already downloaded, reads in data)
  query_maf <- GDCquery_Maf(cancer, 
                            directory = "GDCdata", 
                            pipelines = "mutect2", 
                            save.csv = TRUE)
  # use mutect2: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#pipeline-descriptions

  # extract barcodes from maf query
  barcodes <- sort(query_maf$Tumor_Sample_Barcode) %>%
    unique
  # download clinical data 
  clinical <- GDCquery_clinic(project = paste0("TCGA-", cancer), type = "clinical")
  # extract barcodes from clinical data
  sub_id <- sort(clinical$submitter_id) %>%
    unique()
  clinical <- rename(clinical, Tumor_Sample_Barcode = submitter_id)
  # create maf object with clinical data (defaults to only somatic mutations)
  maf <- read.maf(query_maf, clinicalData = clinical, useAll = FALSE)
  # save maf data to file
  write.mafSummary(maf, basename = paste0("maf/", noquote(cancer)))

  # sample summary
  getSampleSummary(maf)
  # maf summary plot
  jpeg(paste0("data_summary/", noquote(cancer), "_mafSummary.jpg"))
  plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
  dev.off()
  # oncoplot
  jpeg(paste0("data_summary/", noquote(cancer), "_oncoplot.jpg"))
  oncoplot(maf = maf, top = 10, removeNonMutated = TRUE)
  dev.off()
  # plot transitions and transversions
  jpeg(paste0("data_summary/", noquote(cancer), "_titv.jpg"))
  maf_titv <- titv(maf = maf, plot = FALSE, useSyn = TRUE)
  plotTiTv(res = maf_titv)
  dev.off()
}

# test function
create_maf("LUSC")

#### apply function ####

cancers <- c("LUSC", "PAAD", "BRCA", "PRAD", "READ", "COAD", "OV", "KIRC", "HNSC")


